"""
Conservative explicit one-step updater for the barotropic SWE on an Arakawa C-grid.

Matches the MATLAB ``ComputeLeapfrog`` / ``ComputeSpaceU_CS`` / ``ComputeTimeU_FT``
implementation in ``model_handles.m``.

The solver advances the *conservative* (H*u) form:

    1. Continuity RHS from **current** (u, v)  →  eta_new = eta_old + dt * RHSeta
    2. H_new = eta_new + d
    3. Momentum RHS (pressure, advection, diffusion, Coriolis, friction) — all
       terms are in flux form (multiplied by H where appropriate).
    4. u_new = mask_u * (u_old * H_old_u + dt * RHSu) / H_new_u

When called from the leapfrog stepper, ``dt`` is ``2*dt`` (centered leap).
For a plain Euler startup step, ``dt`` is the physical time step and
``eta_old = eta``, ``U_old = U``, etc.
"""

from __future__ import annotations

import logging

import numpy as np

from shel.model.boundary_conditions import get_bc
from shel.model.solvers.common.stencils import avg_x_t_to_u, avg_y_t_to_v
from shel.model.solvers.momentum.diffusion import viscous_tendency
from shel.model.solvers.momentum.friction import bottom_drag_tendency
from shel.model.solvers.momentum.pressure import pressure_gradient
from shel.model.solvers.waterlevel.continuity import continuity_rhs

logger = logging.getLogger(__name__)

Array = np.ndarray


# ---------------------------------------------------------------------------
# Interpolation helpers matching MATLAB fouraverage_u / twoaverage_u
# ---------------------------------------------------------------------------


def _interp_v_to_u(V: Array) -> Array:
    """4-point average of V (ny+1, nx) onto U locations (ny, nx+1).

    Interior columns 1..nx-1; boundary columns 0 and nx are zero.
    Matches MATLAB ``fouraverage_u``.
    """
    ny_p1, nx = V.shape
    ny = ny_p1 - 1
    out = np.zeros((ny, nx + 1))
    out[:, 1:nx] = 0.25 * (
        V[0:ny, 0 : nx - 1]
        + V[1 : ny + 1, 0 : nx - 1]
        + V[0:ny, 1:nx]
        + V[1 : ny + 1, 1:nx]
    )
    return out


def _interp_u_to_v(U: Array) -> Array:
    """4-point average of U (ny, nx+1) onto V locations (ny+1, nx).

    Interior rows 1..ny-1; boundary rows 0 and ny are zero.
    """
    ny, nx_p1 = U.shape
    nx = nx_p1 - 1
    out = np.zeros((ny + 1, nx))
    out[1:ny, :] = 0.25 * (
        U[0 : ny - 1, 0:nx]
        + U[0 : ny - 1, 1 : nx + 1]
        + U[1:ny, 0:nx]
        + U[1:ny, 1 : nx + 1]
    )
    return out


# ---------------------------------------------------------------------------
# Flux-form advection matching MATLAB ComputeSpaceU_CS
# ---------------------------------------------------------------------------


def _advection_u(U: Array, V: Array, H: Array, dx: float, dy: float) -> Array:
    """Centered flux-form advection for U (ny, nx+1).

    MATLAB (1-based i=row=x, j=col=y → Python row=y, col=x):
      Along-x: -0.25/dx * [mask_u(i+1,j)*(u(i+1)+u(i))^2*H(i) - mask_u(i-1,j)*(u(i)+u(i-1))^2*H(i-1)]
      Along-y: -0.0625/dy * [4H*v*(u_j+1 + u_j) - 4H*v*(u_j + u_j-1)] cross terms
    """
    ny, nxp1 = U.shape
    nx = nxp1 - 1
    out = np.zeros_like(U)

    # Along-x advection at interior U faces (cols 1..nx-1)
    # East face: (U[:,c+1] + U[:,c])^2 * H[:,c]
    # West face: (U[:,c] + U[:,c-1])^2 * H[:,c-1]
    c = slice(1, nx)
    east = (U[:, 2:nxp1] + U[:, 1:nx]) ** 2 * H[:, 1:nx]
    west = (U[:, 1:nx] + U[:, 0 : nx - 1]) ** 2 * H[:, 0 : nx - 1]
    out[:, c] = -0.25 / dx * (east - west)

    # Along-y advection at interior U faces (rows 1..ny-2, cols 1..nx-1)
    if ny > 2 and nx > 1:
        r = slice(1, ny - 1)
        # North face terms (row+1 in y-direction)
        # 4-point H average around the north face of U cell:
        H_n = (
            H[2:ny, 1:nx]
            + H[1 : ny - 1, 1:nx]
            + H[2:ny, 0 : nx - 1]
            + H[1 : ny - 1, 0 : nx - 1]
        )
        Vn = V[2:ny, 1:nx] + V[2:ny, 0 : nx - 1]
        Un = U[2:ny, 1:nx] + U[1 : ny - 1, 1:nx]
        # South face terms
        H_s = (
            H[1 : ny - 1, 1:nx]
            + H[0 : ny - 2, 1:nx]
            + H[1 : ny - 1, 0 : nx - 1]
            + H[0 : ny - 2, 0 : nx - 1]
        )
        Vs = V[1 : ny - 1, 1:nx] + V[1 : ny - 1, 0 : nx - 1]
        Us = U[1 : ny - 1, 1:nx] + U[0 : ny - 2, 1:nx]
        out[r, 1:nx] += -0.0625 / dy * (H_n * Vn * Un - H_s * Vs * Us)

    return out


def _advection_v(U: Array, V: Array, H: Array, dx: float, dy: float) -> Array:
    """Centered flux-form advection for V (ny+1, nx).

    Mirrors _advection_u with x↔y transpose logic matching MATLAB's
    ComputeSpaceU_CS called with transposed arguments for V.
    """
    nyp1, nx = V.shape
    ny = nyp1 - 1
    out = np.zeros_like(V)

    # Along-y advection at interior V faces (rows 1..ny-1)
    r = slice(1, ny)
    north = (V[2:nyp1, :] + V[1:ny, :]) ** 2 * H[1:ny, :]
    south = (V[1:ny, :] + V[0 : ny - 1, :]) ** 2 * H[0 : ny - 1, :]
    out[r, :] = -0.25 / dy * (north - south)

    # Along-x advection at interior V faces (rows 1..ny-1, cols 1..nx-2)
    if nx > 2 and ny > 1:
        c = slice(1, nx - 1)
        H_e = (
            H[1:ny, 2:nx]
            + H[1:ny, 1 : nx - 1]
            + H[0 : ny - 1, 2:nx]
            + H[0 : ny - 1, 1 : nx - 1]
        )
        Ue = U[1:ny, 2:nx] + U[0 : ny - 1, 2:nx]
        Ve = V[1:ny, 2:nx] + V[1:ny, 1 : nx - 1]
        H_w = (
            H[1:ny, 1 : nx - 1]
            + H[1:ny, 0 : nx - 2]
            + H[0 : ny - 1, 1 : nx - 1]
            + H[0 : ny - 1, 0 : nx - 2]
        )
        Uw = U[1:ny, 1 : nx - 1] + U[0 : ny - 1, 1 : nx - 1]
        Vw = V[1:ny, 1 : nx - 1] + V[1:ny, 0 : nx - 2]
        out[r, c] += -0.0625 / dx * (H_e * Ue * Ve - H_w * Uw * Vw)

    return out


# ---------------------------------------------------------------------------
# Conservative time step for velocity matching MATLAB ComputeTimeU_FT
# ---------------------------------------------------------------------------


def _time_step_velocity_u(
    U_old: Array,
    RHSu: Array,
    H_old: Array,
    H_new: Array,
    dt: float,
) -> Array:
    """Conservative time step: u_new = (u_old * H_old_u + dt * RHSu) / H_new_u.

    H_old, H_new are on T-grid (ny, nx). We average them to U-faces.
    """
    H_old_u = avg_x_t_to_u(H_old)
    H_new_u = avg_x_t_to_u(H_new)
    # Avoid division by zero at boundaries
    H_new_u = np.where(np.abs(H_new_u) < 1e-12, 1e-12, H_new_u)
    U_new = (U_old * H_old_u + dt * RHSu) / H_new_u
    return U_new


def _time_step_velocity_v(
    V_old: Array,
    RHSv: Array,
    H_old: Array,
    H_new: Array,
    dt: float,
) -> Array:
    """Conservative time step for V (ny+1, nx)."""
    H_old_v = avg_y_t_to_v(H_old)
    H_new_v = avg_y_t_to_v(H_new)
    H_new_v = np.where(np.abs(H_new_v) < 1e-12, 1e-12, H_new_v)
    V_new = (V_old * H_old_v + dt * RHSv) / H_new_v
    return V_new


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def explicit_step(
    eta: Array,
    H: Array,
    U: Array,
    V: Array,
    *,
    dt: float,
    dx: float,
    dy: float,
    g: float = 9.81,
    r: float = 0.0,
    nu: float = 0.0,
    enable_advection: bool = True,
    f: Array | None = None,
    enable_coriolis: bool = False,
    bc_type: str = "closed",
    # --- New conservative parameters (optional for backward compat) ---
    d: Array | None = None,
    eta_old: Array | None = None,
    H_old: Array | None = None,
    U_old: Array | None = None,
    V_old: Array | None = None,
) -> tuple[Array, Array, Array, Array]:
    """Advance (eta, U, V) by one step in conservative flux form.

    When the "old" arrays are not provided (backward-compatible Euler mode),
    they default to the current arrays and the step reduces to a forward-Euler
    update — suitable for a startup step or simple Euler integration.

    For a leapfrog call, pass ``eta_old/U_old/V_old`` as n-1 arrays and
    ``eta/U/V`` as n arrays; set ``dt = 2*physical_dt``.

    Returns (eta_new, U_new, V_new, H_new).
    """
    # Default old-level arrays to current (Euler mode)
    if eta_old is None:
        eta_old = eta
    if H_old is None:
        H_old = H
    if U_old is None:
        U_old = U
    if V_old is None:
        V_old = V
    if d is None:
        # Infer bathymetry from H and eta: d = H - eta
        d = H - eta

    # -------------------------------------------------------------------
    # 1. Continuity: compute RHS from CURRENT (u_n, v_n)
    # -------------------------------------------------------------------
    RHSeta = continuity_rhs(H, U, V, dx, dy)
    eta_new = eta_old + dt * RHSeta

    # -------------------------------------------------------------------
    # 2. Recompute H_new = eta_new + d
    # -------------------------------------------------------------------
    H_new = eta_new + d

    # -------------------------------------------------------------------
    # 3. Momentum tendencies in conservative form (all × H)
    # -------------------------------------------------------------------

    # Pressure gradient: -g * H_u * d_eta/dx  (conservative)
    PG_u, PG_v = pressure_gradient(eta, H, g, dx, dy)
    # Replace NaN at boundary faces with 0
    PG_u = np.where(np.isnan(PG_u), 0.0, PG_u)
    PG_v = np.where(np.isnan(PG_v), 0.0, PG_v)

    # Diffusion (evaluated at old time level, conservative form)
    dU_visc, dV_visc = viscous_tendency(U_old, V_old, H_old, nu, dx, dy)

    # Linear bottom drag (on current velocity, not conservative — acts as damping)
    dU_drag, dV_drag = bottom_drag_tendency(U, V, r)

    # Advection (flux form, conservative)
    if enable_advection:
        Au = _advection_u(U, V, H, dx, dy)
        Av = _advection_v(U, V, H, dx, dy)
    else:
        Au = np.zeros_like(U)
        Av = np.zeros_like(V)

    # Coriolis: f * v_at_u * H_u (conservative form)
    if enable_coriolis and f is not None:
        H_u = avg_x_t_to_u(H)
        H_v = avg_y_t_to_v(H)
        v_at_u = _interp_v_to_u(V)
        u_at_v = _interp_u_to_v(U)
        # Interpolate f from T-grid to U/V faces if it's an array
        if np.ndim(f) == 2:
            f_u = avg_x_t_to_u(f)
            f_v = avg_y_t_to_v(f)
        else:
            f_u = f
            f_v = f
        dU_cor = f_u * v_at_u * H_u
        dV_cor = -f_v * u_at_v * H_v
    else:
        dU_cor = 0.0
        dV_cor = 0.0

    # Combine RHS for momentum
    RHSu = PG_u + Au + dU_visc + dU_cor
    RHSv = PG_v + Av + dV_visc + dV_cor

    # -------------------------------------------------------------------
    # 4. Conservative time step for velocity
    # -------------------------------------------------------------------
    U_new = _time_step_velocity_u(U_old, RHSu, H_old, H_new, dt)
    V_new = _time_step_velocity_v(V_old, RHSv, H_old, H_new, dt)

    # Add drag as a damping correction (not in conservative H*u form)
    # This is a simplification; MATLAB applies it inside the RHS with safeguards.
    if r != 0.0:
        U_new = U_new + dt * dU_drag / np.where(
            np.abs(avg_x_t_to_u(H_new)) < 1e-12, 1e-12, avg_x_t_to_u(H_new)
        )
        V_new = V_new + dt * dV_drag / np.where(
            np.abs(avg_y_t_to_v(H_new)) < 1e-12, 1e-12, avg_y_t_to_v(H_new)
        )

    # -------------------------------------------------------------------
    # 5. Apply boundary conditions
    # -------------------------------------------------------------------
    m_cls, _ = get_bc(bc_type)
    if m_cls is None:
        m_cls, _ = get_bc("closed")
    if m_cls is not None:
        m_cls().apply_uniform(U_new, V_new)

    return eta_new, U_new, V_new, H_new


__all__ = ["explicit_step"]
