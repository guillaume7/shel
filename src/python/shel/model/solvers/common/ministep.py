"""
Minimal explicit one-step updater bundling basic tendencies (Phase 5 prototype).

Included tendencies:
- Pressure gradient (from eta on T to faces U,V)
- Linear bottom drag
- Laplacian viscosity

Excluded for now: advection, Coriolis, BC strategies (we enforce closed box by
zeroing normal-flow boundary faces after update).
"""
from __future__ import annotations

import numpy as np

from shel.model.solvers.momentum.pressure import pressure_gradient
from shel.model.solvers.momentum.friction import bottom_drag_tendency
from shel.model.solvers.momentum.diffusion import viscous_tendency
from shel.model.solvers.waterlevel.continuity import update_free_surface
from shel.model.solvers.momentum.advection import advect_momentum
from shel.model.solvers.common.stencils import avg_x_t_to_u, avg_y_t_to_v
from shel.model.boundary_conditions import get_bc

Array = np.ndarray


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
) -> tuple[Array, Array, Array]:
    """Advance (eta, U, V) by one explicit Euler step with simple physics.

    Returns (eta_next, U_next, V_next). Closed-box boundary enforced by setting
    outer U and V faces to zero after the update.
    """
    # Pressure gradient accelerations on faces
    PG_u, PG_v = pressure_gradient(eta, g, dx, dy)

    # Replace NaNs at boundaries (non-computable faces) with 0 to keep closed box
    if np.isnan(PG_u).any():
        PG_u = np.where(np.isnan(PG_u), 0.0, PG_u)
    if np.isnan(PG_v).any():
        PG_v = np.where(np.isnan(PG_v), 0.0, PG_v)

    # Drag and diffusion tendencies on faces
    dU_drag, dV_drag = bottom_drag_tendency(U, V, r)
    dU_visc, dV_visc = viscous_tendency(U, V, nu, dx, dy)

    # Centered momentum advection (optional)
    if enable_advection:
        Au, Av = advect_momentum(U, V, dx, dy)
    else:
        Au = np.zeros_like(U)
        Av = np.zeros_like(V)

    # Optional Coriolis tendencies: du/dt = f * v, dv/dt = -f * u
    if enable_coriolis and f is not None:
        # Interpolate f to faces
        f_u = avg_x_t_to_u(f)
        f_v = avg_y_t_to_v(f)

        # Interpolate cross velocities to U and V
        # v_at_u: average V at (j,j+1) x (i-1,i)
        ny_p1, nx = V.shape
        ny = ny_p1 - 1
        v_at_u = np.zeros_like(U)
        v_at_u[:, 1:nx] = 0.25 * (
            V[0:ny, 0:nx-1] + V[1:ny+1, 0:nx-1] + V[0:ny, 1:nx] + V[1:ny+1, 1:nx]
        )
        # u_at_v: average U at (j-1,j) x (i,i+1)
        ny_u, nx_p1 = U.shape
        u_at_v = np.zeros_like(V)
        u_at_v[1:ny_u, :] = 0.25 * (
            U[0:ny_u-1, 0:nx_p1-1]
            + U[0:ny_u-1, 1:nx_p1]
            + U[1:ny_u, 0:nx_p1-1]
            + U[1:ny_u, 1:nx_p1]
        )
        dU_cor = f_u * v_at_u
        dV_cor = -f_v * u_at_v
    else:
        dU_cor = 0.0
        dV_cor = 0.0

    # Combine and update velocities
    dU = PG_u + dU_drag + dU_visc + Au + dU_cor
    dV = PG_v + dV_drag + dV_visc + Av + dV_cor
    U_next = U + dt * dU
    V_next = V + dt * dV

    # Apply boundary conditions via strategy layer (default closed)
    m_cls, _ = get_bc(bc_type)
    if m_cls is None:
        m_cls, _ = get_bc("closed")
    if m_cls is not None:
        m_cls().apply_uniform(U_next, V_next)

    # Update free surface using flux divergence of updated velocities
    eta_next = update_free_surface(eta, H, U_next, V_next, dt, dx, dy)

    return eta_next, U_next, V_next


__all__ = ["explicit_step"]
