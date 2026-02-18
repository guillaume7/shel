"""Viscous diffusion tendency for U and V on a C-grid (conservative form).

Matches MATLAB ``ComputeSpaceU_CS`` diffusion:

    ν * d/dx [ H_old * d(u_old)/dx ] + ν * d/dy [ H_old * d(u_old)/dy ]

Evaluated at the **previous** time level (u_old, H_old) for leapfrog stability.
Returns the tendency contribution to the conservative (Hu) equation.
"""

from __future__ import annotations

import numpy as np

from ..common.stencils import avg_x_t_to_u, avg_y_t_to_v

Array = np.ndarray


def viscous_tendency(
    U_old: Array,
    V_old: Array,
    H_old: Array,
    nu: float,
    dx: float,
    dy: float,
) -> tuple[Array, Array]:
    """Compute ν ∇·(H_old ∇ u_old) on U and V grids (conservative form).

    Parameters
    ----------
    U_old : (ny, nx+1)  velocity at previous time level
    V_old : (ny+1, nx)  velocity at previous time level
    H_old : (ny, nx)    water-column height at previous time level
    nu    : kinematic viscosity
    dx, dy : grid spacings
    """
    dU = _conservative_laplacian_u(U_old, H_old, nu, dx, dy)
    dV = _conservative_laplacian_v(V_old, H_old, nu, dx, dy)
    return dU, dV


def _conservative_laplacian_u(
    U: Array, H: Array, nu: float, dx: float, dy: float
) -> Array:
    """Conservative Laplacian for U field (ny, nx+1).

    MATLAB reference (interior 2:M, 2:N-1 in 1-based, x-row, y-col):
        +Nu * mask_u(i+1,j) * H_old(i,j) * (u_old(i+1,j) - u_old(i,j)) / dx
        -Nu * mask_u(i-1,j) * H_old(i-1,j) * (u_old(i,j) - u_old(i-1,j)) / dx
        all / dx

    Python (row=y, col=x): U is (ny, nx+1), H is (ny, nx).
    Interior U columns: 1..nx-1. Interior U rows: 1..ny-2.
    """
    ny, nxp1 = U.shape
    nx = nxp1 - 1
    out = np.zeros_like(U)

    # Diffusion along x (columns of U)
    # At U[:,c] with c in 1..nx-1, neighbours are U[:,c+1] and U[:,c-1].
    # H at T-cell c corresponds to east of U[:,c]: H[:,c], H[:,c-1] is west.
    sl = slice(1, nx)  # interior U columns
    diff_x = np.zeros_like(U)
    # East flux: H[:,c] * (U[:,c+1] - U[:,c]) / dx  for c in 1..nx-1
    # Need c+1 <= nx, so c <= nx-1 is fine since U has nx+1 cols.
    diff_x[:, sl] = (
        nu * H[:, 1:nx] * (U[:, 2:nxp1] - U[:, 1:nx]) / dx
        - nu * H[:, 0 : nx - 1] * (U[:, 1:nx] - U[:, 0 : nx - 1]) / dx
    ) / dx

    # Diffusion along y (rows of U)
    # At U[r,c] neighbours are U[r+1,c] and U[r-1,c].
    diff_y = np.zeros_like(U)
    if ny > 2:
        rsl = slice(1, ny - 1)
        diff_y[rsl, :] = (
            nu
            * H[1 : ny - 1, :].mean(axis=1, keepdims=True)
            * np.ones((1, nxp1))
            * (U[2:ny, :] - U[1 : ny - 1, :])
            / dy
            - nu
            * H[0 : ny - 2, :].mean(axis=1, keepdims=True)
            * np.ones((1, nxp1))
            * (U[1 : ny - 1, :] - U[0 : ny - 2, :])
            / dy
        ) / dy
        # Better: use H averaged to U position along y.
        # MATLAB uses H_old(i,j) and H_old(i, j-1/j+1) for the y-diffusion fluxes.
        # For U at col c (between T-cells c-1 and c):
        # Along-y diffusion uses H_old at the same T-row.
        # Re-do with proper H values:
        diff_y = np.zeros_like(U)
        for c in range(1, nx):
            # H_u at [r, c] ≈ 0.5*(H[r,c-1] + H[r,c])
            H_here = 0.5 * (H[:, c - 1] + H[:, c])
            if ny > 2:
                diff_y[1 : ny - 1, c] = (
                    nu / dy * H_here[1 : ny - 1] * (U[2:ny, c] - U[1 : ny - 1, c])
                    - nu
                    / dy
                    * H_here[0 : ny - 2]
                    * (U[1 : ny - 1, c] - U[0 : ny - 2, c])
                ) / dy

    out = diff_x + diff_y
    return out


def _conservative_laplacian_v(
    V: Array, H: Array, nu: float, dx: float, dy: float
) -> Array:
    """Conservative Laplacian for V field (ny+1, nx).

    Same structure as U but transposed axes:
    V is (ny+1, nx), H is (ny, nx).
    Interior V rows: 1..ny-1. Interior V cols: 1..nx-2.
    """
    nyp1, nx = V.shape
    ny = nyp1 - 1
    out = np.zeros_like(V)

    # Diffusion along y (rows of V)
    rsl = slice(1, ny)  # interior V rows
    diff_y = np.zeros_like(V)
    diff_y[rsl, :] = (
        nu * H[1:ny, :] * (V[2:nyp1, :] - V[1:ny, :]) / dy
        - nu * H[0 : ny - 1, :] * (V[1:ny, :] - V[0 : ny - 1, :]) / dy
    ) / dy

    # Diffusion along x (columns of V)
    diff_x = np.zeros_like(V)
    if nx > 2:
        for r in range(1, ny):
            H_here = 0.5 * (H[r - 1, :] + H[r, :])
            diff_x[r, 1 : nx - 1] = (
                nu / dx * H_here[1 : nx - 1] * (V[r, 2:nx] - V[r, 1 : nx - 1])
                - nu / dx * H_here[0 : nx - 2] * (V[r, 1 : nx - 1] - V[r, 0 : nx - 2])
            ) / dx

    out = diff_x + diff_y
    return out


__all__ = ["viscous_tendency"]
