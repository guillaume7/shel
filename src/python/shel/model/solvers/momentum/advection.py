"""Centered momentum advection for U and V on a C-grid (minimal baseline).

Computes nonlinear advection terms using centered differences:

  A_u = -( u * ∂u/∂x + ṽ * ∂u/∂y ) on U grid
  A_v = -( ũ * ∂v/∂x + v * ∂v/∂y ) on V grid

where ṽ is V interpolated to U locations and ũ is U interpolated to V.
Boundaries are left as zero (no contribution) to match closed-box prototypes.
"""

from __future__ import annotations

import numpy as np

Array = np.ndarray


def _interp_v_to_u(V: Array) -> Array:
    """Interpolate V (ny+1,nx) to U (ny,nx+1) by 4-point average around U faces.

    For U[:, 1:nx], average V at rows (j, j+1) and cols (i-1, i).
    Boundaries (col 0 and nx) remain zero.
    """
    ny_p1, nx = V.shape  # V: (ny+1, nx)
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
    """Interpolate U (ny,nx+1) to V (ny+1,nx) by 4-point average around V nodes.

    For V[1:ny, :], average U at rows (j-1, j) and cols (i, i+1). Top/bottom
    rows (0 and ny) remain zero.
    """
    ny, nx_p1 = U.shape
    out = np.zeros((ny + 1, nx_p1 - 1))
    out[1:ny, :] = 0.25 * (
        U[0 : ny - 1, 0 : nx_p1 - 1]
        + U[0 : ny - 1, 1:nx_p1]
        + U[1:ny, 0 : nx_p1 - 1]
        + U[1:ny, 1:nx_p1]
    )
    return out


def advect_momentum(U: Array, V: Array, dx: float, dy: float) -> tuple[Array, Array]:
    """Return centered advection tendencies (A_u, A_v) on U and V grids."""
    ny, nxp1 = U.shape
    nyp1, nx = V.shape

    Au = np.zeros_like(U)
    Av = np.zeros_like(V)

    # Interpolated velocities to cross locations
    v_at_u = _interp_v_to_u(V)
    u_at_v = _interp_u_to_v(U)

    # ∂u/∂x on U grid (central)
    dudx = np.zeros_like(U)
    dudx[:, 1 : nxp1 - 1] = (U[:, 2:] - U[:, 0 : nxp1 - 2]) / (2.0 * dx)
    # ∂u/∂y on U grid (central)
    dudy = np.zeros_like(U)
    dudy[1 : ny - 1, :] = (U[2:, :] - U[0 : ny - 2, :]) / (2.0 * dy)

    # ∂v/∂x on V grid (central)
    dvdx = np.zeros_like(V)
    dvdx[:, 1 : nx - 1] = (V[:, 2:] - V[:, 0 : nx - 2]) / (2.0 * dx)
    # ∂v/∂y on V grid (central)
    dvdy = np.zeros_like(V)
    dvdy[1 : nyp1 - 1, :] = (V[2:, :] - V[0 : nyp1 - 2, :]) / (2.0 * dy)

    # Nonlinear terms (interiors effectively)
    Au = -(U * dudx + v_at_u * dudy)
    Av = -(u_at_v * dvdx + V * dvdy)

    return Au, Av


__all__ = ["advect_momentum"]
