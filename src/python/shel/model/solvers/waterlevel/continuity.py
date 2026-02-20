"""Continuity / free-surface update on the C-grid.

Discrete continuity in flux form for barotropic SWE:

        ∂η/∂t = -div( H * U, H * V ) |_T

Provides:
    - ``continuity_rhs``: computes -div(HU, HV) at T points
    - ``update_free_surface``: one Euler step  eta_new = eta - dt * div

Staggering:
    - H, eta on T (ny, nx)
    - U on U faces (ny, nx+1), V on V faces (ny+1, nx)
"""

from __future__ import annotations

import numpy as np

from ..common import avg_x_t_to_u, avg_y_t_to_v, div_uv_to_t

Array = np.ndarray


def continuity_rhs(H: Array, U: Array, V: Array, dx: float, dy: float) -> Array:
    """Compute the continuity tendency -div(H*U, H*V) at T points.

    This is the RHS of deta/dt = -div(flux). The caller performs
    the time stepping (Euler or leapfrog).
    """
    H_u = avg_x_t_to_u(H)
    H_v = avg_y_t_to_v(H)
    Hu = H_u * U
    Hv = H_v * V
    return -div_uv_to_t(Hu, Hv, dx, dy)


def update_free_surface(
    eta: Array, H: Array, U: Array, V: Array, dt: float, dx: float, dy: float
) -> Array:
    """Advance eta by one explicit Euler step using flux divergence.

    Kept for backward compatibility. Equivalent to:
        eta_new = eta + dt * continuity_rhs(H, U, V, dx, dy)
    """
    return eta + dt * continuity_rhs(H, U, V, dx, dy)


__all__ = ["continuity_rhs", "update_free_surface"]
