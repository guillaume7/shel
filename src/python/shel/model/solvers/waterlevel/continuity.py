"""Continuity / free-surface update (minimal prototype).

Discrete continuity in flux form on the C-grid for barotropic SWE:

        eta^{n+1} = eta^{n} - dt * div( H * U, H * V ) |_T

Here we provide a single explicit Euler update for initial prototyping and
unit tests. Staggering:
    - H, eta on T (ny, nx)
    - U on U faces (ny, nx+1), V on V faces (ny+1, nx)
"""

from __future__ import annotations

import numpy as np

from ..common import avg_x_t_to_u, avg_y_t_to_v, div_uv_to_t

Array = np.ndarray


def update_free_surface(
    eta: Array, H: Array, U: Array, V: Array, dt: float, dx: float, dy: float
) -> Array:
    """Advance eta by one explicit Euler step using flux divergence.

    Fluxes use H averaged to faces: Hu = H̄_x * U, Hv = H̄_y * V.
    Boundaries are naturally handled by divergence using interior face diffs.
    """
    H_u = avg_x_t_to_u(H)
    H_v = avg_y_t_to_v(H)
    Hu = H_u * U
    Hv = H_v * V
    div = div_uv_to_t(Hu, Hv, dx, dy)
    return eta - dt * div


__all__ = ["update_free_surface"]
