"""Pressure gradient term on an Arakawa C-grid (minimal prototype).

Given free-surface elevation ``eta`` on T points, compute the barotropic
pressure gradient acceleration at U and V faces:

    PG_u = -g * ∂eta/∂x |_U
    PG_v = -g * ∂eta/∂y |_V

Derivatives use centered differences mapped T→U/V. Boundary faces are NaN.
"""

from __future__ import annotations

import numpy as np

from ..common import d_dx_t_to_u, d_dy_t_to_v

Array = np.ndarray


def pressure_gradient(
    eta: Array, g: float, dx: float, dy: float
) -> tuple[Array, Array]:
    """Compute pressure gradient accelerations at faces from eta on centers.

    Parameters
    - eta: (ny, nx) free-surface elevation at T points
    - g: gravitational acceleration (scalar)
    - dx, dy: grid spacing (scalars)

    Returns
    - PG_u: (ny, nx+1) array at U faces
    - PG_v: (ny+1, nx) array at V faces
    """
    d_eta_dx_u = d_dx_t_to_u(eta, dx)
    d_eta_dy_v = d_dy_t_to_v(eta, dy)
    PG_u = -g * d_eta_dx_u
    PG_v = -g * d_eta_dy_v
    return PG_u, PG_v


__all__ = ["pressure_gradient"]
