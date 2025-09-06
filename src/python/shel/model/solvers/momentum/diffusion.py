"""Viscous diffusion tendency for U and V on a C-grid (minimal Laplacian).

Applies ν ∇² U and ν ∇² V using standard 5-point stencil on the respective
staggered grids. Boundaries are left as zeros by default (Dirichlet-like),
and NaN values are preserved if present.
"""

from __future__ import annotations

import numpy as np

Array = np.ndarray


def _laplacian(arr: Array, dx: float, dy: float) -> Array:
    ny, nx = arr.shape
    out = np.zeros_like(arr)
    # interior points
    out[1 : ny - 1, 1 : nx - 1] = (
        arr[1 : ny - 1, 2:]
        - 2 * arr[1 : ny - 1, 1 : nx - 1]
        + arr[1 : ny - 1, 0 : nx - 2]
    ) / (dx * dx) + (
        arr[2:, 1 : nx - 1]
        - 2 * arr[1 : ny - 1, 1 : nx - 1]
        + arr[0 : ny - 2, 1 : nx - 1]
    ) / (
        dy * dy
    )
    # propagate NaNs where any neighbor is NaN
    nan_mask = (
        np.isnan(arr[1 : ny - 1, 2:])
        | np.isnan(arr[1 : ny - 1, 1 : nx - 1])
        | np.isnan(arr[1 : ny - 1, 0 : nx - 2])
        | np.isnan(arr[2:, 1 : nx - 1])
        | np.isnan(arr[0 : ny - 2, 1 : nx - 1])
    )
    out[1 : ny - 1, 1 : nx - 1][nan_mask] = np.nan
    return out


def viscous_tendency(
    U: Array, V: Array, nu: float, dx: float, dy: float
) -> tuple[Array, Array]:
    """Compute ν ∇² U and ν ∇² V on their native staggered grids.

    Parameters
    - U: (ny, nx+1), V: (ny+1, nx)
    - nu: kinematic viscosity
    - dx, dy: spacings
    """
    return nu * _laplacian(U, dx, dy), nu * _laplacian(V, dx, dy)


__all__ = ["viscous_tendency"]
