"""
Common C-grid stencil and interpolation operators for uniform grids.

Conventions (Arakawa C-grid):
- T: scalar at cell centers, shape (ny, nx)
- U: x-velocity at vertical faces, shape (ny, nx+1)
- V: y-velocity at horizontal faces, shape (ny+1, nx)

All arrays are float64. Derivatives use centered differences mapped to the
appropriate staggered locations. Boundaries for T->U and T->V derivatives
that require a neighbor outside the domain are filled with NaN.
"""

from __future__ import annotations

import numpy as np

Array = np.ndarray


def d_dx_t_to_u(T: Array, dx: float) -> Array:
    """Compute ∂T/∂x at U points using centered differences.

    Result shape: (ny, nx+1). Interior columns 1..nx-1 are computed as
    (T[:, i] - T[:, i-1]) / dx. Boundary columns 0 and nx are set to NaN.
    """
    ny, nx = T.shape
    out = np.full((ny, nx + 1), np.nan, dtype=T.dtype)
    # interior U columns correspond to differences between T columns
    out[:, 1:nx] = (T[:, 1:] - T[:, :-1]) / dx
    return out


def d_dy_t_to_v(T: Array, dy: float) -> Array:
    """Compute ∂T/∂y at V points using centered differences.

    Result shape: (ny+1, nx). Interior rows 1..ny-1 are computed as
    (T[j, :] - T[j-1, :]) / dy. Boundary rows 0 and ny are set to NaN.
    """
    ny, nx = T.shape
    out = np.full((ny + 1, nx), np.nan, dtype=T.dtype)
    out[1:ny, :] = (T[1:, :] - T[:-1, :]) / dy
    return out


def avg_x_t_to_u(T: Array) -> Array:
    """Average T to U faces along x (arithmetic mean).

    Result shape: (ny, nx+1). Interior columns are (T[:, i] + T[:, i-1]) / 2.
    Boundaries 0 and nx are NaN (no neighbor beyond domain).
    """
    ny, nx = T.shape
    out = np.empty((ny, nx + 1), dtype=T.dtype)
    # interior arithmetic mean
    out[:, 1:nx] = 0.5 * (T[:, 1:] + T[:, :-1])
    # boundary: nearest-cell copy (matches MATLAB; avoids negative H from extrapolation)
    out[:, 0] = T[:, 0]
    out[:, -1] = T[:, -1]
    return out


def avg_y_t_to_v(T: Array) -> Array:
    """Average T to V faces along y.

    Result shape: (ny+1, nx). Interior rows are (T[j, :] + T[j-1, :]) / 2.
    Boundaries 0 and ny are NaN.
    """
    ny, nx = T.shape
    out = np.empty((ny + 1, nx), dtype=T.dtype)
    out[1:ny, :] = 0.5 * (T[1:, :] + T[:-1, :])
    # boundary: nearest-cell copy (matches MATLAB; avoids negative H from extrapolation)
    out[0, :] = T[0, :]
    out[-1, :] = T[-1, :]
    return out


def avg_x_u_to_t(U: Array) -> Array:
    """Average U faces back to T centers along x.

    Result shape: (ny, nx). Computed as (U[:, 1:] + U[:, :-1]) / 2.
    """
    stacked = np.stack([U[:, 1:], U[:, :-1]], axis=0)
    return np.nanmean(stacked, axis=0)


def avg_y_v_to_t(V: Array) -> Array:
    """Average V faces back to T centers along y.

    Result shape: (ny, nx). Computed as (V[1:, :] + V[:-1, :]) / 2.
    """
    stacked = np.stack([V[1:, :], V[:-1, :]], axis=0)
    return np.nanmean(stacked, axis=0)


def div_uv_to_t(U: Array, V: Array, dx: float, dy: float) -> Array:
    """Compute discrete divergence at T points from U,V on faces.

    div_T = dU/dx + dV/dy mapped onto T centers using face differences.
    Result shape: (ny, nx).
    """
    d_udx = (U[:, 1:] - U[:, :-1]) / dx
    d_vdy = (V[1:, :] - V[:-1, :]) / dy
    return d_udx + d_vdy
