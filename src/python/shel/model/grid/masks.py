"""Grid mask construction utilities (Phase 2).

Faithful Python port of MATLAB mask logic in `model_handles.m` lines ~335-374.

Conventions (Arakawa C-grid):
    T cells (scalar points eta,H,tracers,bathymetry): shape (ny, nx)
    U faces (zonal velocity u): shape (ny, nx+1)
    V faces (meridional velocity v): shape (ny+1, nx)
    Q/F corners / vorticity points: shape (ny+1, nx+1) (future)

Land/Water Encoding:
    1 -> water
    0 -> land (masked / removed)

Functions herein are *pure*: they do not mutate the provided mask arrays.
"""

from __future__ import annotations

from typing import Tuple

import numpy as np
from numpy.typing import NDArray

__all__ = [
    "build_staggered_masks",
    "apply_noslip_flux_masks",
    "build_corner_mask",
    "build_all_masks",
    "mask_velocities",
]


def build_staggered_masks(
    mask_t: NDArray[np.integer],
) -> Tuple[NDArray[np.int_], NDArray[np.int_]]:
    """Build U and V staggered masks from T-cell land mask.

    Replicates the MATLAB loops:
        if mask(i,j) == 0:
            mask_u(i,j) = mask_u(i+1,j) = 0
            mask_v(i,j) = mask_v(i,j+1) = 0

    Vectorized approach:
        For each land T cell, zero the adjacent U/V faces that touch it.

    Parameters
    ----------
    mask_t : (ny, nx) int array
        T-cell mask (1 water, 0 land).

    Returns
    -------
    mask_u : (ny, nx+1) int array
    mask_v : (ny+1, nx) int array
    """
    ny, nx = mask_t.shape
    mask_u = np.ones((ny, nx + 1), dtype=int)
    mask_v = np.ones((ny + 1, nx), dtype=int)

    # Land cell indices
    land = np.where(mask_t == 0)
    if land[0].size == 0:
        return mask_u, mask_v

    i_t = land[0]
    j_t = land[1]

    # U faces: (ny, nx+1). Land cell at (i,j) affects (i,j) and (i,j+1) in MATLAB 1-based.
    # In 0-based Python: affects (i,j) and (i,j+1) provided j+1 <= nx.
    mask_u[i_t, j_t] = 0
    inside = (
        j_t + 1 <= nx
    )  # j_t max is nx-1 so j_t+1 <= nx always True, but keep explicit.
    mask_u[i_t[inside], j_t[inside] + 1] = 0

    # V faces: (ny+1, nx). Land cell at (i,j) affects (i,j) and (i+1,j) in MATLAB 1-based.
    mask_v[i_t, j_t] = 0
    inside_i = i_t + 1 <= ny - 1  # ensure within ny (since v has ny+1)
    mask_v[i_t[inside_i] + 1, j_t[inside_i]] = 0

    return mask_u, mask_v


def apply_noslip_flux_masks(
    mask_t: NDArray[np.integer], mask_u: NDArray[np.int_], mask_v: NDArray[np.int_]
) -> Tuple[NDArray[np.int_], NDArray[np.int_]]:
    """Derive noslip flux masks (mask_u, mask_v) from T mask per MATLAB logic.

    MATLAB (model_handles.m):
        mask_u(i,j) = mask(i  ,j+1) * mask(i  ,j-1) * mask(i-1,j+1) * mask(i-1,j-1)
        mask_v(i,j) = mask(i+1,  j) * mask(i+1,j-1) * mask(i-1,  j) * mask(i-1,j-1)
    for interior ranges (converted to Python 0-based indexing).

    We recompute new arrays (do not mutate inputs) limited to interior indices.

    Parameters
    ----------
    mask_t : (ny, nx) int array
    mask_u, mask_v : precomputed face masks (will be copied)

    Returns
    -------
    noslip_mask_u, noslip_mask_v : int arrays with interior faces zeroed when any adjacent corner T cell is land.
    """
    ny, nx = mask_t.shape
    mu = mask_u.copy()
    mv = mask_v.copy()

    # U faces interior indices: i in [1, ny-1), j in [1, nx)
    # For each U face centered between T cells at (i-1,j-1),(i-1,j),(i,j-1),(i,j)
    # MATLAB expression corresponds to those four T cells.
    if ny > 1 and nx > 1:
        i_idx = np.arange(1, ny)  # 1..ny-1
        j_idx = np.arange(
            1, nx
        )  # 1..nx-1 (since mu has nx+1 columns, interior j faces exclude endpoints)
        # Broadcast create 2D grids
        I, J = np.meshgrid(i_idx, j_idx, indexing="ij")  # shapes (ny-1, nx-1)
        corners = (
            mask_t[I - 1, J] * mask_t[I - 1, J - 1] * mask_t[I, J] * mask_t[I, J - 1]
        )
        mu[I, J] = corners

    # V faces interior indices: i in [1, ny), j in [1, nx-1)
    # V face uses T cells at (i-1,j-1),(i-1,j),(i,j-1),(i,j)
    if ny > 1 and nx > 1:
        i_idx_v = np.arange(1, ny)
        j_idx_v = np.arange(1, nx)
        I, J = np.meshgrid(i_idx_v, j_idx_v, indexing="ij")
        corners_v = (
            mask_t[I - 1, J - 1] * mask_t[I - 1, J] * mask_t[I, J - 1] * mask_t[I, J]
        )
        mv[I, J] = corners_v

    return mu, mv


def build_corner_mask(mask_t: NDArray[np.integer]) -> NDArray[np.int_]:
    """Build corner (Q) mask as product of adjacent T cells.

    A corner (i,j) in Q grid touches up to 4 T cells:
        (i-1,j-1), (i-1,j), (i,j-1), (i,j) in Python 0-based indexing.
    For boundary corners fewer cells exist; we treat missing outside-domain
    cells as water (factor 1) so that edge corners reflect the interior
    water presence without artificially zeroing.
    """
    ny, nx = mask_t.shape
    mq = np.ones((ny + 1, nx + 1), dtype=int)
    # Interior corners depend on 4 adjacent T cells
    mq[1:ny, 1:nx] = (
        mask_t[0 : ny - 1, 0 : nx - 1]
        * mask_t[0 : ny - 1, 1:nx]
        * mask_t[1:ny, 0 : nx - 1]
        * mask_t[1:ny, 1:nx]
    )
    # Edges: use available adjacent T cells (already 1 by default elsewhere)
    # Top edge (i=0): depends on first row of T cells horizontally
    mq[0, 1:nx] = mask_t[0, 0 : nx - 1] * mask_t[0, 1:nx]
    # Bottom edge (i=ny): last row
    mq[ny, 1:nx] = mask_t[ny - 1, 0 : nx - 1] * mask_t[ny - 1, 1:nx]
    # Left edge (j=0)
    mq[1:ny, 0] = mask_t[0 : ny - 1, 0] * mask_t[1:ny, 0]
    # Right edge (j=nx)
    mq[1:ny, nx] = mask_t[0 : ny - 1, nx - 1] * mask_t[1:ny, nx - 1]
    # Corners remain 1 if any adjacent T is water; if all adjacent (1 or 2) T cells are 0 set to 0
    mq[0, 0] = mask_t[0, 0]
    mq[0, nx] = mask_t[0, nx - 1]
    mq[ny, 0] = mask_t[ny - 1, 0]
    mq[ny, nx] = mask_t[ny - 1, nx - 1]
    return mq


def build_all_masks(
    mask_t: NDArray[np.integer], noslip: bool = False
) -> Tuple[NDArray[np.int_], NDArray[np.int_], NDArray[np.int_]]:
    """Convenience wrapper: construct U/V/Q masks with optional noslip processing.

    Returns
    -------
    mask_u, mask_v, mask_q
    """
    mu, mv = build_staggered_masks(mask_t)
    if noslip:
        mu, mv = apply_noslip_flux_masks(mask_t, mu, mv)
    mq = build_corner_mask(mask_t)
    return mu, mv, mq


def mask_velocities(
    u: NDArray, v: NDArray, mask_u: NDArray, mask_v: NDArray
) -> Tuple[NDArray, NDArray]:
    """Apply staggered masks to velocity components (pure function).

    Parameters
    ----------
    u : (ny, nx+1) array
    v : (ny+1, nx) array
    mask_u : (ny, nx+1) int array
    mask_v : (ny+1, nx) int array

    Returns
    -------
    u_masked, v_masked : arrays with land / noslip faces zeroed.
    """
    return u * mask_u, v * mask_v
