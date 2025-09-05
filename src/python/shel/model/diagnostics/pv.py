"""Potential vorticity diagnostic field utilities.

Computes potential vorticity on T cells:

    q = (ζ_T + f) / H

where ζ_T is the relative vorticity interpolated from corner (Q/F) grid
values to T cells (area average of adjoining four corners), f is the
Coriolis parameter on T cells and H the total depth.

Numerical fidelity: mirrors logic in IntegratedDiagnostics._vorticity_t
to avoid duplication and ensure identical discretization of ζ_T.
"""
from __future__ import annotations
from typing import Dict
import numpy as np
from numpy.typing import NDArray
from .integrated import IntegratedDiagnostics
from ..grid import Grid  # type: ignore

def potential_vorticity(u: NDArray, v: NDArray, H: NDArray, f: NDArray, grid: Grid) -> NDArray:
    """Compute Ertel-like barotropic potential vorticity q on T cells.

    Parameters
    ----------
    u, v : ndarray
        Velocity components on C-grid (U and V staggering).
    H : ndarray (ny, nx)
        Total depth (bathymetry + surface elevation).
    f : ndarray (ny, nx)
        Coriolis parameter on T cells.
    grid : Grid
        Grid instance (for dimensions & spacing).

    Returns
    -------
    ndarray (ny, nx)
        Potential vorticity field ( (ζ+f)/H ).
    """
    zeta_t = IntegratedDiagnostics._vorticity_t(u, v, grid)
    return (zeta_t + f) / H

def as_dict(u: NDArray, v: NDArray, H: NDArray, f: NDArray, grid: Grid) -> Dict[str, NDArray]:
    return {"potential_vorticity": potential_vorticity(u, v, H, f, grid)}

__all__ = ["potential_vorticity", "as_dict"]
