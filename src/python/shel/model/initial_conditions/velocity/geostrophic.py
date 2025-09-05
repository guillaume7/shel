"""Geostrophic velocity initial condition (MATLAB parity, Option A).

Mimics the original MATLAB parametric construction rather than computing
finite-difference pressure gradients. In the MATLAB code the geostrophic
"ugeo" and "vgeo" fields are proportional to the displacement from the
domain centre times the elevation field:

    u_T(i,j) =  g/f * (j - j_c) * eta(i,j) / sigma_y^2
    v_T(i,j) = -g/f * (i - i_c) * eta(i,j) / sigma_x^2

with sigma_x = sigma_y = 1 (index units) and i_c, j_c the central T-cell
indices (integer). These T-grid velocities are then interpolated to the
staggered U and V grids. We replicate that logic here for parity and to
keep tests simple. This is a *parametric* geostrophic field, not the
diagnostic one obtained from numerical derivatives of eta.
"""
from __future__ import annotations
import numpy as np
from numpy.typing import NDArray
from ..base import VelocityIC
from ..factory import VELOCITY_REGISTRY


class GeostrophicVelocity(VelocityIC):
    def __init__(self, gravity: float = 9.81):
        self.g = gravity

    def build(self, grid, **params) -> tuple[NDArray, NDArray]:
        eta = params.get("eta")
        coriolis = params.get("coriolis")
        if eta is None or coriolis is None:
            raise ValueError("GeostrophicVelocity requires 'eta' and 'coriolis'")
        ny, nx = eta.shape
        u = np.zeros((ny, nx + 1))
        v = np.zeros((ny + 1, nx))
        f_safe = np.where(coriolis == 0.0, 1e-12, coriolis)
        # Central indices (zero-based) corresponding to MATLAB j_L=i_L=floor(N/2)+1 (1-based)
        j_c = ny // 2
        i_c = nx // 2
        # Index grids (T grid indexing)
        j_idx = np.arange(ny).reshape(ny, 1)
        i_idx = np.arange(nx).reshape(1, nx)
        sigma_x = 1.0
        sigma_y = 1.0
        ug = (self.g / f_safe) * ((j_idx - j_c) / (sigma_y**2)) * eta  # (ny, nx)
        vg = -(self.g / f_safe) * ((i_idx - i_c) / (sigma_x**2)) * eta  # (ny, nx)
        # Interpolate T -> U (average adjacent T columns)
        u[:, 1:-1] = 0.5 * (ug[:, :-1] + ug[:, 1:])
        # Interpolate T -> V (average adjacent T rows)
        v[1:-1, :] = 0.5 * (vg[:-1, :] + vg[1:, :])
        return u, v


VELOCITY_REGISTRY["geostrophic"] = GeostrophicVelocity
