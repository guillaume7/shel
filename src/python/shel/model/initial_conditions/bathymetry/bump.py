"""Gaussian bump bathymetry initial condition."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from ..base import BathymetryIC
from ..factory import BATHYMETRY_REGISTRY


class BumpBathymetry(BathymetryIC):
    def __init__(
        self,
        depth0: float = 1000.0,
        amp: float = 200.0,
        sx: float = 0.2,
        sy: float = 0.2,
        x0: float | None = None,
        y0: float | None = None,
    ):
        self.depth0 = depth0
        self.amp = amp
        self.sx = sx
        self.sy = sy
        self.x0 = x0
        self.y0 = y0

    def build(self, grid, **params) -> NDArray:
        x0 = self.x0 if self.x0 is not None else grid.x_t.mean()
        y0 = self.y0 if self.y0 is not None else grid.y_t.mean()
        X = grid.x_t
        Y = grid.y_t
        rx = (X - x0) / (self.sx * (grid.x_t.max() - grid.x_t.min()))
        ry = (Y - y0) / (self.sy * (grid.y_t.max() - grid.y_t.min()))
        bump = self.amp * np.exp(-(rx**2 + ry**2))
        return np.full(X.shape, self.depth0) - bump


BATHYMETRY_REGISTRY["bump"] = BumpBathymetry
