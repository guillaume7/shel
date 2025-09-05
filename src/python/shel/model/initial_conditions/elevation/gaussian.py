"""Gaussian elevation bump."""
from __future__ import annotations
import numpy as np
from numpy.typing import NDArray
from ..base import ElevationIC
from ..factory import ELEVATION_REGISTRY

class GaussianElevation(ElevationIC):
    def __init__(self, amp: float = 1.0, sx: float = 0.2, sy: float = 0.2, x0: float | None = None, y0: float | None = None):
        self.amp = amp; self.sx = sx; self.sy = sy; self.x0 = x0; self.y0 = y0

    def build(self, grid, **params) -> NDArray:
        X = grid.x_t; Y = grid.y_t
        x0 = self.x0 if self.x0 is not None else X.mean()
        y0 = self.y0 if self.y0 is not None else Y.mean()
        Lx = X.max()-X.min(); Ly = Y.max()-Y.min()
        rx = (X - x0)/(self.sx * Lx); ry = (Y - y0)/(self.sy * Ly)
        return self.amp * np.exp(-(rx**2 + ry**2))

ELEVATION_REGISTRY["gaussian"] = GaussianElevation
