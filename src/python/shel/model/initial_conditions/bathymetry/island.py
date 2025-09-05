"""Circular island (elevated bathymetry -> shallower depth)."""
from __future__ import annotations
import numpy as np
from numpy.typing import NDArray
from ..base import BathymetryIC
from ..factory import BATHYMETRY_REGISTRY

class IslandBathymetry(BathymetryIC):
    def __init__(self, depth0: float = 1000.0, raise_height: float = 900.0, radius_frac: float = 0.1, x0: float | None = None, y0: float | None = None):
        self.depth0 = depth0
        self.raise_height = raise_height
        self.radius_frac = radius_frac
        self.x0 = x0
        self.y0 = y0

    def build(self, grid, **params) -> NDArray:
        X = grid.x_t; Y = grid.y_t
        x0 = self.x0 if self.x0 is not None else X.mean()
        y0 = self.y0 if self.y0 is not None else Y.mean()
        Lx = X.max() - X.min(); Ly = Y.max() - Y.min()
        R = self.radius_frac * min(Lx, Ly)
        r = np.sqrt((X - x0)**2 + (Y - y0)**2)
        island = np.where(r <= R, self.raise_height, 0.0)
        return np.full(X.shape, self.depth0) - island

BATHYMETRY_REGISTRY["island"] = IslandBathymetry
