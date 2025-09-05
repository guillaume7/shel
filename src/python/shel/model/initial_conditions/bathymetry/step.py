"""Step shelf bathymetry IC."""
from __future__ import annotations
import numpy as np
from numpy.typing import NDArray
from ..base import BathymetryIC
from ..factory import BATHYMETRY_REGISTRY

class StepBathymetry(BathymetryIC):
    def __init__(self, depth_shallow: float = 100.0, depth_deep: float = 1000.0, x_frac: float = 0.3):
        self.depth_shallow = depth_shallow
        self.depth_deep = depth_deep
        self.x_frac = x_frac

    def build(self, grid, **params) -> NDArray:
        nx = grid.nx
        cutoff = int(self.x_frac * nx)
        d = np.full((grid.ny, grid.nx), self.depth_deep)
        d[:, :cutoff] = self.depth_shallow
        return d

BATHYMETRY_REGISTRY["step"] = StepBathymetry
