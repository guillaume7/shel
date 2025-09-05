"""Linear shear in x: u = U0 * (y - y_min)/(y_max - y_min); v = 0."""
from __future__ import annotations
import numpy as np
from numpy.typing import NDArray
from ..base import VelocityIC
from ..factory import VELOCITY_REGISTRY

class ShearVelocity(VelocityIC):
    def __init__(self, u0: float = 1.0):
        self.u0 = u0

    def build(self, grid, **params) -> tuple[NDArray, NDArray]:
        y_min = grid.y_u.min(); y_max = grid.y_u.max()
        norm = (grid.y_u - y_min) / (y_max - y_min)
        u = self.u0 * norm
        v = np.zeros((grid.ny + 1, grid.nx))
        return u, v

VELOCITY_REGISTRY["shear"] = ShearVelocity
