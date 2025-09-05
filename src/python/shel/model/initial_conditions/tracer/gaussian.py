"""Gaussian tracer distribution."""
from __future__ import annotations
import numpy as np
from numpy.typing import NDArray
from ..base import TracerIC
from ..factory import TRACER_REGISTRY

class GaussianTracer(TracerIC):
    def __init__(self, c0: float = 1.0, sx: float = 0.2, sy: float = 0.2):
        self.c0 = c0; self.sx = sx; self.sy = sy

    def build(self, grid, **params) -> NDArray:
        X = grid.x_t; Y = grid.y_t
        x0 = X.mean(); y0 = Y.mean()
        Lx = X.max()-X.min(); Ly = Y.max()-Y.min()
        rx = (X - x0)/(self.sx*Lx); ry = (Y - y0)/(self.sy*Ly)
        return self.c0 * np.exp(-(rx**2 + ry**2))

TRACER_REGISTRY["gaussian"] = GaussianTracer
