"""Solid-body rotation velocity: u = -Omega (y - y0), v = Omega (x - x0)."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from ..base import VelocityIC
from ..factory import VELOCITY_REGISTRY


class SolidBodyVelocity(VelocityIC):
    def __init__(self, omega: float = 1e-4):
        self.omega = omega

    def build(self, grid, **params) -> tuple[NDArray, NDArray]:
        x_c = grid.x_t.mean()
        y_c = grid.y_t.mean()
        # Build velocities at their native grids using center coordinates along spans.
        u = -self.omega * (grid.y_u - y_c)
        v = self.omega * (grid.x_v - x_c)
        return u, v


VELOCITY_REGISTRY["solid_body"] = SolidBodyVelocity
