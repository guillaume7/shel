"""Uniform tracer field."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from ..base import TracerIC
from ..factory import TRACER_REGISTRY


class UniformTracer(TracerIC):
    def __init__(self, value: float = 1.0):
        self.value = value

    def build(self, grid, **params) -> NDArray:
        return np.full((grid.ny, grid.nx), self.value)


TRACER_REGISTRY["uniform"] = UniformTracer
