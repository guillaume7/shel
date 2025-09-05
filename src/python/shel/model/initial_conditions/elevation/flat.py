"""Flat zero elevation."""
from __future__ import annotations
import numpy as np
from numpy.typing import NDArray
from ..base import ElevationIC
from ..factory import ELEVATION_REGISTRY

class FlatElevation(ElevationIC):
    def build(self, grid, **params) -> NDArray:  # noqa: D401
        return np.zeros((grid.ny, grid.nx))

ELEVATION_REGISTRY["flat"] = FlatElevation
