"""Surface forcings placeholder (Phase 1) – inert implementation.

Future Phase will move logic into ``forcings/surface/`` submodules.
"""
from __future__ import annotations

from typing import Dict, Any
import numpy as np
from numpy.typing import NDArray


def apply_wind_stress(u: NDArray, v: NDArray, params: Dict[str, Any]) -> tuple[NDArray, NDArray]:
    """Identity wind stress (no-op) for Phase 1 skeleton."""
    return u, v

__all__ = ["apply_wind_stress"]
