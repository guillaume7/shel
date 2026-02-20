"""Water level (continuity) strategy interfaces (Phase 1 skeleton)."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Tuple

from numpy.typing import NDArray


class BaseContinuity(ABC):
    @abstractmethod
    def compute(
        self, eta: NDArray, u: NDArray, v: NDArray, H: NDArray, grid, dt: float
    ) -> Tuple[NDArray, NDArray]: ...


__all__ = ["BaseContinuity"]
