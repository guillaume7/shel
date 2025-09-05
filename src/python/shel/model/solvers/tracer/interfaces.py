"""Tracer solver strategy interfaces (Phase 1 skeleton)."""
from __future__ import annotations
from abc import ABC, abstractmethod
from typing import Tuple
from numpy.typing import NDArray

class BaseTracerAdvection(ABC):
    @abstractmethod
    def compute(self, tracer: NDArray, u, v, H, grid, **kwargs) -> NDArray: ...

__all__ = ["BaseTracerAdvection"]
