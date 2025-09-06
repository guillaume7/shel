"""Tracer solver strategy interfaces (Phase 1 skeleton)."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any

from numpy.typing import NDArray


class BaseTracerAdvection(ABC):
    """Abstract advection interface for tracer transport.

    Implementations should compute the advective tendency of a tracer field.
    """

    @abstractmethod
    def compute(
        self,
        tracer: NDArray,
        u: NDArray,
        v: NDArray,
        H: NDArray,
        grid: Any,
        **kwargs: Any,
    ) -> NDArray:
        """Return advective tendency for tracer given velocities and grid."""


__all__ = ["BaseTracerAdvection"]
