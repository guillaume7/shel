"""Abstract initial condition category interfaces (Phase 1 skeleton)."""

from __future__ import annotations

from abc import ABC, abstractmethod

from numpy.typing import NDArray


class BathymetryIC(ABC):
    @abstractmethod
    def build(self, grid, **params) -> NDArray: ...


class ElevationIC(ABC):
    @abstractmethod
    def build(self, grid, **params) -> NDArray: ...


class VelocityIC(ABC):
    @abstractmethod
    def build(self, grid, **params) -> tuple[NDArray, NDArray]: ...  # u, v


class TracerIC(ABC):
    @abstractmethod
    def build(self, grid, **params) -> NDArray: ...


__all__ = ["BathymetryIC", "ElevationIC", "VelocityIC", "TracerIC"]
