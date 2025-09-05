"""Momentum solver strategy interfaces (Phase 1 skeleton).

Defines abstract base classes for interchangeable momentum-related
numerical schemes (advection, pressure gradient, diffusion, friction).
All `compute` methods return tendencies (du_dt, dv_dt) with no side
effects. Concrete implementations arrive in later phases.
"""
from __future__ import annotations
from abc import ABC, abstractmethod
from typing import Tuple, Protocol
from numpy.typing import NDArray

class MomentumTendency(Protocol):  # lightweight structural type
    def compute(self, u: NDArray, v: NDArray, *args, **kwargs) -> Tuple[NDArray, NDArray]: ...

class BaseMomentumAdvection(ABC):
    @abstractmethod
    def compute(self, u: NDArray, v: NDArray, H: NDArray, grid, **kwargs) -> Tuple[NDArray, NDArray]: ...

class BasePressureGradient(ABC):
    @abstractmethod
    def compute(self, eta: NDArray, H: NDArray, grid, gravity: float) -> Tuple[NDArray, NDArray]: ...

class BaseViscosityOperator(ABC):
    @abstractmethod
    def compute(self, u: NDArray, v: NDArray, grid, nu: float) -> Tuple[NDArray, NDArray]: ...

class BaseFriction(ABC):
    @abstractmethod
    def compute(self, u: NDArray, v: NDArray, H: NDArray, grid, **params) -> Tuple[NDArray, NDArray]: ...

__all__ = [
    "BaseMomentumAdvection",
    "BasePressureGradient",
    "BaseViscosityOperator",
    "BaseFriction",
]
