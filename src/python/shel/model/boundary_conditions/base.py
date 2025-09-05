"""
Functional boundary-condition strategy base classes.

These operate directly on staggered arrays used by the ministep path
and avoid coupling to ModelState. They complement the legacy OO API
in `boundary.py` that works with ModelState and the orchestrated solver.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Protocol

import numpy as np

Array = np.ndarray


class BoundaryCondition(ABC):
    name: str


class MomentumBC(BoundaryCondition, ABC):
    """Momentum (U,V) boundary condition strategy on a C-grid.

    Provides uniform application over the whole domain edges and per-side updates.
    """

    @abstractmethod
    def apply_uniform(self, U: Array, V: Array) -> None:
        """Apply uniformly to all four domain sides in-place."""

    @abstractmethod
    def apply_side(
        self,
        U: Array,
        V: Array,
        side: str,
        *,
        U_old: Array | None = None,
        V_old: Array | None = None,
        H: Array | None = None,
        g: float | None = None,
        dt: float | None = None,
        dx: float | None = None,
        dy: float | None = None,
    eta_old: Array | None = None,
    eta_ext: Array | None = None,
    ) -> None:
        """Apply to a single side in-place. Extra args optional per strategy."""


class EtaBC(BoundaryCondition, ABC):
    """Free-surface (eta) boundary condition strategy."""

    @abstractmethod
    def apply_side_eta(
        self,
        eta_next: Array,
        side: str,
        *,
        eta_old: Array | None = None,
    eta_ext: Array | None = None,
        H: Array | None = None,
        g: float | None = None,
        dt: float | None = None,
        dx: float | None = None,
        dy: float | None = None,
    ) -> None:
        """Apply to a single side in-place. Extra args optional per strategy."""
