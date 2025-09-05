"""Global (time-series) diagnostics accumulator.

Collects scalar integrated diagnostics each model step (or at a chosen
output cadence) into simple Python lists for later serialization or
conversion to DataFrames / Parquet.

Design goals:
 - Stateless computations delegated to IntegratedDiagnostics.
 - Simple append-only structure (can be replaced by ring buffer later).
 - Deterministic ordering of keys for reproducible serialization.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Dict
from numpy.typing import NDArray
from .integrated import IntegratedDiagnostics
from ..grid import Grid  # type: ignore

@dataclass
class GlobalAccumulator:
    times: List[float] = field(default_factory=list)
    kinetic_energy: List[float] = field(default_factory=list)
    potential_energy: List[float] = field(default_factory=list)
    total_energy: List[float] = field(default_factory=list)
    volume: List[float] = field(default_factory=list)
    enstrophy: List[float] = field(default_factory=list)
    potential_enstrophy: List[float] = field(default_factory=list)

    def update(
        self,
        t: float,
        u: NDArray,
        v: NDArray,
        eta: NDArray,
        H: NDArray,
        coriolis: NDArray,
        grid: Grid,
        gravity: float = 9.81,
    ) -> None:
        diags = IntegratedDiagnostics.as_dict(u, v, eta, H, grid, gravity, coriolis)
        self.times.append(t)
        self.kinetic_energy.append(diags["kinetic_energy"])
        self.potential_energy.append(diags["potential_energy"])
        self.total_energy.append(diags["total_energy"])
        self.volume.append(diags["volume"])
        self.enstrophy.append(diags["enstrophy"])
        self.potential_enstrophy.append(diags["potential_enstrophy"])

    def as_dict(self) -> Dict[str, List[float]]:
        return {
            "time": self.times,
            "kinetic_energy": self.kinetic_energy,
            "potential_energy": self.potential_energy,
            "total_energy": self.total_energy,
            "volume": self.volume,
            "enstrophy": self.enstrophy,
            "potential_enstrophy": self.potential_enstrophy,
        }

__all__ = ["GlobalAccumulator"]
