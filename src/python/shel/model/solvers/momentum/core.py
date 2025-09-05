"""Momentum update orchestrator (Phase 1 placeholder).

Future phases: aggregate tendencies from advection, pressure, diffusion,
friction modules returning (du_dt, dv_dt) arrays.
"""
from __future__ import annotations
from numpy.typing import NDArray

def momentum_tendency(*_args, **_kwargs) -> tuple[NDArray, NDArray]:  # type: ignore[name-defined]
    """Placeholder returning NotImplemented to flag unimplemented usage."""
    raise NotImplementedError("Momentum tendency not yet implemented (Phase 1)")

__all__ = ["momentum_tendency"]
