"""Tracer update orchestrator (Phase 1 placeholder)."""

from __future__ import annotations

from typing import Any, Mapping

import numpy as np

from ..common.stepper import apply_tracer_bc_per_side, resolve_bc_sides_from_config

Array = np.ndarray


def tracer_step_minimal(
    C: Array,
    U: Array,
    V: Array,
    *,
    dt: float,
    dx: float,
    dy: float,
    config: Mapping[str, Any] | None = None,
) -> Array:
    """Minimal tracer step that only applies BCs.

    Placeholder for future advection/diffusion; currently returns C after BCs.
    """
    Cn = np.array(C, copy=True)
    bc_sides = resolve_bc_sides_from_config(config)
    apply_tracer_bc_per_side(
        Cn, bc_sides=bc_sides, C_old=C, U=U, V=V, dt=dt, dx=dx, dy=dy
    )
    return Cn


def tracer_tendency(*_args, **_kwargs):  # type: ignore[no-untyped-def]
    """Placeholder tracer tendency function to be implemented in later phases."""
    raise NotImplementedError("Tracer tendency not yet implemented (Phase 1)")


__all__ = ["tracer_tendency", "tracer_step_minimal"]
