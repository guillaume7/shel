"""Deprecated solver factory (class-based solvers removed).

The project now exposes functional time steppers under
``shel.model.solvers.time`` (e.g., ``leapfrog_stepper`` and
``leapfrog_step_with_config``). This module is kept for backward
compatibility only and will raise on use.
"""

from __future__ import annotations

from typing import Any, Dict


class SolverFactory:  # pragma: no cover - deprecated
    """Deprecated factory kept for backward compatibility."""

    @staticmethod
    def create(solver_type: str, config: Dict[str, Any]):  # type: ignore[unused-argument]
        raise RuntimeError(
            "SolverFactory is deprecated. Use functions in shel.model.solvers.time, "
            "such as leapfrog_step_with_config."
        )
