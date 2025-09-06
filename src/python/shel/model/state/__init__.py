"""State subpackage (Phase 1 skeleton).

This subpackage will be split into:
- core.py : CoreState data container (arrays & parameters)
- update.py : Time-stepping orchestration utilities
- diagnostics_buffer.py : Optional rolling accumulators for global diagnostics

Currently re-exports existing ModelState to avoid breaking imports.
Refactor in Phase 5 will migrate logic.
"""

from .model_state import ModelState  # type: ignore  # noqa: F401

__all__ = ["ModelState"]
