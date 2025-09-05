"""Backward compatibility shim for ModelState import path.

Phase 1 relocated the implementation to ``shel.model.state.model_state``.
This module remains so existing imports (``from shel.model.state import ModelState``)
continue to work while new subpackage structure is introduced.
In later phases this file may host only façade helpers or be removed
after a deprecation window.
"""
from .state.model_state import ModelState  # type: ignore  # noqa: F401

__all__ = ["ModelState"]
