"""Factories for continuity (free-surface) strategies (Phase 1 skeleton)."""
from __future__ import annotations
from typing import Dict, Type
from .interfaces import BaseContinuity

CONTINUITY_REGISTRY: Dict[str, Type[BaseContinuity]] = {}

def create_continuity(name: str, **kwargs) -> BaseContinuity:
    try:
        return CONTINUITY_REGISTRY[name.lower()](**kwargs)
    except KeyError as exc:
        raise ValueError(f"Unknown continuity scheme '{name}'") from exc

__all__ = ["create_continuity", "CONTINUITY_REGISTRY"]
