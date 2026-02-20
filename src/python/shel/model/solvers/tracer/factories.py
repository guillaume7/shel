"""Factories for tracer advection strategies (Phase 1 skeleton)."""

from __future__ import annotations

from typing import Dict, Type

from .interfaces import BaseTracerAdvection

TRACER_ADVECTION_REGISTRY: Dict[str, Type[BaseTracerAdvection]] = {}


def create_tracer_advection(name: str, **kwargs) -> BaseTracerAdvection:
    try:
        return TRACER_ADVECTION_REGISTRY[name.lower()](**kwargs)
    except KeyError as exc:
        raise ValueError(f"Unknown tracer advection scheme '{name}'") from exc


__all__ = ["create_tracer_advection", "TRACER_ADVECTION_REGISTRY"]
