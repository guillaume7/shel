"""Factories for momentum strategy components (Phase 1 skeleton)."""

from __future__ import annotations

from typing import Dict, Type

from .interfaces import (
    BaseFriction,
    BaseMomentumAdvection,
    BasePressureGradient,
    BaseViscosityOperator,
)

# Registries (name -> class). Populated in later phases when concrete classes added.
MOMENTUM_ADVECTION_REGISTRY: Dict[str, Type[BaseMomentumAdvection]] = {}
PRESSURE_GRADIENT_REGISTRY: Dict[str, Type[BasePressureGradient]] = {}
VISCOSITY_REGISTRY: Dict[str, Type[BaseViscosityOperator]] = {}
FRICTION_REGISTRY: Dict[str, Type[BaseFriction]] = {}


def create_momentum_advection(name: str, **kwargs) -> BaseMomentumAdvection:
    try:
        return MOMENTUM_ADVECTION_REGISTRY[name.lower()](**kwargs)
    except KeyError as exc:
        raise ValueError(f"Unknown momentum advection scheme '{name}'") from exc


def create_pressure_gradient(name: str, **kwargs) -> BasePressureGradient:
    try:
        return PRESSURE_GRADIENT_REGISTRY[name.lower()](**kwargs)
    except KeyError as exc:
        raise ValueError(f"Unknown pressure gradient scheme '{name}'") from exc


def create_viscosity_operator(name: str, **kwargs) -> BaseViscosityOperator:
    try:
        return VISCOSITY_REGISTRY[name.lower()](**kwargs)
    except KeyError as exc:
        raise ValueError(f"Unknown viscosity operator '{name}'") from exc


def create_friction(name: str, **kwargs) -> BaseFriction:
    try:
        return FRICTION_REGISTRY[name.lower()](**kwargs)
    except KeyError as exc:
        raise ValueError(f"Unknown friction scheme '{name}'") from exc


__all__ = [
    "create_momentum_advection",
    "create_pressure_gradient",
    "create_viscosity_operator",
    "create_friction",
    "MOMENTUM_ADVECTION_REGISTRY",
    "PRESSURE_GRADIENT_REGISTRY",
    "VISCOSITY_REGISTRY",
    "FRICTION_REGISTRY",
]
