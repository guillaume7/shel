from __future__ import annotations

from typing import Dict, Type

from .base import EtaBC, MomentumBC, TracerBC

_MOMENTUM: Dict[str, Type[MomentumBC]] = {}
_ETA: Dict[str, Type[EtaBC]] = {}
_TRACER: Dict[str, Type[TracerBC]] = {}


def register_bc(
    name: str,
    *,
    momentum: Type[MomentumBC] | None = None,
    eta: Type[EtaBC] | None = None,
    tracer: Type[TracerBC] | None = None,
) -> None:
    key = name.lower()
    if momentum is not None:
        _MOMENTUM[key] = momentum
    if eta is not None:
        _ETA[key] = eta
    if tracer is not None:
        _TRACER[key] = tracer


def get_bc(name: str):
    key = name.lower()
    return _MOMENTUM.get(key), _ETA.get(key)


def list_bcs() -> dict:
    return {"momentum": sorted(_MOMENTUM.keys()), "eta": sorted(_ETA.keys())}


def get_tracer_bc(name: str) -> Type[TracerBC] | None:
    return _TRACER.get(name.lower())


def list_tracer_bcs() -> list[str]:
    return sorted(_TRACER.keys())
