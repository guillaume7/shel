from __future__ import annotations

from typing import Dict, Type

from .base import MomentumBC, EtaBC


_MOMENTUM: Dict[str, Type[MomentumBC]] = {}
_ETA: Dict[str, Type[EtaBC]] = {}


def register_bc(name: str, *, momentum: Type[MomentumBC] | None = None, eta: Type[EtaBC] | None = None) -> None:
    key = name.lower()
    if momentum is not None:
        _MOMENTUM[key] = momentum
    if eta is not None:
        _ETA[key] = eta


def get_bc(name: str):
    key = name.lower()
    return _MOMENTUM.get(key), _ETA.get(key)


def list_bcs() -> dict:
    return {"momentum": sorted(_MOMENTUM.keys()), "eta": sorted(_ETA.keys())}
