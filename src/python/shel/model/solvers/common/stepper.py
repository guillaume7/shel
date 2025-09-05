"""
Config-aware convenience stepper for the explicit ministep.

Resolves boundary-condition type from a configuration mapping and forwards
to `explicit_step`.
"""
from __future__ import annotations

from typing import Mapping, Any, Dict

import numpy as np

from .ministep import explicit_step
from shel.model.boundary_conditions import get_bc

Array = np.ndarray


def resolve_bc_type_from_config(config: Mapping[str, Any] | None) -> str:
    """Resolve a uniform bc_type from a config mapping.

    Supported returns: "closed", "freeslip". If boundaries differ or
    unsupported types are found, default to "closed" conservatively.
    """
    if not config:
        return "closed"
    bc_cfg = config.get("boundary_conditions") if isinstance(config, dict) else None
    if not isinstance(bc_cfg, dict):
        return "closed"
    sides = [bc_cfg.get(k) for k in ("west", "east", "south", "north")]
    sides_norm = [str(s).lower() for s in sides if s is not None]
    if not sides_norm:
        return "closed"
    uniq = set(sides_norm)
    if len(uniq) == 1 and list(uniq)[0] in ("closed", "freeslip"):
        return list(uniq)[0]
    # Mixed or unsupported -> fallback to closed for safety
    return "closed"


def resolve_bc_sides_from_config(config: Mapping[str, Any] | None) -> Dict[str, str]:
    """Return a per-side BC mapping: {west,east,south,north} -> type.

    Supported types: "closed", "freeslip", "radiative", "flather". Unknown entries default to
    "closed". Missing config returns all-closed.
    """
    sides = {"west": "closed", "east": "closed", "south": "closed", "north": "closed"}
    if not config or not isinstance(config, dict):
        return sides
    bc_cfg = config.get("boundary_conditions")
    if not isinstance(bc_cfg, dict):
        return sides
    for k in sides.keys():
        val = str(bc_cfg.get(k, "closed")).lower()
        if val in ("closed", "freeslip", "radiative", "flather"):
            sides[k] = val
        else:
            sides[k] = "closed"
    return sides


def explicit_step_with_config(
    eta: Array,
    H: Array,
    U: Array,
    V: Array,
    *,
    dt: float,
    dx: float,
    dy: float,
    g: float = 9.81,
    r: float = 0.0,
    nu: float = 0.0,
    enable_advection: bool = True,
    f: Array | None = None,
    enable_coriolis: bool = False,
    config: Mapping[str, Any] | None = None,
):
    bc_type = resolve_bc_type_from_config(config)
    eta_next, U_next, V_next = explicit_step(
        eta,
        H,
        U,
        V,
        dt=dt,
        dx=dx,
        dy=dy,
        g=g,
        r=r,
        nu=nu,
        enable_advection=enable_advection,
        f=f,
        enable_coriolis=enable_coriolis,
        bc_type=bc_type,
    )
    # Apply per-side overrides if present (e.g., radiative on one boundary)
    bc_sides = resolve_bc_sides_from_config(config)
    # Optional external eta per side for open boundaries (prototype input name)
    eta_ext_map = None
    if isinstance(config, dict):
        eta_ext_map = config.get("boundary_eta_ext")  # expects dict side->2D array

    # Momentum per-side
    for side, bct in bc_sides.items():
        m_cls, _ = get_bc(bct)
        if m_cls is not None:
            m_cls().apply_side(
                U_next,
                V_next,
                side,
                U_old=U,
                V_old=V,
                H=H,
                g=g,
                dt=dt,
                dx=dx,
                dy=dy,
                eta_old=eta,
                eta_ext=None if not isinstance(eta_ext_map, dict) else eta_ext_map.get(side),
            )
    # Eta per-side (radiative only for now)
    for side, bct in bc_sides.items():
        _, e_cls = get_bc(bct)
    if e_cls is not None:
            e_cls().apply_side_eta(
                eta_next,
                side,
                eta_old=eta,
        eta_ext=None if not isinstance(eta_ext_map, dict) else eta_ext_map.get(side),
                H=H,
                g=g,
                dt=dt,
                dx=dx,
                dy=dy,
            )
    return eta_next, U_next, V_next


__all__ = [
    "resolve_bc_type_from_config",
    "resolve_bc_sides_from_config",
    "explicit_step_with_config",
]
