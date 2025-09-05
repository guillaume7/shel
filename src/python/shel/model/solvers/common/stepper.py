"""
Config-aware convenience stepper for the explicit ministep.

Resolves boundary-condition type from a configuration mapping and forwards
to `explicit_step`.
"""
from __future__ import annotations

from typing import Mapping, Any, Dict

import numpy as np

from .ministep import explicit_step
from .boundaries import apply_momentum_per_side, apply_eta_per_side

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

    Supported types: "closed", "freeslip", "radiative". Unknown entries default to
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
        if val in ("closed", "freeslip", "radiative"):
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
    apply_momentum_per_side(U_next, V_next, U, V, H, g, dt, dx, dy, bc_sides)
    # Eta radiative update uses eta before continuity; we approximate using eta from entry
    # Note: For full fidelity, ministep would need to apply eta BCs pre/post continuity consistently.
    apply_eta_per_side(eta_next, eta, H, g, dt, dx, dy, bc_sides)
    return eta_next, U_next, V_next


__all__ = [
    "resolve_bc_type_from_config",
    "resolve_bc_sides_from_config",
    "explicit_step_with_config",
]
