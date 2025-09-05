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
    # Optional external eta per side and eta BC timing/relaxation
    eta_ext_map = None
    eta_bc_stage = "post"
    eta_relax = 1.0
    if isinstance(config, dict):
        eta_ext_map = config.get("boundary_eta_ext")  # dict side->2D array
        eta_bc_stage = str(config.get("eta_bc_stage", "post")).lower()
        try:
            eta_relax = float(config.get("boundary_eta_relax", 1.0))
        except Exception:
            eta_relax = 1.0
    # Momentum relaxation (gamma for Flather)
    mom_relax = 1.0
    if isinstance(config, dict):
        try:
            mom_relax = float(config.get("boundary_momentum_relax", 1.0))
        except Exception:
            mom_relax = 1.0

    # Optionally apply eta BCs before the continuity update by modifying the working eta
    eta_input = eta
    if eta_bc_stage == "pre":
        eta_work = np.array(eta, copy=True)
        bc_sides_pre = resolve_bc_sides_from_config(config)
        for side, bct in bc_sides_pre.items():
            _, e_cls = get_bc(bct)
            if e_cls is not None:
                e_cls().apply_side_eta(
                    eta_work,
                    side,
                    eta_old=eta,
                    eta_ext=None if not isinstance(eta_ext_map, dict) else eta_ext_map.get(side),
                    H=H,
                    g=g,
                    dt=dt,
                    dx=dx,
                    dy=dy,
                    relax=eta_relax,
                )
        eta_input = eta_work
    eta_next, U_next, V_next = explicit_step(
        eta_input,
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
    # Apply per-side overrides if present (e.g., radiative or flather on one boundary)
    bc_sides = resolve_bc_sides_from_config(config)

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
                relax=mom_relax,
            )
    # Eta per-side post continuity (always enforce, even if applied pre)
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
                relax=eta_relax,
            )

    # Optional sponge layer for smoother transition (post-step)
    if isinstance(config, dict):
        sponge = config.get("sponge") or {}
        if isinstance(sponge, dict) and sponge.get("enabled", False):
            width = int(sponge.get("width", 0))
            alpha = float(sponge.get("alpha", 0.0))
            taper = str(sponge.get("taper", "cosine")).lower()
            apply_to = str(sponge.get("apply_to", "both")).lower()
            if width > 0 and alpha > 0.0:
                def w(k: int) -> float:
                    s = min(1.0, max(0.0, k / float(width)))
                    if taper == "linear":
                        return alpha * (1.0 - s)
                    # cosine taper by default
                    return alpha * 0.5 * (1.0 + np.cos(np.pi * s))

                # west/east sponge
                if eta_ext_map and (apply_to in ("eta", "both")):
                    if "west" in bc_sides and isinstance(eta_ext_map, dict) and eta_ext_map.get("west") is not None:
                        ext = eta_ext_map["west"]
                        for k in range(0, width):
                            eta_next[:, k] = (1 - w(k)) * eta_next[:, k] + w(k) * ext[:, 0]
                    if "east" in bc_sides and isinstance(eta_ext_map, dict) and eta_ext_map.get("east") is not None:
                        ext = eta_ext_map["east"]
                        for k in range(0, width):
                            eta_next[:, -1 - k] = (1 - w(k)) * eta_next[:, -1 - k] + w(k) * ext[:, -1]
                if eta_ext_map and (apply_to in ("momentum", "both")):
                    # relax interior normal velocity towards boundary value to smooth gradients
                    if "west" in bc_sides:
                        for k in range(1, width):
                            U_next[:, k] = (1 - w(k)) * U_next[:, k] + w(k) * U_next[:, 0]
                    if "east" in bc_sides:
                        for k in range(1, width):
                            U_next[:, -1 - k] = (1 - w(k)) * U_next[:, -1 - k] + w(k) * U_next[:, -1]
                    if "south" in bc_sides:
                        for k in range(1, width):
                            V_next[k, :] = (1 - w(k)) * V_next[k, :] + w(k) * V_next[0, :]
                    if "north" in bc_sides:
                        for k in range(1, width):
                            V_next[-1 - k, :] = (1 - w(k)) * V_next[-1 - k, :] + w(k) * V_next[-1, :]
    return eta_next, U_next, V_next


__all__ = [
    "resolve_bc_type_from_config",
    "resolve_bc_sides_from_config",
    "explicit_step_with_config",
]
