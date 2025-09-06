"""
Config-aware convenience stepper for the explicit ministep.

Resolves boundary-condition type from a configuration mapping and forwards
to `explicit_step`.
"""

from __future__ import annotations

from typing import Any, Dict, Mapping

import numpy as np

from shel.model.boundary_conditions import get_bc
from shel.model.boundary_conditions.registry import get_tracer_bc

from .ministep import explicit_step

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
    for k in sides:
        val = str(bc_cfg.get(k, "closed")).lower()
        if val in ("closed", "freeslip", "radiative", "flather"):
            sides[k] = val
        else:
            sides[k] = "closed"
    return sides


def _get_eta_bc_params(config: Mapping[str, Any] | None):
    eta_ext_map = None
    eta_bc_stage = "post"
    eta_relax = 1.0
    if isinstance(config, dict):
        eta_ext_map = config.get("boundary_eta_ext")
        eta_bc_stage = str(config.get("eta_bc_stage", "post")).lower()
        try:
            eta_relax = float(config.get("boundary_eta_relax", 1.0))
        except Exception:
            eta_relax = 1.0
    return eta_ext_map, eta_bc_stage, eta_relax


def apply_eta_bc_per_side(
    eta_field: Array,
    *,
    bc_sides: Mapping[str, str],
    eta_old: Array,
    H: Array,
    g: float,
    dt: float,
    dx: float,
    dy: float,
    eta_ext_map: Mapping[str, Array] | None,
    eta_relax: float,
) -> None:
    for side, bct in bc_sides.items():
        _, e_cls = get_bc(bct)
        if e_cls is not None:
            e_cls().apply_side_eta(
                eta_field,
                side,
                eta_old=eta_old,
                eta_ext=(
                    None if not isinstance(eta_ext_map, dict) else eta_ext_map.get(side)
                ),
                H=H,
                g=g,
                dt=dt,
                dx=dx,
                dy=dy,
                relax=eta_relax,
            )


def apply_momentum_bc_per_side(
    U_field: Array,
    V_field: Array,
    *,
    bc_sides: Mapping[str, str],
    U_old: Array,
    V_old: Array,
    eta_old: Array,
    H: Array,
    g: float,
    dt: float,
    dx: float,
    dy: float,
    eta_ext_map: Mapping[str, Array] | None,
    mom_relax: float,
) -> None:
    for side, bct in bc_sides.items():
        m_cls, _ = get_bc(bct)
        if m_cls is not None:
            m_cls().apply_side(
                U_field,
                V_field,
                side,
                U_old=U_old,
                V_old=V_old,
                H=H,
                g=g,
                dt=dt,
                dx=dx,
                dy=dy,
                eta_old=eta_old,
                eta_ext=(
                    None if not isinstance(eta_ext_map, dict) else eta_ext_map.get(side)
                ),
                relax=mom_relax,
            )


def apply_sponge_layer(
    eta_field: Array,
    U_field: Array,
    V_field: Array,
    *,
    bc_sides: Mapping[str, str],
    eta_ext_map: Mapping[str, Array] | None,
    config: Mapping[str, Any] | None,
) -> None:
    if not isinstance(config, dict):
        return
    sponge = config.get("sponge") or {}
    if not isinstance(sponge, dict) or not sponge.get("enabled", False):
        return
    width = int(sponge.get("width", 0))
    alpha = float(sponge.get("alpha", 0.0))
    taper = str(sponge.get("taper", "cosine")).lower()
    apply_to = str(sponge.get("apply_to", "both")).lower()
    if width <= 0 or alpha <= 0.0:
        return

    def w(k: int) -> float:
        s = min(1.0, max(0.0, k / float(width)))
        if taper == "linear":
            return alpha * (1.0 - s)
        return alpha * 0.5 * (1.0 + np.cos(np.pi * s))

    if eta_ext_map and (apply_to in ("eta", "both")):
        if (
            "west" in bc_sides
            and isinstance(eta_ext_map, dict)
            and eta_ext_map.get("west") is not None
        ):
            ext = eta_ext_map["west"]
            for k in range(0, width):
                eta_field[:, k] = (1 - w(k)) * eta_field[:, k] + w(k) * ext[:, 0]
        if (
            "east" in bc_sides
            and isinstance(eta_ext_map, dict)
            and eta_ext_map.get("east") is not None
        ):
            ext = eta_ext_map["east"]
            for k in range(0, width):
                eta_field[:, -1 - k] = (1 - w(k)) * eta_field[:, -1 - k] + w(k) * ext[
                    :, -1
                ]

    if eta_ext_map and (apply_to in ("momentum", "both")):
        if "west" in bc_sides:
            for k in range(1, width):
                U_field[:, k] = (1 - w(k)) * U_field[:, k] + w(k) * U_field[:, 0]
        if "east" in bc_sides:
            for k in range(1, width):
                U_field[:, -1 - k] = (1 - w(k)) * U_field[:, -1 - k] + w(k) * U_field[
                    :, -1
                ]
        if "south" in bc_sides:
            for k in range(1, width):
                V_field[k, :] = (1 - w(k)) * V_field[k, :] + w(k) * V_field[0, :]
        if "north" in bc_sides:
            for k in range(1, width):
                V_field[-1 - k, :] = (1 - w(k)) * V_field[-1 - k, :] + w(k) * V_field[
                    -1, :
                ]


def apply_tracer_bc_per_side(
    C_field: Array,
    *,
    bc_sides: Mapping[str, str],
    C_old: Array | None = None,
    U: Array | None = None,
    V: Array | None = None,
    dt: float | None = None,
    dx: float | None = None,
    dy: float | None = None,
) -> None:
    for side, bct in bc_sides.items():
        t_cls = get_tracer_bc(bct)
        if t_cls is not None:
            t_cls().apply_side_tracer(
                C_field,
                side,
                C_old=C_old,
                U=U,
                V=V,
                dt=dt,
                dx=dx,
                dy=dy,
            )


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
    eta_ext_map, eta_bc_stage, eta_relax = _get_eta_bc_params(config)
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
        apply_eta_bc_per_side(
            eta_work,
            bc_sides=bc_sides_pre,
            eta_old=eta,
            H=H,
            g=g,
            dt=dt,
            dx=dx,
            dy=dy,
            eta_ext_map=eta_ext_map if isinstance(eta_ext_map, dict) else None,
            eta_relax=eta_relax,
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
    apply_momentum_bc_per_side(
        U_next,
        V_next,
        bc_sides=bc_sides,
        U_old=U,
        V_old=V,
        eta_old=eta,
        H=H,
        g=g,
        dt=dt,
        dx=dx,
        dy=dy,
        eta_ext_map=eta_ext_map if isinstance(eta_ext_map, dict) else None,
        mom_relax=mom_relax,
    )
    apply_eta_bc_per_side(
        eta_next,
        bc_sides=bc_sides,
        eta_old=eta,
        H=H,
        g=g,
        dt=dt,
        dx=dx,
        dy=dy,
        eta_ext_map=eta_ext_map if isinstance(eta_ext_map, dict) else None,
        eta_relax=eta_relax,
    )

    # Optional sponge layer for smoother transition (post-step)
    apply_sponge_layer(
        eta_next,
        U_next,
        V_next,
        bc_sides=bc_sides,
        eta_ext_map=eta_ext_map if isinstance(eta_ext_map, dict) else None,
        config=config,
    )
    return eta_next, U_next, V_next


__all__ = [
    "resolve_bc_type_from_config",
    "resolve_bc_sides_from_config",
    "explicit_step_with_config",
    "apply_eta_bc_per_side",
    "apply_momentum_bc_per_side",
    "apply_sponge_layer",
    "apply_tracer_bc_per_side",
]
