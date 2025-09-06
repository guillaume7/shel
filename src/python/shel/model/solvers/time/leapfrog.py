"""
Leapfrog + Asselin filter integrator built on existing tendencies.

Depends on:
- explicit_step: assembles tendencies and applies momentum/eta updates for one Euler step
- asselin_filter: Robert–Asselin filter applied to the middle time level
"""

from __future__ import annotations

from typing import Any, Mapping

import numpy as np

from ..common.ministep import explicit_step
from ..common.stepper import (
    apply_eta_bc_per_side,
    apply_momentum_bc_per_side,
    apply_sponge_layer,
    resolve_bc_sides_from_config,
    resolve_bc_type_from_config,
)
from .asselin import asselin_filter

Array = np.ndarray


def leapfrog_stepper(
    eta_nm1: Array,
    eta_n: Array,
    H: Array,
    U_nm1: Array,
    U_n: Array,
    V_nm1: Array,
    V_n: Array,
    *,
    dt: float,
    dx: float,
    dy: float,
    g: float = 9.81,
    r: float = 0.0,
    nu_visc: float = 0.0,
    enable_advection: bool = True,
    f: Array | None = None,
    enable_coriolis: bool = False,
    bc_type: str = "closed",
    asselin_nu: float = 0.02,
):
    """Advance one leapfrog step using centered time staggering and Asselin filter.

    Returns (eta_np1, U_np1, V_np1, eta_n_f, U_n_f, V_n_f).
    """
    eta_np1, U_np1, V_np1 = explicit_step(
        eta_n,
        H,
        U_n,
        V_n,
        dt=dt,
        dx=dx,
        dy=dy,
        g=g,
        r=r,
        nu=nu_visc,
        enable_advection=enable_advection,
        f=f,
        enable_coriolis=enable_coriolis,
        bc_type=bc_type,
    )
    # Apply Robert–Asselin filter to middle state (n)
    eta_n_f = asselin_filter(eta_nm1, eta_n, eta_np1, asselin_nu)
    U_n_f = asselin_filter(U_nm1, U_n, U_np1, asselin_nu)
    V_n_f = asselin_filter(V_nm1, V_n, V_np1, asselin_nu)
    return eta_np1, U_np1, V_np1, eta_n_f, U_n_f, V_n_f


__all__ = ["leapfrog_stepper", "leapfrog_step_with_config"]


def leapfrog_step_with_config(
    eta_nm1: Array,
    eta_n: Array,
    H: Array,
    U_nm1: Array,
    U_n: Array,
    V_nm1: Array,
    V_n: Array,
    *,
    dt: float,
    dx: float,
    dy: float,
    g: float = 9.81,
    r: float = 0.0,
    nu_visc: float = 0.0,
    enable_advection: bool = True,
    f: Array | None = None,
    enable_coriolis: bool = False,
    config: Mapping[str, Any] | None = None,
    asselin_nu: float = 0.02,
):
    """Leapfrog stepper with config-aware BC and sponge application (DRY with stepper.py)."""
    # Use base step to get n+1
    bc_type = resolve_bc_type_from_config(config)
    eta_np1, U_np1, V_np1 = explicit_step(
        eta_n,
        H,
        U_n,
        V_n,
        dt=dt,
        dx=dx,
        dy=dy,
        g=g,
        r=r,
        nu=nu_visc,
        enable_advection=enable_advection,
        f=f,
        enable_coriolis=enable_coriolis,
        bc_type=bc_type,
    )

    # Apply per-side BCs (momentum then eta, post enforcement)
    eta_ext_map = None
    eta_relax = 1.0
    mom_relax = 1.0
    if isinstance(config, dict):
        eta_ext_map = config.get("boundary_eta_ext")
        try:
            eta_relax = float(config.get("boundary_eta_relax", 1.0))
        except Exception:
            eta_relax = 1.0
        try:
            mom_relax = float(config.get("boundary_momentum_relax", 1.0))
        except Exception:
            mom_relax = 1.0
    bc_sides = resolve_bc_sides_from_config(config)
    apply_momentum_bc_per_side(
        U_np1,
        V_np1,
        bc_sides=bc_sides,
        U_old=U_n,
        V_old=V_n,
        eta_old=eta_n,
        H=H,
        g=g,
        dt=dt,
        dx=dx,
        dy=dy,
        eta_ext_map=eta_ext_map if isinstance(eta_ext_map, dict) else None,
        mom_relax=mom_relax,
    )
    apply_eta_bc_per_side(
        eta_np1,
        bc_sides=bc_sides,
        eta_old=eta_n,
        H=H,
        g=g,
        dt=dt,
        dx=dx,
        dy=dy,
        eta_ext_map=eta_ext_map if isinstance(eta_ext_map, dict) else None,
        eta_relax=eta_relax,
    )

    # Sponge post-step
    apply_sponge_layer(
        eta_np1,
        U_np1,
        V_np1,
        bc_sides=bc_sides,
        eta_ext_map=eta_ext_map if isinstance(eta_ext_map, dict) else None,
        config=config,
    )

    # Asselin filter middle state
    eta_n_f = asselin_filter(eta_nm1, eta_n, eta_np1, asselin_nu)
    U_n_f = asselin_filter(U_nm1, U_n, U_np1, asselin_nu)
    V_n_f = asselin_filter(V_nm1, V_n, V_np1, asselin_nu)
    return eta_np1, U_np1, V_np1, eta_n_f, U_n_f, V_n_f
