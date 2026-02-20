"""
Leapfrog + Asselin filter integrator (conservative form).

Matches MATLAB ``ComputeLeapfrog``:
    1. Compute RHS at time level n from (eta_n, U_n, V_n)
    2. Advance from n-1 using 2*dt (centered leapfrog)
    3. Apply Asselin filter to middle level n

Depends on:
- explicit_step: conservative flux-form step (accepts eta_old, U_old, etc.)
- asselin_filter: Robert–Asselin filter
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
    asselin_nu: float = 0.1,
    d: Array | None = None,
    H_old: Array | None = None,
):
    """Advance one leapfrog step: centered in time (2*dt) with Asselin filter.

    Parameters
    ----------
    eta_nm1, U_nm1, V_nm1 : arrays at time level n-1
    eta_n, U_n, V_n       : arrays at time level n (used for RHS evaluation)
    H                     : water-column height at level n
    dt                    : physical time step (the step internally uses 2*dt)
    d                     : bathymetry depth (constant); inferred from H-eta if None
    H_old                 : water-column height at level n-1; inferred from eta_nm1+d if None

    Returns
    -------
    (eta_np1, U_np1, V_np1, eta_n_f, U_n_f, V_n_f)
    """
    if d is None:
        d = H - eta_n
    if H_old is None:
        H_old = eta_nm1 + d

    # Conservative step: RHS evaluated at n, advance from n-1 with 2*dt
    eta_np1, U_np1, V_np1, H_np1 = explicit_step(
        eta_n,  # current level (for RHS evaluation)
        H,  # current H
        U_n,
        V_n,
        dt=2.0 * dt,  # leapfrog: advance over 2*dt
        dx=dx,
        dy=dy,
        g=g,
        r=r,
        nu=nu_visc,
        enable_advection=enable_advection,
        f=f,
        enable_coriolis=enable_coriolis,
        bc_type=bc_type,
        d=d,
        eta_old=eta_nm1,
        H_old=H_old,
        U_old=U_nm1,
        V_old=V_nm1,
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
    asselin_nu: float = 0.1,
    d: Array | None = None,
    H_old: Array | None = None,
):
    """Leapfrog stepper with config-aware BC and sponge (DRY with stepper.py)."""
    bc_type = resolve_bc_type_from_config(config)
    if d is None:
        d = H - eta_n
    if H_old is None:
        H_old = eta_nm1 + d

    # Conservative step with 2*dt
    eta_np1, U_np1, V_np1, H_np1 = explicit_step(
        eta_n,
        H,
        U_n,
        V_n,
        dt=2.0 * dt,
        dx=dx,
        dy=dy,
        g=g,
        r=r,
        nu=nu_visc,
        enable_advection=enable_advection,
        f=f,
        enable_coriolis=enable_coriolis,
        bc_type=bc_type,
        d=d,
        eta_old=eta_nm1,
        H_old=H_old,
        U_old=U_nm1,
        V_old=V_nm1,
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
