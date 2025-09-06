from __future__ import annotations

import numpy as np

from shel.model.diagnostics import (
    FieldDiagnostics,
    IntegratedDiagnostics,
    potential_vorticity,
)
from shel.model.grid import Grid
from shel.model.initial_conditions.factory import build_initial_state


def make_grid():
    return Grid({"grid": {"nx": 16, "ny": 12, "dx": 1000.0, "dy": 1000.0}})


def test_repeated_diagnostics_idempotent():
    g = make_grid()
    cfg = {
        "bathymetry": {"name": "step"},
        "elevation": {"name": "flat"},
        "velocity": {"name": "solid_body"},
    }
    state = build_initial_state(cfg, g)
    u, v, eta, H, f = (
        state["u"],
        state["v"],
        state["eta"],
        state["H"],
        state["coriolis"],
    )
    diag1 = IntegratedDiagnostics.as_dict(u, v, eta, H, g, 9.81, f)
    diag2 = IntegratedDiagnostics.as_dict(u, v, eta, H, g, 9.81, f)
    for k, v in diag1.items():
        assert np.isclose(v, diag2[k], rtol=0, atol=0)
    fields1 = FieldDiagnostics.as_dict(u, v, g)
    fields2 = FieldDiagnostics.as_dict(u, v, g)
    for k, v in fields1.items():
        assert np.allclose(v, fields2[k], rtol=0, atol=0)
    # PV consistency
    pv1 = potential_vorticity(u, v, H, f, g)
    pv2 = potential_vorticity(u, v, H, f, g)
    assert np.allclose(pv1, pv2)


def test_pv_mean_matches_coriolis_over_depth_constant_case():
    g = make_grid()
    cfg = {"bathymetry": {"name": "step"}, "elevation": {"name": "flat"}}
    state = build_initial_state(cfg, g)
    u, v, H, f = state["u"], state["v"], state["H"], state["coriolis"]
    pv = potential_vorticity(u, v, H, f, g)
    # With zero velocity & flat surface, ζ ≈ 0, so pv ≈ f / H
    # Multiply back by H and compare mean to mean(f)
    mean_f_recon = np.mean(pv * H)
    assert np.isclose(mean_f_recon, np.mean(f), rtol=1e-12, atol=0)
