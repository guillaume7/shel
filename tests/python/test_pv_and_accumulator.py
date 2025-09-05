from __future__ import annotations
import numpy as np
from shel.model.grid import Grid
from shel.model.initial_conditions.factory import build_initial_state
from shel.model.diagnostics import potential_vorticity, GlobalAccumulator, IntegratedDiagnostics


def make_grid(nx=20, ny=16):
    return Grid({"grid": {"nx": nx, "ny": ny, "dx": 1000.0, "dy": 1000.0}})


def test_potential_vorticity_basic_shapes():
    g = make_grid()
    cfg = {
        "bathymetry": {"name": "bump", "params": {"depth0": 100.0, "amp": 5.0}},
        "elevation": {"name": "gaussian", "params": {"amp": 0.5}},
        "velocity": {"name": "solid_body"},
        "coriolis": {"type": "constant", "value": 1e-4},
    }
    state = build_initial_state(cfg, g)
    q = potential_vorticity(state["u"], state["v"], state["H"], state["coriolis"], g)
    assert q.shape == (g.ny, g.nx)
    # PV finite & not all zero
    assert np.isfinite(q).all()
    assert not np.allclose(q, 0.0)


def test_global_accumulator_appends_and_matches_integrated():
    g = make_grid()
    cfg = {"bathymetry": {"name": "step"}, "elevation": {"name": "flat"}}
    state = build_initial_state(cfg, g)
    acc = GlobalAccumulator()
    for step in range(3):
        t = step * 10.0
        acc.update(t, state["u"], state["v"], state["eta"], state["H"], state["coriolis"], g)
    d = acc.as_dict()
    assert d["time"] == [0.0, 10.0, 20.0]
    # Cross-check last entry with direct IntegratedDiagnostics call
    di = IntegratedDiagnostics.as_dict(state["u"], state["v"], state["eta"], state["H"], g, 9.81, state["coriolis"])
    assert np.isclose(d["total_energy"][-1], di["total_energy"])  # unchanged static state
    assert np.isclose(d["volume"][0], d["volume"][1]) and np.isclose(d["volume"][1], d["volume"][2])
