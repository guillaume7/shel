from __future__ import annotations
import numpy as np
from shel.model.grid import Grid
from shel.model.initial_conditions.factory import build_initial_state

def make_grid(nx=20, ny=16):
    cfg = {"grid": {"nx": nx, "ny": ny, "dx": 1000.0, "dy": 1000.0}}
    return Grid(cfg)

def test_build_initial_state_minimal():
    g = make_grid()
    cfg = {
        "bathymetry": {"name": "bump", "params": {"depth0": 1000.0, "amp": 50.0}},
        "elevation": {"name": "gaussian", "params": {"amp": 1.0}},
        "velocity": {"name": "geostrophic"},
        "tracers": [
            {"name": "gaussian", "key": "tracer_main", "params": {"c0": 2.0}},
            {"name": "uniform", "params": {"value": 0.5}}
        ],
        "coriolis": {"type": "constant", "value": 1e-4}
    }
    state = build_initial_state(cfg, g)
    # Shape checks
    ny, nx = g.ny, g.nx
    assert state["h"].shape == (ny, nx)
    assert state["eta"].shape == (ny, nx)
    assert state["H"].shape == (ny, nx)
    assert state["u"].shape == (ny, nx + 1)
    assert state["v"].shape == (ny + 1, nx)
    assert state["coriolis"].shape == (ny, nx)
    assert state["tracer_main"].shape == (ny, nx)
    assert state["tracer_1"].shape == (ny, nx)
    # H consistency
    assert np.allclose(state["H"], state["h"] + state["eta"])
    # Tracer positivity / bounds
    assert state["tracer_main"].max() > state["tracer_main"].mean() > 0
    assert np.isclose(state["tracer_1"].min(), state["tracer_1"].max())
    # Basic mass integral > 0
    mass = state["H"].sum() * g.dx * g.dy
    assert mass > 0

def test_build_initial_state_defaults():
    g = make_grid()
    # Omitting velocity & tracers should yield zeros where expected
    cfg = {"bathymetry": {"name": "step"}, "elevation": {"name": "flat"}}
    state = build_initial_state(cfg, g)
    assert np.allclose(state["eta"], 0.0)
    assert np.allclose(state["u"], 0.0)
    assert np.allclose(state["v"], 0.0)
    assert np.allclose(state["H"], state["h"])  # depth only
