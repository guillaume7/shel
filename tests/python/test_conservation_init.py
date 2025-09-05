from __future__ import annotations
import numpy as np
from shel.model.grid import Grid
from shel.model.initial_conditions.factory import build_initial_state


def make_grid(nx=24, ny=18):
    return Grid({"grid": {"nx": nx, "ny": ny, "dx": 500.0, "dy": 500.0}})


def integrate(field: np.ndarray, dx: float, dy: float) -> float:
    return float(field.sum() * dx * dy)


def test_initial_volume_conservation_multiple_configs():
    g = make_grid()
    cfgs = [
        {"bathymetry": {"name": "bump", "params": {"depth0": 200.0, "amp": 10.0}}, "elevation": {"name": "flat"}},
        {"bathymetry": {"name": "step"}, "elevation": {"name": "gaussian", "params": {"amp": 0.5}}},
        {"bathymetry": {"name": "island"}, "elevation": {"name": "gaussian", "params": {"amp": 0.2}}},
    ]
    for cfg in cfgs:
        state = build_initial_state(cfg, g)
        H = state["H"]; h = state["h"]; eta = state["eta"]
        assert np.allclose(H, h + eta)
    vol_H = integrate(H, g.dx, g.dy)
    vol_parts = integrate(h, g.dx, g.dy) + integrate(eta, g.dx, g.dy)
    # Use relative tolerance due to summation order effects on large magnitudes
    assert np.isclose(vol_H, vol_parts, rtol=1e-12, atol=1e-9)


def test_tracer_mass_integrals():
    g = make_grid()
    cfg = {
        "bathymetry": {"name": "cylinder"},
        "elevation": {"name": "flat"},
        "tracers": [
            {"name": "gaussian", "key": "tr_gauss", "params": {"c0": 3.0}},
            {"name": "uniform", "key": "tr_uni", "params": {"value": 0.5}},
        ],
    }
    state = build_initial_state(cfg, g)
    gauss_mass = integrate(state["tr_gauss"], g.dx, g.dy)
    uni_mass = integrate(state["tr_uni"], g.dx, g.dy)
    # Uniform mass should equal value * domain area
    area = g.nx * g.ny * g.dx * g.dy
    assert np.isclose(uni_mass, 0.5 * area, rtol=0, atol=1e-10)
    # Gaussian positive and less than uniform mass (since localized)
    assert 0 < gauss_mass < uni_mass
