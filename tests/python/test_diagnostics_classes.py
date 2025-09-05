"""Tests for IntegratedDiagnostics and FieldDiagnostics classes."""

import numpy as np

from shel.model.diagnostics import IntegratedDiagnostics, FieldDiagnostics
from shel.model.state import ModelState


def _base_config():
    return {
        "grid": {"nx": 8, "ny": 6, "dx": 500.0, "dy": 500.0},
        "model": {"timestep": 30.0, "gravity": 9.81, "coriolis_parameter": 1e-4},
    }


def test_integrated_diagnostics_basic():
    cfg = _base_config()
    state = ModelState(cfg)
    # Bathymetry and elevation
    depth = np.ones((cfg["grid"]["ny"], cfg["grid"]["nx"])) * 50.0
    state.set_bathymetry(depth)
    eta = np.zeros_like(depth)
    eta[3, 4] = 0.2
    state.set_initial_elevation(eta)
    # Uniform small velocity
    u = np.ones((cfg["grid"]["ny"], cfg["grid"]["nx"] + 1)) * 0.1
    v = np.zeros((cfg["grid"]["ny"] + 1, cfg["grid"]["nx"]))
    state.set_initial_velocities(u, v)

    integ = IntegratedDiagnostics.as_dict(
        state.u,
        state.v,
        state.eta,
        state.H,
        state.grid,
        state.gravity,
        state.coriolis,
    )

    expected_keys = {
        "kinetic_energy",
        "potential_energy",
        "total_energy",
        "volume",
        "enstrophy",
        "potential_enstrophy",
    }
    assert set(integ.keys()) == expected_keys
    assert integ["kinetic_energy"] > 0.0
    assert integ["potential_energy"] > 0.0
    assert np.isclose(
        integ["total_energy"], integ["kinetic_energy"] + integ["potential_energy"], rtol=1e-10
    )
    assert integ["volume"] > 0.0


def test_field_diagnostics_basic():
    cfg = _base_config()
    state = ModelState(cfg)
    depth = np.ones((cfg["grid"]["ny"], cfg["grid"]["nx"])) * 30.0
    state.set_bathymetry(depth)
    # Shear in u
    u = np.zeros((cfg["grid"]["ny"], cfg["grid"]["nx"] + 1))
    for j in range(cfg["grid"]["ny"]):
        u[j, :] = 0.05 * j
    state.set_initial_velocities(u=u)

    fields = FieldDiagnostics.as_dict(state.u, state.v, state.grid)
    # Allow additional diagnostics, but require these two at minimum
    assert {"vorticity", "okubo_weiss"}.issubset(set(fields.keys()))
    vort = fields["vorticity"]
    ow = fields["okubo_weiss"]
    assert vort.shape == (cfg["grid"]["ny"] + 1, cfg["grid"]["nx"] + 1)
    assert ow.shape == (cfg["grid"]["ny"], cfg["grid"]["nx"])
    # With pure shear expect some non-zero vorticity interior
    assert np.any(np.abs(vort[1:-1, 1:-1]) > 0.0)


def test_modelstate_to_dict_structure():
    cfg = _base_config()
    state = ModelState(cfg)
    depth = np.ones((cfg["grid"]["ny"], cfg["grid"]["nx"])) * 10.0
    state.set_bathymetry(depth)
    dct = state.to_dict()

    assert "diagnostics" in dct
    diag = dct["diagnostics"]
    assert set(diag.keys()) == {"integrated", "fields"}
    assert all(k in diag["integrated"] for k in ["kinetic_energy", "volume", "enstrophy"])
    field_arrays = diag["fields"]
    # Ensure lists (serialized) and correct shapes
    assert isinstance(field_arrays["vorticity"], list)
    assert isinstance(field_arrays["okubo_weiss"], list)
    # Spot-check outer list lengths
    assert len(field_arrays["okubo_weiss"]) == cfg["grid"]["ny"]
    assert len(field_arrays["vorticity"]) == cfg["grid"]["ny"] + 1
