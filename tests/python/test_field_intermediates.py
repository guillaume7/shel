import numpy as np
from shel.model.diagnostics import FieldDiagnostics
from shel.model.state import ModelState


def _base_config():
    return {
        "grid": {"nx": 6, "ny": 4, "dx": 500.0, "dy": 400.0},
        "model": {"timestep": 30.0, "gravity": 9.81, "coriolis_parameter": 1e-4},
    }


def test_intermediates_shapes_and_relations():
    cfg = _base_config()
    state = ModelState(cfg)
    depth = np.ones((cfg["grid"]["ny"], cfg["grid"]["nx"])) * 20.0
    state.set_bathymetry(depth)
    # Simple linear shear flow in y for u, zero v
    u = np.zeros((cfg["grid"]["ny"], cfg["grid"]["nx"] + 1))
    for j in range(cfg["grid"]["ny"]):
        u[j, :] = 0.02 * j
    state.set_initial_velocities(u=u)

    inter = FieldDiagnostics.intermediates(state.u, state.v, state.grid)
    ny, nx = cfg["grid"]["ny"], cfg["grid"]["nx"]
    # Check presence of keys
    expected = {"curl_w", "shearrate_w", "strechrate_t", "divergence_t", "okubo_weiss"}
    assert expected.issubset(inter.keys())
    # Shape checks
    assert inter["strechrate_t"].shape == (ny, nx)
    assert inter["divergence_t"].shape == (ny, nx)
    assert inter["okubo_weiss"].shape == (ny, nx)
    # Physical sanity: pure shear in u should give zero divergence (approximately exactly)
    assert np.allclose(inter["divergence_t"], 0.0)
    # Enstrophy_w non-negative
    assert np.min(inter["enstrophy_w"]) >= -1e-12
    # Okubo-Weiss finite
    assert np.isfinite(inter["okubo_weiss"]).all()
