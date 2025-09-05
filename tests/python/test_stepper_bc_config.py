import numpy as np

from shel.model.solvers.common.stepper import resolve_bc_type_from_config, explicit_step_with_config


def test_resolve_bc_type_from_config_uniform_closed():
    cfg = {"boundary_conditions": {"west": "closed", "east": "closed", "south": "closed", "north": "closed"}}
    assert resolve_bc_type_from_config(cfg) == "closed"


def test_resolve_bc_type_from_config_uniform_freeslip():
    cfg = {"boundary_conditions": {"west": "freeslip", "east": "freeslip", "south": "freeslip", "north": "freeslip"}}
    assert resolve_bc_type_from_config(cfg) == "freeslip"


def test_resolve_bc_type_from_config_mixed_defaults_to_closed():
    cfg = {"boundary_conditions": {"west": "closed", "east": "freeslip", "south": "closed", "north": "closed"}}
    assert resolve_bc_type_from_config(cfg) == "closed"


def test_explicit_step_with_config_applies_bc():
    ny, nx = 8, 9
    dx = dy = 1.0
    dt = 0.02

    H = np.ones((ny, nx))
    eta = np.zeros((ny, nx))

    U = 0.1 * np.ones((ny, nx + 1))
    V = -0.1 * np.ones((ny + 1, nx))

    # Freeslip config: should zero normal components at edges
    cfg_fs = {"boundary_conditions": {"west": "freeslip", "east": "freeslip", "south": "freeslip", "north": "freeslip"}}
    eta1, U1, V1 = explicit_step_with_config(
        eta,
        H,
        U,
        V,
        dt=dt,
        dx=dx,
        dy=dy,
        g=9.81,
        r=0.0,
        nu=0.0,
        enable_advection=False,
        f=None,
        enable_coriolis=False,
        config=cfg_fs,
    )
    assert np.allclose(U1[:, 0], 0.0) and np.allclose(U1[:, -1], 0.0)
    assert np.allclose(V1[0, :], 0.0) and np.allclose(V1[-1, :], 0.0)

    # Mixed config: falls back to closed -> same edge zeroing
    cfg_mixed = {"boundary_conditions": {"west": "closed", "east": "freeslip", "south": "closed", "north": "closed"}}
    eta2, U2, V2 = explicit_step_with_config(
        eta,
        H,
        U,
        V,
        dt=dt,
        dx=dx,
        dy=dy,
        g=9.81,
        r=0.0,
        nu=0.0,
        enable_advection=False,
        f=None,
        enable_coriolis=False,
        config=cfg_mixed,
    )
    assert np.allclose(U2[:, 0], 0.0) and np.allclose(U2[:, -1], 0.0)
    assert np.allclose(V2[0, :], 0.0) and np.allclose(V2[-1, :], 0.0)
