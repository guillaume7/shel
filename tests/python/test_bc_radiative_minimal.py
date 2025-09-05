import numpy as np

from shel.model.solvers.common.stepper import explicit_step_with_config


def test_radiative_east_boundary_minimal_outflow():
    # 1D-like pulse traveling east; radiative east boundary should reduce reflection
    ny, nx = 8, 40
    dx = dy = 1.0
    dt = 0.02
    g = 9.81

    H = np.ones((ny, nx)) * 10.0

    x = np.arange(nx)
    X = np.tile(x, (ny, 1))
    # Small right-going slope in eta near the right side
    eta = 0.02 * np.exp(-((X - (nx * 0.7)) ** 2) / (2.0 * 3.0 ** 2))

    U = np.zeros((ny, nx + 1))
    V = np.zeros((ny + 1, nx))

    cfg = {
        "boundary_conditions": {
            "west": "closed",
            "east": "radiative",
            "south": "closed",
            "north": "closed",
        }
    }

    # Step a few times; ensure no explosion and normal east face remains finite
    for _ in range(10):
        eta, U, V = explicit_step_with_config(
            eta,
            H,
            U,
            V,
            dt=dt,
            dx=dx,
            dy=dy,
            g=g,
            r=0.0,
            nu=0.0,
            enable_advection=False,
            f=None,
            enable_coriolis=False,
            config=cfg,
        )
        assert np.isfinite(U[:, -1]).all()
        assert np.isfinite(eta[:, -1]).all()

    # Basic smoke: values near east boundary should remain bounded and not revert to strong reflection pattern
    assert np.abs(U[:, -1]).mean() < 0.5
    assert np.abs(eta[:, -1]).mean() < 0.5
