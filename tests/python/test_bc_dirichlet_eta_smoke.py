import numpy as np

from shel.model.solvers.common.stepper import explicit_step_with_config


def test_dirichlet_eta_east_smoke():
    ny, nx = 12, 30
    dx = dy = 1.0
    dt = 0.05
    g = 9.81

    H = np.ones((ny, nx)) * 10.0

    eta = np.zeros((ny, nx))
    U = np.zeros((ny, nx + 1))
    V = np.zeros((ny + 1, nx))

    eta_ext = np.zeros((ny, nx))
    eta_ext[:, -1] = 0.1

    cfg = {
        "boundary_conditions": {"west": "closed", "east": "dirichlet", "south": "closed", "north": "closed"},
        "boundary_eta_ext": {"east": eta_ext},
        "boundary_eta_relax": 0.5,
    }

    for _ in range(5):
        eta, U, V = explicit_step_with_config(
            eta, H, U, V, dt=dt, dx=dx, dy=dy, g=g, r=0.0, nu=0.0,
            enable_advection=False, f=None, enable_coriolis=False, config=cfg,
        )
        assert np.isfinite(eta).all()
        assert np.isfinite(U).all()
        assert np.isfinite(V).all()

    edge_mean = eta[:, -1].mean()
    assert edge_mean >= 0.0
