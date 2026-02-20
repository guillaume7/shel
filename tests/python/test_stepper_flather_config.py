import numpy as np

from shel.model.solvers.common.stepper import explicit_step_with_config


def test_stepper_flather_east_post_stage_smoke():
    ny, nx = 8, 9
    H = np.full((ny, nx), 10.0)
    eta = np.zeros((ny, nx))
    U = np.zeros((ny, nx + 1))
    V = np.zeros((ny + 1, nx))
    dt = 0.1
    dx = dy = 1.0

    eta_ext = np.zeros_like(eta)
    eta_ext[:, -1] = 0.1

    config = {
        "boundary_conditions": {
            "west": "closed",
            "east": "flather",
            "south": "closed",
            "north": "closed",
        },
        "boundary_eta_ext": {"east": eta_ext},
        "eta_bc_stage": "post",
        "boundary_eta_relax": 1.0,
    }

    eta1, U1, V1, _ = explicit_step_with_config(
        eta, H, U, V, dt=dt, dx=dx, dy=dy, config=config
    )

    # East U should have been nudged from zero (no pressure gradient, so any non-zero must come from BC)
    assert np.any(np.abs(U1[:, -1]) > 0.0)


def test_stepper_eta_pre_stage_applies_dirichlet():
    ny, nx = 6, 6
    H = np.full((ny, nx), 10.0)
    eta = np.zeros((ny, nx))
    U = np.zeros((ny, nx + 1))
    V = np.zeros((ny + 1, nx))
    dt = 0.05
    dx = dy = 1.0

    eta_ext = np.zeros_like(eta)
    eta_ext[:, 0] = 0.2

    config = {
        "boundary_conditions": {
            "west": "flather",
            "east": "closed",
            "south": "closed",
            "north": "closed",
        },
        "boundary_eta_ext": {"west": eta_ext},
        "eta_bc_stage": "pre",
        "boundary_eta_relax": 1.0,
    }

    eta1, U1, V1, _ = explicit_step_with_config(
        eta, H, U, V, dt=dt, dx=dx, dy=dy, config=config
    )

    # Since we applied dirichlet pre-continuity, eta after step should reflect boundary value
    assert np.allclose(eta1[:, 0], eta_ext[:, 0])
