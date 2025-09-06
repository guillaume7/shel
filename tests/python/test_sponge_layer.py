import numpy as np

from shel.model.solvers.common.stepper import explicit_step_with_config


def test_sponge_eta_blending_west_linear():
    ny, nx = 6, 8
    H = np.full((ny, nx), 10.0)
    eta = np.zeros((ny, nx))
    U = np.zeros((ny, nx + 1))
    V = np.zeros((ny + 1, nx))

    eta_ext = np.zeros_like(eta)
    eta_ext[:, 0] = 1.0

    cfg = {
        "boundary_conditions": {
            "west": "flather",
            "east": "closed",
            "south": "closed",
            "north": "closed",
        },
        "boundary_eta_ext": {"west": eta_ext},
        "eta_bc_stage": "post",
        "boundary_eta_relax": 1.0,
        "sponge": {
            "enabled": True,
            "width": 3,
            "alpha": 0.5,
            "taper": "linear",
            "apply_to": "eta",
        },
    }

    eta1, U1, V1 = explicit_step_with_config(
        eta, H, U, V, dt=0.1, dx=1.0, dy=1.0, config=cfg
    )

    # boundary should be 1.0 (Dirichlet), next column should be blended toward 1.0
    assert np.allclose(eta1[:, 0], 1.0)
    assert np.all(eta1[:, 1] > 0.0)
    assert np.all(eta1[:, 2] >= 0.0)


def test_sponge_momentum_relaxes_toward_boundary_u_east():
    ny, nx = 6, 8
    H = np.full((ny, nx), 10.0)
    eta = np.zeros((ny, nx))
    U = np.zeros((ny, nx + 1))
    V = np.zeros((ny + 1, nx))

    eta_ext = np.zeros_like(eta)
    eta_ext[:, -1] = 0.2

    cfg = {
        "boundary_conditions": {
            "west": "closed",
            "east": "flather",
            "south": "closed",
            "north": "closed",
        },
        "boundary_eta_ext": {"east": eta_ext},
        "boundary_momentum_relax": 1.0,
        "sponge": {
            "enabled": True,
            "width": 3,
            "alpha": 0.5,
            "taper": "cosine",
            "apply_to": "momentum",
        },
    }

    eta1, U1, V1 = explicit_step_with_config(
        eta, H, U, V, dt=0.1, dx=1.0, dy=1.0, config=cfg
    )

    # Boundary U at east should be non-zero, and interior adjacent column should be between 0 and boundary value
    ub = np.abs(U1[:, -1])
    assert np.any(ub > 0)
    u_inner = np.abs(U1[:, -2])
    assert np.all(u_inner <= ub + 1e-12)
    assert np.any(u_inner > 0)
