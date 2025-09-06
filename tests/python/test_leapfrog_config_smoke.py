import numpy as np

from shel.model.solvers.time import leapfrog_step_with_config


def build_flather_east_with_sponge_config(ny: int, nx: int,
                                           *,
                                           eta_east_value: float = 0.05,
                                           eta_relax: float = 0.5,
                                           momentum_relax: float = 0.3,
                                           sponge_width: int = 3,
                                           sponge_alpha: float = 0.2):
    eta_ext = np.zeros((ny, nx))
    eta_ext[:, -1] = eta_east_value
    return {
        "boundary_conditions": {
            "west": "closed",
            "east": "flather",
            "south": "closed",
            "north": "closed",
        },
        "boundary_eta_ext": {"east": eta_ext},
        "eta_bc_stage": "post",
        "boundary_eta_relax": float(eta_relax),
        "boundary_momentum_relax": float(momentum_relax),
        "sponge": {
            "enabled": True,
            "width": int(sponge_width),
            "alpha": float(sponge_alpha),
            "taper": "cosine",
            "apply_to": "both",
        },
    }


def test_leapfrog_with_config_flather_east_and_sponge_smoke():
    ny, nx = 40, 60
    dx = dy = 1.0
    dt = 0.05
    g = 9.81

    # Constant depth, no rotation, no viscosity/drag/advection
    H = np.ones((ny, nx)) * 10.0

    eta0 = np.zeros((ny, nx))
    U0 = np.zeros((ny, nx + 1))
    V0 = np.zeros((ny + 1, nx))

    cfg = build_flather_east_with_sponge_config(ny, nx,
                                                eta_east_value=0.04,
                                                eta_relax=0.6,
                                                momentum_relax=0.4,
                                                sponge_width=4,
                                                sponge_alpha=0.2)

    # Startup: duplicate initial to emulate n-1 and n states
    eta_nm1, U_nm1, V_nm1 = eta0.copy(), U0.copy(), V0.copy()
    eta_n, U_n, V_n = eta0.copy(), U0.copy(), V0.copy()

    # Run a handful of steps; boundary/sponge will impose east-side eta
    steps = 60
    for _ in range(steps):
        eta_np1, U_np1, V_np1, eta_n_f, U_n_f, V_n_f = leapfrog_step_with_config(
            eta_nm1, eta_n, H, U_nm1, U_n, V_nm1, V_n,
            dt=dt, dx=dx, dy=dy, g=g, r=0.0, nu_visc=0.0,
            enable_advection=False, f=None, enable_coriolis=False,
            config=cfg, asselin_nu=0.02,
        )
        # advance time levels
        eta_nm1, U_nm1, V_nm1 = eta_n, U_n, V_n
        eta_n, U_n, V_n = eta_np1, U_np1, V_np1

    # Sanity: no NaNs and east column influenced by boundary > interior average
    assert np.isfinite(eta_n).all()
    assert np.isfinite(U_n).all()
    assert np.isfinite(V_n).all()

    interior_mean = eta_n[:, 2:-2].mean()
    east_edge_mean = eta_n[:, -2:].mean()

    assert east_edge_mean > interior_mean - 1e-12
    # Some small positive influence expected toward external value
    assert east_edge_mean >= 0.0
