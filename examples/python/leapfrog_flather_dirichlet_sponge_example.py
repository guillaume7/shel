import numpy as np

from shel.model.solvers.time import leapfrog_step_with_config


def build_config(ny, nx):
    eta_ext = np.zeros((ny, nx))
    eta_ext[:, -1] = 0.05  # east external elevation
    return {
        "boundary_conditions": {
            "west": "closed",
            "east": "flather",  # momentum + eta
            "south": "closed",
            "north": "dirichlet",  # eta blend at north edge
        },
        "boundary_eta_ext": {"east": eta_ext, "north": np.tile(0.02, (ny, nx))},
        "boundary_eta_relax": 0.5,  # conservative eta blending
        "boundary_momentum_relax": 0.4,  # soften Flather correction
        "eta_bc_stage": "post",
        "sponge": {
            "enabled": True,
            "width": 4,
            "alpha": 0.2,
            "taper": "cosine",
            "apply_to": "both",
        },
    }


def main():
    ny, nx = 40, 60
    dx = dy = 1.0
    dt = 0.05
    g = 9.81

    H = np.ones((ny, nx)) * 10.0

    eta_nm1 = np.zeros((ny, nx))
    U_nm1 = np.zeros((ny, nx + 1))
    V_nm1 = np.zeros((ny + 1, nx))

    eta_n = eta_nm1.copy()
    U_n = U_nm1.copy()
    V_n = V_nm1.copy()

    cfg = build_config(ny, nx)

    steps = 100
    for _ in range(steps):
        eta_np1, U_np1, V_np1, eta_n_f, U_n_f, V_n_f = leapfrog_step_with_config(
            eta_nm1,
            eta_n,
            H,
            U_nm1,
            U_n,
            V_nm1,
            V_n,
            dt=dt,
            dx=dx,
            dy=dy,
            g=g,
            r=0.0,
            nu_visc=0.0,
            enable_advection=False,
            f=None,
            enable_coriolis=False,
            config=cfg,
            asselin_nu=0.02,
        )
        eta_nm1, U_nm1, V_nm1 = eta_n, U_n, V_n
        eta_n, U_n, V_n = eta_np1, U_np1, V_np1

    print(
        "Final eta stats: min=",
        float(eta_n.min()),
        "max=",
        float(eta_n.max()),
        "mean=",
        float(eta_n.mean()),
    )


if __name__ == "__main__":
    main()
