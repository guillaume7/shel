import numpy as np

from shel.model.solvers.common.ministep import explicit_step


def test_ministep_coriolis_inertial_oscillation_small_dt():
    # Small local inertial response: u accelerates v and vice versa
    ny, nx = 6, 7
    dx = dy = 1.0
    dt = 0.01
    g = 9.81

    H = np.ones((ny, nx)) * 10.0
    eta = np.zeros((ny, nx))

    # Start with small U, zero V
    U = np.zeros((ny, nx + 1))
    V = np.zeros((ny + 1, nx))
    U[:, 1:-1] = 0.2

    # Constant f-plane Coriolis parameter on T points
    f = np.ones((ny, nx)) * 1e-3

    eta1, U1, V1, H1 = explicit_step(
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
        f=f,
        enable_coriolis=True,
    )

    # Expect V to have grown slightly due to coriolis from U
    # and U to have changed slightly (but not exploded)
    assert np.isfinite(U1).all() and np.isfinite(V1).all()
    assert np.abs(V1[2:-2, 2:-2]).mean() > 0.0
    # Volume still conserved (eta initially zero and no flux divergence with constant U at faces)
    vol0 = float(H.sum())
    vol1 = float(H1.sum())
    assert np.isclose(vol0, vol1, rtol=1e-12, atol=1e-12)
