import numpy as np

from shel.model.solvers.common.ministep import explicit_step


def test_closed_bc_zeroes_normal_faces_in_ministep():
    ny, nx = 12, 14
    dx = dy = 1.0
    dt = 0.03

    H = np.ones((ny, nx)) * 10.0
    x = np.arange(nx)
    y = np.arange(ny)
    X, Y = np.meshgrid(x, y)
    eta = 0.05 * np.exp(-(((X - nx / 2) ** 2 + (Y - ny / 2) ** 2) / (2.0 * 3.0**2)))

    rng = np.random.default_rng(42)
    U = 0.01 * rng.standard_normal((ny, nx + 1))
    V = 0.01 * rng.standard_normal((ny + 1, nx))

    # Deliberately set non-zero boundary values to check they get zeroed
    U[:, 0] = 0.1
    U[:, -1] = -0.1
    V[0, :] = 0.2
    V[-1, :] = -0.2

    eta1, U1, V1, _ = explicit_step(
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
        enable_advection=True,
        f=np.zeros_like(H),
        enable_coriolis=False,
        bc_type="closed",
    )

    assert np.allclose(U1[:, 0], 0.0)
    assert np.allclose(U1[:, -1], 0.0)
    assert np.allclose(V1[0, :], 0.0)
    assert np.allclose(V1[-1, :], 0.0)

    # Sanity check no NaNs
    for arr in (eta1, U1, V1):
        assert np.isfinite(arr).all()
