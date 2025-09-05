import numpy as np

from shel.model.solvers.common.ministep import explicit_step


def test_freeslip_bc_zero_normal_and_zero_tangential_gradient():
    ny, nx = 10, 11
    dx = dy = 1.0
    dt = 0.02

    H = np.ones((ny, nx)) * 8.0
    eta = np.zeros((ny, nx))

    # Seed interior velocities with a smooth field
    y = np.arange(ny)
    x = np.arange(nx)
    uy = 0.01 * np.sin(2 * np.pi * y / max(1, ny - 1))[:, None]
    vx = 0.01 * np.cos(2 * np.pi * x / max(1, nx - 1))[None, :]
    U = uy * np.ones((ny, nx + 1))
    V = vx * np.ones((ny + 1, nx))

    U[:, 0] = 0.05  # non-zero to check enforcement
    U[:, -1] = -0.03
    V[0, :] = 0.07
    V[-1, :] = -0.04

    eta1, U1, V1 = explicit_step(
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
        bc_type="freeslip",
    )

    # Normal components zero at boundaries
    assert np.allclose(U1[:, 0], 0.0)
    assert np.allclose(U1[:, -1], 0.0)
    assert np.allclose(V1[0, :], 0.0)
    assert np.allclose(V1[-1, :], 0.0)

    # Tangential components have zero normal gradient (copied from interior)
    # West/East walls: tangential is V; check equality with adjacent interior column
    if nx > 2 and ny > 2:
        assert np.allclose(V1[1:-1, 0], V1[1:-1, 1])
        assert np.allclose(V1[1:-1, -1], V1[1:-1, -2])
        # South/North walls: tangential is U; check equality with adjacent interior row
        assert np.allclose(U1[0, 1:-1], U1[1, 1:-1])
        assert np.allclose(U1[-1, 1:-1], U1[-2, 1:-1])

    # Sanity: no NaNs
    for arr in (eta1, U1, V1):
        assert np.isfinite(arr).all()
