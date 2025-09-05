import numpy as np

from shel.model.solvers.common.ministep import explicit_step


def domain_volume(field, dx, dy):
    return float(field.sum() * dx * dy)


def test_ministep_closed_box_volume_conservation_over_few_steps():
    ny, nx = 20, 24
    dx, dy = 1.0, 1.0
    dt = 0.05
    g = 9.81
    r = 0.01
    nu = 0.001

    # Flat bottom depth and small gaussian pulse in eta
    H0 = 10.0
    H = np.ones((ny, nx)) * H0
    x = np.arange(nx)
    y = np.arange(ny)
    X, Y = np.meshgrid(x, y)
    cx, cy = nx / 2.0, ny / 2.0
    eta = 0.1 * np.exp(-(((X - cx) ** 2 + (Y - cy) ** 2) / (2.0 * 4.0 ** 2)))

    U = np.zeros((ny, nx + 1))
    V = np.zeros((ny + 1, nx))

    vol0 = domain_volume(eta + H, dx, dy)

    steps = 5
    for _ in range(steps):
        eta, U, V = explicit_step(eta, H, U, V, dt=dt, dx=dx, dy=dy, g=g, r=r, nu=nu)

    volN = domain_volume(eta + H, dx, dy)

    assert np.isclose(vol0, volN, rtol=1e-12, atol=1e-10)


def test_ministep_stability_small_dt():
    # Ensure no NaNs/Infs arise for small dt and reasonable coeffs
    ny, nx = 10, 12
    dx, dy = 1.0, 1.0
    dt = 0.02

    H = np.ones((ny, nx)) * 5.0
    eta = np.zeros((ny, nx))

    rng = np.random.default_rng(0)
    U = 0.01 * rng.standard_normal((ny, nx + 1))
    V = 0.01 * rng.standard_normal((ny + 1, nx))

    eta1, U1, V1 = explicit_step(eta, H, U, V, dt=dt, dx=dx, dy=dy, g=9.81, r=0.0, nu=0.0)

    for arr in (eta1, U1, V1):
        assert np.isfinite(arr).all()
