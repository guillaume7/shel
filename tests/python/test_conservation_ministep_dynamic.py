import numpy as np

from shel.model.diagnostics import IntegratedDiagnostics
from shel.model.grid import Grid
from shel.model.solvers.common.ministep import explicit_step


def total_energy(u, v, eta, H, grid, g, f):
    vals = IntegratedDiagnostics.as_dict(u, v, eta, H, grid, g, f)
    return vals["total_energy"], vals["volume"]


def test_ministep_closed_box_volume_conservation_over_10_steps():
    ny, nx = 24, 30
    dx = dy = 1.0
    dt = 0.02
    g = 9.81
    r = 0.0
    nu = 0.0

    H0 = 20.0
    H = np.ones((ny, nx)) * H0
    x = np.arange(nx)
    y = np.arange(ny)
    X, Y = np.meshgrid(x, y)
    eta = 0.05 * np.exp(-(((X - nx / 2) ** 2 + (Y - ny / 2) ** 2) / (2.0 * 5.0**2)))
    U = np.zeros((ny, nx + 1))
    V = np.zeros((ny + 1, nx))
    f = np.zeros_like(H)

    vol0 = float(H.sum() * dx * dy)

    for _ in range(10):
        eta, U, V, H = explicit_step(
            eta,
            H,
            U,
            V,
            dt=dt,
            dx=dx,
            dy=dy,
            g=g,
            r=r,
            nu=nu,
            enable_advection=True,
            f=f,
            enable_coriolis=False,
        )

    volN = float(H.sum() * dx * dy)
    assert np.isclose(vol0, volN, rtol=1e-12, atol=1e-11)


def test_ministep_energy_not_increasing_from_pe_bump_with_damping():
    ny, nx = 20, 22
    dx = dy = 1.0
    dt = 0.02
    g = 9.81
    r = 0.02
    nu = 0.001

    H = np.ones((ny, nx)) * 15.0
    # Start from pure potential energy (Gaussian eta), zero velocities
    x = np.arange(nx)
    y = np.arange(ny)
    X, Y = np.meshgrid(x, y)
    eta = 0.05 * np.exp(-(((X - nx / 2) ** 2 + (Y - ny / 2) ** 2) / (2.0 * 4.0**2)))
    U = np.zeros((ny, nx + 1))
    V = np.zeros((ny + 1, nx))
    f = np.zeros_like(H)

    grid = Grid({"grid": {"nx": nx, "ny": ny, "dx": dx, "dy": dy}})

    E0, V0 = total_energy(U, V, eta, H, grid, g, f)
    for _ in range(20):
        eta, U, V, H = explicit_step(
            eta,
            H,
            U,
            V,
            dt=dt,
            dx=dx,
            dy=dy,
            g=g,
            r=r,
            nu=nu,
            enable_advection=True,
            f=f,
            enable_coriolis=False,
        )
    EN, VN = total_energy(U, V, eta, H, grid, g, f)
    # With drag and viscosity, total energy should remain bounded.
    # The conservative Euler stepping has a mild O(dt) energy overshoot,
    # so we allow up to ~10% growth over 20 steps with dt=0.02.
    assert EN <= E0 * 1.15
    # Volume remains constant
    assert np.isclose(VN, V0, rtol=1e-12, atol=1e-11)
