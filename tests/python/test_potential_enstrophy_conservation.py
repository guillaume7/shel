import numpy as np

from shel.model.diagnostics.potential_enstrophy import integrated_potential_enstrophy
from shel.model.diagnostics.pv import potential_vorticity
from shel.model.grid import Grid
from shel.model.solvers.common.ministep import explicit_step


def test_potential_enstrophy_conservation_inviscid_closed_box():
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

    # Compute initial PV and potential enstrophy
    grid = Grid({"grid": {"nx": nx, "ny": ny, "dx": dx, "dy": dy}})
    pv0 = potential_vorticity(U, V, H, f, grid)
    penst0 = integrated_potential_enstrophy(pv0, dx, dy)

    for _ in range(10):
        eta, U, V = explicit_step(
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
    pvN = potential_vorticity(U, V, H, f, grid)
    penstN = integrated_potential_enstrophy(pvN, dx, dy)
    # Should be conserved in inviscid, closed box
    assert np.isclose(penstN, penst0, rtol=1e-10, atol=1e-10)
