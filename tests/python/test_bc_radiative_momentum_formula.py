import numpy as np

from shel.model.boundary_conditions.registry import get_bc


def test_sommerfeld_momentum_formula_east_west():
    ny, nx = 6, 8
    dx = dy = 1.0
    dt = 0.1
    g = 9.0
    H = np.ones((ny, nx)) * 4.0  # c = sqrt(g*H) = sqrt(36)=6

    U = np.zeros((ny, nx + 1))
    V = np.zeros((ny + 1, nx))

    # Linear profile in x so (U_old[:,1]-U_old[:,0]) = const
    U_old = np.tile(np.linspace(0.0, 1.0, nx + 1), (ny, 1))
    V_old = V.copy()

    m_cls, _ = get_bc("radiative")
    assert m_cls is not None
    bc = m_cls()

    # WEST
    bc.apply_side(U, V, "west", U_old=U_old, V_old=V_old, H=H, g=g, dt=dt, dx=dx, dy=dy)
    c = np.sqrt(g * H[:, 0].mean())
    r = c * dt / dx
    expected_west = U_old[:, 0] - r * (U_old[:, 1] - U_old[:, 0])
    np.testing.assert_allclose(U[:, 0], expected_west, rtol=0, atol=1e-12)

    # EAST
    U.fill(0.0)
    bc.apply_side(U, V, "east", U_old=U_old, V_old=V_old, H=H, g=g, dt=dt, dx=dx, dy=dy)
    expected_east = U_old[:, -1] - r * (U_old[:, -1] - U_old[:, -2])
    np.testing.assert_allclose(U[:, -1], expected_east, rtol=0, atol=1e-12)


def test_sommerfeld_momentum_formula_south_north():
    ny, nx = 6, 8
    dx = dy = 1.0
    dt = 0.1
    g = 9.0
    H = np.ones((ny, nx)) * 4.0  # c = 6

    U = np.zeros((ny, nx + 1))
    V = np.zeros((ny + 1, nx))

    U_old = U.copy()
    # Linear profile in y so (V_old[1,:]-V_old[0,:]) = const
    vline = np.linspace(0.0, 1.0, ny + 1)[:, None]
    V_old = np.tile(vline, (1, nx))

    m_cls, _ = get_bc("radiative")
    assert m_cls is not None
    bc = m_cls()

    # SOUTH
    bc.apply_side(U, V, "south", U_old=U_old, V_old=V_old, H=H, g=g, dt=dt, dx=dx, dy=dy)
    c = np.sqrt(g * H[0, :].mean())
    r = c * dt / dy
    expected_south = V_old[0, :] - r * (V_old[1, :] - V_old[0, :])
    np.testing.assert_allclose(V[0, :], expected_south, rtol=0, atol=1e-12)

    # NORTH
    V.fill(0.0)
    bc.apply_side(U, V, "north", U_old=U_old, V_old=V_old, H=H, g=g, dt=dt, dx=dx, dy=dy)
    expected_north = V_old[-1, :] - r * (V_old[-1, :] - V_old[-2, :])
    np.testing.assert_allclose(V[-1, :], expected_north, rtol=0, atol=1e-12)
