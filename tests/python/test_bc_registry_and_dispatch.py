import numpy as np

from shel.model.boundary_conditions import get_bc, list_bcs


def test_bc_registry_lists_known_bcs():
    bcs = list_bcs()
    assert "closed" in bcs["momentum"]
    assert "freeslip" in bcs["momentum"]
    assert "radiative" in bcs["momentum"] and "radiative" in bcs["eta"]


def test_apply_helpers_bridge_to_strategies():
    ny, nx = 6, 7
    U = np.ones((ny, nx + 1))
    V = -np.ones((ny + 1, nx))

    # closed uniform via strategy
    m_cls, _ = get_bc("closed")
    assert m_cls is not None
    m_cls().apply_uniform(U, V)
    assert np.allclose(U[:, 0], 0.0) and np.allclose(U[:, -1], 0.0)
    assert np.allclose(V[0, :], 0.0) and np.allclose(V[-1, :], 0.0)

    # freeslip side
    U2 = np.zeros_like(U)
    V2 = np.zeros_like(V)
    V2[1:-1, 1] = 2.0
    m_cls, _ = get_bc("freeslip")
    assert m_cls is not None
    m_cls().apply_side(U2, V2, "west")
    assert np.allclose(U2[:, 0], 0.0)
    if ny > 2:
        assert np.allclose(V2[1:-1, 0], V2[1:-1, 1])


def test_radiative_eta_and_momentum_side_do_not_crash():
    ny, nx = 5, 9
    H = np.ones((ny, nx)) * 10.0
    eta = np.zeros((ny, nx))
    U = np.zeros((ny, nx + 1))
    V = np.zeros((ny + 1, nx))
    dt = 0.03
    dx = dy = 1.0
    g = 9.81

    # Sanity: side application runs
    m_cls, e_cls = get_bc("radiative")
    assert m_cls is not None and e_cls is not None
    m_cls().apply_side(
        U, V, "east", U_old=U.copy(), V_old=V.copy(), H=H, g=g, dt=dt, dx=dx, dy=dy
    )
    e_cls().apply_side_eta(
        eta, "east", eta_old=eta.copy(), H=H, g=g, dt=dt, dx=dx, dy=dy
    )
