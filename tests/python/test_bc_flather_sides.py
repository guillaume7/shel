import numpy as np

from shel.model.boundary_conditions import get_bc


def test_flather_west_and_south_sides():
    ny, nx = 6, 7
    H = np.full((ny, nx), 10.0)
    g = 9.81
    U_old = np.zeros((ny, nx + 1))
    V_old = np.zeros((ny + 1, nx))
    U = U_old.copy()
    V = V_old.copy()
    eta_old = np.zeros((ny, nx))

    m_cls, _ = get_bc("flather")
    # Narrow optional type for Pylance (get_bc may return None for unknown names)
    assert m_cls is not None
    bc = m_cls()

    # West boundary
    eta_ext_w = np.zeros_like(eta_old)
    eta_ext_w[:, 0] = 0.05
    bc.apply_side(
        U,
        V,
        "west",
        U_old=U_old,
        V_old=V_old,
        H=H,
        g=g,
        dt=1.0,
        dx=1.0,
        dy=1.0,
        eta_old=eta_old,
        eta_ext=eta_ext_w,
    )
    c = float(np.sqrt(g * np.mean(H[:, 0])))
    Hm = float(np.maximum(1e-12, np.mean(H[:, 0])))
    expected_w = (c / Hm) * (eta_ext_w[:, 0] - eta_old[:, 0])
    np.testing.assert_allclose(U[:, 0], expected_w)

    # South boundary
    U.fill(0.0)
    V.fill(0.0)
    eta_ext_s = np.zeros_like(eta_old)
    eta_ext_s[0, :] = 0.02
    bc.apply_side(
        U,
        V,
        "south",
        U_old=U_old,
        V_old=V_old,
        H=H,
        g=g,
        dt=1.0,
        dx=1.0,
        dy=1.0,
        eta_old=eta_old,
        eta_ext=eta_ext_s,
    )
    c = float(np.sqrt(g * np.mean(H[0, :])))
    Hm = float(np.maximum(1e-12, np.mean(H[0, :])))
    expected_s = (c / Hm) * (eta_ext_s[0, :] - eta_old[0, :])
    np.testing.assert_allclose(V[0, :], expected_s)
