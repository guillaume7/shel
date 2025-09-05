import numpy as np
import pytest

from shel.model.boundary_conditions import get_bc


def test_flather_updates_u_on_east_with_external_eta():
    ny, nx = 4, 5
    H = np.full((ny, nx), 10.0)
    g = 9.81
    U_old = np.zeros((ny, nx + 1))
    V_old = np.zeros((ny + 1, nx))
    U = U_old.copy()
    V = V_old.copy()
    eta_old = np.zeros((ny, nx))
    # impose external elevated water level at east boundary
    eta_ext = np.zeros_like(eta_old)
    eta_ext[:, -1] = 0.1  # 10 cm higher outside

    m_cls, _ = get_bc("flather")
    assert m_cls is not None
    bc = m_cls()
    bc.apply_side(
        U,
        V,
        "east",
        U_old=U_old,
        V_old=V_old,
        H=H,
        g=g,
        dt=1.0,
        dx=1.0,
        dy=1.0,
        eta_old=eta_old,
        eta_ext=eta_ext,
    )
    # expected correction ~ (c/H_mean) * (eta_ext-eta_old)
    c = float(np.sqrt(g * np.mean(H)))
    H_mean = float(np.maximum(1e-12, np.mean(H)))
    expected = (c / H_mean) * (eta_ext[:, -1] - eta_old[:, -1])
    np.testing.assert_allclose(U[:, -1], expected, rtol=1e-12, atol=1e-12)
    # Other edges unchanged
    assert np.all(U[:, 0] == 0.0)
    assert np.all(V[0, :] == 0.0)
    assert np.all(V[-1, :] == 0.0)
