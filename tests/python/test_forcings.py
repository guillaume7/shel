import numpy as np

from shel.model.forcings.bottom import linear_drag, quadratic_drag
from shel.model.forcings.surface import surface_pressure, wind_stress


def test_wind_stress_scalar():
    tau = wind_stress(10.0, 2.0)
    assert np.isclose(tau, 1.225 * 1.3e-3 * 8.0 * 8.0)


def test_wind_stress_vector():
    U_air = np.array([10.0, 0.0])
    U_surface = np.array([2.0, 0.0])
    tau = wind_stress(U_air, U_surface)
    assert np.allclose(tau, [1.225 * 1.3e-3 * 8.0 * 8.0, 0.0])


def test_surface_pressure():
    assert np.isclose(surface_pressure(101325.0), 0.0)
    assert np.isclose(surface_pressure(101425.0), 100.0)


def test_linear_drag():
    U = np.array([1.0, -2.0])
    r = 0.01
    drag = linear_drag(U, r)
    assert np.allclose(drag, [-0.01, 0.02])


def test_quadratic_drag():
    U = np.array([2.0, -3.0])
    Cd = 0.0025
    H = 10.0
    drag = quadratic_drag(U, Cd, H)
    expected = -Cd * np.abs(U) * U / H
    assert np.allclose(drag, expected)
