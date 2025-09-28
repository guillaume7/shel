import numpy as np

from shel.model.forcings.bottom.quadratic_drag import quadratic_drag


def test_quadratic_drag_zero_velocity():
    U = np.zeros(5)
    Cd = 0.0025
    H = 10.0
    drag = quadratic_drag(U, Cd, H)
    assert np.allclose(drag, 0.0)


def test_quadratic_drag_sign_and_scaling():
    U = np.array([1.0, -2.0, 3.0])
    Cd = 0.0025
    H = 10.0
    drag = quadratic_drag(U, Cd, H)
    # Drag should oppose velocity
    assert np.all(np.sign(drag) == -np.sign(U))
    # Magnitude scaling
    expected = -Cd * np.abs(U) * U / H
    assert np.allclose(drag, expected)
