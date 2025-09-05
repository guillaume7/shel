import numpy as np

from shel.model.solvers.momentum.friction import bottom_drag_tendency
from shel.model.solvers.momentum.diffusion import viscous_tendency


def test_bottom_drag_simple():
    ny, nx = 6, 8
    U = np.ones((ny, nx + 1)) * 3.0
    V = np.ones((ny + 1, nx)) * -2.0
    r = 0.5
    dU, dV = bottom_drag_tendency(U, V, r)
    assert dU.shape == U.shape and dV.shape == V.shape
    assert np.allclose(dU, -r * U)
    assert np.allclose(dV, -r * V)


def test_viscous_tendency_linear_field_zero():
    # Laplacian of linear fields is zero
    ny, nx = 10, 12
    dx, dy = 2.0, 3.0
    nu = 1.0
    x = np.arange(nx + 1)
    y = np.arange(ny)
    X, Y = np.meshgrid(x, y)
    U = 2.0 * X + 1.0 * Y
    x2 = np.arange(nx)
    y2 = np.arange(ny + 1)
    X2, Y2 = np.meshgrid(x2, y2)
    V = -0.5 * X2 + 0.25 * Y2

    dU, dV = viscous_tendency(U, V, nu, dx, dy)

    # Interior should be exactly zero; boundaries left at zero by implementation
    assert np.allclose(dU[1:-1, 1:-1], 0.0)
    assert np.allclose(dV[1:-1, 1:-1], 0.0)


def test_viscous_tendency_quadratic_field_constant():
    # Laplacian of ax^2 + by^2 is 2a + 2b; check constant interior values
    ny, nx = 14, 15
    dx, dy = 1.0, 1.5
    nu = 0.3

    x = np.arange(nx + 1)
    y = np.arange(ny)
    X, Y = np.meshgrid(x, y)
    a, b = 0.2, -0.1
    U = a * X**2 + b * Y**2

    x2 = np.arange(nx)
    y2 = np.arange(ny + 1)
    X2, Y2 = np.meshgrid(x2, y2)
    V = a * X2**2 + b * Y2**2

    dU, dV = viscous_tendency(U, V, nu, dx, dy)

    expected = nu * (2 * a / (dx * dx) + 2 * b / (dy * dy))
    assert np.allclose(dU[2:-2, 2:-2], expected)
    assert np.allclose(dV[2:-2, 2:-2], expected)
