import numpy as np

from shel.model.solvers.momentum.advection import advect_momentum


def test_advection_uniform_velocity_zero():
    # Constant velocity field should yield zero advection
    ny, nx = 16, 18
    dx = dy = 1.0
    U = np.ones((ny, nx + 1)) * 0.3
    V = np.ones((ny + 1, nx)) * -0.2
    Au, Av = advect_momentum(U, V, dx, dy)
    assert np.allclose(Au, 0.0)
    assert np.allclose(Av, 0.0)


def test_advection_pure_shear_u_of_y_zero():
    # u = s*y, v = 0 -> ∂u/∂x = 0 and v=0 so Au=0; Av=0
    ny, nx = 20, 22
    dx = dy = 1.0
    s = 0.05
    y = np.arange(ny)
    U = s * y[:, None] * np.ones((ny, nx + 1))
    V = np.zeros((ny + 1, nx))
    Au, Av = advect_momentum(U, V, dx, dy)
    assert np.allclose(Au[2:-2, 2:-2], 0.0)
    assert np.allclose(Av[2:-2, 2:-2], 0.0)


def test_advection_pure_shear_v_of_x_zero():
    # u = 0, v = s*x -> ∂v/∂y = 0 and u=0 so Av=0; Au=0
    ny, nx = 18, 19
    dx = dy = 1.0
    s = -0.07
    x = np.arange(nx)
    V = s * x[None, :] * np.ones((ny + 1, nx))
    U = np.zeros((ny, nx + 1))
    Au, Av = advect_momentum(U, V, dx, dy)
    assert np.allclose(Au[2:-2, 2:-2], 0.0)
    assert np.allclose(Av[2:-2, 2:-2], 0.0)
