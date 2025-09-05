import numpy as np
import pytest

from shel.model.solvers.common import (
    d_dx_t_to_u,
    d_dy_t_to_v,
    avg_x_t_to_u,
    avg_y_t_to_v,
    avg_x_u_to_t,
    avg_y_v_to_t,
    div_uv_to_t,
)


def make_linear_T(nx=8, ny=6, ax=2.0, ay=3.0, c=1.0):
    x = np.arange(nx)
    y = np.arange(ny)
    X, Y = np.meshgrid(x, y)
    return ax * X + ay * Y + c


def test_t_to_u_v_derivatives_linear_field():
    nx, ny = 10, 7
    dx, dy = 2.0, 3.0
    T = make_linear_T(nx, ny, ax=4.0, ay=-2.0, c=0.5)

    dTx = d_dx_t_to_u(T, dx)
    dTy = d_dy_t_to_v(T, dy)

    # Check shapes
    assert dTx.shape == (ny, nx + 1)
    assert dTy.shape == (ny + 1, nx)

    # Interior should be constant and equal to the analytic derivatives
    assert np.allclose(dTx[:, 1:nx], 4.0 / dx)
    assert np.allclose(dTy[1:ny, :], -2.0 / dy)

    # Boundaries should be NaN
    assert np.all(np.isnan(dTx[:, 0]))
    assert np.all(np.isnan(dTx[:, -1]))
    assert np.all(np.isnan(dTy[0, :]))
    assert np.all(np.isnan(dTy[-1, :]))


def test_t_u_v_averaging_round_trip_centering():
    nx, ny = 12, 5
    T = make_linear_T(nx, ny, ax=1.0, ay=2.0, c=0.0)

    Tu = avg_x_t_to_u(T)
    Tv = avg_y_t_to_v(T)

    # Shapes
    assert Tu.shape == (ny, nx + 1)
    assert Tv.shape == (ny + 1, nx)

    # Back to T centers
    T_from_u = avg_x_u_to_t(Tu)
    T_from_v = avg_y_v_to_t(Tv)

    # For a linear field, averaging to faces and back should reproduce T exactly
    assert np.allclose(T_from_u, T)
    assert np.allclose(T_from_v, T)


def test_divergence_zero_for_solid_body():
    nx, ny = 16, 16
    dx, dy = 1.0, 1.0
    # Solid body rotation u = -omega*y at U faces, v = omega*x at V faces
    omega = 0.1
    # Build U (ny, nx+1) with u at face centers aligned with y index
    y = np.arange(ny)
    U = np.zeros((ny, nx + 1))
    U[:] = -omega * y[:, None]
    # Build V (ny+1, nx) with v proportional to x index
    x = np.arange(nx)
    V = np.zeros((ny + 1, nx))
    V[:] = omega * x[None, :]

    div = div_uv_to_t(U, V, dx, dy)
    assert div.shape == (ny, nx)
    assert np.allclose(div, 0.0)
