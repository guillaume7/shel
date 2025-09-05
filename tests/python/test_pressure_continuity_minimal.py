import numpy as np

from shel.model.solvers.momentum.pressure import pressure_gradient
from shel.model.solvers.waterlevel.continuity import update_free_surface


def test_pressure_gradient_linear_eta():
    ny, nx = 8, 10
    dx, dy = 2.0, 3.0
    g = 9.81
    # eta = ax*x + ay*y
    ax, ay = 0.1, -0.05
    x = np.arange(nx)
    y = np.arange(ny)
    X, Y = np.meshgrid(x, y)
    eta = ax * X + ay * Y

    PG_u, PG_v = pressure_gradient(eta, g, dx, dy)

    # Interior faces should be constant equal to -g * gradient components / spacing
    assert np.allclose(PG_u[:, 1:nx], -(g * ax) / dx)
    assert np.allclose(PG_v[1:ny, :], -(g * ay) / dy)

    # Boundaries NaN due to one-sided lack in minimal operator
    assert np.all(np.isnan(PG_u[:, 0])) and np.all(np.isnan(PG_u[:, -1]))
    assert np.all(np.isnan(PG_v[0, :])) and np.all(np.isnan(PG_v[-1, :]))


ess = 1e-12

def test_continuity_conserves_volume_closed_box():
    ny, nx = 12, 14
    dx, dy = 1.0, 1.0
    dt = 0.1

    # Flat H and zero-normal-flow boundaries realized by zero velocity at domain outer faces
    H = np.ones((ny, nx)) * 10.0
    eta = np.zeros((ny, nx))

    # Build a divergence-free interior velocity field (solid body) with zeros at physical boundaries
    omega = 0.2
    x = np.arange(nx)
    y = np.arange(ny)
    # U (ny, nx+1)
    U = np.zeros((ny, nx + 1))
    U[:, 1:-1] = -omega * y[:, None]
    # V (ny+1, nx)
    V = np.zeros((ny + 1, nx))
    V[1:-1, :] = omega * x[None, :]

    def vol(field):
        return float(field.sum() * dx * dy)

    vol_before = vol(eta + H)
    eta_next = update_free_surface(eta, H, U, V, dt, dx, dy)
    vol_after = vol(eta_next + H)

    # Volume should be unchanged to numerical precision
    assert np.isclose(vol_before, vol_after, rtol=1e-12, atol=1e-12)


def test_continuity_eta_response_simple_flux():
    # Simple 1D-like test: nonzero U flux in a row should change eta accordingly
    ny, nx = 4, 6
    dx, dy = 1.0, 1.0
    dt = 0.5
    H = np.ones((ny, nx))
    eta = np.zeros((ny, nx))

    U = np.zeros((ny, nx + 1))
    V = np.zeros((ny + 1, nx))

    # Create a uniform eastward U in the interior faces of row j=2
    j = 2
    U[j, 1:-1] = 1.0

    eta_next = update_free_surface(eta, H, U, V, dt, dx, dy)

    # Divergence pattern: dU/dx positive at columns where U decreases west->east
    # Here U is 0,1,1,1,0 -> differences at T centers [1,0,0,-1] along the row
    expected = np.zeros_like(eta)
    expected[j, 0] = -1.0  # inflow from west
    expected[j, -1] = 1.0  # outflow to east

    # The interior should be unchanged; edge T cells respond with opposite signs
    # eta_next = eta - dt * div(HU, HV) with H=1, V=0
    assert np.isclose(eta_next[j, 0], -dt)
    assert np.isclose(eta_next[j, -1], dt)
    assert np.allclose(eta_next[j, 1:-1], 0.0)
