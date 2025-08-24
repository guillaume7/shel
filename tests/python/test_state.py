"""
Tests for the SHEL model state.
"""

import pytest
import numpy as np

from shel.model.state import ModelState


def test_model_state_initialization():
    """Test that model state is correctly initialized."""
    config = {
        "grid": {
            "nx": 50,
            "ny": 40,
            "dx": 1000.0,
            "dy": 1000.0,
            "x_origin": 0.0,
            "y_origin": 0.0,
        },
        "model": {
            "timestep": 60.0,
            "num_steps": 100,
            "gravity": 9.81,
            "viscosity": 10.0,
            "bottom_drag_coef": 0.0025,
            "coriolis_parameter": 1e-4,
        },
    }

    state = ModelState(config)

    # Check grid and time properties
    assert state.grid.nx == 50
    assert state.grid.ny == 40
    assert state.time == 0.0
    assert state.timestep == 60.0
    assert state.step == 0

    # Check physical parameters
    assert state.gravity == 9.81
    assert state.viscosity == 10.0
    assert state.bottom_drag_coef == 0.0025

    # Check that fields are initialized to zeros with the correct shape
    assert state.eta.shape == (40, 50)
    assert state.u.shape == (40, 51)
    assert state.v.shape == (41, 50)
    assert state.d.shape == (40, 50)
    assert state.H.shape == (40, 50)
    assert np.allclose(state.eta, 0.0)
    assert np.allclose(state.u, 0.0)
    assert np.allclose(state.v, 0.0)
    assert np.allclose(state.d, 0.0)


def test_set_bathymetry():
    """Test setting bathymetry."""
    config = {
        "grid": {"nx": 10, "ny": 8, "dx": 1000.0, "dy": 1000.0},
        "model": {"timestep": 60.0},
    }

    state = ModelState(config)

    # Create a bathymetry field
    bathymetry = np.ones((8, 10)) * 1000.0

    # Set bathymetry
    state.set_bathymetry(bathymetry)

    # Check that bathymetry was set correctly
    assert np.allclose(state.d, 1000.0)

    # Check that total depth was updated (H = eta + d)
    assert np.allclose(state.H, 1000.0)  # eta is 0, so H = d

    # Test with invalid shape
    with pytest.raises(ValueError):
        state.set_bathymetry(np.ones((5, 5)))


def test_set_initial_elevation():
    """Test setting initial elevation."""
    config = {
        "grid": {"nx": 10, "ny": 8, "dx": 1000.0, "dy": 1000.0},
        "model": {"timestep": 60.0},
    }

    state = ModelState(config)

    # Set bathymetry first (required for H calculation)
    bathymetry = np.ones((8, 10)) * 1000.0
    state.set_bathymetry(bathymetry)

    # Create an elevation field with a Gaussian bump
    x = np.linspace(-5, 5, 10)
    y = np.linspace(-4, 4, 8)
    X, Y = np.meshgrid(x, y)
    elevation = np.exp(-(X**2 + Y**2) / 2)

    # Set initial elevation
    state.set_initial_elevation(elevation)

    # Check that elevation was set correctly
    assert np.allclose(state.eta, elevation)
    assert np.allclose(state.eta_old, elevation)
    assert np.allclose(state.eta_new, elevation)

    # Check that total depth was updated (H = eta + d)
    assert np.allclose(state.H, 1000.0 + elevation)

    # Test with invalid shape
    with pytest.raises(ValueError):
        state.set_initial_elevation(np.ones((5, 5)))


def test_energy_calculation():
    """Test energy calculation."""
    config = {
        "grid": {"nx": 10, "ny": 8, "dx": 1000.0, "dy": 1000.0},
        "model": {"timestep": 60.0, "gravity": 9.81},
    }

    state = ModelState(config)

    # Set up a simple test case
    bathymetry = np.ones((8, 10)) * 1000.0
    state.set_bathymetry(bathymetry)

    # Set up elevation
    elevation = np.zeros((8, 10))
    elevation[4, 5] = 1.0  # 1m bump in the middle
    state.set_initial_elevation(elevation)

    # Set up velocities
    u = np.zeros((8, 11))
    u[:, 5:] = 1.0  # 1 m/s flow to the right in half the domain
    v = np.zeros((9, 10))

    state.set_initial_velocities(u, v)

    # Calculate energies
    ke = state.compute_kinetic_energy()
    pe = state.compute_potential_energy()
    te = state.compute_total_energy()

    # Check that energies are positive
    assert ke > 0
    assert pe > 0
    assert np.isclose(te, ke + pe)

    # Check volume conservation
    initial_volume = state.compute_volume()
    assert initial_volume > 0


def test_time_stepping():
    """Test time stepping and level swapping."""
    config = {
        "grid": {"nx": 10, "ny": 8, "dx": 1000.0, "dy": 1000.0},
        "model": {"timestep": 60.0},
    }

    state = ModelState(config)

    # Set initial state
    bathymetry = np.ones((8, 10)) * 1000.0
    state.set_bathymetry(bathymetry)

    elevation = np.zeros((8, 10))
    elevation[4, 5] = 1.0
    state.set_initial_elevation(elevation)

    # Initial time and step
    assert state.time == 0.0
    assert state.step == 0

    # Advance time
    state.advance_time()
    assert state.time == 60.0
    assert state.step == 1

    # Modify state at new time level
    state.eta_new[4, 5] = 0.9  # Decreased amplitude

    # Swap time levels
    state.swap_time_levels()

    # Check that new became current, current became old
    assert state.eta[4, 5] == 0.9
    assert state.eta_old[4, 5] == 1.0
