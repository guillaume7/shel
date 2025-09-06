"""
Pytest configuration file.
"""

import os
import sys

import numpy as np
import pytest

# Add the src directory to the path so we can import the package
sys.path.insert(
    0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src/python"))
)


@pytest.fixture
def sample_grid():
    """
    Return a sample grid configuration for testing.
    """
    return {
        "nx": 100,
        "ny": 100,
        "dx": 1000.0,
        "dy": 1000.0,
        "x_origin": 0.0,
        "y_origin": 0.0,
    }


@pytest.fixture
def sample_bathymetry(sample_grid):  # pylint: disable=redefined-outer-name
    """
    Return a sample bathymetry field for testing.
    """
    nx, ny = sample_grid["nx"], sample_grid["ny"]
    depth = 1000.0 * np.ones((ny, nx))

    # Add a bump in the middle
    x_center = nx // 2
    y_center = ny // 2
    radius = min(nx, ny) // 10

    for i in range(ny):
        for j in range(nx):
            dist = np.sqrt((i - y_center) ** 2 + (j - x_center) ** 2)
            if dist < radius:
                depth[i, j] = 1000.0 - 500.0 * (1.0 - dist / radius)

    return depth


@pytest.fixture
def sample_config():
    """
    Return a sample configuration for testing.
    """
    return {
        "model": {
            "timestep": 60.0,  # seconds
            "num_steps": 1000,
            "output_interval": 10,
            "coriolis_parameter": 1e-4,  # f-plane approximation
            "gravity": 9.81,  # m/s^2
            "viscosity": 10.0,  # m^2/s
            "bottom_drag_coef": 0.0025,  # dimensionless
        },
        "boundary_conditions": {
            "north": "closed",
            "south": "closed",
            "east": "closed",
            "west": "closed",
        },
        "initial_conditions": {
            "type": "gaussian_bump",
            "amplitude": 1.0,  # meters
            "sigma": 10000.0,  # meters
            "x_center": 50000.0,  # meters
            "y_center": 50000.0,  # meters
        },
    }
