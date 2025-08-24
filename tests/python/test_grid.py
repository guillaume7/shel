"""
Tests for the SHEL model grid implementation.
"""
import pytest
import numpy as np

from shel.model.grid import Grid


def test_grid_initialization():
    """Test that grid is correctly initialized."""
    config = {
        'grid': {
            'nx': 100,
            'ny': 80,
            'dx': 1000.0,
            'dy': 1000.0,
            'x_origin': 0.0,
            'y_origin': 0.0
        }
    }

    grid = Grid(config)

    assert grid.nx == 100
    assert grid.ny == 80
    assert grid.dx == 1000.0
    assert grid.dy == 1000.0
    assert grid.x_origin == 0.0
    assert grid.y_origin == 0.0

    # Check that coordinate arrays have correct shape
    assert grid.x_t.shape == (80, 100)
    assert grid.y_t.shape == (80, 100)
    assert grid.x_u.shape == (80, 101)
    assert grid.y_u.shape == (80, 101)
    assert grid.x_v.shape == (81, 100)
    assert grid.y_v.shape == (81, 100)
    assert grid.x_q.shape == (81, 101)
    assert grid.y_q.shape == (81, 101)


def test_grid_coordinates():
    """Test that grid coordinates are correctly calculated."""
    config = {
        'grid': {
            'nx': 10,
            'ny': 8,
            'dx': 100.0,
            'dy': 200.0,
            'x_origin': 1000.0,
            'y_origin': 2000.0
        }
    }

    grid = Grid(config)

    # Check T-grid (cell centers)
    assert np.isclose(grid.x_t[0, 0], 1000.0 + 100.0/2)  # x_origin + dx/2
    assert np.isclose(grid.y_t[0, 0], 2000.0 + 200.0/2)  # y_origin + dy/2
    assert np.isclose(grid.x_t[7, 9], 1000.0 + 9*100.0 + 100.0/2)  # Last cell center
    assert np.isclose(grid.y_t[7, 9], 2000.0 + 7*200.0 + 200.0/2)  # Last cell center

    # Check U-grid (east/west cell faces)
    assert np.isclose(grid.x_u[0, 0], 1000.0)  # x_origin
    assert np.isclose(grid.y_u[0, 0], 2000.0 + 200.0/2)  # y_origin + dy/2
    assert np.isclose(grid.x_u[7, 10], 1000.0 + 10*100.0)  # Last U-point
    assert np.isclose(grid.y_u[7, 10], 2000.0 + 7*200.0 + 200.0/2)  # Last U-point

    # Check V-grid (north/south cell faces)
    assert np.isclose(grid.x_v[0, 0], 1000.0 + 100.0/2)  # x_origin + dx/2
    assert np.isclose(grid.y_v[0, 0], 2000.0)  # y_origin
    assert np.isclose(grid.x_v[8, 9], 1000.0 + 9*100.0 + 100.0/2)  # Last V-point
    assert np.isclose(grid.y_v[8, 9], 2000.0 + 8*200.0)  # Last V-point


def test_grid_coriolis():
    """Test calculation of Coriolis parameter."""
    config = {
        'grid': {
            'nx': 10,
            'ny': 10,
            'dx': 10000.0,
            'dy': 10000.0,
            'x_origin': 0.0,
            'y_origin': 0.0
        }
    }

    grid = Grid(config)

    # Test f-plane (constant f)
    f0 = 1e-4
    f = grid.compute_coriolis(f0)
    assert f.shape == (10, 10)
    assert np.allclose(f, f0)

    # Test beta-plane (f varies with y)
    f0 = 1e-4
    beta = 1e-11
    f = grid.compute_coriolis(f0, beta)
    assert f.shape == (10, 10)

    # Check that f increases with y
    assert np.all(np.diff(f, axis=0) > 0)

    # Check center value
    center_y = grid.y_t[5, 5]
    y_ref = 0.0 + 10 * 10000.0 / 2  # Middle of domain
    assert np.isclose(f[5, 5], f0 + beta * (center_y - y_ref))


def test_land_mask():
    """Test setting land mask."""
    config = {
        'grid': {
            'nx': 10,
            'ny': 8,
            'dx': 100.0,
            'dy': 100.0
        }
    }

    grid = Grid(config)

    # Create a land mask with land in the corners
    mask = np.zeros((8, 10))
    mask[0, 0] = 1
    mask[0, 9] = 1
    mask[7, 0] = 1
    mask[7, 9] = 1

    grid.set_land_mask(mask)

    assert grid.mask.shape == (8, 10)
    assert grid.mask[0, 0] == 1
    assert grid.mask[0, 9] == 1
    assert grid.mask[7, 0] == 1
    assert grid.mask[7, 9] == 1
    assert grid.mask[4, 5] == 0

    # Test with invalid shape
    with pytest.raises(ValueError):
        grid.set_land_mask(np.zeros((5, 5)))
