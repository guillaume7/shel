"""
Tests for the SHEL model grid implementation.
"""

import numpy as np
import pytest

from shel.model.grid import (
    Grid,
    apply_noslip_flux_masks,
    build_all_masks,
    build_corner_mask,
    build_staggered_masks,
)
from shel.model.state import ModelState


def test_grid_initialization():
    """Test that grid is correctly initialized."""
    config = {
        "grid": {
            "nx": 100,
            "ny": 80,
            "dx": 1000.0,
            "dy": 1000.0,
            "x_origin": 0.0,
            "y_origin": 0.0,
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
        "grid": {
            "nx": 10,
            "ny": 8,
            "dx": 100.0,
            "dy": 200.0,
            "x_origin": 1000.0,
            "y_origin": 2000.0,
        }
    }

    grid = Grid(config)

    # Check T-grid (cell centers)
    assert np.isclose(grid.x_t[0, 0], 1000.0 + 100.0 / 2)  # x_origin + dx/2
    assert np.isclose(grid.y_t[0, 0], 2000.0 + 200.0 / 2)  # y_origin + dy/2
    assert np.isclose(
        grid.x_t[7, 9], 1000.0 + 9 * 100.0 + 100.0 / 2
    )  # Last cell center
    assert np.isclose(
        grid.y_t[7, 9], 2000.0 + 7 * 200.0 + 200.0 / 2
    )  # Last cell center

    # Check U-grid (east/west cell faces)
    assert np.isclose(grid.x_u[0, 0], 1000.0)  # x_origin
    assert np.isclose(grid.y_u[0, 0], 2000.0 + 200.0 / 2)  # y_origin + dy/2
    assert np.isclose(grid.x_u[7, 10], 1000.0 + 10 * 100.0)  # Last U-point
    assert np.isclose(grid.y_u[7, 10], 2000.0 + 7 * 200.0 + 200.0 / 2)  # Last U-point

    # Check V-grid (north/south cell faces)
    assert np.isclose(grid.x_v[0, 0], 1000.0 + 100.0 / 2)  # x_origin + dx/2
    assert np.isclose(grid.y_v[0, 0], 2000.0)  # y_origin
    assert np.isclose(grid.x_v[8, 9], 1000.0 + 9 * 100.0 + 100.0 / 2)  # Last V-point
    assert np.isclose(grid.y_v[8, 9], 2000.0 + 8 * 200.0)  # Last V-point


def test_grid_coriolis():
    """Test calculation of Coriolis parameter."""
    config = {
        "grid": {
            "nx": 10,
            "ny": 10,
            "dx": 10000.0,
            "dy": 10000.0,
            "x_origin": 0.0,
            "y_origin": 0.0,
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
    config = {"grid": {"nx": 10, "ny": 8, "dx": 100.0, "dy": 100.0}}

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


def test_staggered_masks_all_water():
    config = {"grid": {"nx": 4, "ny": 3, "dx": 1.0, "dy": 1.0}}
    grid = Grid(config)
    mask_t = np.ones((3, 4), dtype=int)
    mu, mv = build_staggered_masks(mask_t)
    assert mu.shape == (3, 5)
    assert mv.shape == (4, 4)
    assert np.all(mu == 1)
    assert np.all(mv == 1)


def test_staggered_masks_with_land_cell():
    # Single land cell should zero adjacent faces.
    config = {"grid": {"nx": 3, "ny": 3, "dx": 1.0, "dy": 1.0}}
    grid = Grid(config)
    mask_t = np.ones((3, 3), dtype=int)
    mask_t[1, 1] = 0  # center land
    mu, mv = build_staggered_masks(mask_t)
    # Adjacent U faces: (1,1) and (1,2)
    assert mu[1, 1] == 0 and mu[1, 2] == 0
    # Adjacent V faces: (1,1) and (2,1)
    assert mv[1, 1] == 0 and mv[2, 1] == 0
    # Non-adjacent faces remain 1
    assert mu[0, 0] == 1 and mv[0, 0] == 1


def test_noslip_flux_masks():
    # Build masks for a plus-shaped water area to ensure interior reductions.
    config = {"grid": {"nx": 5, "ny": 5, "dx": 1.0, "dy": 1.0}}
    grid = Grid(config)
    mask_t = np.zeros((5, 5), dtype=int)
    mask_t[2, 2] = 1
    mask_t[2, 1] = 1
    mask_t[2, 3] = 1
    mask_t[1, 2] = 1
    mask_t[3, 2] = 1
    mu, mv = build_staggered_masks(mask_t)
    mu_ns, mv_ns = apply_noslip_flux_masks(mask_t, mu, mv)
    # Faces around isolated water arms should be zeroed except those fully surrounded by water.
    # The central cross has no 2x2 all-water block except at the center which lacks corners -> resulting interior masks zero.
    interior_u = mu_ns[
        1 : 5 - 0 - 1, 1:5
    ]  # subset representing interior region; simple checks
    assert np.any(interior_u == 0)
    interior_v = mv_ns[1:5, 1 : 5 - 0 - 1]
    assert np.any(interior_v == 0)


def test_corner_mask():
    mask_t = np.array([[1, 1, 1], [1, 0, 1], [1, 1, 1]])
    mq = build_corner_mask(mask_t)
    # Corner mask shape
    assert mq.shape == (mask_t.shape[0] + 1, mask_t.shape[1] + 1)
    # Corner above the land cell (between land and waters) should be 0 only if all contributing T cells 0 -> here mixed so expect product zero? Land at (1,1) zeros interior corners sharing it.
    # Interior corners: (1,1),(1,2),(2,1),(2,2) each include land cell -> zero
    assert mq[1, 1] == 0 and mq[1, 2] == 0 and mq[2, 1] == 0 and mq[2, 2] == 0
    # Outer corner far from land remains 1
    assert mq[0, 0] == 1


def test_model_state_mask_update():
    config = {
        "grid": {"nx": 3, "ny": 3, "dx": 1.0, "dy": 1.0},
        "model": {"timestep": 1.0},
    }
    state = ModelState(config)
    land_mask = np.ones((3, 3), dtype=int)
    land_mask[1, 1] = 0
    state.set_land_mask(land_mask)
    assert (
        state.mask_u is not None
        and state.mask_v is not None
        and state.mask_q is not None
    )
    # Faces adjacent to land zero
    assert state.mask_u[1, 1] == 0 and state.mask_u[1, 2] == 0
    assert state.mask_v[1, 1] == 0 and state.mask_v[2, 1] == 0
    # Corner masks near land zero
    assert state.mask_q[1, 1] == 0
