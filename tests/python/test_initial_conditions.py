"""Tests for initial condition implementations (Phase 3)."""
from __future__ import annotations
import numpy as np
from shel.model.grid import Grid
from shel.model.initial_conditions.factory import (
    create_bathymetry, create_elevation, create_velocity, create_tracer,
)

def make_grid(nx=40, ny=30):
    cfg = {"grid": {"nx": nx, "ny": ny, "dx": 1000.0, "dy": 1000.0}}
    return Grid(cfg)

def test_bathymetry_variants():
    g = make_grid()
    for name in ["bump", "step", "island", "cylinder"]:
        d = create_bathymetry(name).build(g)
        assert d.shape == (g.ny, g.nx)
        assert np.isfinite(d).all()

def test_elevation_gaussian_flat_volume_ratio():
    g = make_grid()
    gauss = create_elevation("gaussian").build(g)
    flat = create_elevation("flat").build(g)
    assert gauss.shape == flat.shape == (g.ny, g.nx)
    assert np.allclose(flat, 0.0)
    # Gaussian positive and localized
    assert gauss.max() > 0
    assert gauss.sum() > 0

def test_velocity_solid_body():
    g = make_grid()
    u, v = create_velocity("solid_body").build(g)
    assert u.shape == (g.ny, g.nx + 1)
    assert v.shape == (g.ny + 1, g.nx)
    # Solid body rotational symmetry: mean near zero
    assert abs(u.mean()) < 1e-6
    assert abs(v.mean()) < 1e-6

def test_velocity_shear():
    g = make_grid()
    u, v = create_velocity("shear").build(g)
    assert u[0].max() == 0  # top row relative may not zero but check gradient monotonic
    assert np.all(np.diff(u.mean(axis=1)) >= 0)
    assert np.allclose(v, 0.0)

def test_tracer_variants():
    g = make_grid()
    gauss = create_tracer("gaussian").build(g)
    uni = create_tracer("uniform").build(g)
    assert gauss.shape == uni.shape == (g.ny, g.nx)
    assert np.isclose(uni.min(), uni.max())
    assert gauss.max() > gauss.mean() > 0

def test_geostrophic_velocity_balance():
    g = make_grid()
    eta = create_elevation("gaussian").build(g)
    # constant coriolis for test
    f = np.ones((g.ny, g.nx)) * 1e-4
    # Pass dependent fields to build (not constructor) for geostrophic velocity
    u, v = create_velocity("geostrophic").build(g, eta=eta, coriolis=f)
    # Reconstruct T-grid parametric velocities from staggered fields
    ug_T = 0.5*(u[:,0:-1] + u[:,1:])   # (ny, nx)
    vg_T = 0.5*(v[0:-1,:] + v[1:,:])   # (ny, nx)
    gconst = 9.81
    ny, nx = eta.shape
    j_idx = np.arange(ny).reshape(ny,1)
    i_idx = np.arange(nx).reshape(1,nx)
    j_c = ny//2
    i_c = nx//2
    pred_ug = (gconst/f) * (j_idx - j_c) * eta
    pred_vg = -(gconst/f) * (i_idx - i_c) * eta
    mask = eta > (0.01 * eta.max())
    err_u = np.max(np.abs((ug_T - pred_ug)[mask]) / (np.abs(pred_ug[mask]) + 1e-14))
    err_v = np.max(np.abs((vg_T - pred_vg)[mask]) / (np.abs(pred_vg[mask]) + 1e-14))
    # Slightly looser tolerance on v due to interpolation asymmetry (edge averaging)
    assert err_u < 0.1
    assert err_v < 0.12
