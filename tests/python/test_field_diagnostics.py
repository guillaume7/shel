from __future__ import annotations
import numpy as np
from shel.model.grid import Grid
from shel.model.initial_conditions.factory import create_velocity
from shel.model.diagnostics.fields import FieldDiagnostics

def make_grid(nx=32, ny=24):
    cfg = {"grid": {"nx": nx, "ny": ny, "dx": 1000.0, "dy": 1000.0}}
    return Grid(cfg)

def analytic_linear_field(grid: Grid):
    # Construct a linear velocity field u = ax, v = by so divergence = a + b (constant)
    a, b = 1e-5, -2e-5
    # approximate x at U points and y at V points via grid coordinates already provided
    u = a * (grid.x_u - grid.x_u.min())
    v = b * (grid.y_v - grid.y_v.min())
    return u, v, a, b

def test_divergence_constant_linear():
    g = make_grid()
    u, v, a, b = analytic_linear_field(g)
    div = FieldDiagnostics.divergence(u, v, g)
    # interior mean close to a+b
    assert np.isclose(div.mean(), a + b, rtol=0, atol=5e-7)


def test_solid_body_vorticity_and_divergence():
    g = make_grid()
    u_sb, v_sb = create_velocity("solid_body").build(g)
    div = FieldDiagnostics.divergence(u_sb, v_sb, g)
    # Solid body rotation is (nearly) non-divergent
    assert np.allclose(div, 0.0, atol=1e-10)


def test_shapes_and_nonneg_quadratics():
    g = make_grid()
    # Use shear velocity IC for non-trivial gradient
    u_sh, v_sh = create_velocity("shear").build(g)
    shear = FieldDiagnostics.shear_rate(u_sh, v_sh, g)
    stretch = FieldDiagnostics.stretch_rate(u_sh, v_sh, g)
    div = FieldDiagnostics.divergence(u_sh, v_sh, g)
    assert shear.shape == (g.ny + 1, g.nx + 1)
    assert stretch.shape == (g.ny, g.nx)
    assert div.shape == (g.ny, g.nx)
    # Shear rate squared should be non-negative (sanity)
    assert (shear**2 >= 0).all()
