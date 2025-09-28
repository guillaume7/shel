import numpy as np

from shel.model.solvers.momentum.advection_upwind import upwind_advection
from shel.model.solvers.momentum.biharmonic_diffusion import biharmonic_diffusion


def test_upwind_advection_shape_and_sign():
    field = np.arange(10)
    velocity = np.ones(10)
    dx = 1.0
    tendency = upwind_advection(field, velocity, dx, axis=0)
    assert tendency.shape == field.shape
    # For positive velocity, tendency should be positive for increasing field
    assert np.all(tendency[1:] > 0)


def test_biharmonic_diffusion_zero_field():
    field = np.zeros((5, 5))
    nu4 = 1e-6
    dx = dy = 1.0
    tendency = biharmonic_diffusion(field, nu4, dx, dy)
    assert np.allclose(tendency, 0.0)


def test_biharmonic_diffusion_quadratic():
    # For a quadratic field, biharmonic should be zero in the interior
    x = np.arange(5)
    y = np.arange(5)
    X, Y = np.meshgrid(x, y)
    field = X**2 + Y**2
    nu4 = 1e-6
    dx = dy = 1.0
    tendency = biharmonic_diffusion(field, nu4, dx, dy)
    # Check interior points (excluding boundaries)
    interior = tendency[1:-1, 1:-1]
    assert np.allclose(interior, 0.0, atol=1e-4)
