import numpy as np

from shel.model.solvers.tracer.advection_centered import centered_tracer_advection


def test_centered_tracer_advection_shape():
    tracer = np.ones((5, 5))
    u = np.ones((5, 5))
    v = np.ones((5, 5))
    dx = dy = 1.0
    tendency = centered_tracer_advection(tracer, u, v, dx, dy)
    assert tendency.shape == tracer.shape


def test_centered_tracer_advection_mask():
    tracer = np.random.rand(5, 5)
    u = np.ones((5, 5))
    v = np.ones((5, 5))
    dx = dy = 1.0
    mask = np.zeros((5, 5), dtype=bool)
    mask[2, 2] = True
    tendency = centered_tracer_advection(tracer, u, v, dx, dy, mask=mask)
    # Only masked location should be nonzero
    assert np.count_nonzero(tendency) == 1


def test_centered_tracer_advection_conservation():
    tracer = np.random.rand(5, 5)
    u = np.zeros((5, 5))
    v = np.zeros((5, 5))
    dx = dy = 1.0
    tendency = centered_tracer_advection(tracer, u, v, dx, dy)
    # No velocity: tendency should be zero
    assert np.allclose(tendency, 0.0)
