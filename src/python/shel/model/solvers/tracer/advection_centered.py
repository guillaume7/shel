import numpy as np


def centered_tracer_advection(tracer, u, v, dx, dy, mask=None):
    """
    Centered advection for tracer field on C-grid.
    Args:
        tracer: array, tracer field (T points)
        u, v: arrays, velocity fields (U, V points)
        dx, dy: grid spacing
        mask: optional boolean mask (True=wet)
    Returns:
        tendency: array, same shape as tracer
    """
    # Compute fluxes (simple central difference)
    flux_x = (np.roll(tracer, -1, axis=1) - np.roll(tracer, 1, axis=1)) / (2 * dx) * u
    flux_y = (np.roll(tracer, -1, axis=0) - np.roll(tracer, 1, axis=0)) / (2 * dy) * v
    tendency = -(flux_x + flux_y)
    if mask is not None:
        tendency = np.where(mask, tendency, 0.0)
    return tendency
