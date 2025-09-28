import numpy as np


def upwind_advection(field, velocity, dx, axis=0):
    """
    2nd order upwind advection for a scalar field.
    Args:
        field: array, scalar field to advect
        velocity: array, velocity field (same shape)
        dx: grid spacing
        axis: axis along which to compute advection (0=y, 1=x)
    Returns:
        tendency: array, same shape as field
    """
    f = np.asarray(field)
    u = np.asarray(velocity)
    # Shifted fields for upwind
    f_p1 = np.roll(f, -1, axis=axis)
    f_m1 = np.roll(f, 1, axis=axis)
    # Upwind scheme
    tendency = np.where(
        u > 0,
        (3 * f - 4 * f_m1 + f_m1) / (2 * dx) * u,
        (-3 * f + 4 * f_p1 - f_p1) / (2 * dx) * u,
    )
    return tendency
