import numpy as np


def linear_drag(U, r):
    """
    Linear bottom drag tendency: -r * U
    Args:
        U: array-like, velocity at bottom (m/s)
        r: array-like or float, drag coefficient (s^-1)
    Returns:
        drag tendency (same shape as U)
    """
    return -np.asarray(r) * np.asarray(U)
