import numpy as np


def quadratic_drag(U, Cd, H):
    """
    Quadratic bottom drag tendency: -Cd * |U| * U / H
    Args:
        U: array-like, velocity at bottom (m/s)
        Cd: array-like or float, drag coefficient (dimensionless)
        H: array-like or float, water depth (m)
    Returns:
        drag tendency (same shape as U)
    """
    mag = np.abs(U)
    return -np.asarray(Cd) * mag * np.asarray(U) / np.asarray(H)
