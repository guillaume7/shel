__all__ = ["linear_drag", "quadratic_drag"]


import numpy as np

"""
Bottom drag coefficient utilities (SHEL best practice).
Supports linear and quadratic drag, spatially variable.
"""


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
