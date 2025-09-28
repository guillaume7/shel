import numpy as np


def biharmonic_diffusion(field, nu4, dx, dy):
    """
    2D biharmonic diffusion operator (∇⁴) for a scalar field.
    Args:
        field: array, scalar field
        nu4: biharmonic viscosity coefficient
        dx, dy: grid spacing
    Returns:
        tendency: array, same shape as field
    """
    f = np.asarray(field)
    # Laplacian
    lap = (
        -4 * f
        + np.roll(f, 1, axis=0)
        + np.roll(f, -1, axis=0)
        + np.roll(f, 1, axis=1)
        + np.roll(f, -1, axis=1)
    ) / (dx * dy)
    # Biharmonic: Laplacian of Laplacian
    lap2 = (
        -4 * lap
        + np.roll(lap, 1, axis=0)
        + np.roll(lap, -1, axis=0)
        + np.roll(lap, 1, axis=1)
        + np.roll(lap, -1, axis=1)
    ) / (dx * dy)
    return nu4 * lap2
