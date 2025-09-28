import numpy as np


def enstrophy_field(vorticity):
    """
    Pointwise enstrophy: 0.5 * vorticity**2
    Args:
        vorticity: array
    Returns:
        enstrophy: array
    """
    return 0.5 * np.square(vorticity)


def integrated_enstrophy(vorticity, dx, dy, mask=None):
    """
    Domain-integrated enstrophy.
    Args:
        vorticity: array
        dx, dy: grid spacing
        mask: optional boolean mask (True=wet)
    Returns:
        scalar enstrophy
    """
    enst = enstrophy_field(vorticity)
    if mask is not None:
        enst = np.where(mask, enst, 0.0)
    return np.sum(enst) * dx * dy
