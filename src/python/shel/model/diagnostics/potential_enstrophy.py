import numpy as np


def potential_enstrophy_field(pv):
    """
    Pointwise potential enstrophy: 0.5 * PV**2
    Args:
        pv: array (potential vorticity)
    Returns:
        potential enstrophy: array
    """
    return 0.5 * np.square(pv)


def integrated_potential_enstrophy(pv, dx, dy, mask=None):
    """
    Domain-integrated potential enstrophy.
    Args:
        pv: array
        dx, dy: grid spacing
        mask: optional boolean mask (True=wet)
    Returns:
        scalar potential enstrophy
    """
    penst = potential_enstrophy_field(pv)
    if mask is not None:
        penst = np.where(mask, penst, 0.0)
    return np.sum(penst) * dx * dy
