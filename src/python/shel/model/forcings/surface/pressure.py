import numpy as np


def surface_pressure(P_atm, P_ref=101325.0):
    """
    Compute surface pressure anomaly (Pa).
    Args:
        P_atm: array-like, atmospheric pressure (Pa)
        P_ref: reference pressure (Pa)
    Returns:
        pressure anomaly (Pa)
    """
    return np.asarray(P_atm) - P_ref
