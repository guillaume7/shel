__all__ = ["wind_stress", "surface_pressure"]


import numpy as np

"""
Surface forcing module: wind stress and pressure (SHEL best practice).
"""


def wind_stress(U_air, U_surface, rho_air=1.225, Cd=1.3e-3):
    """
    Compute wind stress (vectorized) at the surface.
    tau = rho_air * Cd * |U_air - U_surface| * (U_air - U_surface)
    Args:
        U_air: array-like, wind velocity at 10m (m/s)
        U_surface: array-like, surface velocity (m/s)
        rho_air: air density (kg/m^3)
        Cd: drag coefficient (dimensionless)
    Returns:
        tau: wind stress (N/m^2), same shape as U_air
    """
    rel = np.asarray(U_air) - np.asarray(U_surface)
    mag = np.linalg.norm(rel, axis=0) if rel.ndim > 1 else np.abs(rel)
    tau = rho_air * Cd * mag * rel
    return tau


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
