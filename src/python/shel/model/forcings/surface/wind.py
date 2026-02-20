import numpy as np


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
