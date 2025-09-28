import numpy as np


def kinetic_energy(u, v, H):
    """Total kinetic energy: 0.5 * H * (u^2 + v^2)"""
    return 0.5 * np.sum(H * (u**2 + v**2))


def potential_energy(eta, H, g=9.81):
    """Total potential energy: 0.5 * g * (eta^2)"""
    return 0.5 * g * np.sum(eta**2)


def wind_work(tau, u, v):
    """Work rate from wind stress: sum(tau * (u, v))"""
    # tau, u, v are arrays; sum over domain
    if isinstance(tau, tuple) or isinstance(tau, list):
        tau_u, tau_v = tau
        return np.sum(tau_u * u + tau_v * v)
    return np.sum(tau * u)  # scalar wind stress


def drag_dissipation(drag, u, v):
    """Dissipation rate from drag: sum(drag * (u, v))"""
    if isinstance(drag, tuple) or isinstance(drag, list):
        drag_u, drag_v = drag
        return np.sum(drag_u * u + drag_v * v)
    return np.sum(drag * u)
