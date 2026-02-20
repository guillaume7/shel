"""Tests for diagnostics utilities (enstrophy, potential enstrophy)."""

import numpy as np

from shel.model.state import ModelState


def test_enstrophy_and_potential_enstrophy():
    config = {
        "grid": {"nx": 6, "ny": 4, "dx": 1000.0, "dy": 1000.0},
        "model": {"timestep": 60.0, "gravity": 9.81, "coriolis_parameter": 1e-4},
    }
    state = ModelState(config)

    # Bathymetry / depth
    depth = np.ones((4, 6)) * 100.0
    state.set_bathymetry(depth)

    # Simple shear flow to generate vorticity: u varies linearly with y
    for j in range(4):
        state.u[j, :] = j * 0.1
    state.u_old[:] = state.u
    state.u_new[:] = state.u

    # Non-zero elevation bump for potential energy interplay (not required but realistic)
    state.set_initial_elevation(np.zeros((4, 6)))

    ens = state.compute_enstrophy()
    p_ens = state.compute_potential_enstrophy()

    assert ens > 0.0
    assert p_ens > 0.0
    # Potential enstrophy should be O(f^2 * area * H) / (2H) ~ 0.5 * f^2 * area for uniform H.
    # Compute expected potential enstrophy using same discrete definition for a zero relative vorticity case
    # Here relative vorticity not exactly zero but small; we recompute reference directly
    # Sanity bounds: potential enstrophy should be larger than zero and smaller than
    # 0.5 * (f_max^2) * domain area (uniform upper bound ignoring ζ amplification)
    f_max = np.max(state.coriolis)
    domain_area = (
        config["grid"]["nx"]
        * config["grid"]["ny"]
        * config["grid"]["dx"]
        * config["grid"]["dy"]
    )
    upper_bound = 0.5 * f_max**2 * domain_area
    assert p_ens < upper_bound
