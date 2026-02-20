"""Pressure gradient term on an Arakawa C-grid.

Conservative flux form matching MATLAB ``ComputeSpaceU_CS``:

    PG_u = -g * H_u * ∂eta/∂x |_U
    PG_v = -g * H_v * ∂eta/∂y |_V

where H_u is the water-column height averaged to U-faces and similarly for V.
This is a tendency in the conservative (Hu) equation, not a bare acceleration.
"""

from __future__ import annotations

import numpy as np

from ..common import avg_x_t_to_u, avg_y_t_to_v, d_dx_t_to_u, d_dy_t_to_v

Array = np.ndarray


def pressure_gradient(
    eta: Array, H: Array, g: float, dx: float, dy: float
) -> tuple[Array, Array]:
    """Compute pressure gradient tendencies at faces (conservative form).

    Parameters
    ----------
    eta : (ny, nx) free-surface elevation at T points
    H   : (ny, nx) total water-column height (eta + d) at T points
    g   : gravitational acceleration (scalar)
    dx, dy : grid spacing (scalars)

    Returns
    -------
    PG_u : (ny, nx+1) pressure gradient tendency at U faces
    PG_v : (ny+1, nx) pressure gradient tendency at V faces
    """
    H_u = avg_x_t_to_u(H)
    H_v = avg_y_t_to_v(H)
    d_eta_dx_u = d_dx_t_to_u(eta, dx)
    d_eta_dy_v = d_dy_t_to_v(eta, dy)
    PG_u = -g * H_u * d_eta_dx_u
    PG_v = -g * H_v * d_eta_dy_v
    return PG_u, PG_v


__all__ = ["pressure_gradient"]
