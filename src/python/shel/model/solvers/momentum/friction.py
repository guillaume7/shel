"""Linear bottom friction (drag) tendency for barotropic momentum (minimal).

Implements dU/dt = -r * U and dV/dt = -r * V, where r is a scalar drag
coefficient. Shapes are preserved and NaNs (e.g., at boundaries) are propagated
by multiplication.
"""

from __future__ import annotations

import numpy as np

Array = np.ndarray


def bottom_drag_tendency(U: Array, V: Array, r: float) -> tuple[Array, Array]:
    """Compute linear bottom drag tendencies for U and V faces.

    Parameters
    - U: (ny, nx+1) velocity at U faces
    - V: (ny+1, nx) velocity at V faces
    - r: linear drag coefficient (s^-1)
    """
    return -r * U, -r * V


__all__ = ["bottom_drag_tendency"]
