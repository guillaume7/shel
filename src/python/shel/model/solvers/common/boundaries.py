"""
Functional boundary-condition helpers for solver steps (ministep path).

These operate directly on staggered velocity arrays (U at west/east faces,
V at south/north faces) and enforce simple strategies without depending on
ModelState.
"""
from __future__ import annotations

from typing import Dict

import numpy as np

Array = np.ndarray


def apply_closed(U: Array, V: Array) -> None:
    """Enforce closed (no normal flow) boundaries in-place.

    Sets U at west/east faces and V at south/north faces to zero.
    """
    U[:, 0] = 0.0
    U[:, -1] = 0.0
    V[0, :] = 0.0
    V[-1, :] = 0.0


def apply_freeslip(U: Array, V: Array) -> None:
        """Enforce free-slip boundaries in-place.

        - Zero normal flow at domain edges (same as closed)
        - Zero tangential shear at walls via zero-gradient on tangential component
            (copy nearest interior value to boundary-adjacent locations)
        """
        # West/East: normal component is U at faces -> zero it
        U[:, 0] = 0.0
        U[:, -1] = 0.0
        # Tangential along west/east walls is V; zero-gradient across wall
        if V.shape[0] > 2:
                V[1:-1, 0] = V[1:-1, 1]
                V[1:-1, -1] = V[1:-1, -2]

        # South/North: normal component is V at faces -> zero it
        V[0, :] = 0.0
        V[-1, :] = 0.0
        # Tangential along south/north walls is U; zero-gradient across wall
        if U.shape[0] > 2:
                U[0, 1:-1] = U[1, 1:-1]
                U[-1, 1:-1] = U[-2, 1:-1]


def _mean_c_along_side(H: Array, g: float, side: str) -> float:
        if side == "west":
            edge = H[:, 0]
        elif side == "east":
            edge = H[:, -1]
        elif side == "south":
            edge = H[0, :]
        elif side == "north":
            edge = H[-1, :]
        else:
            edge = H
        return float(np.sqrt(g * float(np.mean(edge))))


def apply_closed_side(U: Array, V: Array, side: str) -> None:
        if side == "west":
            U[:, 0] = 0.0
        elif side == "east":
            U[:, -1] = 0.0
        elif side == "south":
            V[0, :] = 0.0
        elif side == "north":
            V[-1, :] = 0.0


def apply_freeslip_side(U: Array, V: Array, side: str) -> None:
        if side == "west":
            U[:, 0] = 0.0
            if V.shape[0] > 2:
                V[1:-1, 0] = V[1:-1, 1]
        elif side == "east":
            U[:, -1] = 0.0
            if V.shape[0] > 2:
                V[1:-1, -1] = V[1:-1, -2]
        elif side == "south":
            V[0, :] = 0.0
            if U.shape[0] > 2:
                U[0, 1:-1] = U[1, 1:-1]
        elif side == "north":
            V[-1, :] = 0.0
            if U.shape[0] > 2:
                U[-1, 1:-1] = U[-2, 1:-1]


def apply_radiative_momentum_side(
        U_next: Array,
        V_next: Array,
        U_old: Array,
        V_old: Array,
        H: Array,
        g: float,
        dt: float,
        dx: float,
        dy: float,
        side: str,
    ) -> None:
        """Sommerfeld-type update for normal velocity at boundary using old fields.

        u_b^{n+1} = u_b^{n} - r (u_i^{n} - u_b^{n}), r = c dt / Δ
        """
        c = _mean_c_along_side(H, g, side)
        if side == "west":
            r = c * dt / dx
            U_next[:, 0] = U_old[:, 0] - r * (U_old[:, 1] - U_old[:, 0])
        elif side == "east":
            r = c * dt / dx
            U_next[:, -1] = U_old[:, -1] - r * (U_old[:, -1] - U_old[:, -2])
        elif side == "south":
            r = c * dt / dy
            V_next[0, :] = V_old[0, :] - r * (V_old[1, :] - V_old[0, :])
        elif side == "north":
            r = c * dt / dy
            V_next[-1, :] = V_old[-1, :] - r * (V_old[-1, :] - V_old[-2, :])


def apply_radiative_eta_side(
        eta_next: Array,
        eta_old: Array,
        H: Array,
        g: float,
        dt: float,
        dx: float,
        dy: float,
        side: str,
    ) -> None:
        c = _mean_c_along_side(H, g, side)
        if side == "west":
            r = c * dt / dx
            eta_next[:, 0] = eta_old[:, 0] - r * (eta_old[:, 1] - eta_old[:, 0])
        elif side == "east":
            r = c * dt / dx
            eta_next[:, -1] = eta_old[:, -1] - r * (eta_old[:, -1] - eta_old[:, -2])
        elif side == "south":
            r = c * dt / dy
            eta_next[0, :] = eta_old[0, :] - r * (eta_old[1, :] - eta_old[0, :])
        elif side == "north":
            r = c * dt / dy
            eta_next[-1, :] = eta_old[-1, :] - r * (eta_old[-1, :] - eta_old[-2, :])


def apply_momentum_per_side(
        U_next: Array,
        V_next: Array,
        U_old: Array,
        V_old: Array,
        H: Array,
        g: float,
        dt: float,
        dx: float,
        dy: float,
        bc_sides: Dict[str, str],
    ) -> None:
        for side, bct in bc_sides.items():
            if bct == "closed":
                apply_closed_side(U_next, V_next, side)
            elif bct == "freeslip":
                apply_freeslip_side(U_next, V_next, side)
            elif bct == "radiative":
                apply_radiative_momentum_side(U_next, V_next, U_old, V_old, H, g, dt, dx, dy, side)
            else:
                apply_closed_side(U_next, V_next, side)


def apply_eta_per_side(
        eta_next: Array,
        eta_old: Array,
        H: Array,
        g: float,
        dt: float,
        dx: float,
        dy: float,
        bc_sides: Dict[str, str],
    ) -> None:
        for side, bct in bc_sides.items():
            if bct == "radiative":
                apply_radiative_eta_side(eta_next, eta_old, H, g, dt, dx, dy, side)
            # closed/freeslip for eta are implicitly handled via velocities


__all__ = [
        "apply_closed",
        "apply_freeslip",
        "apply_closed_side",
        "apply_freeslip_side",
        "apply_radiative_momentum_side",
        "apply_radiative_eta_side",
        "apply_momentum_per_side",
        "apply_eta_per_side",
    ]
