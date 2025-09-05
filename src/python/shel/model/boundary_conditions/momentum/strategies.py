from __future__ import annotations

import numpy as np

from ..base import MomentumBC, Array
from ..common.utils import mean_c_along_side, mean_H_along_side


class ClosedBC(MomentumBC):
    name = "closed"

    def apply_uniform(self, U: Array, V: Array) -> None:
        U[:, 0] = 0.0
        U[:, -1] = 0.0
        V[0, :] = 0.0
        V[-1, :] = 0.0

    def apply_side(self, U: Array, V: Array, side: str, **_: object) -> None:
        if side == "west":
            U[:, 0] = 0.0
        elif side == "east":
            U[:, -1] = 0.0
        elif side == "south":
            V[0, :] = 0.0
        elif side == "north":
            V[-1, :] = 0.0


class FreeSlipBC(MomentumBC):
    name = "freeslip"

    def apply_uniform(self, U: Array, V: Array) -> None:
        U[:, 0] = 0.0
        U[:, -1] = 0.0
        V[0, :] = 0.0
        V[-1, :] = 0.0
        if V.shape[0] > 2:
            V[1:-1, 0] = V[1:-1, 1]
            V[1:-1, -1] = V[1:-1, -2]
        if U.shape[0] > 2:
            U[0, 1:-1] = U[1, 1:-1]
            U[-1, 1:-1] = U[-2, 1:-1]

    def apply_side(self, U: Array, V: Array, side: str, **_: object) -> None:
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


class RadiativeSommerfeldBC(MomentumBC):
    name = "radiative"

    def apply_uniform(self, U: Array, V: Array) -> None:
        # No-op for uniform; use side-specific application
        return None

    def apply_side(
        self,
        U: Array,
        V: Array,
        side: str,
        *,
        U_old: Array | None = None,
        V_old: Array | None = None,
        H: Array | None = None,
        g: float | None = None,
        dt: float | None = None,
        dx: float | None = None,
    dy: float | None = None,
    eta_old: Array | None = None,
    eta_ext: Array | None = None,
    relax: float | None = None,
    ) -> None:
        assert U_old is not None and V_old is not None and H is not None
        assert g is not None and dt is not None and dx is not None and dy is not None
        c = mean_c_along_side(H, g, side)
        if side in ("west", "east"):
            r = c * dt / dx
            if side == "west":
                U[:, 0] = U_old[:, 0] - r * (U_old[:, 1] - U_old[:, 0])
            else:
                U[:, -1] = U_old[:, -1] - r * (U_old[:, -1] - U_old[:, -2])
        else:
            r = c * dt / dy
            if side == "south":
                V[0, :] = V_old[0, :] - r * (V_old[1, :] - V_old[0, :])
            else:
                V[-1, :] = V_old[-1, :] - r * (V_old[-1, :] - V_old[-2, :])


class FlatherBC(MomentumBC):
    name = "flather"

    def apply_uniform(self, U: Array, V: Array) -> None:
        # No uniform behavior; use side-specific application only
        return None

    def apply_side(
        self,
        U: Array,
        V: Array,
        side: str,
        *,
        U_old: Array | None = None,
        V_old: Array | None = None,
        H: Array | None = None,
        g: float | None = None,
        dt: float | None = None,
        dx: float | None = None,
        dy: float | None = None,
        eta_old: Array | None = None,
        eta_ext: Array | None = None,
        relax: float | None = None,
    ) -> None:
        # Minimal Flather: if external eta is not provided, no-op.
        if H is None or g is None or eta_old is None or eta_ext is None:
            return None
        c = mean_c_along_side(H, g, side)
        # Use mean H along the selected side to scale
        H_mean = mean_H_along_side(H, side)
        gamma = 1.0 if relax is None else float(relax)
        if gamma < 0.0:
            gamma = 0.0
        if gamma > 1.0:
            gamma = 1.0
        if side in ("west", "east"):
            corr = (
                gamma * (c / H_mean) * (eta_ext[:, 0] - eta_old[:, 0])
                if side == "west"
                else gamma * (c / H_mean) * (eta_ext[:, -1] - eta_old[:, -1])
            )
            if U_old is not None:
                if side == "west":
                    U[:, 0] = U_old[:, 0] + corr
                else:
                    U[:, -1] = U_old[:, -1] + corr
        else:
            corr = (
                gamma * (c / H_mean) * (eta_ext[0, :] - eta_old[0, :])
                if side == "south"
                else gamma * (c / H_mean) * (eta_ext[-1, :] - eta_old[-1, :])
            )
            if V_old is not None:
                if side == "south":
                    V[0, :] = V_old[0, :] + corr
                else:
                    V[-1, :] = V_old[-1, :] + corr
