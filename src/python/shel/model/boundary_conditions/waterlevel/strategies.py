from __future__ import annotations

from ..base import EtaBC, Array
from ..common.utils import mean_c_along_side


class RadiativeSommerfeldEtaBC(EtaBC):
    name = "radiative"

    def apply_side_eta(
        self,
        eta_next: Array,
        side: str,
        *,
        eta_old: Array | None = None,
        eta_ext: Array | None = None,
        H: Array | None = None,
        g: float | None = None,
        dt: float | None = None,
        dx: float | None = None,
        dy: float | None = None,
    ) -> None:
        assert eta_old is not None and H is not None
        assert g is not None and dt is not None and dx is not None and dy is not None
        c = mean_c_along_side(H, g, side)
        if side in ("west", "east"):
            r = c * dt / dx
            if side == "west":
                eta_next[:, 0] = eta_old[:, 0] - r * (eta_old[:, 1] - eta_old[:, 0])
            else:
                eta_next[:, -1] = eta_old[:, -1] - r * (eta_old[:, -1] - eta_old[:, -2])
        else:
            r = c * dt / dy
            if side == "south":
                eta_next[0, :] = eta_old[0, :] - r * (eta_old[1, :] - eta_old[0, :])
            else:
                eta_next[-1, :] = eta_old[-1, :] - r * (eta_old[-1, :] - eta_old[-2, :])


class FlatherEtaBC(EtaBC):
    name = "flather"

    def apply_side_eta(
        self,
        eta_next: Array,
        side: str,
        *,
        eta_old: Array | None = None,
        eta_ext: Array | None = None,
        H: Array | None = None,
        g: float | None = None,
        dt: float | None = None,
        dx: float | None = None,
        dy: float | None = None,
    ) -> None:
        # Simple Dirichlet: if external eta provided, set boundary to it; otherwise no-op
        if eta_ext is None:
            return None
        if side == "west":
            eta_next[:, 0] = eta_ext[:, 0]
        elif side == "east":
            eta_next[:, -1] = eta_ext[:, -1]
        elif side == "south":
            eta_next[0, :] = eta_ext[0, :]
        elif side == "north":
            eta_next[-1, :] = eta_ext[-1, :]
