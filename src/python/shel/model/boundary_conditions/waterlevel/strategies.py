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
        relax: float | None = None,
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
        relax: float | None = None,
    ) -> None:
        # Simple Dirichlet: if external eta provided, set boundary to it; otherwise no-op
        if eta_ext is None:
            return None
        alpha = 1.0 if relax is None else float(relax)
        alpha = 0.0 if alpha < 0.0 else (1.0 if alpha > 1.0 else alpha)
        if side == "west":
            eta_next[:, 0] = (1 - alpha) * eta_next[:, 0] + alpha * eta_ext[:, 0]
        elif side == "east":
            eta_next[:, -1] = (1 - alpha) * eta_next[:, -1] + alpha * eta_ext[:, -1]
        elif side == "south":
            eta_next[0, :] = (1 - alpha) * eta_next[0, :] + alpha * eta_ext[0, :]
        elif side == "north":
            eta_next[-1, :] = (1 - alpha) * eta_next[-1, :] + alpha * eta_ext[-1, :]


class DirichletEtaBC(EtaBC):
    name = "dirichlet"

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
        relax: float | None = None,
    ) -> None:
        if eta_ext is None:
            return None
        alpha = 1.0 if relax is None else float(relax)
        alpha = 0.0 if alpha < 0.0 else (1.0 if alpha > 1.0 else alpha)
        if side == "west":
            eta_next[:, 0] = (1 - alpha) * eta_next[:, 0] + alpha * eta_ext[:, 0]
        elif side == "east":
            eta_next[:, -1] = (1 - alpha) * eta_next[:, -1] + alpha * eta_ext[:, -1]
        elif side == "south":
            eta_next[0, :] = (1 - alpha) * eta_next[0, :] + alpha * eta_ext[0, :]
        elif side == "north":
            eta_next[-1, :] = (1 - alpha) * eta_next[-1, :] + alpha * eta_ext[-1, :]
