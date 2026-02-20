from __future__ import annotations

from ..base import Array, TracerBC


class TracerClosedBC(TracerBC):
    name = "closed"

    def apply_side_tracer(
        self,
        C: Array,
        side: str,
        *,
        C_old: Array | None = None,
        U: Array | None = None,
        V: Array | None = None,
        dt: float | None = None,
        dx: float | None = None,
        dy: float | None = None,
    ) -> None:
        # Zero-gradient (copy interior) to avoid artificial fluxes
        if side == "west" and C.shape[1] > 1:
            C[:, 0] = C[:, 1]
        elif side == "east" and C.shape[1] > 1:
            C[:, -1] = C[:, -2]
        elif side == "south" and C.shape[0] > 1:
            C[0, :] = C[1, :]
        elif side == "north" and C.shape[0] > 1:
            C[-1, :] = C[-2, :]


class TracerRadiativeBC(TracerBC):
    name = "radiative"

    def apply_side_tracer(
        self,
        C: Array,
        side: str,
        *,
        C_old: Array | None = None,
        U: Array | None = None,
        V: Array | None = None,
        dt: float | None = None,
        dx: float | None = None,
        dy: float | None = None,
    ) -> None:
        # Simple outward one-sided copy consistent with advection outflow
        if side == "west" and C.shape[1] > 1:
            C[:, 0] = C[:, 1]
        elif side == "east" and C.shape[1] > 1:
            C[:, -1] = C[:, -2]
        elif side == "south" and C.shape[0] > 1:
            C[0, :] = C[1, :]
        elif side == "north" and C.shape[0] > 1:
            C[-1, :] = C[-2, :]


__all__ = ["TracerClosedBC", "TracerRadiativeBC"]
