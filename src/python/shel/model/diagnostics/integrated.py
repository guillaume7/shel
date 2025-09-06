"""Domain-integrated diagnostics for the SHEL model.

This module defines :class:`IntegratedDiagnostics` providing static
methods that compute scalar, domain-integrated physical quantities
useful for monitoring numerical stability, conservation and physical
fidelity (energies, enstrophy, volume, etc.).

The design deliberately follows a functional style (stateless static
methods) emulating the original MATLAB functional utilities while
ensuring clear separation between integrated (scalar) metrics and
pointwise diagnostic fields (implemented in ``fields.py``).
"""

from __future__ import annotations

from typing import Dict, Tuple

import numpy as np
from numpy.typing import NDArray

from shel.model.grid import Grid

RHO0: float = 1000.0  # Reference density (kg/m^3)


class IntegratedDiagnostics:
    """Namespace for integrated (domain aggregated) diagnostics.

    All methods are static and expect arrays on an Arakawa C-grid.
    """

    # ------------------------------ helpers ------------------------------
    @staticmethod
    def _interpolate_H_to_u(H: NDArray, grid: Grid) -> NDArray:
        ny, nx = grid.ny, grid.nx
        out = np.zeros((ny, nx + 1), dtype=H.dtype)
        out[:, 1:nx] = 0.5 * (H[:, 0 : nx - 1] + H[:, 1:nx])
        out[:, 0] = H[:, 0]
        out[:, nx] = H[:, nx - 1]
        return out

    @staticmethod
    def _interpolate_H_to_v(H: NDArray, grid: Grid) -> NDArray:
        ny, nx = grid.ny, grid.nx
        out = np.zeros((ny + 1, nx), dtype=H.dtype)
        out[1:ny, :] = 0.5 * (H[0 : ny - 1, :] + H[1:ny, :])
        out[0, :] = H[0, :]
        out[ny, :] = H[ny - 1, :]
        return out

    @staticmethod
    def _vorticity_q(u: NDArray, v: NDArray, grid: Grid) -> NDArray:
        nx, ny = grid.nx, grid.ny
        dx, dy = grid.dx, grid.dy
        zeta = np.zeros((ny + 1, nx + 1), dtype=u.dtype)
        zeta[1:ny, 1:nx] = (v[1:ny, 1:nx] - v[1:ny, 0 : nx - 1]) / dx - (
            u[1:ny, 1:nx] - u[0 : ny - 1, 1:nx]
        ) / dy
        return zeta

    @staticmethod
    def _vorticity_t(u: NDArray, v: NDArray, grid: Grid) -> NDArray:
        zeta_q = IntegratedDiagnostics._vorticity_q(u, v, grid)
        return 0.25 * (
            zeta_q[0:-1, 0:-1] + zeta_q[1:, 0:-1] + zeta_q[0:-1, 1:] + zeta_q[1:, 1:]
        )

    # ------------------------------ energies ------------------------------
    @staticmethod
    def kinetic_energy(
        u: NDArray, v: NDArray, H: NDArray, grid: Grid, rho: float = RHO0
    ) -> float:
        area_t, area_u, area_v = grid.compute_areas()
        H_u = IntegratedDiagnostics._interpolate_H_to_u(H, grid)
        H_v = IntegratedDiagnostics._interpolate_H_to_v(H, grid)
        ke_u = 0.5 * rho * np.sum(H_u[:, 1:-1] * u[:, 1:-1] ** 2 * area_u[:, 1:-1])
        ke_v = 0.5 * rho * np.sum(H_v[1:-1, :] * v[1:-1, :] ** 2 * area_v[1:-1, :])
        return float(ke_u + ke_v)

    @staticmethod
    def potential_energy(
        eta: NDArray, grid: Grid, gravity: float, rho: float = RHO0
    ) -> float:
        area_t, _, _ = grid.compute_areas()
        return float(0.5 * rho * gravity * np.sum(eta**2 * area_t))

    @staticmethod
    def total_energy(
        u: NDArray,
        v: NDArray,
        eta: NDArray,
        H: NDArray,
        grid: Grid,
        gravity: float,
        rho: float = RHO0,
    ) -> float:
        return IntegratedDiagnostics.kinetic_energy(
            u, v, H, grid, rho
        ) + IntegratedDiagnostics.potential_energy(eta, grid, gravity, rho)

    @staticmethod
    def volume(H: NDArray, grid: Grid) -> float:
        area_t, _, _ = grid.compute_areas()
        return float(np.sum(H * area_t))

    # ------------------------------ enstrophy ------------------------------
    @staticmethod
    def enstrophy(u: NDArray, v: NDArray, grid: Grid) -> float:
        zeta = IntegratedDiagnostics._vorticity_q(u, v, grid)
        zeta_int = zeta[1:-1, 1:-1]
        area = grid.dx * grid.dy
        return float(0.5 * np.sum(zeta_int**2) * area)

    @staticmethod
    def potential_enstrophy(
        u: NDArray, v: NDArray, H: NDArray, f: NDArray, grid: Grid
    ) -> float:
        zeta_t = IntegratedDiagnostics._vorticity_t(u, v, grid)
        q = (zeta_t + f) / H
        area_t, _, _ = grid.compute_areas()
        return float(0.5 * np.sum(H * q**2 * area_t))

    # ------------------------------ collections ------------------------------
    @staticmethod
    def energy_partition(
        u: NDArray,
        v: NDArray,
        eta: NDArray,
        H: NDArray,
        grid: Grid,
        gravity: float,
        rho: float = RHO0,
    ) -> Tuple[float, float, float]:
        ke = IntegratedDiagnostics.kinetic_energy(u, v, H, grid, rho)
        pe = IntegratedDiagnostics.potential_energy(eta, grid, gravity, rho)
        return ke, pe, ke + pe

    @staticmethod
    def as_dict(
        u: NDArray,
        v: NDArray,
        eta: NDArray,
        H: NDArray,
        grid: Grid,
        gravity: float,
        coriolis: NDArray,
    ) -> Dict[str, float]:
        ke, pe, te = IntegratedDiagnostics.energy_partition(u, v, eta, H, grid, gravity)
        return {
            "kinetic_energy": ke,
            "potential_energy": pe,
            "total_energy": te,
            "volume": IntegratedDiagnostics.volume(H, grid),
            "enstrophy": IntegratedDiagnostics.enstrophy(u, v, grid),
            "potential_enstrophy": IntegratedDiagnostics.potential_enstrophy(
                u, v, H, coriolis, grid
            ),
        }


__all__ = ["IntegratedDiagnostics", "RHO0"]
