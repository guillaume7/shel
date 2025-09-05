"""Local field diagnostics for the SHEL model.

This module defines :class:`FieldDiagnostics` which provides static
methods that compute spatially varying (pointwise) diagnostic fields
(e.g. vorticity, Okubo–Weiss) separate from the domain-integrated
scalars in ``integrated.py``.
"""
from __future__ import annotations

from typing import Dict

import numpy as np
from numpy.typing import NDArray

from shel.model.grid import Grid


class FieldDiagnostics:
    """Namespace for local diagnostic fields (static methods)."""

    @staticmethod
    def vorticity(u: NDArray, v: NDArray, grid: Grid) -> NDArray:
        nx, ny = grid.nx, grid.ny
        dx, dy = grid.dx, grid.dy
        zeta = np.zeros((ny + 1, nx + 1), dtype=u.dtype)
        zeta[1:ny, 1:nx] = (
            (v[1:ny, 1:nx] - v[1:ny, 0 : nx - 1]) / dx
            - (u[1:ny, 1:nx] - u[0 : ny - 1, 1:nx]) / dy
        )
        return zeta

    @staticmethod
    def okubo_weiss(u: NDArray, v: NDArray, grid: Grid) -> NDArray:
        # Use the faithful MATLAB-style pipeline implemented in intermediates()
        return FieldDiagnostics.intermediates(u, v, grid)["okubo_weiss"]

    @staticmethod
    def divergence(u: NDArray, v: NDArray, grid: Grid) -> NDArray:
        """Horizontal divergence on T cells (ny, nx)."""
        return FieldDiagnostics.intermediates(u, v, grid)["divergence_t"]

    @staticmethod
    def shear_rate(u: NDArray, v: NDArray, grid: Grid) -> NDArray:
        """Shear rate on corner/Q grid (ny+1, nx+1)."""
        return FieldDiagnostics.intermediates(u, v, grid)["shearrate_w"]

    @staticmethod
    def stretch_rate(u: NDArray, v: NDArray, grid: Grid) -> NDArray:
        """Stretch rate on T cells (ny, nx)."""
        return FieldDiagnostics.intermediates(u, v, grid)["strechrate_t"]  # MATLAB-spelled key

    # ------------------------------------------------------------------
    # MATLAB-faithful intermediate diagnostics (curl, shear, stretch, etc.)
    # We mimic the original model_handles.m logic by working in a transposed
    # orientation so indices map cleanly onto the MATLAB (M,N) ordering. After
    # computing we transpose back to the native (ny, nx) orientation.
    # ------------------------------------------------------------------
    @staticmethod
    def intermediates(u: NDArray, v: NDArray, grid: Grid) -> Dict[str, NDArray]:
        ny, nx = grid.ny, grid.nx
        dx, dy = grid.dx, grid.dy
        # Build water mask (MATLAB mask==1 water, here grid.mask==0 water)
        water_mask_native = (grid.mask == 0).astype(u.dtype)  # (ny, nx)
        # Transpose to MATLAB orientation where M=nx (x extent), N=ny (y extent)
        # so arrays become U(M+1,N), V(M,N+1)
        U = u.T  # (nx+1, ny)
        V = v.T  # (nx, ny+1)
        mask_T = water_mask_native.T  # (nx, ny)
        M, N = nx, ny
        # Flux masks
        mask_u = np.ones((M + 1, N), dtype=u.dtype)
        mask_v = np.ones((M, N + 1), dtype=u.dtype)
        # Where T cell is land (mask_T==0) zero adjacent flux masks
        land_i, land_j = np.where(mask_T == 0)
        for ii, jj in zip(land_i, land_j):
            mask_u[ii, jj] = 0.0
            mask_u[ii + 1, jj] = 0.0
            mask_v[ii, jj] = 0.0
            mask_v[ii, jj + 1] = 0.0
        dA = dx * dy
        # ---------------- curl_w (W/F cell – corners) -----------------
        curl_w = np.zeros((M + 1, N + 1), dtype=u.dtype)
        # interior 2..M, 2..N in MATLAB => indices 1:M-1,1:N-1 zero-based
        if M > 0 and N > 0:
            curl_w[1:M, 1:N] = (
                (mask_v[1:M, 1:N] * V[1:M, 1:N] - mask_v[0:M-1, 1:N] * V[0:M-1, 1:N]) * dy
                - (mask_u[1:M, 1:N] * U[1:M, 1:N] - mask_u[1:M, 0:N-1] * U[1:M, 0:N-1]) * dx
            ) / dA
        # corners
        curl_w[0, 0] = (mask_v[0, 0] * V[0, 0] * dy - mask_u[0, 0] * U[0, 0] * dx) / dA
        curl_w[0, N] = (mask_v[0, N] * V[0, N] * dy - (-mask_u[0, N - 1] * U[0, N - 1]) * dx) / dA
        curl_w[M, 0] = ((-mask_v[M - 1, 0] * V[M - 1, 0]) * dy - (mask_u[M, 0] * U[M, 0]) * dx) / dA
        curl_w[M, N] = ((-mask_v[M - 1, N] * V[M - 1, N]) * dy - (-mask_u[M, N - 1] * U[M, N - 1]) * dx) / dA
        # western boundary (excluding corners)
        if M > 1:
            curl_w[1:M, 0] = (
                (mask_v[1:M, 0] * V[1:M, 0] - mask_v[0:M-1, 0] * V[0:M-1, 0]) * dy
                - (mask_u[1:M, 0] * U[1:M, 0]) * dx
            ) / dA
            curl_w[1:M, N] = (
                (mask_v[1:M, N] * V[1:M, N] - mask_v[0:M-1, N] * V[0:M-1, N]) * dy
                - (-mask_u[1:M, N - 1] * U[1:M, N - 1]) * dx
            ) / dA
        if N > 1:
            curl_w[0, 1:N] = (
                (mask_v[0, 1:N] * V[0, 1:N]) * dy
                - (mask_u[0, 1:N] * U[0, 1:N] - mask_u[0, 0:N-1] * U[0, 0:N-1]) * dx
            ) / dA
            curl_w[M, 1:N] = (
                (-mask_v[M - 1, 1:N] * V[M - 1, 1:N]) * dy
                - (mask_u[M, 1:N] * U[M, 1:N] - mask_u[M, 0:N-1] * U[M, 0:N-1]) * dx
            ) / dA
        # ---------------- shear rate (shearrate_w) -------------------
        shearrate_w = np.zeros_like(curl_w)
        if M > 0 and N > 0:
            shearrate_w[1:M, 1:N] = (
                (V[1:M, 1:N] - V[0:M-1, 1:N]) * dy + (U[1:M, 1:N] - U[1:M, 0:N-1]) * dx
            ) / dA
        shearrate_w[0, 0] = (V[0, 0] * dy + U[0, 0] * dx) / dA
        shearrate_w[0, N] = (V[0, N] * dy + (-U[0, N - 1]) * dx) / dA
        shearrate_w[M, 0] = ((-V[M - 1, 0]) * dy + U[M, 0] * dx) / dA
        shearrate_w[M, N] = ((-V[M - 1, N]) * dy + (-U[M, N - 1]) * dx) / dA
        if M > 1:
            shearrate_w[1:M, 0] = ((V[1:M, 0] - V[0:M-1, 0]) * dy + U[1:M, 0] * dx) / dA
            shearrate_w[1:M, N] = ((V[1:M, N] - V[0:M-1, N]) * dy + (-U[1:M, N - 1]) * dx) / dA
        if N > 1:
            shearrate_w[0, 1:N] = (V[0, 1:N] * dy + (U[0, 1:N] - U[0, 0:N-1]) * dx) / dA
            shearrate_w[M, 1:N] = ((-V[M - 1, 1:N]) * dy + (U[M, 1:N] - U[M, 0:N-1]) * dx) / dA
        # ---------------- stretch & divergence (T-cells) --------------
        strechrate_t = np.zeros((M, N), dtype=u.dtype)
        divergence_t = np.zeros((M, N), dtype=u.dtype)
        if M > 0 and N > 0:
            strechrate_t[:, :] = mask_T * (
                ((U[1:M + 1, :] - U[0:M, :]) * dy - (V[:, 1:N + 1] - V[:, 0:N]) * dx) / dA
            )
            divergence_t[:, :] = mask_T * (
                ((U[1:M + 1, 0:N] - U[0:M, 0:N]) * dy + (V[0:M, 1:N + 1] - V[0:M, 0:N]) * dx) / dA
            )
        # ---------------- quadratic quantities -----------------------
        enstrophy_w = 0.5 * curl_w**2
        sqshearrate_w = 0.5 * shearrate_w**2
        sqshearrate_t = 0.25 * (
            sqshearrate_w[0:M, 0:N]
            + sqshearrate_w[1:M + 1, 0:N]
            + sqshearrate_w[0:M, 1:N + 1]
            + sqshearrate_w[1:M + 1, 1:N + 1]
        )
        sqstrechrate_t = 0.5 * strechrate_t**2 * mask_T
        sqdivergence_t = 0.5 * divergence_t**2 * mask_T
        # Okubo–Weiss on F then interpolate F->T then add/subtract terms
        okuboweiss_w = sqshearrate_w - enstrophy_w
        # interpolate corner to T
        okuboweiss_t = 0.25 * (
            okuboweiss_w[0:M, 0:N]
            + okuboweiss_w[1:M + 1, 0:N]
            + okuboweiss_w[0:M, 1:N + 1]
            + okuboweiss_w[1:M + 1, 1:N + 1]
        )
        okuboweiss_t = okuboweiss_t + sqstrechrate_t - sqdivergence_t
        # Transpose back to native orientation (ny,nx)
        def tb(a: NDArray) -> NDArray:
            return a.T
        intermediates_native = {
            "curl_w": tb(curl_w),
            "shearrate_w": tb(shearrate_w),
            "strechrate_t": tb(strechrate_t),
            "divergence_t": tb(divergence_t),
            "enstrophy_w": tb(enstrophy_w),
            "sqshearrate_w": tb(sqshearrate_w),
            "sqstrechrate_t": tb(sqstrechrate_t),
            "sqdivergence_t": tb(sqdivergence_t),
            "okubo_weiss": tb(okuboweiss_t),
        }
        return intermediates_native

    @staticmethod
    def as_dict(u: NDArray, v: NDArray, grid: Grid) -> Dict[str, NDArray]:
        return {
            "vorticity": FieldDiagnostics.vorticity(u, v, grid),
            "okubo_weiss": FieldDiagnostics.okubo_weiss(u, v, grid),
            "divergence": FieldDiagnostics.divergence(u, v, grid),
            "shear_rate": FieldDiagnostics.shear_rate(u, v, grid),
            "stretch_rate": FieldDiagnostics.stretch_rate(u, v, grid),
        }


__all__ = ["FieldDiagnostics"]
