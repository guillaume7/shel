"""
Leapfrog solver implementation for SHEL.

This module implements the leapfrog time-stepping scheme for
solving the shallow water equations.
"""

import logging
from typing import Dict, Any

import numpy as np
from numpy.typing import NDArray

from shel.model.solvers.base import Solver
from shel.model.state import ModelState

logger = logging.getLogger(__name__)


class LeapfrogSolver(Solver):
    """
    Leapfrog solver for shallow water equations.

    The leapfrog scheme is a second-order accurate time-stepping scheme
    that uses central differences in space and time. It requires three
    time levels (old, current, new) for each variable.
    """

    def __init__(self, config: Dict[str, Any]):
        """
        Initialize the leapfrog solver.

        Args:
            config: Configuration dictionary
        """
        super().__init__(config)

        # Robert-Asselin filter coefficient
        self.ra_filter_coef = config["model"].get("ra_filter_coef", 0.1)

        logger.info(
            f"Leapfrog solver initialized: dt={self.timestep}, "
            f"ra_filter_coef={self.ra_filter_coef}"
        )

    def step(self, state: ModelState) -> None:
        """
        Advance the model state by one time step using the leapfrog scheme.

        Args:
            state: Current model state
        """
        # Store the grid parameters for convenience
        dx = state.grid.dx
        dy = state.grid.dy

        # 1. Update water elevation (eta)
        self._update_elevation(state, dx, dy)

        # 2. Update velocities (u, v)
        self._update_velocities(state, dx, dy)

        # 3. Apply boundary conditions
        self.apply_boundary_conditions(state)

        # 4. Update total water depth (H = eta + d)
        state.H_new = state.eta_new + state.d

        # 5. Apply Robert-Asselin filter to dampen the computational mode
        if self.ra_filter_coef > 0:
            self._apply_robert_asselin_filter(state)

        # 6. Advance time and swap time levels
        state.advance_time()
        state.swap_time_levels()

    def _update_elevation(self, state: ModelState, dx: float, dy: float) -> None:
        """
        Update the water elevation using the continuity equation.

        The continuity equation in flux form is:
        dη/dt = -∇·(H*u)

        Args:
            state: Current model state
            dx: Grid spacing in x-direction
            dy: Grid spacing in y-direction
        """
        # Compute divergence of momentum flux (H*u, H*v)
        div_Hu = np.zeros_like(state.eta)

        # Interior points
        for j in range(state.grid.ny):
            for i in range(state.grid.nx):
                # Interpolate H to u-points and v-points
                H_u_east = 0.5 * (
                    state.H[j, i]
                    + (state.H[j, i + 1] if i < state.grid.nx - 1 else state.H[j, i])
                )
                H_u_west = 0.5 * (
                    state.H[j, i] + (state.H[j, i - 1] if i > 0 else state.H[j, i])
                )
                H_v_north = 0.5 * (
                    state.H[j, i]
                    + (state.H[j + 1, i] if j < state.grid.ny - 1 else state.H[j, i])
                )
                H_v_south = 0.5 * (
                    state.H[j, i] + (state.H[j - 1, i] if j > 0 else state.H[j, i])
                )

                # Compute flux divergence
                div_Hu[j, i] = (
                    H_u_east * state.u[j, i + 1] - H_u_west * state.u[j, i]
                ) / dx + (
                    H_v_north * state.v[j + 1, i] - H_v_south * state.v[j, i]
                ) / dy

        # Update eta using leapfrog scheme
        state.eta_new = state.eta_old - 2 * self.timestep * div_Hu

    def _update_velocities(self, state: ModelState, dx: float, dy: float) -> None:
        """
        Update the velocity fields using the momentum equations.

        The momentum equations are:
        du/dt = -u·∇u - g·∂η/∂x + f·v + ν·∇²u - Cd·|u|·u/H
        dv/dt = -v·∇v - g·∂η/∂y - f·u + ν·∇²v - Cd·|v|·v/H

        Args:
            state: Current model state
            dx: Grid spacing in x-direction
            dy: Grid spacing in y-direction
        """
        # Update u-velocity
        self._update_u_velocity(state, dx, dy)

        # Update v-velocity
        self._update_v_velocity(state, dx, dy)

    def _update_u_velocity(self, state: ModelState, dx: float, dy: float) -> None:
        """
        Update the u-velocity field.

        Args:
            state: Current model state
            dx: Grid spacing in x-direction
            dy: Grid spacing in y-direction
        """
        # Compute pressure gradient term
        pg_u = np.zeros_like(state.u)
        for j in range(state.grid.ny):
            for i in range(1, state.grid.nx):
                pg_u[j, i] = (
                    -self.gravity * (state.eta[j, i] - state.eta[j, i - 1]) / dx
                )

        # Compute advection term
        adv_u = np.zeros_like(state.u)
        for j in range(1, state.grid.ny - 1):
            for i in range(1, state.grid.nx):
                # Compute advection in x-direction
                if state.u[j, i] > 0:
                    adv_x = state.u[j, i] * (state.u[j, i] - state.u[j, i - 1]) / dx
                else:
                    adv_x = state.u[j, i] * (state.u[j, i + 1] - state.u[j, i]) / dx

                # Compute advection in y-direction
                v_avg = 0.25 * (
                    state.v[j, i - 1]
                    + state.v[j, i]
                    + state.v[j + 1, i - 1]
                    + state.v[j + 1, i]
                )
                if v_avg > 0:
                    adv_y = v_avg * (state.u[j, i] - state.u[j - 1, i]) / dy
                else:
                    adv_y = v_avg * (state.u[j + 1, i] - state.u[j, i]) / dy

                adv_u[j, i] = -(adv_x + adv_y)

        # Compute Coriolis term
        cor_u = self.compute_coriolis_term_u(state)

        # Compute viscosity term
        vis_u = self.compute_viscosity_term(state.u, dx, dy)

        # Compute bottom drag term
        drag_u = self.compute_bottom_drag_u(state)

        # Combine all terms and update u using leapfrog scheme
        du_dt = pg_u + adv_u + cor_u + vis_u + drag_u
        state.u_new = state.u_old + 2 * self.timestep * du_dt

    def _update_v_velocity(self, state: ModelState, dx: float, dy: float) -> None:
        """
        Update the v-velocity field.

        Args:
            state: Current model state
            dx: Grid spacing in x-direction
            dy: Grid spacing in y-direction
        """
        # Compute pressure gradient term
        pg_v = np.zeros_like(state.v)
        for j in range(1, state.grid.ny):
            for i in range(state.grid.nx):
                pg_v[j, i] = (
                    -self.gravity * (state.eta[j, i] - state.eta[j - 1, i]) / dy
                )

        # Compute advection term
        adv_v = np.zeros_like(state.v)
        for j in range(1, state.grid.ny):
            for i in range(1, state.grid.nx - 1):
                # Compute advection in x-direction
                u_avg = 0.25 * (
                    state.u[j - 1, i]
                    + state.u[j - 1, i + 1]
                    + state.u[j, i]
                    + state.u[j, i + 1]
                )
                if u_avg > 0:
                    adv_x = u_avg * (state.v[j, i] - state.v[j, i - 1]) / dx
                else:
                    adv_x = u_avg * (state.v[j, i + 1] - state.v[j, i]) / dx

                # Compute advection in y-direction
                if state.v[j, i] > 0:
                    adv_y = state.v[j, i] * (state.v[j, i] - state.v[j - 1, i]) / dy
                else:
                    adv_y = state.v[j, i] * (state.v[j + 1, i] - state.v[j, i]) / dy

                adv_v[j, i] = -(adv_x + adv_y)

        # Compute Coriolis term
        cor_v = self.compute_coriolis_term_v(state)

        # Compute viscosity term
        vis_v = self.compute_viscosity_term(state.v, dx, dy)

        # Compute bottom drag term
        drag_v = self.compute_bottom_drag_v(state)

        # Combine all terms and update v using leapfrog scheme
        dv_dt = pg_v + adv_v + cor_v + vis_v + drag_v
        state.v_new = state.v_old + 2 * self.timestep * dv_dt

    def _apply_robert_asselin_filter(self, state: ModelState) -> None:
        """
        Apply the Robert-Asselin filter to dampen the computational mode.

        The Robert-Asselin filter is applied to the current time level:
        φ_t = φ_t + 0.5 * ν * (φ_t+1 - 2 * φ_t + φ_t-1)

        Args:
            state: Current model state
        """
        # Apply filter to eta
        state.eta = state.eta + 0.5 * self.ra_filter_coef * (
            state.eta_new - 2 * state.eta + state.eta_old
        )

        # Apply filter to u
        state.u = state.u + 0.5 * self.ra_filter_coef * (
            state.u_new - 2 * state.u + state.u_old
        )

        # Apply filter to v
        state.v = state.v + 0.5 * self.ra_filter_coef * (
            state.v_new - 2 * state.v + state.v_old
        )

        # Update H after filtering eta
        state.H = state.eta + state.d

    def apply_boundary_conditions(self, state: ModelState) -> None:
        """
        Apply boundary conditions to the model state.

        Args:
            state: Current model state
        """
        # Apply closed boundary conditions (no-slip)
        self._apply_closed_boundaries(state)

        # Other boundary conditions can be added here (open, radiative, etc.)

    def _apply_closed_boundaries(self, state: ModelState) -> None:
        """
        Apply closed (no-slip) boundary conditions.

        Args:
            state: Current model state
        """
        # No normal flow at boundaries
        # Western boundary
        state.u_new[:, 0] = 0

        # Eastern boundary
        state.u_new[:, -1] = 0

        # Southern boundary
        state.v_new[0, :] = 0

        # Northern boundary
        state.v_new[-1, :] = 0
