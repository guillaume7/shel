"""
Base solver class for SHEL.

This module provides the base class for numerical solvers.
"""

import logging
from abc import ABC, abstractmethod
from typing import Dict, Any

from shel.model.state import ModelState

logger = logging.getLogger(__name__)


class Solver(ABC):
    """
    Abstract base class for numerical solvers.

    This class defines the interface that all numerical solvers must implement.
    Concrete solver implementations should inherit from this class and implement
    the step method.
    """

    def __init__(self, config: Dict[str, Any]):
        """
        Initialize the solver with configuration parameters.

        Args:
            config: Configuration dictionary
        """
        self.config = config
        self.timestep = config["model"]["timestep"]

        # Physical parameters
        self.gravity = config["model"].get("gravity", 9.81)  # m/s^2
        self.viscosity = config["model"].get("viscosity", 0.0)  # m^2/s
        self.bottom_drag_coef = config["model"].get(
            "bottom_drag_coef", 0.0
        )  # dimensionless

    @abstractmethod
    def step(self, state: ModelState) -> None:
        """
        Advance the model state by one time step.

        Args:
            state: Current model state
        """
        pass

    def apply_boundary_conditions(self, state: ModelState) -> None:
        """
        Apply boundary conditions to the model state.

        Args:
            state: Current model state
        """
        # Default implementation does nothing - concrete solvers should override this
        pass

    def compute_viscosity_term(self, field, dx: float, dy: float) -> Any:
        """
        Compute the viscosity term for a field using the Laplacian operator.

        Args:
            field: Field to compute viscosity term for
            dx: Grid spacing in x-direction
            dy: Grid spacing in y-direction

        Returns:
            Viscosity term
        """
        # Default implementation - concrete solvers may override this
        import numpy as np

        # Compute second derivatives
        d2f_dx2 = np.zeros_like(field)
        d2f_dy2 = np.zeros_like(field)

        # Interior points
        d2f_dx2[1:-1, 1:-1] = (
            field[1:-1, 2:] - 2 * field[1:-1, 1:-1] + field[1:-1, :-2]
        ) / (dx**2)
        d2f_dy2[1:-1, 1:-1] = (
            field[2:, 1:-1] - 2 * field[1:-1, 1:-1] + field[:-2, 1:-1]
        ) / (dy**2)

        # Combine to form Laplacian
        return self.viscosity * (d2f_dx2 + d2f_dy2)

    def compute_coriolis_term_u(self, state: ModelState) -> Any:
        """
        Compute the Coriolis term for the u-velocity.

        The Coriolis term for u is: f * v

        Args:
            state: Current model state

        Returns:
            Coriolis term for u-velocity
        """
        import numpy as np

        # Interpolate Coriolis parameter to u-points
        f_u = np.zeros_like(state.u)
        for i in range(state.grid.nx + 1):
            if i == 0:
                f_u[:, i] = state.coriolis[:, 0]
            elif i == state.grid.nx:
                f_u[:, i] = state.coriolis[:, -1]
            else:
                f_u[:, i] = 0.5 * (state.coriolis[:, i - 1] + state.coriolis[:, i])

        # Interpolate v to u-points
        v_u = np.zeros_like(state.u)
        for j in range(state.grid.ny):
            for i in range(state.grid.nx + 1):
                if i == 0:
                    v_avg = 0.25 * (state.v[j, 0] + state.v[j + 1, 0])
                elif i == state.grid.nx:
                    v_avg = 0.25 * (state.v[j, -1] + state.v[j + 1, -1])
                else:
                    v_avg = 0.25 * (
                        state.v[j, i - 1]
                        + state.v[j + 1, i - 1]
                        + state.v[j, i]
                        + state.v[j + 1, i]
                    )
                v_u[j, i] = v_avg

        # Compute Coriolis term
        return f_u * v_u

    def compute_coriolis_term_v(self, state: ModelState) -> Any:
        """
        Compute the Coriolis term for the v-velocity.

        The Coriolis term for v is: -f * u

        Args:
            state: Current model state

        Returns:
            Coriolis term for v-velocity
        """
        import numpy as np

        # Interpolate Coriolis parameter to v-points
        f_v = np.zeros_like(state.v)
        for j in range(state.grid.ny + 1):
            if j == 0:
                f_v[j, :] = state.coriolis[0, :]
            elif j == state.grid.ny:
                f_v[j, :] = state.coriolis[-1, :]
            else:
                f_v[j, :] = 0.5 * (state.coriolis[j - 1, :] + state.coriolis[j, :])

        # Interpolate u to v-points
        u_v = np.zeros_like(state.v)
        for j in range(state.grid.ny + 1):
            for i in range(state.grid.nx):
                if j == 0:
                    u_avg = 0.25 * (state.u[0, i] + state.u[0, i + 1])
                elif j == state.grid.ny:
                    u_avg = 0.25 * (state.u[-1, i] + state.u[-1, i + 1])
                else:
                    u_avg = 0.25 * (
                        state.u[j - 1, i]
                        + state.u[j - 1, i + 1]
                        + state.u[j, i]
                        + state.u[j, i + 1]
                    )
                u_v[j, i] = u_avg

        # Compute Coriolis term (negative because of coordinate system)
        return -f_v * u_v

    def compute_bottom_drag_u(self, state: ModelState) -> Any:
        """
        Compute the bottom drag term for the u-velocity.

        The bottom drag term for u is: -Cd * |u| * u / H

        Args:
            state: Current model state

        Returns:
            Bottom drag term for u-velocity
        """
        import numpy as np

        # Interpolate H to u-points
        H_u = np.zeros_like(state.u)
        for i in range(state.grid.nx + 1):
            if i == 0:
                H_u[:, i] = state.H[:, 0]
            elif i == state.grid.nx:
                H_u[:, i] = state.H[:, -1]
            else:
                H_u[:, i] = 0.5 * (state.H[:, i - 1] + state.H[:, i])

        # Compute velocity magnitude at u-points
        # For simplicity, we just use |u| instead of sqrt(u^2 + v^2)
        u_mag = np.abs(state.u)

        # Compute bottom drag term
        return -self.bottom_drag_coef * u_mag * state.u / np.maximum(H_u, 0.1)

    def compute_bottom_drag_v(self, state: ModelState) -> Any:
        """
        Compute the bottom drag term for the v-velocity.

        The bottom drag term for v is: -Cd * |v| * v / H

        Args:
            state: Current model state

        Returns:
            Bottom drag term for v-velocity
        """
        import numpy as np

        # Interpolate H to v-points
        H_v = np.zeros_like(state.v)
        for j in range(state.grid.ny + 1):
            if j == 0:
                H_v[j, :] = state.H[0, :]
            elif j == state.grid.ny:
                H_v[j, :] = state.H[-1, :]
            else:
                H_v[j, :] = 0.5 * (state.H[j - 1, :] + state.H[j, :])

        # Compute velocity magnitude at v-points
        # For simplicity, we just use |v| instead of sqrt(u^2 + v^2)
        v_mag = np.abs(state.v)

        # Compute bottom drag term
        return -self.bottom_drag_coef * v_mag * state.v / np.maximum(H_v, 0.1)
