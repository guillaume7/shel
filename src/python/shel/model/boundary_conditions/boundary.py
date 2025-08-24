"""
Boundary conditions for the shallow water model.

This module provides classes for implementing different types of boundary conditions.
"""

import logging
from abc import ABC, abstractmethod
from typing import Dict, Any

import numpy as np

from shel.model.state import ModelState

logger = logging.getLogger(__name__)


class BoundaryCondition(ABC):
    """
    Abstract base class for boundary conditions.

    This class defines the interface that all boundary condition implementations
    must follow.
    """

    def __init__(self, config: Dict[str, Any]):
        """
        Initialize the boundary condition.

        Args:
            config: Configuration dictionary
        """
        self.config = config

    @abstractmethod
    def apply(self, state: ModelState) -> None:
        """
        Apply the boundary condition to the model state.

        Args:
            state: Current model state
        """
        pass


class ClosedBoundaryCondition(BoundaryCondition):
    """
    Closed (no-slip) boundary condition.

    This boundary condition enforces no normal flow at the boundary
    (u=0 at east/west boundaries, v=0 at north/south boundaries).
    """

    def apply(self, state: ModelState) -> None:
        """
        Apply closed boundary conditions.

        Args:
            state: Current model state
        """
        # Western boundary
        state.u_new[:, 0] = 0

        # Eastern boundary
        state.u_new[:, -1] = 0

        # Southern boundary
        state.v_new[0, :] = 0

        # Northern boundary
        state.v_new[-1, :] = 0


class FreeslipBoundaryCondition(BoundaryCondition):
    """
    Free-slip boundary condition.

    This boundary condition enforces no normal flow at the boundary,
    but allows tangential flow (no stress at the boundary).
    """

    def apply(self, state: ModelState) -> None:
        """
        Apply free-slip boundary conditions.

        Args:
            state: Current model state
        """
        # Western boundary
        state.u_new[:, 0] = 0
        # Free-slip: zero gradient in tangential velocity
        state.v_new[1:-1, 0] = state.v_new[1:-1, 1]

        # Eastern boundary
        state.u_new[:, -1] = 0
        # Free-slip: zero gradient in tangential velocity
        state.v_new[1:-1, -1] = state.v_new[1:-1, -2]

        # Southern boundary
        state.v_new[0, :] = 0
        # Free-slip: zero gradient in tangential velocity
        state.u_new[0, 1:-1] = state.u_new[1, 1:-1]

        # Northern boundary
        state.v_new[-1, :] = 0
        # Free-slip: zero gradient in tangential velocity
        state.u_new[-1, 1:-1] = state.u_new[-2, 1:-1]


class RadiativeBoundaryCondition(BoundaryCondition):
    """
    Radiative (Sommerfeld) boundary condition.

    This boundary condition allows waves to exit the domain with
    minimal reflection, based on the Sommerfeld radiation condition.
    """

    def __init__(self, config: Dict[str, Any]):
        """
        Initialize the radiative boundary condition.

        Args:
            config: Configuration dictionary
        """
        super().__init__(config)

        # Wave speed for radiation condition
        self.gravity = config["model"].get("gravity", 9.81)
        self.wave_speed = np.sqrt(
            self.gravity * config.get("bathymetry", {}).get("depth", 1000.0)
        )

        # Boundaries to apply radiative condition
        self.boundaries = config.get("boundary_conditions", {}).get(
            "radiative", ["east"]
        )

        logger.info(
            f"Radiative boundary condition initialized: "
            f"wave_speed={self.wave_speed}, boundaries={self.boundaries}"
        )

    def apply(self, state: ModelState) -> None:
        """
        Apply radiative boundary conditions.

        Args:
            state: Current model state
        """
        dt = state.timestep
        dx = state.grid.dx
        dy = state.grid.dy

        # Western boundary
        if "west" in self.boundaries:
            c = self.wave_speed
            r = c * dt / dx

            # Apply Sommerfeld condition to u
            state.u_new[:, 0] = state.u[:, 0] - r * (state.u[:, 1] - state.u[:, 0])

            # Apply Sommerfeld condition to eta
            state.eta_new[:, 0] = state.eta[:, 0] - r * (
                state.eta[:, 1] - state.eta[:, 0]
            )

        # Eastern boundary
        if "east" in self.boundaries:
            c = self.wave_speed
            r = c * dt / dx

            # Apply Sommerfeld condition to u
            state.u_new[:, -1] = state.u[:, -1] - r * (state.u[:, -1] - state.u[:, -2])

            # Apply Sommerfeld condition to eta
            state.eta_new[:, -1] = state.eta[:, -1] - r * (
                state.eta[:, -1] - state.eta[:, -2]
            )

        # Southern boundary
        if "south" in self.boundaries:
            c = self.wave_speed
            r = c * dt / dy

            # Apply Sommerfeld condition to v
            state.v_new[0, :] = state.v[0, :] - r * (state.v[1, :] - state.v[0, :])

            # Apply Sommerfeld condition to eta
            state.eta_new[0, :] = state.eta[0, :] - r * (
                state.eta[1, :] - state.eta[0, :]
            )

        # Northern boundary
        if "north" in self.boundaries:
            c = self.wave_speed
            r = c * dt / dy

            # Apply Sommerfeld condition to v
            state.v_new[-1, :] = state.v[-1, :] - r * (state.v[-1, :] - state.v[-2, :])

            # Apply Sommerfeld condition to eta
            state.eta_new[-1, :] = state.eta[-1, :] - r * (
                state.eta[-1, :] - state.eta[-2, :]
            )


class BoundaryConditionFactory:
    """Factory class for creating boundary condition instances."""

    @staticmethod
    def create(boundary_type: str, config: Dict[str, Any]) -> BoundaryCondition:
        """
        Create a boundary condition instance of the specified type.

        Args:
            boundary_type: Type of boundary condition to create
            config: Configuration dictionary

        Returns:
            BoundaryCondition instance

        Raises:
            ValueError: If the boundary condition type is not supported
        """
        if boundary_type.lower() == "closed":
            logger.info("Creating closed boundary condition")
            return ClosedBoundaryCondition(config)
        elif boundary_type.lower() == "freeslip":
            logger.info("Creating free-slip boundary condition")
            return FreeslipBoundaryCondition(config)
        elif boundary_type.lower() == "radiative":
            logger.info("Creating radiative boundary condition")
            return RadiativeBoundaryCondition(config)
        else:
            raise ValueError(f"Unsupported boundary condition type: {boundary_type}")
