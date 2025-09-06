"""
Initial conditions for water elevation.

This module provides factory classes and implementations for
different types of initial water elevation conditions.
"""

import logging
import math
from typing import Any, Dict

import numpy as np
from numpy.typing import NDArray

logger = logging.getLogger(__name__)


class WaterlevelInitialCondition:
    """Base class for water elevation initial conditions."""

    @staticmethod
    def create(name: str, ny: int, nx: int, params: Dict[str, Any]) -> NDArray:
        """
        Factory method to create initial water elevation field.

        Args:
            name: Type of initial condition
            ny: Number of grid cells in y-direction
            nx: Number of grid cells in x-direction
            params: Parameters for the initial condition

        Returns:
            Initial water elevation field

        Raises:
            ValueError: If the initial condition type is not supported
        """
        if name.lower() == "flat":
            return WaterlevelInitialCondition.flat(ny, nx, params)
        elif name.lower() == "gaussian_bump":
            return WaterlevelInitialCondition.gaussian_bump(ny, nx, params)
        elif name.lower() == "gaussian_depression":
            return WaterlevelInitialCondition.gaussian_depression(ny, nx, params)
        elif name.lower() == "sinusoidal":
            return WaterlevelInitialCondition.sinusoidal(ny, nx, params)
        else:
            raise ValueError(f"Unsupported water elevation initial condition: {name}")

    @staticmethod
    def flat(ny: int, nx: int, params: Dict[str, Any]) -> NDArray:
        """
        Create a flat water surface.

        Args:
            ny: Number of grid cells in y-direction
            nx: Number of grid cells in x-direction
            params: Parameters for the initial condition

        Returns:
            Initial water elevation field (all zeros)
        """
        eta0 = params.get("initial_elevation", 0.0)
        return np.ones((ny, nx)) * eta0

    @staticmethod
    def gaussian_bump(ny: int, nx: int, params: Dict[str, Any]) -> NDArray:
        """
        Create a Gaussian bump in water elevation.

        Args:
            ny: Number of grid cells in y-direction
            nx: Number of grid cells in x-direction
            params: Parameters for the initial condition

        Returns:
            Initial water elevation field with Gaussian bump
        """
        dx = params["grid"]["dx"]
        dy = params["grid"]["dy"]
        x_origin = params["grid"].get("x_origin", 0.0)
        y_origin = params["grid"].get("y_origin", 0.0)

        # Gaussian bump parameters
        amplitude = params["initial_conditions"].get("amplitude", 1.0)
        sigma = params["initial_conditions"].get("sigma", 10.0)
        x_center = params["initial_conditions"].get("x_center", (nx * dx) / 2)
        y_center = params["initial_conditions"].get("y_center", (ny * dy) / 2)

        # Create coordinate arrays
        x = x_origin + np.arange(nx) * dx + dx / 2
        y = y_origin + np.arange(ny) * dy + dy / 2
        X, Y = np.meshgrid(x, y)

        # Calculate distance from center
        distance_squared = (X - x_center) ** 2 + (Y - y_center) ** 2

        # Create Gaussian bump
        eta = amplitude * np.exp(-distance_squared / (2 * sigma**2))

        logger.info(f"Created Gaussian bump: amplitude={amplitude}, sigma={sigma}")
        return eta

    @staticmethod
    def gaussian_depression(ny: int, nx: int, params: Dict[str, Any]) -> NDArray:
        """
        Create a Gaussian depression in water elevation.

        Args:
            ny: Number of grid cells in y-direction
            nx: Number of grid cells in x-direction
            params: Parameters for the initial condition

        Returns:
            Initial water elevation field with Gaussian depression
        """
        # Just negate the Gaussian bump
        bump = WaterlevelInitialCondition.gaussian_bump(ny, nx, params)
        return -bump

    @staticmethod
    def sinusoidal(ny: int, nx: int, params: Dict[str, Any]) -> NDArray:
        """
        Create a sinusoidal water elevation field.

        Args:
            ny: Number of grid cells in y-direction
            nx: Number of grid cells in x-direction
            params: Parameters for the initial condition

        Returns:
            Initial water elevation field with sinusoidal pattern
        """
        dx = params["grid"]["dx"]
        dy = params["grid"]["dy"]
        x_origin = params["grid"].get("x_origin", 0.0)
        y_origin = params["grid"].get("y_origin", 0.0)

        # Sinusoidal parameters
        amplitude = params["initial_conditions"].get("amplitude", 1.0)
        wavelength_x = params["initial_conditions"].get("wavelength_x", nx * dx)
        wavelength_y = params["initial_conditions"].get("wavelength_y", ny * dy)

        # Create coordinate arrays
        x = x_origin + np.arange(nx) * dx + dx / 2
        y = y_origin + np.arange(ny) * dy + dy / 2
        X, Y = np.meshgrid(x, y)

        # Calculate wavenumbers
        k_x = 2 * math.pi / wavelength_x if wavelength_x > 0 else 0
        k_y = 2 * math.pi / wavelength_y if wavelength_y > 0 else 0

        # Create sinusoidal field
        eta = amplitude * np.sin(k_x * X) * np.sin(k_y * Y)

        logger.info(
            f"Created sinusoidal elevation: amplitude={amplitude}, "
            f"wavelength_x={wavelength_x}, wavelength_y={wavelength_y}"
        )
        return eta
