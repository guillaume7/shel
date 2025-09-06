"""
Initial conditions for bottom bathymetry.

This module provides factory classes and implementations for
different types of bottom bathymetry configurations.
"""

import logging
import math
from typing import Any, Dict

import numpy as np
from numpy.typing import NDArray

logger = logging.getLogger(__name__)


class BathymetryInitialCondition:
    """Base class for bathymetry initial conditions."""

    @staticmethod
    def create(name: str, ny: int, nx: int, params: Dict[str, Any]) -> NDArray:
        """
        Factory method to create initial bathymetry field.

        Args:
            name: Type of initial condition
            ny: Number of grid cells in y-direction
            nx: Number of grid cells in x-direction
            params: Parameters for the initial condition

        Returns:
            Initial bathymetry field

        Raises:
            ValueError: If the initial condition type is not supported
        """
        if name.lower() == "flat":
            return BathymetryInitialCondition.flat(ny, nx, params)
        elif name.lower() == "sloping":
            return BathymetryInitialCondition.sloping(ny, nx, params)
        elif name.lower() == "gaussian_bump":
            return BathymetryInitialCondition.gaussian_bump(ny, nx, params)
        elif name.lower() == "gaussian_depression":
            return BathymetryInitialCondition.gaussian_depression(ny, nx, params)
        elif name.lower() == "channel":
            return BathymetryInitialCondition.channel(ny, nx, params)
        elif name.lower() == "from_file":
            return BathymetryInitialCondition.from_file(ny, nx, params)
        else:
            raise ValueError(f"Unsupported bathymetry initial condition: {name}")

    @staticmethod
    def flat(ny: int, nx: int, params: Dict[str, Any]) -> NDArray:
        """
        Create a flat bottom bathymetry.

        Args:
            ny: Number of grid cells in y-direction
            nx: Number of grid cells in x-direction
            params: Parameters for the initial condition

        Returns:
            Initial bathymetry field (constant depth)
        """
        depth = params.get("bathymetry", {}).get("depth", 1000.0)
        return np.ones((ny, nx)) * depth

    @staticmethod
    def sloping(ny: int, nx: int, params: Dict[str, Any]) -> NDArray:
        """
        Create a sloping bottom bathymetry.

        Args:
            ny: Number of grid cells in y-direction
            nx: Number of grid cells in x-direction
            params: Parameters for the initial condition

        Returns:
            Initial bathymetry field with linear slope
        """
        dx = params["grid"]["dx"]
        dy = params["grid"]["dy"]
        x_origin = params["grid"].get("x_origin", 0.0)
        y_origin = params["grid"].get("y_origin", 0.0)

        # Slope parameters
        depth_min = params.get("bathymetry", {}).get("depth_min", 100.0)
        depth_max = params.get("bathymetry", {}).get("depth_max", 1000.0)
        direction = params.get("bathymetry", {}).get("slope_direction", "x")

        # Create coordinate arrays
        x = x_origin + np.arange(nx) * dx + dx / 2
        y = y_origin + np.arange(ny) * dy + dy / 2
        X, Y = np.meshgrid(x, y)

        # Create bathymetry field with slope
        if direction.lower() == "x":
            # Slope in x-direction
            x_norm = (X - x_origin) / (nx * dx)
            d = depth_min + (depth_max - depth_min) * x_norm
        elif direction.lower() == "y":
            # Slope in y-direction
            y_norm = (Y - y_origin) / (ny * dy)
            d = depth_min + (depth_max - depth_min) * y_norm
        else:
            raise ValueError(f"Unsupported slope direction: {direction}")

        logger.info(
            f"Created sloping bathymetry: depth_min={depth_min}, "
            f"depth_max={depth_max}, direction={direction}"
        )
        return d

    @staticmethod
    def gaussian_bump(ny: int, nx: int, params: Dict[str, Any]) -> NDArray:
        """
        Create a Gaussian bump in bathymetry.

        Args:
            ny: Number of grid cells in y-direction
            nx: Number of grid cells in x-direction
            params: Parameters for the initial condition

        Returns:
            Initial bathymetry field with Gaussian bump
        """
        dx = params["grid"]["dx"]
        dy = params["grid"]["dy"]
        x_origin = params["grid"].get("x_origin", 0.0)
        y_origin = params["grid"].get("y_origin", 0.0)

        # Gaussian bump parameters
        depth = params.get("bathymetry", {}).get("depth", 1000.0)
        amplitude = params.get("bathymetry", {}).get("amplitude", 500.0)
        sigma = params.get("bathymetry", {}).get("sigma", 10.0)
        x_center = params.get("bathymetry", {}).get("x_center", (nx * dx) / 2)
        y_center = params.get("bathymetry", {}).get("y_center", (ny * dy) / 2)

        # Create coordinate arrays
        x = x_origin + np.arange(nx) * dx + dx / 2
        y = y_origin + np.arange(ny) * dy + dy / 2
        X, Y = np.meshgrid(x, y)

        # Calculate distance from center
        distance_squared = (X - x_center) ** 2 + (Y - y_center) ** 2

        # Create flat bathymetry with Gaussian bump
        # Note: bathymetry is positive downward, so the bump is a shallower region
        d = depth - amplitude * np.exp(-distance_squared / (2 * sigma**2))

        # Ensure minimum depth
        d = np.maximum(d, 1.0)  # Minimum depth of 1 meter

        logger.info(
            f"Created Gaussian bump in bathymetry: depth={depth}, "
            f"amplitude={amplitude}, sigma={sigma}"
        )
        return d

    @staticmethod
    def gaussian_depression(ny: int, nx: int, params: Dict[str, Any]) -> NDArray:
        """
        Create a Gaussian depression in bathymetry.

        Args:
            ny: Number of grid cells in y-direction
            nx: Number of grid cells in x-direction
            params: Parameters for the initial condition

        Returns:
            Initial bathymetry field with Gaussian depression
        """
        dx = params["grid"]["dx"]
        dy = params["grid"]["dy"]
        x_origin = params["grid"].get("x_origin", 0.0)
        y_origin = params["grid"].get("y_origin", 0.0)

        # Gaussian depression parameters
        depth = params.get("bathymetry", {}).get("depth", 1000.0)
        amplitude = params.get("bathymetry", {}).get("amplitude", 500.0)
        sigma = params.get("bathymetry", {}).get("sigma", 10.0)
        x_center = params.get("bathymetry", {}).get("x_center", (nx * dx) / 2)
        y_center = params.get("bathymetry", {}).get("y_center", (ny * dy) / 2)

        # Create coordinate arrays
        x = x_origin + np.arange(nx) * dx + dx / 2
        y = y_origin + np.arange(ny) * dy + dy / 2
        X, Y = np.meshgrid(x, y)

        # Calculate distance from center
        distance_squared = (X - x_center) ** 2 + (Y - y_center) ** 2

        # Create flat bathymetry with Gaussian depression
        # Note: bathymetry is positive downward, so the depression is a deeper region
        d = depth + amplitude * np.exp(-distance_squared / (2 * sigma**2))

        logger.info(
            f"Created Gaussian depression in bathymetry: depth={depth}, "
            f"amplitude={amplitude}, sigma={sigma}"
        )
        return d

    @staticmethod
    def channel(ny: int, nx: int, params: Dict[str, Any]) -> NDArray:
        """
        Create a channel bathymetry.

        Args:
            ny: Number of grid cells in y-direction
            nx: Number of grid cells in x-direction
            params: Parameters for the initial condition

        Returns:
            Initial bathymetry field with channel
        """
        dx = params["grid"]["dx"]
        dy = params["grid"]["dy"]
        x_origin = params["grid"].get("x_origin", 0.0)
        y_origin = params["grid"].get("y_origin", 0.0)

        # Channel parameters
        depth_shallow = params.get("bathymetry", {}).get("depth_shallow", 10.0)
        depth_deep = params.get("bathymetry", {}).get("depth_deep", 1000.0)
        width = params.get("bathymetry", {}).get("channel_width", nx * dx / 3)
        direction = params.get("bathymetry", {}).get("channel_direction", "x")
        center = params.get("bathymetry", {}).get("channel_center", None)

        if center is None:
            if direction.lower() == "x":
                center = ny * dy / 2
            else:
                center = nx * dx / 2

        # Create coordinate arrays
        x = x_origin + np.arange(nx) * dx + dx / 2
        y = y_origin + np.arange(ny) * dy + dy / 2
        X, Y = np.meshgrid(x, y)

        # Create bathymetry field with channel
        if direction.lower() == "x":
            # Channel along x-direction
            distance = np.abs(Y - center)
            d = np.where(distance < width / 2, depth_deep, depth_shallow)
        elif direction.lower() == "y":
            # Channel along y-direction
            distance = np.abs(X - center)
            d = np.where(distance < width / 2, depth_deep, depth_shallow)
        else:
            raise ValueError(f"Unsupported channel direction: {direction}")

        logger.info(
            f"Created channel bathymetry: depth_shallow={depth_shallow}, "
            f"depth_deep={depth_deep}, width={width}, direction={direction}"
        )
        return d

    @staticmethod
    def from_file(ny: int, nx: int, params: Dict[str, Any]) -> NDArray:
        """
        Load bathymetry from a file.

        Args:
            ny: Number of grid cells in y-direction
            nx: Number of grid cells in x-direction
            params: Parameters for the initial condition

        Returns:
            Initial bathymetry field loaded from file

        Raises:
            ValueError: If the file format is not supported
            FileNotFoundError: If the file doesn't exist
        """
        import os

        from shel.io import netcdf_reader

        file_path = params.get("bathymetry", {}).get("file_path")
        if file_path is None:
            raise ValueError("No file_path specified for bathymetry")

        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Bathymetry file not found: {file_path}")

        # Determine file type from extension
        _, ext = os.path.splitext(file_path)

        if ext.lower() == ".nc":
            # NetCDF file
            variable_name = params.get("bathymetry", {}).get("variable_name", "depth")
            d = netcdf_reader.read_variable(file_path, variable_name)
        else:
            raise ValueError(f"Unsupported bathymetry file format: {ext}")

        # Check dimensions
        if d.shape != (ny, nx):
            logger.warning(
                f"Bathymetry dimensions {d.shape} don't match grid dimensions ({ny}, {nx})"
            )
            # Resize if necessary
            from scipy.interpolate import griddata

            # Create target grid
            x_target = np.linspace(0, 1, nx)
            y_target = np.linspace(0, 1, ny)
            X_target, Y_target = np.meshgrid(x_target, y_target)

            # Create source grid
            x_source = np.linspace(0, 1, d.shape[1])
            y_source = np.linspace(0, 1, d.shape[0])
            X_source, Y_source = np.meshgrid(x_source, y_source)

            # Flatten source grid and values
            points = np.column_stack((X_source.flatten(), Y_source.flatten()))
            values = d.flatten()

            # Interpolate to target grid
            d_resized = griddata(points, values, (X_target, Y_target), method="linear")

            logger.info(f"Resized bathymetry from {d.shape} to {d_resized.shape}")
            d = d_resized

        logger.info(f"Loaded bathymetry from file: {file_path}")
        return d
