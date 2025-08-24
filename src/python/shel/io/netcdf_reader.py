"""
NetCDF I/O utilities for SHEL.

This module provides functions for reading and writing NetCDF files.
"""

import logging
import os
from typing import Dict, Any, List, Optional, Tuple, Union

import numpy as np
import xarray as xr
from numpy.typing import NDArray

logger = logging.getLogger(__name__)


def read_variable(file_path: str, variable_name: str) -> NDArray:
    """
    Read a variable from a NetCDF file.

    Args:
        file_path: Path to the NetCDF file
        variable_name: Name of the variable to read

    Returns:
        Array containing the variable data

    Raises:
        FileNotFoundError: If the file doesn't exist
        KeyError: If the variable doesn't exist in the file
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"NetCDF file not found: {file_path}")

    with xr.open_dataset(file_path) as ds:
        if variable_name not in ds:
            raise KeyError(f"Variable '{variable_name}' not found in file {file_path}")

        data = ds[variable_name].values

    logger.info(f"Read variable '{variable_name}' from {file_path}")
    return data


def read_grid(file_path: str) -> Tuple[NDArray, NDArray, Dict[str, Any]]:
    """
    Read grid coordinates from a NetCDF file.

    Args:
        file_path: Path to the NetCDF file

    Returns:
        Tuple containing x-coordinates, y-coordinates, and grid attributes

    Raises:
        FileNotFoundError: If the file doesn't exist
        KeyError: If required coordinate variables are not found
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"NetCDF file not found: {file_path}")

    with xr.open_dataset(file_path) as ds:
        # Look for common coordinate variable names
        x_vars = ["x", "lon", "longitude", "x_rho", "x_t"]
        y_vars = ["y", "lat", "latitude", "y_rho", "y_t"]

        # Find x-coordinate
        x_var = next((var for var in x_vars if var in ds), None)
        if x_var is None:
            raise KeyError(f"No x-coordinate variable found in {file_path}")

        # Find y-coordinate
        y_var = next((var for var in y_vars if var in ds), None)
        if y_var is None:
            raise KeyError(f"No y-coordinate variable found in {file_path}")

        # Get coordinate values
        x = ds[x_var].values
        y = ds[y_var].values

        # Get grid attributes
        grid_attrs = {
            "dx": float(np.diff(x).mean()) if len(x) > 1 else 1.0,
            "dy": float(np.diff(y).mean()) if len(y) > 1 else 1.0,
            "nx": len(x),
            "ny": len(y),
            "x_origin": float(x.min()),
            "y_origin": float(y.min()),
        }

    logger.info(
        f"Read grid from {file_path}: {grid_attrs['nx']}x{grid_attrs['ny']} cells"
    )
    return x, y, grid_attrs


def read_bathymetry(file_path: str, variable_name: Optional[str] = None) -> NDArray:
    """
    Read bathymetry data from a NetCDF file.

    Args:
        file_path: Path to the NetCDF file
        variable_name: Name of the bathymetry variable (if None, tries common names)

    Returns:
        Bathymetry array

    Raises:
        FileNotFoundError: If the file doesn't exist
        KeyError: If the bathymetry variable cannot be found
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"NetCDF file not found: {file_path}")

    with xr.open_dataset(file_path) as ds:
        # If variable name not provided, try common bathymetry variable names
        if variable_name is None:
            common_names = [
                "h",
                "depth",
                "bathy",
                "bathymetry",
                "topo",
                "topography",
                "d",
            ]
            for name in common_names:
                if name in ds:
                    variable_name = name
                    break

        if variable_name is None or variable_name not in ds:
            raise KeyError(f"Bathymetry variable not found in {file_path}")

        bathy = ds[variable_name].values

    logger.info(f"Read bathymetry from {file_path} using variable '{variable_name}'")
    return bathy


def write_model_state(state: Dict[str, Any], file_path: str) -> None:
    """
    Write model state to a NetCDF file.

    Args:
        state: Model state dictionary
        file_path: Path to save the NetCDF file
    """
    # Create coordinate arrays
    nx = state["grid"]["nx"]
    ny = state["grid"]["ny"]
    dx = state["grid"]["dx"]
    dy = state["grid"]["dy"]
    x_origin = state["grid"]["x_origin"]
    y_origin = state["grid"]["y_origin"]

    x = x_origin + np.arange(nx) * dx + dx / 2
    y = y_origin + np.arange(ny) * dy + dy / 2
    x_u = x_origin + np.arange(nx + 1) * dx
    y_v = y_origin + np.arange(ny + 1) * dy

    # Create dataset
    ds = xr.Dataset(
        data_vars={
            "eta": (["y", "x"], np.array(state["fields"]["eta"])),
            "u": (["y", "x_u"], np.array(state["fields"]["u"])),
            "v": (["y_v", "x"], np.array(state["fields"]["v"])),
            "d": (["y", "x"], np.array(state["fields"]["d"])),
            "H": (["y", "x"], np.array(state["fields"]["H"])),
        },
        coords={
            "x": ("x", x),
            "y": ("y", y),
            "x_u": ("x_u", x_u),
            "y_v": ("y_v", y_v),
        },
        attrs={
            "title": "SHEL model output",
            "time": state["time"],
            "step": state["step"],
            "gravity": state["parameters"]["gravity"],
            "viscosity": state["parameters"]["viscosity"],
            "bottom_drag_coef": state["parameters"]["bottom_drag_coef"],
            "kinetic_energy": state["diagnostics"]["kinetic_energy"],
            "potential_energy": state["diagnostics"]["potential_energy"],
            "total_energy": state["diagnostics"]["total_energy"],
            "volume": state["diagnostics"]["volume"],
        },
    )

    # Save to file
    ds.to_netcdf(file_path)
    logger.info(f"Wrote model state to {file_path}")


def append_timeseries(timeseries_data: Dict[str, Any], file_path: str) -> None:
    """
    Append timeseries data to a NetCDF file.

    Args:
        timeseries_data: Dictionary containing timeseries data
        file_path: Path to the NetCDF file
    """
    # Check if file exists
    if os.path.exists(file_path):
        # Append to existing file
        with xr.open_dataset(file_path) as ds:
            # Get current time dimension length
            time_dim = len(ds.time)

            # Create new dataset with the new data point
            new_ds = xr.Dataset(
                data_vars={
                    var: (["time"], [value])
                    for var, value in timeseries_data.items()
                    if var != "time"
                },
                coords={"time": (["time"], [timeseries_data["time"]])},
            )

            # Combine datasets and write to file
            combined_ds = xr.concat([ds, new_ds], dim="time")
            combined_ds.to_netcdf(file_path + ".tmp")

        # Replace original file with updated file
        os.replace(file_path + ".tmp", file_path)
    else:
        # Create new file
        ds = xr.Dataset(
            data_vars={
                var: (["time"], [value])
                for var, value in timeseries_data.items()
                if var != "time"
            },
            coords={"time": (["time"], [timeseries_data["time"]])},
        )
        ds.to_netcdf(file_path)

    logger.info(f"Appended timeseries data to {file_path}")
