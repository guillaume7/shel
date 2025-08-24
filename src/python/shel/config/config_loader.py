"""
Configuration loading and validation for SHEL.
"""

import logging
import os
from typing import Dict, Any

import yaml

logger = logging.getLogger(__name__)


def load_config(config_path: str) -> Dict[str, Any]:
    """
    Load and validate a configuration from a YAML file.

    Args:
        config_path: Path to the YAML configuration file

    Returns:
        Configuration dictionary

    Raises:
        FileNotFoundError: If the config file doesn't exist
        ValueError: If the config is invalid
    """
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Configuration file not found: {config_path}")

    logger.info(f"Loading configuration from {config_path}")

    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    # Validate the configuration
    validate_config(config)

    return config


def validate_config(config: Dict[str, Any]) -> None:
    """
    Validate a configuration dictionary.

    Args:
        config: Configuration dictionary

    Raises:
        ValueError: If the configuration is invalid
    """
    # Required top-level sections
    required_sections = ["model", "grid", "initial_conditions", "boundary_conditions"]

    for section in required_sections:
        if section not in config:
            raise ValueError(f"Missing required configuration section: {section}")

    # Validate model section
    required_model_params = ["timestep", "num_steps", "gravity"]
    for param in required_model_params:
        if param not in config["model"]:
            raise ValueError(f"Missing required model parameter: {param}")

    # Validate grid section
    required_grid_params = ["nx", "ny", "dx", "dy"]
    for param in required_grid_params:
        if param not in config["grid"]:
            raise ValueError(f"Missing required grid parameter: {param}")

    # Validate initial conditions section
    if "type" not in config["initial_conditions"]:
        raise ValueError("Missing 'type' in initial_conditions section")

    # Validate boundary conditions section
    required_boundaries = ["north", "south", "east", "west"]
    for boundary in required_boundaries:
        if boundary not in config["boundary_conditions"]:
            raise ValueError(f"Missing {boundary} boundary condition")


def create_default_config() -> Dict[str, Any]:
    """
    Create a default configuration dictionary.

    Returns:
        Default configuration
    """
    return {
        "model": {
            "timestep": 60.0,  # seconds
            "num_steps": 1000,
            "output_interval": 10,
            "coriolis_parameter": 1e-4,  # f-plane approximation
            "gravity": 9.81,  # m/s^2
            "viscosity": 10.0,  # m^2/s
            "bottom_drag_coef": 0.0025,  # dimensionless
        },
        "grid": {
            "nx": 100,
            "ny": 100,
            "dx": 1000.0,  # meters
            "dy": 1000.0,  # meters
            "x_origin": 0.0,
            "y_origin": 0.0,
        },
        "boundary_conditions": {
            "north": "closed",
            "south": "closed",
            "east": "closed",
            "west": "closed",
        },
        "initial_conditions": {
            "type": "gaussian_bump",
            "amplitude": 1.0,  # meters
            "sigma": 10000.0,  # meters
            "x_center": 50000.0,  # meters
            "y_center": 50000.0,  # meters
        },
    }


def save_config(config: Dict[str, Any], config_path: str) -> None:
    """
    Save a configuration to a YAML file.

    Args:
        config: Configuration dictionary
        config_path: Path to save the YAML file
    """
    with open(config_path, "w") as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)

    logger.info(f"Configuration saved to {config_path}")
