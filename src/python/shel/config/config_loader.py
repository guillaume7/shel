"""
Configuration loading and validation for SHEL.
"""

import logging
import os
from typing import Any, Dict

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

    logger.info("Loading configuration from %s", config_path)

    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    # Validate and normalize configuration
    validate_config(config)
    config = normalize_config(config)

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


def normalize_config(config: Dict[str, Any]) -> Dict[str, Any]:
    """
    Normalize configuration values and structures.

    This bridges differences between the GUI-emitted schema and the
    model's back-end expectations.

    - initial_conditions.type: map 'gaussian' -> 'gaussian_bump'
      and flatten nested gaussian params into top-level keys expected by
      WaterlevelInitialCondition (amplitude, sigma, x_center, y_center).
      If 'custom' is requested (not yet supported), fallback to 'flat' with a warning.

    - bathymetry.type: map aliases 'linear_slope' -> 'sloping',
      'bump' -> 'gaussian_bump', 'custom' -> 'from_file'. Also flatten
      nested parameter blocks into the keys expected by BathymetryInitialCondition.

    Args:
        config: Raw configuration dictionary

    Returns:
        Normalized configuration dictionary (same object mutated for convenience)
    """
    logger = logging.getLogger(__name__)

    # --- Initial conditions (water elevation) ---
    ic = config.get("initial_conditions", {})
    ic_type = str(ic.get("type", "flat")).lower()
    if ic_type == "gaussian":
        ic_type = "gaussian_bump"
        logger.info("Normalizing initial_conditions.type: gaussian -> gaussian_bump")
    elif ic_type == "custom":
        # Not supported yet in water elevation path; fallback to flat to avoid crash
        logger.warning(
            "initial_conditions.type 'custom' not supported; falling back to 'flat'"
        )
        ic_type = "flat"
    ic["type"] = ic_type

    # Flatten gaussian parameters if provided under a nested block from GUI
    gauss_block = ic.get("gaussian")
    if gauss_block and ic_type == "gaussian_bump":
        # Map GUI keys (x0,y0) to model keys (x_center,y_center)
        ic.setdefault("amplitude", gauss_block.get("amplitude"))
        ic.setdefault("sigma", gauss_block.get("sigma"))
        if "x_center" not in ic and "x0" in gauss_block:
            ic["x_center"] = gauss_block.get("x0")
        if "y_center" not in ic and "y0" in gauss_block:
            ic["y_center"] = gauss_block.get("y0")

    config["initial_conditions"] = ic

    # --- Bathymetry ---
    bathy = config.get("bathymetry", {})
    btype = str(bathy.get("type", "flat")).lower()
    alias_map = {
        "linear_slope": "sloping",
        "bump": "gaussian_bump",
        "custom": "from_file",
    }
    if btype in alias_map:
        logger.info("Normalizing bathymetry.type: %s -> %s", btype, alias_map[btype])
        btype = alias_map[btype]
    bathy["type"] = btype

    # Flatten nested parameter blocks into expected keys
    # Flat depth
    if "flat" in bathy and isinstance(bathy["flat"], dict):
        if "depth" in bathy["flat"] and "depth" not in bathy:
            bathy["depth"] = bathy["flat"]["depth"]

    # Linear slope params
    if "linear_slope" in bathy and isinstance(bathy["linear_slope"], dict):
        ls = bathy["linear_slope"]
        if "min_depth" in ls and "depth_min" not in bathy:
            bathy["depth_min"] = ls["min_depth"]
        if "max_depth" in ls and "depth_max" not in bathy:
            bathy["depth_max"] = ls["max_depth"]

    # Bump params
    if "bump" in bathy and isinstance(bathy["bump"], dict):
        bb = bathy["bump"]
        if "amplitude" in bb and "amplitude" not in bathy:
            bathy["amplitude"] = bb["amplitude"]
        if "sigma" in bb and "sigma" not in bathy:
            bathy["sigma"] = bb["sigma"]
        if "x0" in bb and "x_center" not in bathy:
            bathy["x_center"] = bb["x0"]
        if "y0" in bb and "y_center" not in bathy:
            bathy["y_center"] = bb["y0"]
        if "base_depth" in bb and "depth" not in bathy:
            bathy["depth"] = bb["base_depth"]

    # Custom bathy file
    if "custom" in bathy and isinstance(bathy["custom"], dict):
        cb = bathy["custom"]
        if "file" in cb and "file_path" not in bathy:
            bathy["file_path"] = cb["file"]

    config["bathymetry"] = bathy

    return config


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

    logger.info("Configuration saved to %s", config_path)
