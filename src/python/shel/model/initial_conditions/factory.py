"""Initial condition factory (Phase 1 skeleton).

Provides registry-based creation of category-specific initial condition
strategies. Concrete classes self-register in later phases.
"""
from __future__ import annotations
from typing import Dict, Type, Any
import numpy as np
from .base import BathymetryIC, ElevationIC, VelocityIC, TracerIC

BATHYMETRY_REGISTRY: Dict[str, Type[BathymetryIC]] = {}
ELEVATION_REGISTRY: Dict[str, Type[ElevationIC]] = {}
VELOCITY_REGISTRY: Dict[str, Type[VelocityIC]] = {}
TRACER_REGISTRY: Dict[str, Type[TracerIC]] = {}

def create_bathymetry(name: str, **kwargs) -> BathymetryIC:
    try: return BATHYMETRY_REGISTRY[name.lower()](**kwargs)
    except KeyError as exc: raise ValueError(f"Unknown bathymetry IC '{name}'") from exc

def create_elevation(name: str, **kwargs) -> ElevationIC:
    try: return ELEVATION_REGISTRY[name.lower()](**kwargs)
    except KeyError as exc: raise ValueError(f"Unknown elevation IC '{name}'") from exc

def create_velocity(name: str, **kwargs) -> VelocityIC:
    try: return VELOCITY_REGISTRY[name.lower()](**kwargs)
    except KeyError as exc: raise ValueError(f"Unknown velocity IC '{name}'") from exc

def create_tracer(name: str, **kwargs) -> TracerIC:
    try: return TRACER_REGISTRY[name.lower()](**kwargs)
    except KeyError as exc: raise ValueError(f"Unknown tracer IC '{name}'") from exc

__all__ = [
    "create_bathymetry","create_elevation","create_velocity","create_tracer",
    "BATHYMETRY_REGISTRY","ELEVATION_REGISTRY","VELOCITY_REGISTRY","TRACER_REGISTRY",
    "build_initial_state"
]


def build_initial_state(cfg: Dict[str, Any], grid):
    """Composite initial state builder.

    Parameters
    ----------
    cfg : dict
        Configuration dictionary with optional sections:
          bathymetry: {name: str, params: {...}}
          elevation: {name: str, params: {...}}
          velocity: {name: str, params: {...}}  (may depend on elevation, coriolis)
          tracers: [ {name: str, key: str (optional), params: {...}}, ... ]
          coriolis: {type: 'constant', value: float} (future: beta-plane / map)
    grid : Grid
        Grid instance.

    Returns
    -------
    dict
        Keys: H (total depth), h (bathymetry), eta (free surface), u, v, coriolis,
        plus tracer_<i> or tracer:<key> entries for tracers.
    """
    state: Dict[str, Any] = {}

    # Bathymetry (static depth below datum)
    bathy_cfg = cfg.get("bathymetry", {"name": "bump", "params": {}})
    b_name = bathy_cfg.get("name", "bump")
    b_params = bathy_cfg.get("params", {})
    h = create_bathymetry(b_name, **b_params).build(grid)
    state["h"] = h

    # Elevation (free-surface displacement)
    elev_cfg = cfg.get("elevation", {"name": "flat", "params": {}})
    e_name = elev_cfg.get("name", "flat")
    e_params = elev_cfg.get("params", {})
    eta = create_elevation(e_name, **e_params).build(grid)
    state["eta"] = eta

    # Total depth H = h + eta (assuming h already positive depth; if h represents bathymetric depth)
    H = h + eta
    state["H"] = H

    # Coriolis parameter (allow constant for now)
    coriolis_cfg = cfg.get("coriolis", {"type": "constant", "value": 1e-4})
    if coriolis_cfg.get("type", "constant") == "constant":
        f_val = float(coriolis_cfg.get("value", 1e-4))
        coriolis = np.full_like(eta, f_val, dtype=float)
    else:
        raise ValueError("Unsupported coriolis type")
    state["coriolis"] = coriolis

    # Velocity (may need eta & coriolis)
    vel_cfg = cfg.get("velocity", None)
    if vel_cfg is not None:
        v_name = vel_cfg.get("name", "solid_body")
        v_params = vel_cfg.get("params", {})
        # Build with dependent fields if requested
        if v_name.lower() == "geostrophic":
            u, v = create_velocity(v_name, **v_params).build(grid, eta=eta, coriolis=coriolis)
        else:
            u, v = create_velocity(v_name, **v_params).build(grid)
    else:
        # default zero velocity
        ny, nx = eta.shape
        u = np.zeros((ny, nx + 1))
        v = np.zeros((ny + 1, nx))
    state["u"], state["v"] = u, v

    # Tracers list
    tracer_list = cfg.get("tracers", [])
    for idx, t_spec in enumerate(tracer_list):
        t_name = t_spec.get("name")
        if t_name is None:
            raise ValueError("Tracer spec missing 'name'")
        t_params = t_spec.get("params", {})
        field = create_tracer(t_name, **t_params).build(grid)
        key = t_spec.get("key") or f"tracer_{idx}"
        state[key] = field

    # Basic validation
    ny, nx = eta.shape
    assert h.shape == (ny, nx)
    assert H.shape == (ny, nx)
    assert u.shape == (ny, nx + 1)
    assert v.shape == (ny + 1, nx)
    # Conservation check placeholder: mass = sum(H) * dx*dy (caller can compute)
    return state
