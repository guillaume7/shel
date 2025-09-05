"""Initial conditions package for SHEL.

This package uses self-registration: each concrete initial condition
module inserts its class into a central registry dictionary located in
``factory.py``. For that to work when users import only the factory
symbols (as in tests), we must eagerly import the concrete modules so
their side-effect registration executes.

The explicit imports below ensure that ``from shel.model.initial_conditions.factory
import create_bathymetry`` (which implicitly imports this package's
``__init__`` first) will populate all registries before any factory
creation calls are made.
"""

# Bathymetry ICs
from .bathymetry import bump, step, island, cylinder  # noqa: F401
# Elevation ICs
from .elevation import gaussian as elevation_gaussian, flat  # noqa: F401
# Velocity ICs
from .velocity import geostrophic, solid_body, shear  # noqa: F401
# Tracer ICs
from .tracer import gaussian as tracer_gaussian, uniform  # noqa: F401

__all__ = []  # Registries & factories are exposed via factory.py
