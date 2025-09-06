"""Boundary condition package.

Exports legacy ModelState-based API (boundary.py) and new functional strategy
interfaces used by the ministep path.
"""

# Legacy OO API kept temporarily for backward compatibility during refactor.
# Can be removed once all orchestrated paths are migrated.
from .base import BoundaryCondition, EtaBC, MomentumBC
from .momentum.strategies import FlatherBC
from .registry import get_bc, list_bcs, register_bc
from .strategies import ClosedBC, FreeSlipBC, RadiativeSommerfeldBC
from .waterlevel.strategies import RadiativeSommerfeldEtaBC as RadiativeEtaBC

__all__ = [
    # Strategy API
    "BoundaryCondition",
    "MomentumBC",
    "EtaBC",
    "get_bc",
    "register_bc",
    "list_bcs",
    "ClosedBC",
    "FreeSlipBC",
    "RadiativeSommerfeldBC",
    "FlatherBC",
    "RadiativeEtaBC",
]
