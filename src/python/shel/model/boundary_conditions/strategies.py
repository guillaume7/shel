from __future__ import annotations

from .registry import register_bc, list_tracer_bcs
from .momentum.strategies import (
    ClosedBC,
    FreeSlipBC,
    RadiativeSommerfeldBC as RadiativeMomentumBC,
    FlatherBC,
)
from .waterlevel.strategies import RadiativeSommerfeldEtaBC as RadiativeEtaBC
from .waterlevel.strategies import FlatherEtaBC, DirichletEtaBC
from .tracer import TracerClosedBC, TracerRadiativeBC


# Register domain-specific strategies under common names
register_bc("closed", momentum=ClosedBC)
register_bc("freeslip", momentum=FreeSlipBC)
register_bc("radiative", momentum=RadiativeMomentumBC, eta=RadiativeEtaBC)
register_bc("flather", momentum=FlatherBC, eta=FlatherEtaBC)
register_bc("dirichlet", eta=DirichletEtaBC)
register_bc("closed", tracer=TracerClosedBC)
register_bc("radiative", tracer=TracerRadiativeBC)

# Re-export canonical names expected by callers
RadiativeSommerfeldBC = RadiativeMomentumBC


