from __future__ import annotations

from .registry import register_bc
from .momentum.strategies import (
    ClosedBC,
    FreeSlipBC,
    RadiativeSommerfeldBC as RadiativeMomentumBC,
    FlatherBC,
)
from .waterlevel.strategies import RadiativeSommerfeldEtaBC as RadiativeEtaBC
from .waterlevel.strategies import FlatherEtaBC


# Register domain-specific strategies under common names
register_bc("closed", momentum=ClosedBC)
register_bc("freeslip", momentum=FreeSlipBC)
register_bc("radiative", momentum=RadiativeMomentumBC, eta=RadiativeEtaBC)
register_bc("flather", momentum=FlatherBC, eta=FlatherEtaBC)

# Re-export canonical names expected by callers
RadiativeSommerfeldBC = RadiativeMomentumBC


