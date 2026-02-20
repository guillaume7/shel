from __future__ import annotations

from .momentum.strategies import ClosedBC, FlatherBC, FreeSlipBC
from .momentum.strategies import RadiativeSommerfeldBC as RadiativeMomentumBC
from .registry import list_tracer_bcs, register_bc
from .tracer import TracerClosedBC, TracerRadiativeBC
from .waterlevel.strategies import DirichletEtaBC, FlatherEtaBC
from .waterlevel.strategies import RadiativeSommerfeldEtaBC as RadiativeEtaBC

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
