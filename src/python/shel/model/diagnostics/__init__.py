"""Diagnostics package for SHEL model.

Namespaces:
	- :mod:`integrated` with :class:`IntegratedDiagnostics` for scalar
		domain‑integrated quantities
	- :mod:`fields` with :class:`FieldDiagnostics` for pointwise local
		diagnostic fields (vorticity, Okubo–Weiss, etc.)

Legacy monolithic ``Diagnostics`` class has been removed after refactor.
"""

from .integrated import IntegratedDiagnostics, RHO0  # noqa: F401
from .fields import FieldDiagnostics  # noqa: F401
from .pv import potential_vorticity  # noqa: F401
from .time_series import GlobalAccumulator  # noqa: F401

__all__ = [
	"IntegratedDiagnostics",
	"FieldDiagnostics",
	"potential_vorticity",
	"GlobalAccumulator",
	"RHO0",
]

