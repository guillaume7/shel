"""Boundary condition package (Phase 1 skeleton).

Defines abstract base & factory registries in `boundary.py` (legacy
implementation). Future refactor will split domain-specific strategies
into subpackages matching prompt blueprint.
"""

from .boundary import (
	BoundaryCondition,
	ClosedBoundaryCondition,
	FreeslipBoundaryCondition,
	RadiativeBoundaryCondition,
	BoundaryConditionFactory,
)

__all__ = [
	"BoundaryCondition",
	"ClosedBoundaryCondition",
	"FreeslipBoundaryCondition",
	"RadiativeBoundaryCondition",
	"BoundaryConditionFactory",
]

