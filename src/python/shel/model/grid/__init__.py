"""Grid package (Phase 2 in progress).

Exports:
	Grid : core geometry & coordinate arrays
	build_staggered_masks / apply_noslip_flux_masks : mask utilities
"""
from .core import Grid  # type: ignore
from . import masks as _masks
from .masks import (
	build_staggered_masks,
	apply_noslip_flux_masks,
	build_corner_mask,
	build_all_masks,
	mask_velocities,
)

__all__ = [
	"Grid",
	"build_staggered_masks",
	"apply_noslip_flux_masks",
	"build_corner_mask",
	"build_all_masks",
	"mask_velocities",
]

