"""
Schemes bundle for orchestrated solver configuration (Phase 9).
Encapsulates operator selection, BC strategies, sponge, and time stepper.
"""

from __future__ import annotations

from typing import Any, Dict


class Schemes:
    def __init__(self, config: Dict[str, Any]):
        self.config = config
        self.bc_sides = self._resolve_bc_sides()
        self.time_stepper = config.get("time_stepper", "explicit")
        self.solver = config.get("solver", "explicit_step")
        self.sponge = config.get("sponge", None)
        # Add more operator selection as needed

    def _resolve_bc_sides(self) -> Dict[str, str]:
        bc_cfg = self.config.get("boundary_conditions", {})
        sides = {
            "west": "closed",
            "east": "closed",
            "south": "closed",
            "north": "closed",
        }
        for k in sides:
            val = str(bc_cfg.get(k, "closed")).lower()
            if val in ("closed", "freeslip", "radiative", "flather"):
                sides[k] = val
        return sides

    def as_dict(self) -> Dict[str, Any]:
        return {
            "bc_sides": self.bc_sides,
            "time_stepper": self.time_stepper,
            "solver": self.solver,
            "sponge": self.sponge,
        }
