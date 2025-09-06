"""Temporary bridge module wrapping existing ModelState class.

Phase 1 simply relocates the original implementation (formerly
`shel.model.state`) under the state subpackage without *yet* splitting
into core/update/diagnostics components. All behavior remains identical.

Import path compatibility: `from shel.model.state import ModelState`
continues to function because `state.__init__` re-exports this symbol.

Subsequent phases will extract smaller, testable units.
"""

from __future__ import annotations

import logging
from typing import Any, Dict, Optional

import numpy as np
from numpy.typing import NDArray

from shel.model.diagnostics import FieldDiagnostics, IntegratedDiagnostics
from shel.model.grid import Grid, build_all_masks

# Re-import everything from the legacy module (which we rename shortly)
# During Phase 1 we copy the code directly below (no logic changes).


logger = logging.getLogger(__name__)


class ModelState:
    """(Unmodified) model state container.

    Original docstring shortened for brevity; full description remains
    in legacy commit history. No behavioral modifications in Phase 1.
    """

    def __init__(self, config: Dict[str, Any]):
        self.grid = Grid(config)
        self.time = 0.0
        self.timestep = config["model"]["timestep"]
        self.step = 0
        self.gravity = config["model"].get("gravity", 9.81)
        self.viscosity = config["model"].get("viscosity", 0.0)
        self.bottom_drag_coef = config["model"].get("bottom_drag_coef", 0.0)
        f0 = config["model"].get("coriolis_parameter", 0.0)
        beta = config["model"].get("beta", 0.0)
        self.coriolis = self.grid.compute_coriolis(f0, beta)
        self._initialize_fields()
        # Staggered masks (initialized to all-water; updated when bathymetry or explicit mask set)
        self.mask_u = None
        self.mask_v = None
        self.mask_q = None
        logger.info("Model state initialized (Phase 2 masks integrated)")

    def _initialize_fields(self) -> None:
        nx, ny = self.grid.nx, self.grid.ny
        self.d = np.zeros((ny, nx))
        self.eta = np.zeros((ny, nx))
        self.eta_old = self.eta.copy()
        self.eta_new = self.eta.copy()
        self.H = np.zeros((ny, nx))
        self.H_old = self.H.copy()
        self.H_new = self.H.copy()
        self.u = np.zeros((ny, nx + 1))
        self.u_old = self.u.copy()
        self.u_new = self.u.copy()
        self.v = np.zeros((ny + 1, nx))
        self.v_old = self.v.copy()
        self.v_new = self.v.copy()
        self.tracers: Dict[str, NDArray] = {}

    def set_bathymetry(self, d: NDArray) -> None:
        if d.shape != (self.grid.ny, self.grid.nx):
            raise ValueError("Bathymetry shape mismatch")
        self.d = np.maximum(d, 0.1)
        self.H = self.eta + self.d
        self.H_old = self.eta_old + self.d
        self.H_new = self.eta_new + self.d
        # Infer land mask from negative depths (MATLAB logic used d<0 as land). Here treat d<=0.1 threshold? Keep explicit: user may set mask separately.
        # If user encodes land via negative (legacy), build land mask.
        land_mask = (self.d > 0).astype(int)
        self._update_staggered_masks(land_mask)

    def set_land_mask(self, mask: NDArray) -> None:
        """Explicitly set land mask (1 water / 0 land) and rebuild staggered masks."""
        if mask.shape != (self.grid.ny, self.grid.nx):
            raise ValueError("Land mask shape mismatch")
        self._update_staggered_masks(mask.astype(int))

    def set_initial_elevation(self, eta: NDArray) -> None:
        if eta.shape != (self.grid.ny, self.grid.nx):
            raise ValueError("Elevation shape mismatch")
        self.eta = eta.copy()
        self.eta_old = eta.copy()
        self.eta_new = eta.copy()
        self.H = self.eta + self.d
        self.H_old = self.eta_old + self.d
        self.H_new = self.eta_new + self.d

    def set_initial_velocities(
        self, u: Optional[NDArray] = None, v: Optional[NDArray] = None
    ) -> None:
        if u is not None:
            if u.shape != (self.grid.ny, self.grid.nx + 1):
                raise ValueError("U shape mismatch")
            self.u = u.copy()
            self.u_old = u.copy()
            self.u_new = u.copy()
        if v is not None:
            if v.shape != (self.grid.ny + 1, self.grid.nx):
                raise ValueError("V shape mismatch")
            self.v = v.copy()
            self.v_old = v.copy()
            self.v_new = v.copy()

    def add_tracer(self, name: str, initial_value: NDArray) -> None:
        if initial_value.shape != (self.grid.ny, self.grid.nx):
            raise ValueError("Tracer shape mismatch")
        self.tracers[name] = initial_value.copy()

    # Diagnostics wrappers
    def compute_vorticity(self) -> NDArray:
        return FieldDiagnostics.vorticity(self.u, self.v, self.grid)

    def compute_kinetic_energy(self) -> float:
        u, v = self._masked_uv()
        return IntegratedDiagnostics.kinetic_energy(u, v, self.H, self.grid)

    def compute_potential_energy(self) -> float:
        return IntegratedDiagnostics.potential_energy(self.eta, self.grid, self.gravity)

    def compute_total_energy(self) -> float:
        return IntegratedDiagnostics.total_energy(
            self.u, self.v, self.eta, self.H, self.grid, self.gravity
        )

    def compute_volume(self) -> float:
        return IntegratedDiagnostics.volume(self.H, self.grid)

    def compute_okubo_weiss(self) -> NDArray:
        return FieldDiagnostics.okubo_weiss(self.u, self.v, self.grid)

    def compute_enstrophy(self) -> float:
        return IntegratedDiagnostics.enstrophy(self.u, self.v, self.grid)

    def compute_potential_enstrophy(self) -> float:
        return IntegratedDiagnostics.potential_enstrophy(
            self.u, self.v, self.H, self.coriolis, self.grid
        )

    # Time management
    def advance_time(self) -> None:
        self.time += self.timestep
        self.step += 1

    def swap_time_levels(self) -> None:
        self.eta_old = self.eta.copy()
        self.eta = self.eta_new.copy()
        self.H_old = self.H.copy()
        self.H = self.H_new.copy()
        self.u_old = self.u.copy()
        self.u = self.u_new.copy()
        self.v_old = self.v.copy()
        self.v = self.v_new.copy()

    # Internal helpers
    def _update_staggered_masks(self, mask_t: NDArray) -> None:
        mu, mv, mq = build_all_masks(mask_t, noslip=True)
        self.mask_t = mask_t
        self.mask_u = mu
        self.mask_v = mv
        self.mask_q = mq

    def _masked_uv(self):
        if self.mask_u is None or self.mask_v is None:
            return self.u, self.v
        return self.u * self.mask_u, self.v * self.mask_v

    def _serialize_masks(self):
        def opt(a):
            return a.tolist() if a is not None else None

        return {
            "t": opt(getattr(self, "mask_t", None)),
            "u": opt(self.mask_u),
            "v": opt(self.mask_v),
            "q": opt(self.mask_q),
        }

    def to_dict(self) -> Dict[str, Any]:
        integ = IntegratedDiagnostics.as_dict(
            self.u, self.v, self.eta, self.H, self.grid, self.gravity, self.coriolis
        )
        fields = FieldDiagnostics.as_dict(self.u, self.v, self.grid)
        return {
            "time": self.time,
            "step": self.step,
            "grid": {
                "nx": self.grid.nx,
                "ny": self.grid.ny,
                "dx": self.grid.dx,
                "dy": self.grid.dy,
                "x_origin": self.grid.x_origin,
                "y_origin": self.grid.y_origin,
            },
            "fields": {
                "eta": self.eta.tolist(),
                "u": self.u.tolist(),
                "v": self.v.tolist(),
                "d": self.d.tolist(),
                "H": self.H.tolist(),
            },
            "masks": self._serialize_masks(),
            "parameters": {
                "gravity": self.gravity,
                "viscosity": self.viscosity,
                "bottom_drag_coef": self.bottom_drag_coef,
            },
            "diagnostics": {
                "integrated": integ,
                "fields": {k: v.tolist() for k, v in fields.items()},
            },
        }


__all__ = ["ModelState"]
