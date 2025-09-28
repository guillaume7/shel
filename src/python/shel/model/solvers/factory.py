"""
SolverFactory: Unified interface for time integration schemes.

Provides a factory for instantiating solver step functions or objects
based on a string key and config. Supports leapfrog, explicit Euler (ministep),
and can be extended for other schemes. Returns a callable with .step(state) interface.
"""

from typing import Any, Callable, Dict

from shel.model.solvers.common.ministep import explicit_step
from shel.model.solvers.time.leapfrog import leapfrog_step_with_config


class SolverStepWrapper:
    """Wraps a step function to provide a .step(state) interface."""

    def __init__(self, step_func: Callable, config: Dict[str, Any]):
        self.step_func = step_func
        self.config = config

    def step(self, state):
        # For leapfrog, call with previous/current arrays and config
        if self.step_func is leapfrog_step_with_config:
            eta_np1, U_np1, V_np1, *_ = self.step_func(
                eta_nm1=state.eta_old,
                eta_n=state.eta,
                H=state.H,
                U_nm1=state.u_old,
                U_n=state.u,
                V_nm1=state.v_old,
                V_n=state.v,
                dx=state.grid.dx,
                dy=state.grid.dy,
                dt=state.timestep,
                g=state.gravity,
                nu_visc=state.viscosity,
                r=state.bottom_drag_coef,
                f=getattr(state, "coriolis", None),
                enable_coriolis=(getattr(state, "coriolis", 0.0) != 0.0),
                config=self.config,
            )
            state.eta_old = state.eta.copy()
            state.u_old = state.u.copy()
            state.v_old = state.v.copy()
            state.eta = eta_np1
            state.u = U_np1
            state.v = V_np1
        elif self.step_func is explicit_step:
            eta_np1, U_np1, V_np1 = self.step_func(
                eta=state.eta,
                H=state.H,
                U=state.u,
                V=state.v,
                dx=state.grid.dx,
                dy=state.grid.dy,
                dt=state.timestep,
                g=state.gravity,
                nu=self.config.get("viscosity", 0.0),
                r=self.config.get("bottom_drag_coef", 0.0),
                f=getattr(state, "coriolis", None),
                enable_coriolis=(getattr(state, "coriolis", 0.0) != 0.0),
                bc_type=self.config.get("boundary_conditions", {}).get(
                    "type", "closed"
                ),
            )
            state.eta = eta_np1
            state.u = U_np1
            state.v = V_np1
        else:
            raise NotImplementedError("Unknown step function")


class SolverFactory:
    """Factory for solver step wrappers."""

    @staticmethod
    def create(solver_type: str, config: Dict[str, Any]) -> SolverStepWrapper:
        if solver_type == "leapfrog":
            return SolverStepWrapper(leapfrog_step_with_config, config)
        elif solver_type in ("explicit_euler", "ministep", "euler"):
            return SolverStepWrapper(explicit_step, config)
        else:
            raise ValueError(f"Unknown solver type: {solver_type}")
