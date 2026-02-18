"""
Solver lifecycle manager for SHEL web UI.

Manages the solver execution, configuration, and state updates via WebSocket.
"""

import asyncio
import logging
from enum import Enum
from typing import Any, Dict, Optional

from shel.model.model_runner import ModelRunner
from shel.web.publisher import WebSocketPublisher

logger = logging.getLogger(__name__)


class SolverState(Enum):
    """Solver execution state."""

    IDLE = "idle"
    RUNNING = "running"
    PAUSED = "paused"
    STOPPED = "stopped"


class SolverManager:
    """
    Manages the SHEL solver lifecycle for the web UI.

    Coordinates solver execution, configuration updates, and WebSocket publishing.
    """

    def __init__(
        self, connection_manager, default_config: Optional[Dict[str, Any]] = None
    ):
        """
        Initialize solver manager.

        Args:
            connection_manager: WebSocket connection manager
            default_config: Default solver configuration
        """
        self.connection_manager = connection_manager
        self.publisher = WebSocketPublisher(connection_manager)
        self.config = default_config or self._get_default_config()
        self.state = SolverState.IDLE
        self.runner: Optional[ModelRunner] = None
        self.task: Optional[asyncio.Task] = None
        self.current_step = 0
        self.total_steps = 0
        self.pause_event = asyncio.Event()
        self.pause_event.set()  # Not paused initially

    def _get_default_config(self) -> Dict[str, Any]:
        """Get default solver configuration."""
        return {
            "grid": {
                "nx": 100,
                "ny": 100,
                "dx": 1000.0,
                "dy": 1000.0,
            },
            "model": {
                "timestep": 0.5,
                "num_steps": 2000,
                "output_interval": 10,
                "gravity": 9.81,
                "viscosity": 0.0,
                "coriolis_parameter": 1.01e-4,
                "h0": 10.0,
            },
            "physics": {},
            "initial_conditions": {
                "type": "gaussian_bump",
                "amplitude": 0.1,
                "sigma": 20000.0,
            },
            "boundary_conditions": {
                "waterlevel": {
                    "north": "radiation",
                    "south": "radiation",
                    "east": "radiation",
                    "west": "radiation",
                },
                "momentum": {
                    "north": "radiation",
                    "south": "radiation",
                    "east": "radiation",
                    "west": "radiation",
                },
            },
        }

    def update_config(self, updates: Dict[str, Any]):
        """
        Update solver configuration.

        Args:
            updates: Configuration updates (partial)
        """

        # Deep merge updates into config
        def deep_update(d, u):
            for k, v in u.items():
                if isinstance(v, dict):
                    d[k] = deep_update(d.get(k, {}), v)
                else:
                    d[k] = v
            return d

        deep_update(self.config, updates)
        logger.info("Configuration updated: %s", updates)

    async def start(self):
        """Start the solver."""
        if self.state == SolverState.RUNNING:
            raise RuntimeError("Solver is already running")

        if self.state == SolverState.PAUSED:
            # Resume from pause
            self.pause_event.set()
            self.state = SolverState.RUNNING
            logger.info("Solver resumed")
            return

        # Start new simulation
        self.state = SolverState.RUNNING
        self.current_step = 0
        self.total_steps = self.config["model"]["num_steps"]
        self.pause_event.set()

        # Create and run solver task
        try:
            self.task = asyncio.create_task(self._run_solver())
            logger.info("Solver task created")
        except Exception as e:
            logger.error("Failed to create solver task: %s", e)
            self.state = SolverState.IDLE
            raise
        logger.info("Solver started")

    async def pause(self):
        """Pause the solver."""
        if self.state != SolverState.RUNNING:
            raise RuntimeError("Solver is not running")

        self.state = SolverState.PAUSED
        self.pause_event.clear()
        logger.info("Solver paused")

    async def stop(self):
        """Stop the solver without resetting state."""
        if self.state == SolverState.IDLE or self.state == SolverState.STOPPED:
            return

        self.state = SolverState.STOPPED
        self.pause_event.set()  # Unblock if paused

        if self.task:
            self.task.cancel()
            try:
                await self.task
            except asyncio.CancelledError:
                pass
            self.task = None

        logger.info("Solver stopped (state preserved)")

    async def reset(self):
        """Reset the solver to initial state."""
        logger.info("Resetting solver from state: %s", self.state)
        await self.stop()
        self.runner = None
        self.current_step = 0
        self.total_steps = 0
        self.state = SolverState.IDLE
        self.pause_event.set()  # Ensure new start is not blocked
        logger.info("Solver reset complete")

    async def _run_solver(self):
        """Run the solver in the background."""
        try:
            # Set event loop for publisher
            loop = asyncio.get_running_loop()
            self.publisher.set_event_loop(loop)

            # Create model runner
            logger.info("Creating ModelRunner with config: %s", self.config)
            self.runner = ModelRunner(self.config)
            self.runner._publisher = self.publisher  # Replace ZeroMQ publisher

            # Initialize model
            self.runner.initialize()
            logger.info("Model initialized")

            # Publish initial state
            self._publish_state(0)

            # Time stepping loop
            num_steps = self.config["model"]["num_steps"]
            output_interval = self.config["model"].get("output_interval", 10)

            logger.info("Starting simulation loop for %d steps", num_steps)
            for step in range(1, num_steps + 1):
                # Check for pause
                await self.pause_event.wait()

                # Check for stop
                if self.state == SolverState.STOPPED:
                    logger.info("Simulation loop stopped by user")
                    break

                # Advance one time step
                if step % output_interval == 0 or step == 1 or step == num_steps:
                    logger.debug("Advancing step %d", step)
                self.runner._advance_timestep()

                # Update simulation time
                self.runner.state.time += self.runner.dt

                self.current_step = step

                # Publish state at intervals
                if step % output_interval == 0 or step == num_steps:
                    self._publish_state(step)

            # Simulation complete
            if self.state == SolverState.RUNNING:
                self.state = SolverState.IDLE
                logger.info("Simulation complete")

        except asyncio.CancelledError:
            logger.info("Solver task cancelled")
            raise
        except Exception as e:
            logger.error("Solver error: %s", e, exc_info=True)
            self.state = SolverState.IDLE
            raise

    def _publish_state(self, step: int):
        """
        Publish current solver state.

        Args:
            step: Current time step
        """
        if not self.runner:
            return

        state = self.runner.state
        t = state.time

        if (
            not hasattr(state, "eta")
            or state.eta is None
            or not hasattr(state, "u")
            or state.u is None
            or not hasattr(state, "v")
            or state.v is None
        ):
            logger.warning("Simulation state incomplete, skipping publish")
            return

        # Publish eta
        self.publisher.publish_state_eta(t, state.eta)

        # Publish velocity
        self.publisher.publish_state_velocity(t, state.u, state.v)

        # Publish diagnostics
        try:
            diagnostics = {
                "total_energy": float(state.compute_total_energy()),
                "kinetic_energy": float(state.compute_kinetic_energy()),
                "potential_energy": float(state.compute_potential_energy()),
                "volume": float(state.compute_volume()),
                "enstrophy": (
                    float(state.compute_enstrophy())
                    if hasattr(state, "compute_enstrophy")
                    else 0.0
                ),
            }
            self.publisher.publish_diag_global(t, diagnostics)
        except Exception as e:
            logger.error("Error computing diagnostics: %s", e)

        # Publish progress
        percent = (step / self.total_steps) * 100 if self.total_steps > 0 else 0
        self.publisher.publish_progress(
            step, t, f"Step {step}/{self.total_steps}", percent
        )

    def get_status(self) -> Dict[str, Any]:
        """Get current solver status."""
        try:
            return {
                "running": self.state == SolverState.RUNNING,
                "step": self.current_step,
                "time": (
                    self.runner.state.time
                    if (self.runner and self.runner.state)
                    else 0.0
                ),
                "progress_percent": (
                    (self.current_step / self.total_steps * 100)
                    if self.total_steps > 0
                    else 0.0
                ),
            }
        except Exception as e:
            logger.error("Error in get_status: %s", e)
            return {"running": False, "step": 0, "time": 0.0, "progress_percent": 0.0}
