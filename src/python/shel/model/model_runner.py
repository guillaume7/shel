"""
Model runner implementation for SHEL.

This module provides the ModelRunner class which coordinates the execution
of the shallow water model, including initialization, time stepping, output,
and communication with the GUI.
"""

import logging
import os
import time
from typing import Dict, Any, List, Union

import numpy as np
import zmq
from numpy.typing import NDArray

from shel.model.state import ModelState
from shel.model.solvers.factory import SolverFactory
from shel.model.boundary_conditions import get_bc
from shel.model.initial_conditions.waterlevel import WaterlevelInitialCondition
from shel.model.initial_conditions.bottom import BathymetryInitialCondition
from shel.io import netcdf_reader, parquet_reader

logger = logging.getLogger(__name__)


class ModelRunner:
    """
    Class for running the shallow water model.

    This class coordinates the execution of the model, including:
    - Initialization of the model state
    - Setting up solvers and boundary conditions
    - Time stepping
    - Output to files
    - Communication with the GUI via ZeroMQ

    Attributes:
        config (Dict[str, Any]): Configuration dictionary
        state (ModelState): Model state
        solver: Numerical solver
        boundary_conditions (List): List of boundary conditions
        zmq_context (zmq.Context): ZeroMQ context
        zmq_publisher (zmq.Socket): ZeroMQ publisher socket
    """

    def __init__(self, config: Dict[str, Any]):
        """
        Initialize the model runner.

        Args:
            config: Configuration dictionary
        """
        self.config = config

        # Initialize the model state
        self.state = ModelState(config)

        # Set up the solver
        solver_type = config["model"].get("solver", "leapfrog")
        self.solver = SolverFactory.create(solver_type, config)

        # Prepare boundary condition strategies per side using the new registry
        self.bc_sides = self._resolve_bc_sides_from_config(config)
        # Instantiate one strategy per type (stateless strategies can be reused)
        self._bc_momentum_by_type = {}
        self._bc_eta_by_type = {}
        for bct in set(self.bc_sides.values()):
            m_cls, e_cls = get_bc(bct)
            if m_cls is not None:
                self._bc_momentum_by_type[bct] = m_cls()
            if e_cls is not None:
                self._bc_eta_by_type[bct] = e_cls()

        # Set up ZeroMQ for real-time updates
        self.zmq_context = None
        self.zmq_publisher = None
        if config.get("communication", {}).get("enable_zmq", True):
            self._setup_zmq()

        logger.info("Model runner initialized")

    def _setup_zmq(self) -> None:
        """Set up ZeroMQ publisher socket for real-time updates."""
        port = self.config.get("communication", {}).get("zmq_pub_port", 5556)

        self.zmq_context = zmq.Context()
        self.zmq_publisher = self.zmq_context.socket(zmq.PUB)
        self.zmq_publisher.bind(f"tcp://*:{port}")

        logger.info(f"ZeroMQ publisher started on port {port}")

    def initialize(self) -> None:
        """Initialize the model state with initial conditions."""
        # Set up bathymetry
        bathymetry_type = self.config.get("bathymetry", {}).get("type", "flat")
        bathymetry = BathymetryInitialCondition.create(
            bathymetry_type, self.state.grid.ny, self.state.grid.nx, self.config
        )
        self.state.set_bathymetry(bathymetry)

        # Set up initial water elevation
        elevation_type = self.config.get("initial_conditions", {}).get("type", "flat")
        elevation = WaterlevelInitialCondition.create(
            elevation_type, self.state.grid.ny, self.state.grid.nx, self.config
        )
        self.state.set_initial_elevation(elevation)

        # Set up initial velocities (zero by default)
        self.state.set_initial_velocities()

        logger.info("Model initialized with initial conditions")

    def run(self) -> None:
        """Run the model for the specified number of time steps."""
        # Get configuration parameters
        num_steps = self.config["model"]["num_steps"]
        output_interval = self.config["model"].get("output_interval", 10)

        # Set up output directory
        output_dir = self.config.get("output", {}).get("directory", "./output")
        os.makedirs(output_dir, exist_ok=True)

        # Initialize timers
        start_time = time.time()
        last_output_time = start_time

        # Run time steps
        logger.info(f"Starting model run: {num_steps} steps")
        for step in range(num_steps):
            # Advance the model by one time step
            self.solver.step(self.state)

            # Apply boundary conditions per side using strategy layer
            self._apply_boundary_conditions_strategies()

            # Output at specified intervals
            if step % output_interval == 0:
                step_start_time = time.time()

                # Write output files
                if self.config.get("output", {}).get("enabled", True):
                    self._write_output(step, output_dir)

                # Publish state update via ZeroMQ
                if self.zmq_publisher is not None:
                    self._publish_state_update(step)

                # Log progress
                current_time = time.time()
                steps_per_second = output_interval / (current_time - last_output_time)
                last_output_time = current_time

                elapsed = current_time - start_time
                estimated_total = elapsed * num_steps / (step + 1)
                remaining = estimated_total - elapsed

                logger.info(
                    f"Step {step+1}/{num_steps} ({(step+1)/num_steps*100:.1f}%) "
                    f"- {steps_per_second:.1f} steps/s "
                    f"- Elapsed: {elapsed:.1f}s "
                    f"- Remaining: {remaining:.1f}s"
                )

        # Final output
        if self.config.get("output", {}).get("enabled", True):
            self._write_output(num_steps, output_dir, is_final=True)

        # Clean up ZeroMQ
        if self.zmq_publisher is not None and self.zmq_context is not None:
            self._publish_final_update(num_steps)
            self.zmq_publisher.close()
            self.zmq_context.term()

        # Log completion
        total_time = time.time() - start_time
        logger.info(
            f"Model run completed: {num_steps} steps in {total_time:.1f}s "
            f"({num_steps/total_time:.1f} steps/s)"
        )

    @staticmethod
    def _resolve_bc_sides_from_config(config: Dict[str, Any]) -> Dict[str, str]:
        sides = {"west": "closed", "east": "closed", "south": "closed", "north": "closed"}
        bc_cfg = config.get("boundary_conditions", {}) if isinstance(config, dict) else {}
        if not isinstance(bc_cfg, dict):
            return sides
        for k in sides.keys():
            val = str(bc_cfg.get(k, "closed")).lower()
            if val in ("closed", "freeslip", "radiative"):
                sides[k] = val
            else:
                sides[k] = "closed"
        return sides

    def _apply_boundary_conditions_strategies(self) -> None:
        s = self.state
        dt = s.timestep
        dx, dy = s.grid.dx, s.grid.dy
        g = s.gravity
        H = s.H
        # Momentum per-side
        for side, bct in self.bc_sides.items():
            m = self._bc_momentum_by_type.get(bct)
            if m is not None:
                m.apply_side(
                    s.u_new,
                    s.v_new,
                    side,
                    U_old=s.u,
                    V_old=s.v,
                    H=H,
                    g=g,
                    dt=dt,
                    dx=dx,
                    dy=dy,
                )
        # Eta per-side
        for side, bct in self.bc_sides.items():
            e = self._bc_eta_by_type.get(bct)
            if e is not None:
                e.apply_side_eta(
                    s.eta_new,
                    side,
                    eta_old=s.eta,
                    H=H,
                    g=g,
                    dt=dt,
                    dx=dx,
                    dy=dy,
                )

    def _write_output(self, step: int, output_dir: str, is_final: bool = False) -> None:
        """
        Write output files for the current state.

        Args:
            step: Current step number
            output_dir: Directory to write output files to
            is_final: Whether this is the final output
        """
        # Create state dictionary for serialization
        state_dict = self.state.to_dict()

        # Write NetCDF output
        nc_file = os.path.join(output_dir, f"state_{step:06d}.nc")
        netcdf_reader.write_model_state(state_dict, nc_file)

        # Write timeseries data to Parquet file
        timeseries_file = os.path.join(output_dir, "timeseries.parquet")
        timeseries_data = {
            "time": self.state.time,
            "step": step,
            "kinetic_energy": self.state.compute_kinetic_energy(),
            "potential_energy": self.state.compute_potential_energy(),
            "total_energy": self.state.compute_total_energy(),
            "volume": self.state.compute_volume(),
            "max_elevation": float(np.max(self.state.eta)),
            "min_elevation": float(np.min(self.state.eta)),
            "max_velocity": float(
                max(np.max(np.abs(self.state.u)), np.max(np.abs(self.state.v)))
            ),
        }
        parquet_reader.append_timeseries(timeseries_data, timeseries_file)

        logger.debug(f"Output written for step {step}")

    def _publish_state_update(self, step: int) -> None:
        """
        Publish a state update via ZeroMQ.

        Args:
            step: Current step number
        """
        # Create a subset of the state to publish
        # Only send the most important fields to reduce message size
        message = {
            "time": self.state.time,
            "step": step,
            "status": "running",
            "max_elevation": float(np.max(self.state.eta)),
            "min_elevation": float(np.min(self.state.eta)),
            "max_velocity": float(
                max(np.max(np.abs(self.state.u)), np.max(np.abs(self.state.v)))
            ),
            "kinetic_energy": float(self.state.compute_kinetic_energy()),
            "potential_energy": float(self.state.compute_potential_energy()),
            "total_energy": float(self.state.compute_total_energy()),
            "volume": float(self.state.compute_volume()),
            # Send downsampled fields for visualization
            "fields": self._downsample_fields_for_message(),
        }

        if self.zmq_publisher is None:
            logging.warning(
                "ZeroMQ publisher is not set up. Cannot publish state update."
            )
            return

        self.zmq_publisher.send_json(message)

    def _publish_final_update(self, num_steps: int) -> None:
        """
        Publish a final state update via ZeroMQ.

        Args:
            num_steps: Total number of steps
        """
        message = {
            "time": self.state.time,
            "step": num_steps,
            "status": "complete",
            "max_elevation": float(np.max(self.state.eta)),
            "min_elevation": float(np.min(self.state.eta)),
            "max_velocity": float(
                max(np.max(np.abs(self.state.u)), np.max(np.abs(self.state.v)))
            ),
            "kinetic_energy": float(self.state.compute_kinetic_energy()),
            "potential_energy": float(self.state.compute_potential_energy()),
            "total_energy": float(self.state.compute_total_energy()),
            "volume": float(self.state.compute_volume()),
            "fields": self._downsample_fields_for_message(),
        }

        if self.zmq_publisher is None:
            logging.warning(
                "ZeroMQ publisher is not set up. Cannot publish final update."
            )
            return

        self.zmq_publisher.send_json(message)

    def _downsample_fields_for_message(self) -> Dict[str, Union[List[float], int]]:
        """
        Downsample model fields for inclusion in ZeroMQ messages.

        To reduce message size, we downsample the fields to a maximum size.

        Returns:
            Dictionary with downsampled fields
        """
        # Maximum dimensions for message fields
        max_dim = 50

        # Get current dimensions
        ny, nx = self.state.grid.ny, self.state.grid.nx

        # Calculate stride for downsampling
        stride_y = max(1, ny // max_dim)
        stride_x = max(1, nx // max_dim)

        # Downsample fields
        eta_ds = self.state.eta[::stride_y, ::stride_x].tolist()
        u_ds = self.state.u[::stride_y, ::stride_x].tolist()
        v_ds = self.state.v[::stride_y, ::stride_x].tolist()

        return {
            "eta": eta_ds,
            "u": u_ds,
            "v": v_ds,
            "stride_x": stride_x,
            "stride_y": stride_y,
        }


def run(config: Dict[str, Any]) -> None:
    """
    Run the model with the given configuration.

    This is the main entry point for running the model.

    Args:
        config: Configuration dictionary
    """
    # Create and initialize the model runner
    runner = ModelRunner(config)
    runner.initialize()

    # Run the model
    runner.run()
