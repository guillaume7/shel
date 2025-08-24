"""
Model state management for SHEL.

This module provides classes for managing the state of the SHEL model,
including water elevation, velocities, tracers, and diagnostics.
"""

import logging
from typing import Dict, Any, List, Optional

import numpy as np
from numpy.typing import NDArray

from shel.model.grid import Grid

logger = logging.getLogger(__name__)


class ModelState:
    """
    Class to manage the state of the shallow water model.

    This class stores all the state variables of the model, including:
    - Water elevation (eta)
    - Velocities (u, v)
    - Tracers
    - Total water depth (H = eta + d)
    - Bathymetry (d)

    It also provides methods to compute derived quantities like:
    - Vorticity
    - Kinetic energy
    - Potential energy
    - Total energy

    Attributes:
        grid (Grid): The grid on which the model is defined
        time (float): Current model time (seconds)
        timestep (float): Model timestep (seconds)
        step (int): Current step number
    """

    def __init__(self, config: Dict[str, Any]):
        """
        Initialize the model state from configuration.

        Args:
            config: Configuration dictionary
        """
        # Set up the grid
        self.grid = Grid(config)

        # Current time info
        self.time = 0.0
        self.timestep = config["model"]["timestep"]
        self.step = 0

        # Physical parameters
        self.gravity = config["model"].get("gravity", 9.81)  # m/s^2
        self.viscosity = config["model"].get("viscosity", 0.0)  # m^2/s
        self.bottom_drag_coef = config["model"].get(
            "bottom_drag_coef", 0.0
        )  # dimensionless

        # Coriolis parameter
        f0 = config["model"].get("coriolis_parameter", 0.0)
        beta = config["model"].get("beta", 0.0)
        self.coriolis = self.grid.compute_coriolis(f0, beta)

        # Initialize fields
        self._initialize_fields()

        logger.info("Model state initialized")

    def _initialize_fields(self) -> None:
        """Initialize all model state fields with zeros."""
        nx, ny = self.grid.nx, self.grid.ny

        # Bathymetry (bottom depth, positive downward)
        self.d = np.zeros((ny, nx))

        # Water elevation (positive upward)
        self.eta = np.zeros((ny, nx))
        self.eta_old = np.zeros((ny, nx))
        self.eta_new = np.zeros((ny, nx))

        # Total water depth (H = eta + d)
        self.H = np.zeros((ny, nx))
        self.H_old = np.zeros((ny, nx))
        self.H_new = np.zeros((ny, nx))

        # Velocities
        # U-velocity at east/west cell faces
        self.u = np.zeros((ny, nx + 1))
        self.u_old = np.zeros((ny, nx + 1))
        self.u_new = np.zeros((ny, nx + 1))

        # V-velocity at north/south cell faces
        self.v = np.zeros((ny + 1, nx))
        self.v_old = np.zeros((ny + 1, nx))
        self.v_new = np.zeros((ny + 1, nx))

        # Tracers
        self.tracers: Dict[str, NDArray] = {}

    def set_bathymetry(self, d: NDArray) -> None:
        """
        Set the bathymetry (bottom depth).

        Args:
            d: Bathymetry array (positive downward)

        Raises:
            ValueError: If dimensions don't match grid dimensions
        """
        if d.shape != (self.grid.ny, self.grid.nx):
            raise ValueError(
                f"Bathymetry shape {d.shape} doesn't match grid dimensions "
                f"({self.grid.ny}, {self.grid.nx})"
            )

        # Ensure minimum depth
        min_depth = 0.1  # minimum depth to prevent division by zero
        self.d = np.maximum(d, min_depth)

        # Update total depth
        self.H = self.eta + self.d
        self.H_old = self.eta_old + self.d
        self.H_new = self.eta_new + self.d

        logger.info(f"Bathymetry set: min={self.d.min():.2f}, max={self.d.max():.2f}")

    def set_initial_elevation(self, eta: NDArray) -> None:
        """
        Set the initial water elevation.

        Args:
            eta: Water elevation array (positive upward)

        Raises:
            ValueError: If dimensions don't match grid dimensions
        """
        if eta.shape != (self.grid.ny, self.grid.nx):
            raise ValueError(
                f"Elevation shape {eta.shape} doesn't match grid dimensions "
                f"({self.grid.ny}, {self.grid.nx})"
            )

        self.eta = eta.copy()
        self.eta_old = eta.copy()
        self.eta_new = eta.copy()

        # Update total depth
        self.H = self.eta + self.d
        self.H_old = self.eta_old + self.d
        self.H_new = self.eta_new + self.d

        logger.info(
            f"Initial elevation set: min={self.eta.min():.2f}, max={self.eta.max():.2f}"
        )

    def set_initial_velocities(
        self, u: Optional[NDArray] = None, v: Optional[NDArray] = None
    ) -> None:
        """
        Set the initial velocity fields.

        Args:
            u: U-velocity array (east-west) at U-points
            v: V-velocity array (north-south) at V-points

        Raises:
            ValueError: If dimensions don't match grid dimensions
        """
        if u is not None:
            if u.shape != (self.grid.ny, self.grid.nx + 1):
                raise ValueError(
                    f"U-velocity shape {u.shape} doesn't match required dimensions "
                    f"({self.grid.ny}, {self.grid.nx + 1})"
                )

            self.u = u.copy()
            self.u_old = u.copy()
            self.u_new = u.copy()
            logger.info(
                f"Initial U-velocity set: min={self.u.min():.2f}, max={self.u.max():.2f}"
            )

        if v is not None:
            if v.shape != (self.grid.ny + 1, self.grid.nx):
                raise ValueError(
                    f"V-velocity shape {v.shape} doesn't match required dimensions "
                    f"({self.grid.ny + 1}, {self.grid.nx})"
                )

            self.v = v.copy()
            self.v_old = v.copy()
            self.v_new = v.copy()
            logger.info(
                f"Initial V-velocity set: min={self.v.min():.2f}, max={self.v.max():.2f}"
            )

    def add_tracer(self, name: str, initial_value: NDArray) -> None:
        """
        Add a tracer field to the model.

        Args:
            name: Name of the tracer
            initial_value: Initial values of the tracer

        Raises:
            ValueError: If dimensions don't match grid dimensions
        """
        if initial_value.shape != (self.grid.ny, self.grid.nx):
            raise ValueError(
                f"Tracer shape {initial_value.shape} doesn't match grid dimensions "
                f"({self.grid.ny}, {self.grid.nx})"
            )

        self.tracers[name] = initial_value.copy()
        logger.info(
            f"Tracer {name} added: min={initial_value.min():.2f}, max={initial_value.max():.2f}"
        )

    def compute_vorticity(self) -> NDArray:
        """
        Compute the vorticity field (curl of velocity).

        Returns:
            Vorticity field at Q-points
        """
        nx, ny = self.grid.nx, self.grid.ny
        dx, dy = self.grid.dx, self.grid.dy

        # Initialize vorticity array (at Q-points, which are cell corners)
        vorticity = np.zeros((ny + 1, nx + 1))

        # Interior points
        for j in range(1, ny):
            for i in range(1, nx):
                # Compute curl(u,v) = dv/dx - du/dy at each cell corner
                dv_dx = (self.v[j, i] - self.v[j, i - 1]) / dx
                du_dy = (self.u[j, i] - self.u[j - 1, i]) / dy
                vorticity[j, i] = dv_dx - du_dy

        return vorticity

    def compute_kinetic_energy(self) -> float:
        """
        Compute the total kinetic energy of the system.

        Returns:
            Total kinetic energy (J)
        """
        # Compute cell areas
        area_t, area_u, area_v = self.grid.compute_areas()

        # Interpolate H to U and V points
        H_u = np.zeros_like(self.u)
        for i in range(self.grid.nx + 1):
            if i == 0:
                H_u[:, i] = self.H[:, 0]
            elif i == self.grid.nx:
                H_u[:, i] = self.H[:, self.grid.nx - 1]
            else:
                H_u[:, i] = 0.5 * (self.H[:, i - 1] + self.H[:, i])

        H_v = np.zeros_like(self.v)
        for j in range(self.grid.ny + 1):
            if j == 0:
                H_v[j, :] = self.H[0, :]
            elif j == self.grid.ny:
                H_v[j, :] = self.H[self.grid.ny - 1, :]
            else:
                H_v[j, :] = 0.5 * (self.H[j - 1, :] + self.H[j, :])

        # Compute KE = 0.5 * ρ * ∫ H * (u² + v²) dA
        # For simplicity, we use ρ = 1000 kg/m³
        rho = 1000.0

        # U-contribution (excluding boundary points)
        ke_u = 0.5 * rho * np.sum(H_u[:, 1:-1] * self.u[:, 1:-1] ** 2 * area_u[:, 1:-1])

        # V-contribution (excluding boundary points)
        ke_v = 0.5 * rho * np.sum(H_v[1:-1, :] * self.v[1:-1, :] ** 2 * area_v[1:-1, :])

        return ke_u + ke_v

    def compute_potential_energy(self) -> float:
        """
        Compute the total potential energy of the system.

        Returns:
            Total potential energy (J)
        """
        # Compute cell areas
        area_t, _, _ = self.grid.compute_areas()

        # Compute PE = 0.5 * ρ * g * ∫ η² dA
        # For simplicity, we use ρ = 1000 kg/m³
        rho = 1000.0

        # PE contribution (excluding boundary points if they are not physical)
        pe = 0.5 * rho * self.gravity * np.sum(self.eta**2 * area_t)

        return pe

    def compute_total_energy(self) -> float:
        """
        Compute the total energy of the system (kinetic + potential).

        Returns:
            Total energy (J)
        """
        return self.compute_kinetic_energy() + self.compute_potential_energy()

    def compute_volume(self) -> float:
        """
        Compute the total water volume in the domain.

        Returns:
            Total water volume (m³)
        """
        # Compute cell areas
        area_t, _, _ = self.grid.compute_areas()

        # Compute V = ∫ H dA
        volume = np.sum(self.H * area_t)

        return volume

    def compute_okubo_weiss(self) -> NDArray:
        """
        Compute the Okubo-Weiss parameter.

        The Okubo-Weiss parameter is used to identify vortices in the flow.
        It is defined as Q = s_n² + s_s² - ω², where s_n is the normal strain,
        s_s is the shear strain, and ω is the vorticity.

        Returns:
            Okubo-Weiss parameter field at T-points
        """
        nx, ny = self.grid.nx, self.grid.ny
        dx, dy = self.grid.dx, self.grid.dy

        # Initialize arrays
        okubo_weiss = np.zeros((ny, nx))

        # Compute velocity derivatives at T-points
        du_dx = np.zeros((ny, nx))
        du_dy = np.zeros((ny, nx))
        dv_dx = np.zeros((ny, nx))
        dv_dy = np.zeros((ny, nx))

        # du/dx at T-points
        for i in range(nx):
            du_dx[:, i] = (self.u[:, i + 1] - self.u[:, i]) / dx

        # dv/dy at T-points
        for j in range(ny):
            dv_dy[j, :] = (self.v[j + 1, :] - self.v[j, :]) / dy

        # du/dy at T-points
        for j in range(ny):
            if j == 0:
                du_dy[j, :] = (self.u[j + 1, 1:-1] - self.u[j, 1:-1]) / dy
            elif j == ny - 1:
                du_dy[j, :] = (self.u[j, 1:-1] - self.u[j - 1, 1:-1]) / dy
            else:
                du_dy[j, :] = (self.u[j + 1, 1:-1] - self.u[j - 1, 1:-1]) / (2 * dy)

        # dv/dx at T-points
        for i in range(nx):
            if i == 0:
                dv_dx[:, i] = (self.v[1:-1, i + 1] - self.v[1:-1, i]) / dx
            elif i == nx - 1:
                dv_dx[:, i] = (self.v[1:-1, i] - self.v[1:-1, i - 1]) / dx
            else:
                dv_dx[:, i] = (self.v[1:-1, i + 1] - self.v[1:-1, i - 1]) / (2 * dx)

        # Compute strain and vorticity
        normal_strain = du_dx - dv_dy
        shear_strain = du_dy + dv_dx
        vorticity = dv_dx - du_dy

        # Compute Okubo-Weiss parameter
        okubo_weiss = normal_strain**2 + shear_strain**2 - vorticity**2

        return okubo_weiss

    def advance_time(self) -> None:
        """Advance the model time by one timestep."""
        self.time += self.timestep
        self.step += 1

    def swap_time_levels(self) -> None:
        """
        Swap time levels for leap-frog scheme.

        For a leap-frog time stepping scheme, we need to keep track of
        variables at multiple time levels. This method updates the time
        levels after a time step.
        """
        # Swap eta time levels
        self.eta_old = self.eta.copy()
        self.eta = self.eta_new.copy()

        # Swap H time levels
        self.H_old = self.H.copy()
        self.H = self.H_new.copy()

        # Swap u time levels
        self.u_old = self.u.copy()
        self.u = self.u_new.copy()

        # Swap v time levels
        self.v_old = self.v.copy()
        self.v = self.v_new.copy()

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert the model state to a dictionary for serialization.

        Returns:
            Dictionary representation of the model state
        """
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
            "parameters": {
                "gravity": self.gravity,
                "viscosity": self.viscosity,
                "bottom_drag_coef": self.bottom_drag_coef,
            },
            "diagnostics": {
                "kinetic_energy": self.compute_kinetic_energy(),
                "potential_energy": self.compute_potential_energy(),
                "total_energy": self.compute_total_energy(),
                "volume": self.compute_volume(),
            },
        }
