"""
Grid implementation for SHEL (Arakawa C-grid).

This module implements the staggered grid system used in SHEL, which follows
the Arakawa C-grid approach commonly used in ocean modeling.
"""

import logging
from typing import Any, Dict, Tuple

import numpy as np
from numpy.typing import NDArray

logger = logging.getLogger(__name__)


class Grid:
    """
    Arakawa C-grid implementation for shallow water equations.

    The Arakawa C-grid is a staggered grid where different variables are stored at different
    points to improve numerical stability:
    - Elevation (eta) and tracers are stored at cell centers (T-points)
    - U-velocity is stored at the middle of east/west cell faces (U-points)
    - V-velocity is stored at the middle of north/south cell faces (V-points)
    - Vorticity is computed at cell corners (Q-points)

    Attributes:
        nx (int): Number of grid cells in x-direction
        ny (int): Number of grid cells in y-direction
        dx (float): Grid spacing in x-direction (meters)
        dy (float): Grid spacing in y-direction (meters)
        x_origin (float): x-coordinate of the bottom-left corner
        y_origin (float): y-coordinate of the bottom-left corner
    """

    def __init__(self, config: Dict[str, Any]):
        """
        Initialize the grid from configuration.

        Args:
            config: Configuration dictionary containing grid parameters
        """
        # Grid dimensions
        self.nx = config["grid"]["nx"]
        self.ny = config["grid"]["ny"]
        self.dx = config["grid"]["dx"]
        self.dy = config["grid"]["dy"]
        self.x_origin = config["grid"].get("x_origin", 0.0)
        self.y_origin = config["grid"].get("y_origin", 0.0)

        # Initialize grid coordinate arrays
        self._initialize_coordinates()

        # Land mask (0 for water, 1 for land)
        self.mask = np.zeros((self.ny, self.nx), dtype=int)

        logger.info(
            "Grid initialized: %sx%s cells, %sx%s m resolution",
            self.nx,
            self.ny,
            self.dx,
            self.dy,
        )

    def _initialize_coordinates(self) -> None:
        """Initialize all grid coordinate arrays."""
        # T-grid (cell centers) coordinates
        self.x_t, self.y_t = self._create_t_grid()

        # U-grid (east/west cell faces) coordinates
        self.x_u, self.y_u = self._create_u_grid()

        # V-grid (north/south cell faces) coordinates
        self.x_v, self.y_v = self._create_v_grid()

        # Q-grid (cell corners) coordinates
        self.x_q, self.y_q = self._create_q_grid()

    def _create_t_grid(self) -> Tuple[NDArray, NDArray]:
        """
        Create T-grid (cell centers) coordinate arrays.

        Returns:
            Tuple of x and y coordinate arrays for T-grid
        """
        # 1D coordinate arrays
        x = self.x_origin + np.arange(self.nx) * self.dx + self.dx / 2
        y = self.y_origin + np.arange(self.ny) * self.dy + self.dy / 2

        # 2D meshgrid
        x_t, y_t = np.meshgrid(x, y)
        return x_t, y_t

    def _create_u_grid(self) -> Tuple[NDArray, NDArray]:
        """
        Create U-grid (east/west cell faces) coordinate arrays.

        Returns:
            Tuple of x and y coordinate arrays for U-grid
        """
        # 1D coordinate arrays
        x = self.x_origin + np.arange(self.nx + 1) * self.dx
        y = self.y_origin + np.arange(self.ny) * self.dy + self.dy / 2

        # 2D meshgrid
        x_u, y_u = np.meshgrid(x, y)
        return x_u, y_u

    def _create_v_grid(self) -> Tuple[NDArray, NDArray]:
        """
        Create V-grid (north/south cell faces) coordinate arrays.

        Returns:
            Tuple of x and y coordinate arrays for V-grid
        """
        # 1D coordinate arrays
        x = self.x_origin + np.arange(self.nx) * self.dx + self.dx / 2
        y = self.y_origin + np.arange(self.ny + 1) * self.dy

        # 2D meshgrid
        x_v, y_v = np.meshgrid(x, y)
        return x_v, y_v

    def _create_q_grid(self) -> Tuple[NDArray, NDArray]:
        """
        Create Q-grid (cell corners) coordinate arrays.

        Returns:
            Tuple of x and y coordinate arrays for Q-grid
        """
        # 1D coordinate arrays
        x = self.x_origin + np.arange(self.nx + 1) * self.dx
        y = self.y_origin + np.arange(self.ny + 1) * self.dy

        # 2D meshgrid
        x_q, y_q = np.meshgrid(x, y)
        return x_q, y_q

    def set_land_mask(self, mask: NDArray) -> None:
        """
        Set the land mask for the grid.

        Args:
            mask: Land mask array (0 for water, 1 for land)

        Raises:
            ValueError: If mask dimensions don't match grid dimensions
        """
        if mask.shape != (self.ny, self.nx):
            raise ValueError(
                "Mask shape %s doesn't match grid dimensions (%s, %s)",
                mask.shape,
                self.ny,
                self.nx,
            )

        self.mask = mask.astype(int)
        logger.info("Land mask set: %s land cells", np.sum(self.mask))

    def compute_areas(self) -> Tuple[NDArray, NDArray, NDArray]:
        """
        Compute the areas of T, U, and V cells.

        Returns:
            Tuple of T, U, and V cell area arrays
        """
        # T-cell areas (all equal for uniform grid)
        area_t = np.ones((self.ny, self.nx)) * self.dx * self.dy

        # U-cell areas
        area_u = np.ones((self.ny, self.nx + 1)) * self.dy * self.dx

        # V-cell areas
        area_v = np.ones((self.ny + 1, self.nx)) * self.dx * self.dy

        return area_t, area_u, area_v

    def compute_coriolis(self, f0: float, beta: float = 0.0) -> NDArray:
        """
        Compute Coriolis parameter at each grid point.

        Args:
            f0: Coriolis parameter at reference latitude (1/s)
            beta: Beta-plane parameter (variation of f with y) (1/m/s)

        Returns:
            Array of Coriolis parameter at each T-point
        """
        # For f-plane, beta=0 and f is constant
        if beta == 0:
            return np.ones((self.ny, self.nx)) * f0

        # For beta-plane, f varies with y
        y_ref = self.y_origin + self.ny * self.dy / 2  # Middle of domain
        y_rel = self.y_t - y_ref
        return f0 + beta * y_rel
