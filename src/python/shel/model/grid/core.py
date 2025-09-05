"""Core grid implementation (relocated from grid.py in Phase 1)."""
from __future__ import annotations
import logging
from typing import Dict, Any, Tuple
import numpy as np
from numpy.typing import NDArray

logger = logging.getLogger(__name__)

class Grid:
    """Arakawa C-grid implementation for shallow water equations (unchanged)."""
    def __init__(self, config: Dict[str, Any]):
        self.nx = config["grid"]["nx"]
        self.ny = config["grid"]["ny"]
        self.dx = config["grid"]["dx"]
        self.dy = config["grid"]["dy"]
        self.x_origin = config["grid"].get("x_origin", 0.0)
        self.y_origin = config["grid"].get("y_origin", 0.0)
        self._initialize_coordinates()
        self.mask = np.zeros((self.ny, self.nx), dtype=int)
        logger.info(
            f"Grid initialized: {self.nx}x{self.ny} cells, {self.dx}x{self.dy} m resolution"
        )

    def _initialize_coordinates(self) -> None:
        self.x_t, self.y_t = self._create_t_grid()
        self.x_u, self.y_u = self._create_u_grid()
        self.x_v, self.y_v = self._create_v_grid()
        self.x_q, self.y_q = self._create_q_grid()

    def _create_t_grid(self) -> Tuple[NDArray, NDArray]:
        x = self.x_origin + np.arange(self.nx) * self.dx + self.dx / 2
        y = self.y_origin + np.arange(self.ny) * self.dy + self.dy / 2
        x_t, y_t = np.meshgrid(x, y)
        return x_t, y_t

    def _create_u_grid(self) -> Tuple[NDArray, NDArray]:
        x = self.x_origin + np.arange(self.nx + 1) * self.dx
        y = self.y_origin + np.arange(self.ny) * self.dy + self.dy / 2
        x_u, y_u = np.meshgrid(x, y)
        return x_u, y_u

    def _create_v_grid(self) -> Tuple[NDArray, NDArray]:
        x = self.x_origin + np.arange(self.nx) * self.dx + self.dx / 2
        y = self.y_origin + np.arange(self.ny + 1) * self.dy
        x_v, y_v = np.meshgrid(x, y)
        return x_v, y_v

    def _create_q_grid(self) -> Tuple[NDArray, NDArray]:
        x = self.x_origin + np.arange(self.nx + 1) * self.dx
        y = self.y_origin + np.arange(self.ny + 1) * self.dy
        x_q, y_q = np.meshgrid(x, y)
        return x_q, y_q

    def set_land_mask(self, mask: NDArray) -> None:
        if mask.shape != (self.ny, self.nx):
            raise ValueError(
                f"Mask shape {mask.shape} doesn't match grid dimensions ({self.ny}, {self.nx})"
            )
        self.mask = mask.astype(int)
        logger.info(f"Land mask set: {np.sum(self.mask)} land cells")

    def compute_areas(self):
        area_t = np.ones((self.ny, self.nx)) * self.dx * self.dy
        area_u = np.ones((self.ny, self.nx + 1)) * self.dy * self.dx
        area_v = np.ones((self.ny + 1, self.nx)) * self.dx * self.dy
        return area_t, area_u, area_v

    def compute_coriolis(self, f0: float, beta: float = 0.0) -> NDArray:
        if beta == 0:
            return np.ones((self.ny, self.nx)) * f0
        y_ref = self.y_origin + self.ny * self.dy / 2
        y_rel = self.y_t - y_ref
        return f0 + beta * y_rel

__all__ = ["Grid"]
