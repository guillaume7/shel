"""
Solver factory for SHEL.

This module provides a factory for creating solver instances.
"""

import logging
from typing import Dict, Any

from shel.model.solvers.base import Solver
from shel.model.solvers.leapfrog import LeapfrogSolver

logger = logging.getLogger(__name__)


class SolverFactory:
    """Factory class for creating solver instances."""

    @staticmethod
    def create(solver_type: str, config: Dict[str, Any]) -> Solver:
        """
        Create a solver instance of the specified type.

        Args:
            solver_type: Type of solver to create
            config: Configuration dictionary

        Returns:
            Solver instance

        Raises:
            ValueError: If the solver type is not supported
        """
        if solver_type.lower() == "leapfrog":
            logger.info("Creating leapfrog solver")
            return LeapfrogSolver(config)
        else:
            raise ValueError(f"Unsupported solver type: {solver_type}")
