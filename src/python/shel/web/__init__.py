"""Web server package for SHEL browser-based UI."""

from shel.web.publisher import WebSocketPublisher
from shel.web.server import (
    app,
    get_connection_manager,
    publish_diag_global,
    publish_progress,
    publish_state_eta,
    publish_state_velocity,
)
from shel.web.solver_manager import SolverManager

__all__ = [
    "app",
    "get_connection_manager",
    "publish_state_eta",
    "publish_state_velocity",
    "publish_diag_global",
    "publish_progress",
    "WebSocketPublisher",
    "SolverManager",
]
