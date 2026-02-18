"""
WebSocket publisher for SHEL solver state updates.

Replaces ZeroMQ PUB/SUB with WebSocket for browser-based visualization.
"""

import asyncio
import base64
import logging
from typing import Any, Dict, Optional

import numpy as np

logger = logging.getLogger(__name__)


def encode_array_to_base64(arr: np.ndarray) -> str:
    """
    Encode a numpy array to base64 string.

    Args:
        arr: Numpy array to encode

    Returns:
        Base64-encoded string
    """
    return base64.b64encode(arr.tobytes()).decode("utf-8")


def create_array_message(
    arr: np.ndarray, units: Optional[str] = None
) -> Dict[str, Any]:
    """
    Create a message dictionary for a numpy array.

    Args:
        arr: Numpy array
        units: Optional units string

    Returns:
        Dictionary with shape, dtype, data (base64), and optional units
    """
    msg = {
        "shape": list(arr.shape),
        "dtype": str(arr.dtype),
        "data": encode_array_to_base64(arr),
    }
    if units:
        msg["units"] = units
    return msg


class WebSocketPublisher:
    """
    Publisher that sends solver state updates via WebSocket.

    Provides the same interface as ZeroMQ publisher but uses WebSocket
    for browser compatibility.
    """

    def __init__(self, connection_manager):
        """
        Initialize WebSocket publisher.

        Args:
            connection_manager: ConnectionManager instance from web server
        """
        self.manager = connection_manager
        self.loop = None

    def set_event_loop(self, loop: asyncio.AbstractEventLoop):
        """Set the event loop for async operations."""
        self.loop = loop

    def publish_state_eta(self, t: float, eta: np.ndarray):
        """
        Publish water elevation state.

        Args:
            t: Simulation time
            eta: Water elevation array (ny, nx)
        """
        message = {
            "topic": "state.eta",
            "t": t,
            **create_array_message(eta, units="meters"),
        }
        self._publish_async(message)

    def publish_state_velocity(self, t: float, U: np.ndarray, V: np.ndarray):
        """
        Publish velocity state.

        Args:
            t: Simulation time
            U: U-velocity array
            V: V-velocity array
        """
        message = {
            "topic": "state.velocity",
            "t": t,
            "U": create_array_message(U),
            "V": create_array_message(V),
            "units": "m/s",
        }
        self._publish_async(message)

    def publish_diag_global(self, t: float, diagnostics: Dict[str, float]):
        """
        Publish global diagnostics.

        Args:
            t: Simulation time
            diagnostics: Dictionary with energy, volume, enstrophy, etc.
        """
        message = {
            "topic": "diag.global",
            "t": t,
            **diagnostics,
        }
        self._publish_async(message)

    def publish_progress(self, step: int, t: float, message_text: str, percent: float):
        """
        Publish simulation progress.

        Args:
            step: Current time step
            t: Simulation time
            message_text: Status message
            percent: Progress percentage (0-100)
        """
        message = {
            "topic": "event.progress",
            "step": step,
            "t": t,
            "message": message_text,
            "percent": percent,
        }
        self._publish_async(message)

    def _publish_async(self, message: Dict[str, Any]):
        """
        Publish a message asynchronously.

        Args:
            message: Message dictionary to publish
        """
        if self.loop and self.loop.is_running():
            # Schedule coroutine in the event loop
            asyncio.run_coroutine_threadsafe(self.manager.broadcast(message), self.loop)
        else:
            logger.warning("Event loop not set or not running, message not published")

    def close(self):
        """Close the publisher (no-op for WebSocket)."""
        pass
