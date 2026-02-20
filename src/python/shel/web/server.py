"""
WebSocket server for SHEL web UI.

Replaces ZeroMQ PUB/SUB with WebSocket for browser-based real-time visualization.
"""

import json
import logging
from contextlib import asynccontextmanager
from typing import Set

from fastapi import FastAPI, Request, WebSocket, WebSocketDisconnect
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse

from shel.web.api import router as api_router

logger = logging.getLogger(__name__)


class ConnectionManager:
    """Manages WebSocket connections and broadcasts messages to clients."""

    def __init__(self):
        self.active_connections: Set[WebSocket] = set()

    async def connect(self, websocket: WebSocket):
        """Accept and register a new WebSocket connection."""
        await websocket.accept()
        self.active_connections.add(websocket)
        logger.info(
            "Client connected. Total connections: %d", len(self.active_connections)
        )

    def disconnect(self, websocket: WebSocket):
        """Remove a WebSocket connection."""
        self.active_connections.discard(websocket)
        logger.info(
            "Client disconnected. Total connections: %d", len(self.active_connections)
        )

    async def broadcast(self, message: dict):
        """Broadcast a message to all connected clients."""
        if not self.active_connections:
            return

        message_json = json.dumps(message)
        disconnected = set()

        for connection in self.active_connections:
            try:
                await connection.send_text(message_json)
            except Exception as e:
                logger.warning("Failed to send to client: %s", e)
                disconnected.add(connection)

        # Clean up disconnected clients
        for connection in disconnected:
            self.disconnect(connection)


# Global connection manager
manager = ConnectionManager()


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Lifespan context manager for FastAPI app."""
    from shel.web.api import set_solver_manager
    from shel.web.solver_manager import SolverManager

    logger.info("Starting SHEL WebSocket server")

    # Initialize solver manager
    solver_manager = SolverManager(manager)
    set_solver_manager(solver_manager)
    logger.info("Solver manager initialized")

    yield

    # Cleanup
    await solver_manager.stop()
    logger.info("Shutting down SHEL WebSocket server")


# Create FastAPI app
app = FastAPI(
    title="SHEL Web Server",
    description="WebSocket server for SHEL shallow water model visualization",
    version="1.0.0",
    lifespan=lifespan,
)


@app.exception_handler(Exception)
async def global_exception_handler(request: Request, exc: Exception):
    logger.error("Unhandled exception in %s: %s", request.url.path, exc, exc_info=True)
    return JSONResponse(
        status_code=500,
        content={"detail": str(exc)},
    )


# Enable CORS for local development
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # In production, specify exact origins
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Include API router
app.include_router(api_router)


@app.websocket("/ws")
async def websocket_endpoint(websocket: WebSocket):
    """
    WebSocket endpoint for real-time solver updates.

    Clients connect here to receive state and diagnostic updates.
    """
    await manager.connect(websocket)
    try:
        while True:
            # Keep connection alive and handle incoming messages
            data = await websocket.receive_text()
            # Echo back for now (can add client→server commands later)
            logger.debug("Received from client: %s", data)
    except WebSocketDisconnect:
        manager.disconnect(websocket)


@app.get("/api/health")
async def health_check():
    """Health check endpoint."""
    return {
        "status": "ok",
        "connections": len(manager.active_connections),
    }


async def publish_state_eta(t: float, eta: dict):
    """
    Publish water elevation state update.

    Args:
        t: Simulation time
        eta: Dictionary with 'shape', 'dtype', 'data' (base64), 'units'
    """
    message = {
        "topic": "state.eta",
        "t": t,
        **eta,
    }
    await manager.broadcast(message)


async def publish_state_velocity(t: float, U: dict, V: dict):
    """
    Publish velocity state update.

    Args:
        t: Simulation time
        U: Dictionary with 'shape', 'dtype', 'data' (base64)
        V: Dictionary with 'shape', 'dtype', 'data' (base64)
    """
    message = {
        "topic": "state.velocity",
        "t": t,
        "U": U,
        "V": V,
        "units": "m/s",
    }
    await manager.broadcast(message)


async def publish_diag_global(t: float, diagnostics: dict):
    """
    Publish global diagnostics update.

    Args:
        t: Simulation time
        diagnostics: Dictionary with energy, enstrophy, volume, etc.
    """
    message = {
        "topic": "diag.global",
        "t": t,
        **diagnostics,
    }
    await manager.broadcast(message)


async def publish_progress(step: int, t: float, message_text: str, percent: float):
    """
    Publish simulation progress event.

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
    await manager.broadcast(message)


def get_connection_manager() -> ConnectionManager:
    """Get the global connection manager instance."""
    return manager
