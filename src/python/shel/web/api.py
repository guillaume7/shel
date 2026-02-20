"""
REST API endpoints for SHEL web UI.

Provides configuration, control, and status endpoints.
"""

import logging
from typing import Optional

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field

from shel.web.solver_manager import SolverState

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api", tags=["api"])

# Global solver manager (will be set by server on startup)
_solver_manager = None


def set_solver_manager(manager):
    """Set the global solver manager instance."""
    global _solver_manager
    _solver_manager = manager


# --- Request/Response Models ---


class ConfigUpdate(BaseModel):
    """Configuration update request."""

    grid_nx: Optional[int] = Field(None, ge=10, le=1000, description="Grid points in x")
    grid_ny: Optional[int] = Field(None, ge=10, le=1000, description="Grid points in y")
    dt: Optional[float] = Field(None, gt=0, description="Time step (seconds)")
    duration: Optional[float] = Field(
        None, gt=0, description="Simulation duration (seconds)"
    )
    viscosity: Optional[float] = Field(
        None, ge=0, description="Horizontal viscosity (m²/s)"
    )
    coriolis_f: Optional[float] = Field(None, description="Coriolis parameter (rad/s)")


class SolverStatus(BaseModel):
    """Solver status response."""

    running: bool
    step: int
    time: float
    progress_percent: float


# --- Endpoints ---


@router.get("/config")
async def get_config():
    """Get current solver configuration."""
    if not _solver_manager:
        raise HTTPException(status_code=503, detail="Solver manager not initialized")
    return _solver_manager.config


@router.post("/config")
async def update_config(config: ConfigUpdate):
    """Update solver configuration."""
    if not _solver_manager:
        raise HTTPException(status_code=503, detail="Solver manager not initialized")

    # Convert to nested dict structure matching ModelRunner expectations
    updates = {}
    if config.grid_nx is not None:
        updates.setdefault("grid", {})["nx"] = config.grid_nx
    if config.grid_ny is not None:
        updates.setdefault("grid", {})["ny"] = config.grid_ny
    if config.dt is not None:
        updates.setdefault("model", {})["timestep"] = config.dt
    if config.duration is not None and (config.dt or 0.5) > 0:
        dt = config.dt if config.dt is not None else 0.5
        updates.setdefault("model", {})["num_steps"] = int(config.duration / dt)
    if config.viscosity is not None:
        updates.setdefault("model", {})["viscosity"] = config.viscosity
    if config.coriolis_f is not None:
        updates.setdefault("model", {})["coriolis_parameter"] = config.coriolis_f

    _solver_manager.update_config(updates)
    logger.info("Configuration updated: %s", updates)
    return {"status": "ok", "config": _solver_manager.config}


@router.get("/status")
async def get_status() -> SolverStatus:
    """Get current solver status."""
    if not _solver_manager:
        raise HTTPException(status_code=503, detail="Solver manager not initialized")
    return SolverStatus(**_solver_manager.get_status())


@router.post("/control/start")
async def start_simulation():
    """Start the simulation."""
    if not _solver_manager:
        raise HTTPException(status_code=503, detail="Solver manager not initialized")

    try:
        await _solver_manager.start()
        logger.info("Simulation started")
        return {"status": "started"}
    except RuntimeError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.post("/control/pause")
async def pause_simulation():
    """Pause the simulation."""
    if not _solver_manager:
        raise HTTPException(status_code=503, detail="Solver manager not initialized")

    try:
        await _solver_manager.pause()
        logger.info("Simulation paused")
        return {"status": "paused"}
    except RuntimeError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.post("/control/stop")
async def stop_simulation():
    """Stop/Pause the simulation."""
    if not _solver_manager:
        raise HTTPException(status_code=503, detail="Solver manager not initialized")

    try:
        if _solver_manager.state == SolverState.RUNNING:
            await _solver_manager.pause()
        logger.info("Simulation stopped (paused/resumable)")
        return {"status": "stopped"}
    except Exception as e:
        logger.error("Error in stop: %s", e)
        # Even if it fails, we return 200 to keep the UI happy
        return {"status": "stopped", "error": str(e)}


@router.post("/control/reset")
async def reset_simulation():
    """Reset simulation to initial state."""
    if not _solver_manager:
        raise HTTPException(status_code=503, detail="Solver manager not initialized")

    try:
        await _solver_manager.reset()
        logger.info("Simulation reset")
        return {"status": "reset"}
    except Exception as e:
        logger.error("Error in reset: %s", e, exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/debug/health")
async def health_check():
    """Health check endpoint."""
    if not _solver_manager:
        return {"status": "error", "message": "Solver manager not initialized"}
    return {
        "status": "ok",
        "state": _solver_manager.state.value,
        "has_runner": _solver_manager.runner is not None,
        "has_task": _solver_manager.task is not None,
    }
