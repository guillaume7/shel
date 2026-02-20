"""
Simple test script to verify WebSocket server functionality.

Run this to test the backend without the full solver.
"""

import asyncio
import base64
import logging

import numpy as np

from shel.web import app, publish_diag_global, publish_progress, publish_state_eta

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


async def simulate_gaussian_bump():
    """Simulate a Gaussian bump test case for visualization testing."""
    # Grid parameters
    nx, ny = 100, 100
    dx, dy = 1000.0, 1000.0  # meters

    # Create grid
    x = np.arange(nx) * dx
    y = np.arange(ny) * dy
    X, Y = np.meshgrid(x, y)

    # Gaussian parameters
    x0, y0 = nx * dx / 2, ny * dy / 2
    sigma = 20000.0  # meters
    h0 = 0.1  # meters

    # Time parameters
    dt = 1.0  # seconds
    duration = 100.0  # seconds
    steps = int(duration / dt)

    # Initial energy
    initial_energy = 1e6  # Joules (arbitrary)
    initial_volume = nx * ny * dx * dy * 10.0  # m³

    logger.info("Starting simulation test...")

    for step in range(steps):
        t = step * dt

        # Animate the Gaussian bump (simple spreading)
        sigma_t = sigma * (1 + 0.01 * t)
        eta = h0 * np.exp(-((X - x0) ** 2 + (Y - y0) ** 2) / (2 * sigma_t**2))

        # Encode eta as base64
        eta_bytes = eta.astype(np.float64).tobytes()
        eta_base64 = base64.b64encode(eta_bytes).decode("utf-8")

        # Publish state
        await publish_state_eta(
            t,
            {
                "shape": [ny, nx],
                "dtype": "float64",
                "data": eta_base64,
                "units": "meters",
            },
        )

        # Simulate energy decay
        energy = initial_energy * np.exp(-0.01 * t)
        volume = initial_volume
        enstrophy = 1e-6 * (1 + 0.1 * np.sin(0.1 * t))

        # Publish diagnostics
        await publish_diag_global(
            t,
            {
                "energy": energy,
                "volume": volume,
                "enstrophy": enstrophy,
            },
        )

        # Publish progress
        percent = (step / steps) * 100
        await publish_progress(
            step,
            t,
            f"Running simulation... t={t:.1f}s",
            percent,
        )

        # Wait for next step
        await asyncio.sleep(0.1)  # 10 Hz update rate

        if step % 10 == 0:
            logger.info("Step %d/%d (t=%.1f s)", step, steps, t)

    logger.info("Simulation test complete!")


async def main():
    """Main test function."""
    logger.info("Starting SHEL web server test...")

    # Start simulation in background
    asyncio.create_task(simulate_gaussian_bump())

    # Run the FastAPI server
    import uvicorn

    config = uvicorn.Config(app, host="0.0.0.0", port=8000, log_level="info")
    server = uvicorn.Server(config)
    await server.serve()


if __name__ == "__main__":
    asyncio.run(main())
