#!/usr/bin/env python3
"""
Startup script for SHEL web server.

Run with: python -m shel.web.run_server
Or: uvicorn shel.web.server:app --reload
"""

import logging
import sys

import uvicorn

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)

logger = logging.getLogger(__name__)


def main():
    """Start the SHEL web server."""
    logger.info("Starting SHEL web server...")

    uvicorn.run(
        "shel.web.server:app",
        host="0.0.0.0",
        port=8000,
        reload=True,  # Auto-reload on code changes during development
        log_level="info",
    )


if __name__ == "__main__":
    main()
