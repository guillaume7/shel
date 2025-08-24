#!/usr/bin/env python
"""
Command-line interface for the SHEL model.
"""
import argparse
import logging
import os
import sys
from typing import Dict, Any, Optional

import yaml

from shel.config import config_loader
from shel.model import model_runner


def setup_logging(verbosity: int) -> None:
    """
    Set up logging based on verbosity level.

    Args:
        verbosity: Verbosity level (0-3)
    """
    log_levels = {
        0: logging.WARNING,
        1: logging.INFO,
        2: logging.DEBUG,
        3: logging.DEBUG,
    }
    log_level = log_levels.get(verbosity, logging.INFO)

    # Configure root logger
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout),
        ],
    )

    # If verbosity is at the highest level, enable debug logging for all modules
    if verbosity >= 3:
        logging.getLogger().setLevel(logging.DEBUG)


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns:
        Parsed arguments
    """
    parser = argparse.ArgumentParser(description="SHEL - SHallow-water numerical modEL")

    parser.add_argument(
        "-c",
        "--config",
        type=str,
        required=True,
        help="Path to the YAML configuration file",
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        type=str,
        default="./output",
        help="Directory for output files (default: ./output)",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase verbosity (can be specified multiple times)",
    )

    parser.add_argument(
        "-p",
        "--pub-port",
        type=int,
        default=5556,
        help="Port for ZeroMQ publisher (default: 5556)",
    )

    parser.add_argument(
        "--no-output",
        action="store_true",
        help="Do not write output files (useful for benchmarking)",
    )

    return parser.parse_args()


def main() -> int:
    """
    Main entry point for the CLI.

    Returns:
        Exit code
    """
    args = parse_args()
    setup_logging(args.verbose)

    logger = logging.getLogger(__name__)
    logger.info("Starting SHEL model")

    # Ensure output directory exists
    if not args.no_output:
        os.makedirs(args.output_dir, exist_ok=True)

    # Load configuration
    try:
        config = config_loader.load_config(args.config)
    except Exception as e:
        logger.error(f"Failed to load configuration: {e}")
        return 1

    # Add CLI arguments to config
    config["output"] = {
        "directory": args.output_dir,
        "enabled": not args.no_output,
    }

    config["communication"] = {
        "zmq_pub_port": args.pub_port,
    }

    # Run the model
    try:
        model_runner.run(config)
    except Exception as e:
        logger.error(f"Model run failed: {e}")
        if args.verbose >= 2:
            logger.exception("Stack trace:")
        return 2

    logger.info("Model run completed successfully")
    return 0


if __name__ == "__main__":
    sys.exit(main())
