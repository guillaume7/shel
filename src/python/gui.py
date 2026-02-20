#!/usr/bin/env python
"""
Entry point for the SHEL GUI.

This script is the entry point for the SHEL GUI, allowing users to
interact with the model through a graphical interface.
"""
import logging
import os
import sys

from shel.gui.main_window import main

if __name__ == "__main__":
    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    # Run the GUI
    sys.exit(main())
