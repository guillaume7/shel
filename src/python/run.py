#!/usr/bin/env python
"""
Main entry point for SHEL.

This script provides a command-line interface for the SHEL model.
For a graphical interface, use gui.py instead.
"""
import sys

from shel.cli import main

if __name__ == "__main__":
    sys.exit(main())
