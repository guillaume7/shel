"""
Minimal Parquet writers for SHEL outputs.

These functions wrap the existing utilities under `shel.io.parquet_reader`
to provide a stable outputs-domain API used by the solver/manager.
"""

from __future__ import annotations

import logging
from typing import Any, Dict

from ....io import parquet_reader

logger = logging.getLogger(__name__)


def write_timeseries(timeseries_data: Dict[str, Any], file_path: str) -> None:
    """Write full timeseries to Parquet file (overwrites)."""
    logger.info("Writing Parquet timeseries to %s", file_path)
    parquet_reader.write_timeseries(timeseries_data, file_path)


def append_timeseries(new_data: Dict[str, Any], file_path: str) -> None:
    """Append a single row to an existing Parquet timeseries file (or create)."""
    logger.info("Appending Parquet timeseries to %s", file_path)
    parquet_reader.append_timeseries(new_data, file_path)
