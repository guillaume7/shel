"""
Minimal NetCDF writers for SHEL outputs.

These functions wrap the existing I/O utilities under `shel.io.netcdf_reader`
to provide a stable path under the model outputs domain.
"""

from __future__ import annotations

import logging
from typing import Any, Dict

from ....io import netcdf_reader

logger = logging.getLogger(__name__)


def write_snapshot(state: Dict[str, Any], file_path: str) -> None:
    """Write a full model snapshot to a NetCDF file.

    Contract:
    - `state` dict contains keys: `grid`, `fields`, `parameters`, `diagnostics`, `time`, `step`.
    - Arrays follow SHEL stagger conventions and shapes.
    - File path is overwritten if exists.
    """
    logger.info("Writing NetCDF snapshot to %s", file_path)
    netcdf_reader.write_model_state(state, file_path)


def append_timeseries(timeseries_data: Dict[str, Any], file_path: str) -> None:
    """Append a single timeseries row to NetCDF `time` dimension file.

    Contract:
    - `timeseries_data` includes a `time` key (float/int) and any number of scalar metrics.
    - Creates a new file if it does not exist, otherwise appends along `time`.
    """
    logger.info("Appending NetCDF timeseries to %s", file_path)
    netcdf_reader.append_timeseries(timeseries_data, file_path)
