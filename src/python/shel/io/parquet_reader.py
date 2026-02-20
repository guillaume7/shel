"""
Parquet I/O utilities for SHEL.

This module provides functions for reading and writing Parquet files,
which are particularly efficient for timeseries data.
"""

import logging
import os
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

logger = logging.getLogger(__name__)


def write_timeseries(timeseries_data: Dict[str, Any], file_path: str) -> None:
    """
    Write timeseries data to a Parquet file.

    Args:
        timeseries_data: Dictionary with keys as column names and values as lists of data
        file_path: Path to save the Parquet file
    """
    # Convert to pandas DataFrame
    df = pd.DataFrame(timeseries_data)

    # Convert to pyarrow Table
    table = pa.Table.from_pandas(df)

    # Write to Parquet file
    pq.write_table(table, file_path)

    logger.info("Wrote timeseries data to %s", file_path)


def append_timeseries(new_data: Dict[str, Any], file_path: str) -> None:
    """
    Append a new row of timeseries data to an existing Parquet file.

    Args:
        new_data: Dictionary with keys as column names and values as single data points
        file_path: Path to the Parquet file
    """
    # Check if file exists
    if os.path.exists(file_path):
        # Read existing data
        df = pd.read_parquet(file_path)

        # Convert new data to DataFrame (single row)
        new_df = pd.DataFrame([new_data])

        # Append new data
        df = pd.concat([df, new_df], ignore_index=True)
    else:
        # Create new DataFrame
        df = pd.DataFrame([new_data])

    # Write back to file
    df.to_parquet(file_path)

    logger.info("Appended timeseries data to %s", file_path)


def read_timeseries(file_path: str) -> pd.DataFrame:
    """
    Read timeseries data from a Parquet file.

    Args:
        file_path: Path to the Parquet file

    Returns:
        DataFrame containing the timeseries data

    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Parquet file not found: {file_path}")

    # Read Parquet file
    df = pd.read_parquet(file_path)

    logger.info(
        "Read timeseries data from %s: %d rows, %d columns",
        file_path,
        len(df),
        len(df.columns),
    )
    return df


def get_timeseries_variables(file_path: str) -> List[str]:
    """
    Get the list of variables in a timeseries Parquet file.

    Args:
        file_path: Path to the Parquet file

    Returns:
        List of variable names

    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Parquet file not found: {file_path}")

    # Read Parquet file metadata
    metadata = pq.read_metadata(file_path)
    schema = metadata.schema

    # Get column names from schema
    variables = [schema.names[i] for i in range(schema.num_fields)]

    logger.info("Found %s variables in %s", len(variables), file_path)
    return variables


def get_timeseries_range(file_path: str, variable: str) -> Dict[str, float]:
    """
    Get the range of values for a variable in a timeseries Parquet file.

    Args:
        file_path: Path to the Parquet file
        variable: Name of the variable

    Returns:
        Dictionary with min and max values

    Raises:
        FileNotFoundError: If the file doesn't exist
        KeyError: If the variable doesn't exist in the file
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Parquet file not found: {file_path}")

    # Read Parquet file
    df = pd.read_parquet(file_path, columns=[variable])

    if variable not in df.columns:
        raise KeyError(f"Variable '{variable}' not found in file {file_path}")

    # Calculate min and max
    min_val = df[variable].min()
    max_val = df[variable].max()

    logger.info(
        "Variable '%s' range in %s: %s to %s",
        variable,
        file_path,
        min_val,
        max_val,
    )
    return {"min": min_val, "max": max_val}
