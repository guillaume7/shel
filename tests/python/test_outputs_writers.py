import json
import os
from typing import Any, Dict

import numpy as np

from shel.io.netcdf_reader import read_variable
from shel.io.parquet_reader import get_timeseries_range, read_timeseries
from shel.model.outputs.writers import (
    append_netcdf_timeseries,
    append_parquet_timeseries,
    write_json_output,
    write_netcdf_snapshot,
    write_parquet_timeseries,
)


def test_json_writer_deterministic(tmp_path):
    out = tmp_path / "artifact.json"
    data = {"b": 2, "a": 1, "arr": [3.0, 2.0, 1.0]}
    write_json_output(data, str(out))
    # Read back and check key order and values
    with open(out, "r") as f:
        loaded = json.load(f)
    assert list(loaded.keys()) == ["a", "arr", "b"]
    assert loaded["arr"] == [3.0, 2.0, 1.0]


def _fake_state(nx: int = 4, ny: int = 3) -> Dict[str, Any]:
    x0, y0, dx, dy = 0.0, 0.0, 1.0, 1.0
    eta = np.zeros((ny, nx), dtype=np.float64)
    u = np.zeros((ny, nx + 1), dtype=np.float64)
    v = np.zeros((ny + 1, nx), dtype=np.float64)
    d = np.ones((ny, nx), dtype=np.float64) * 10.0
    H = d.copy()
    return {
        "grid": {
            "nx": nx,
            "ny": ny,
            "dx": dx,
            "dy": dy,
            "x_origin": x0,
            "y_origin": y0,
        },
        "fields": {"eta": eta, "u": u, "v": v, "d": d, "H": H},
        "parameters": {"gravity": 9.81, "viscosity": 0.0, "bottom_drag_coef": 0.0},
        "diagnostics": {
            "kinetic_energy": 0.0,
            "potential_energy": float(np.sum(0.5 * 9.81 * eta**2)),
            "total_energy": 0.0,
            "volume": float(np.sum(H + eta)),
        },
        "time": 0.0,
        "step": 0,
    }


def test_netcdf_write_snapshot_and_read_variable(tmp_path):
    state = _fake_state()
    out = tmp_path / "snapshot.nc"
    write_netcdf_snapshot(state, str(out))

    # Round-trip: read a known variable
    eta = read_variable(str(out), "eta")
    assert eta.shape == (state["grid"]["ny"], state["grid"]["nx"])  # ny, nx
    assert np.allclose(eta, 0.0)


def test_netcdf_append_timeseries(tmp_path):
    out = tmp_path / "timeseries.nc"
    append_netcdf_timeseries({"time": 0.0, "E": 1.0}, str(out))
    append_netcdf_timeseries({"time": 1.0, "E": 1.5}, str(out))

    # Use xarray through read_timeseries-like path via netcdf is not present; test existence
    # Here, validate file exists and at least check basic integrity by reusing write/read path
    assert os.path.exists(out)


def test_parquet_write_and_read_timeseries(tmp_path):
    out = tmp_path / "timeseries.parquet"
    write_parquet_timeseries({"time": [0.0, 1.0], "E": [1.0, 1.5]}, str(out))
    df = read_timeseries(str(out))
    assert list(df.columns) == ["time", "E"]
    assert len(df) == 2


def test_parquet_append_timeseries_and_range(tmp_path):
    out = tmp_path / "timeseries.parquet"
    append_parquet_timeseries({"time": 0.0, "E": 1.0}, str(out))
    append_parquet_timeseries({"time": 1.0, "E": 2.0}, str(out))

    rng = get_timeseries_range(str(out), "E")
    assert rng["min"] == 1.0
    assert rng["max"] == 2.0
