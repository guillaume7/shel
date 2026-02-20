import os
import shutil
import tempfile

from shel.model.outputs.manager import OutputManager
from shel.model.outputs.writers import (
    write_json_output,
    write_netcdf_snapshot,
    write_parquet_timeseries,
)


def test_output_manager_scheduling():
    tmpdir = tempfile.mkdtemp()
    try:
        mgr = OutputManager(
            tmpdir,
            {"snapshot_interval": 2, "timeseries_interval": 1, "json_interval": 3},
        )
        state = {
            "grid": {
                "nx": 2,
                "ny": 2,
                "dx": 1.0,
                "dy": 1.0,
                "x_origin": 0.0,
                "y_origin": 0.0,
            },
            "fields": {
                "eta": [[0, 0], [0, 0]],
                "u": [[0, 0, 0], [0, 0, 0]],
                "v": [[0, 0], [0, 0], [0, 0]],
                "d": [[1, 1], [1, 1]],
                "H": [[1, 1], [1, 1]],
            },
            "parameters": {"gravity": 9.81, "viscosity": 0.0, "bottom_drag_coef": 0.0},
            "diagnostics": {
                "kinetic_energy": 0.0,
                "potential_energy": 0.0,
                "total_energy": 0.0,
                "volume": 0.0,
            },
            "time": 0.0,
            "step": 0,
        }
        # Use lists for timeseries columns
        timeseries = {"time": [0.0], "E": [1.0]}
        diagnostics = {"step": 0, "energy": 1.0}
        # Step 0: all outputs
        mgr.maybe_write_snapshot(state, 0, write_netcdf_snapshot)
        mgr.maybe_write_timeseries(timeseries, 0, write_parquet_timeseries)
        mgr.maybe_write_json(diagnostics, 0, write_json_output)
        # Step 1: timeseries only
        mgr.maybe_write_snapshot(state, 1, write_netcdf_snapshot)
        mgr.maybe_write_timeseries(timeseries, 1, write_parquet_timeseries)
        mgr.maybe_write_json(diagnostics, 1, write_json_output)
        # Step 2: snapshot and timeseries
        mgr.maybe_write_snapshot(state, 2, write_netcdf_snapshot)
        mgr.maybe_write_timeseries(timeseries, 2, write_parquet_timeseries)
        mgr.maybe_write_json(diagnostics, 2, write_json_output)
        # Step 3: timeseries and json
        mgr.maybe_write_snapshot(state, 3, write_netcdf_snapshot)
        mgr.maybe_write_timeseries(timeseries, 3, write_parquet_timeseries)
        mgr.maybe_write_json(diagnostics, 3, write_json_output)
        # Check files
        files = os.listdir(tmpdir)
        assert "snapshot_step0.nc" in files
        assert "snapshot_step2.nc" in files
        assert "timeseries.parquet" in files
        assert "diagnostics_step0.json" in files
        assert "diagnostics_step3.json" in files
    finally:
        shutil.rmtree(tmpdir)
