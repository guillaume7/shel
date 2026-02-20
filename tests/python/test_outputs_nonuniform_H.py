import os
import tempfile

import numpy as np

from shel.io import netcdf_reader, parquet_reader
from shel.model.outputs.writers.json_writer import write_json
from shel.model.outputs.writers.netcdf import write_snapshot
from shel.model.outputs.writers.parquet import write_timeseries


def test_netcdf_nonuniform_H():
    H = np.array([[1.0, 2.0], [3.0, 4.0]])
    eta = np.array([[0.1, 0.2], [0.3, 0.4]])
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
            "H": H,
            "eta": eta,
            "u": np.zeros((2, 3)),
            "v": np.zeros((3, 2)),
            "d": np.zeros((2, 2)),
        },
        "parameters": {"gravity": 9.81, "viscosity": 1e-6, "bottom_drag_coef": 0.0025},
        "diagnostics": {
            "kinetic_energy": 0.0,
            "potential_energy": 0.0,
            "total_energy": 0.0,
            "volume": 0.0,
        },
        "time": 0.0,
        "step": 0,
    }
    with tempfile.TemporaryDirectory() as tmpdir:
        path = os.path.join(tmpdir, "snap.nc")
        write_snapshot(state, path)
        H_out = netcdf_reader.read_variable(path, "H")
        eta_out = netcdf_reader.read_variable(path, "eta")
        assert np.allclose(H_out, H)
        assert np.allclose(eta_out, eta)


def test_parquet_nonuniform_H():
    H = np.array([[1.0, 2.0], [3.0, 4.0]])
    eta = np.array([[0.1, 0.2], [0.3, 0.4]])
    timeseries = {"step": [0], "H": [H.tolist()], "eta": [eta.tolist()]}
    with tempfile.TemporaryDirectory() as tmpdir:
        path = os.path.join(tmpdir, "ts.parquet")
        write_timeseries(timeseries, path)
        df = parquet_reader.read_timeseries(path)
    H_out = np.stack(df["H"][0])
    eta_out = np.stack(df["eta"][0])
    assert np.allclose(H_out, H)
    assert np.allclose(eta_out, eta)


def test_json_nonuniform_H():
    import json

    H = np.array([[1.0, 2.0], [3.0, 4.0]])
    eta = np.array([[0.1, 0.2], [0.3, 0.4]])
    with tempfile.TemporaryDirectory() as tmpdir:
        path = os.path.join(tmpdir, "snap.json")
        write_json({"H": H.tolist(), "eta": eta.tolist()}, path)
        with open(path, "r") as f:
            out = json.load(f)
        assert np.allclose(out["H"], H)
        assert np.allclose(out["eta"], eta)
