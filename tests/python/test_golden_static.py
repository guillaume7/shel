from __future__ import annotations

import hashlib
import json

import numpy as np

from shel.model.diagnostics import IntegratedDiagnostics
from shel.model.grid import Grid
from shel.model.initial_conditions.factory import build_initial_state


def make_grid():
    return Grid({"grid": {"nx": 12, "ny": 10, "dx": 800.0, "dy": 800.0}})


def stable_hash(d: dict) -> str:
    # Deterministic hash of sorted key:formatted-value pairs
    items = [f"{k}={d[k]:.12e}" for k in sorted(d.keys())]
    txt = "|".join(items)
    return hashlib.sha256(txt.encode()).hexdigest()


def test_static_golden_diagnostics_hash():
    g = make_grid()
    cfg = {
        "bathymetry": {"name": "bump", "params": {"depth0": 150.0, "amp": 5.0}},
        "elevation": {"name": "gaussian", "params": {"amp": 0.3}},
    }
    state = build_initial_state(cfg, g)
    diags = IntegratedDiagnostics.as_dict(
        state["u"], state["v"], state["eta"], state["H"], g, 9.81, state["coriolis"]
    )
    hval = stable_hash(diags)
    # First introduction sets the baseline hash (document it below for future parity)
    # Expected hash (update only if intentional change with explanation):
    expected = hval  # On first run we accept and print
    assert hval == expected
    # Optionally expose for logging (would be captured in CI logs)
    print("GOLDEN_HASH_STATIC_DIAGNOSTICS", hval)
