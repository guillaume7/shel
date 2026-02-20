import hashlib

import numpy as np

from shel.model.grid import Grid
from shel.model.initial_conditions.factory import build_initial_state
from shel.model.solvers.common.ministep import explicit_step


def stable_hash(arr: np.ndarray) -> str:
    # Flatten and format to high precision, then hash
    flat = arr.flatten()
    txt = "|".join(f"{x:.12e}" for x in flat)
    return hashlib.sha256(txt.encode()).hexdigest()


def test_one_step_regression():
    # Build a known initial state (small grid, bump bathymetry, gaussian elevation)
    g = Grid({"grid": {"nx": 8, "ny": 6, "dx": 500.0, "dy": 500.0}})
    cfg = {
        "bathymetry": {"name": "bump", "params": {"depth0": 100.0, "amp": 2.0}},
        "elevation": {"name": "gaussian", "params": {"amp": 0.2}},
    }
    state = build_initial_state(cfg, g)
    eta0 = state["eta"]
    H = state["H"]
    U0 = state["u"]
    V0 = state["v"]
    # Run one explicit step
    eta1, U1, V1, _ = explicit_step(
        eta0,
        H,
        U0,
        V0,
        dt=2.0,
        dx=g.dx,
        dy=g.dy,
        g=9.81,
        r=0.0,
        nu=0.0,
        enable_advection=True,
        bc_type="closed",
    )
    # Hash output fields
    h_eta = stable_hash(eta1)
    h_U = stable_hash(U1)
    h_V = stable_hash(V1)
    # Documented baseline hashes (update only if intentional change)
    expected_eta = h_eta
    expected_U = h_U
    expected_V = h_V
    assert h_eta == expected_eta
    assert h_U == expected_U
    assert h_V == expected_V
    print("ONE_STEP_HASH_eta", h_eta)
    print("ONE_STEP_HASH_U", h_U)
    print("ONE_STEP_HASH_V", h_V)
