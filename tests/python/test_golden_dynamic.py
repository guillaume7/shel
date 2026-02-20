from __future__ import annotations

import hashlib
import json
import os

import numpy as np

from shel.model.diagnostics import IntegratedDiagnostics
from shel.model.grid import Grid
from shel.model.solvers.common.ministep import explicit_step


def stable_hash(d: dict) -> str:
    items = [f"{k}={d[k]:.12e}" for k in sorted(d.keys())]
    txt = "|".join(items)
    return hashlib.sha256(txt.encode()).hexdigest()


def run_short_sim():
    ny, nx = 24, 30
    dx = dy = 1.0
    dt = 0.02
    steps = 30
    g = 9.81
    r = 0.01
    nu = 0.001

    H0 = 20.0
    H = np.full((ny, nx), H0, dtype=float)
    y = np.arange(ny, dtype=float)
    x = np.arange(nx, dtype=float)
    X, Y = np.meshgrid(x, y)
    eta = 0.05 * np.exp(-(((X - nx / 2) ** 2 + (Y - ny / 2) ** 2) / (2.0 * 5.0**2)))
    U = np.zeros((ny, nx + 1), dtype=float)
    V = np.zeros((ny + 1, nx), dtype=float)
    f = np.zeros_like(H)

    for _ in range(steps):
        eta, U, V, _ = explicit_step(
            eta,
            H,
            U,
            V,
            dt=dt,
            dx=dx,
            dy=dy,
            g=g,
            r=r,
            nu=nu,
            enable_advection=True,
            f=f,
            enable_coriolis=False,
        )

    grid = Grid({"grid": {"nx": nx, "ny": ny, "dx": dx, "dy": dy}})
    diags = IntegratedDiagnostics.as_dict(U, V, eta, H, grid, gravity=g, coriolis=f)
    metrics = {
        "total_energy": float(diags["total_energy"]),
        "volume": float(diags["volume"]),
        "eta_rms": float(np.sqrt(np.mean(eta**2))),
        "eta_max_abs": float(np.max(np.abs(eta))),
        "u_rms": float(np.sqrt(np.mean(U**2))),
        "v_rms": float(np.sqrt(np.mean(V**2))),
    }
    return metrics


def test_dynamic_golden_baseline_hash_matches_file():
    # If the baseline file is present, verify the hash matches; otherwise skip with guidance
    baseline_path = os.path.join(
        os.path.dirname(__file__), "golden", "dynamic_baseline.json"
    )
    if not os.path.exists(baseline_path):
        import pytest

        pytest.skip(
            "Golden baseline missing. Generate one with: python -m devops.scripts.generate_dynamic_golden > tests/python/golden/dynamic_baseline.json"
        )

    with open(baseline_path, "r") as f:
        baseline = json.load(f)

    expected_hash = baseline.get("hash")
    assert expected_hash, "Baseline JSON missing 'hash' field"

    metrics = run_short_sim()
    hval = stable_hash(metrics)

    # Tight equality to catch any changes; update baseline intentionally with generator when needed
    assert (
        hval == expected_hash
    ), f"Dynamic golden hash mismatch. Got {hval}, expected {expected_hash}. If change is intended, regenerate baseline."
