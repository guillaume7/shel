"""
Generate a multi-step dynamic golden baseline for regression tests.

This script runs a short, deterministic shallow-water integration using the
explicit ministep and prints a JSON blob with parameters, metrics, and a
stable hash. Pipe it to tests/python/golden/dynamic_baseline.json to enable
the regression test.

Usage (from repo root):
  python -m devops.scripts.generate_dynamic_golden > tests/python/golden/dynamic_baseline.json
"""

from __future__ import annotations

import json
import hashlib
import numpy as np

from shel.model.solvers.common.ministep import explicit_step
from shel.model.diagnostics import IntegratedDiagnostics
from shel.model.grid import Grid


def stable_hash(d: dict) -> str:
    items = [f"{k}={d[k]:.12e}" for k in sorted(d.keys())]
    txt = "|".join(items)
    return hashlib.sha256(txt.encode()).hexdigest()


def run_short_sim(params: dict) -> dict:
    ny = params.get("ny", 24)
    nx = params.get("nx", 30)
    dx = params.get("dx", 1.0)
    dy = params.get("dy", 1.0)
    dt = params.get("dt", 0.02)
    steps = params.get("steps", 30)
    g = params.get("g", 9.81)
    r = params.get("r", 0.01)
    nu = params.get("nu", 0.001)
    H0 = params.get("H0", 20.0)
    eta_amp = params.get("eta_amp", 0.05)
    sigma = params.get("sigma", 5.0)

    # Fields
    H = np.full((ny, nx), H0, dtype=float)
    y = np.arange(ny, dtype=float)
    x = np.arange(nx, dtype=float)
    X, Y = np.meshgrid(x, y)
    eta = eta_amp * np.exp(-(((X - nx / 2) ** 2 + (Y - ny / 2) ** 2) / (2.0 * sigma ** 2)))
    U = np.zeros((ny, nx + 1), dtype=float)
    V = np.zeros((ny + 1, nx), dtype=float)
    f = np.zeros_like(H)

    # Run
    for _ in range(steps):
        eta, U, V = explicit_step(
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

    # Metrics
    grid = Grid({"grid": {"nx": nx, "ny": ny, "dx": dx, "dy": dy}})
    diags = IntegratedDiagnostics.as_dict(U, V, eta, H, grid=grid, gravity=g, coriolis=f)
    # Compose a compact set of metrics for hashing and human-inspection
    metrics = {
        "total_energy": float(diags["total_energy"]),
        "volume": float(diags["volume"]),
        "eta_rms": float(np.sqrt(np.mean(eta**2))),
        "eta_max_abs": float(np.max(np.abs(eta))),
        "u_rms": float(np.sqrt(np.mean(U**2))),
        "v_rms": float(np.sqrt(np.mean(V**2))),
    }
    params_out = {
        "ny": ny,
        "nx": nx,
        "dx": dx,
        "dy": dy,
        "dt": dt,
        "steps": steps,
        "g": g,
        "r": r,
        "nu": nu,
        "H0": H0,
        "eta_amp": eta_amp,
        "sigma": sigma,
    }

    hval = stable_hash(metrics)
    return {"params": params_out, "metrics": metrics, "hash": hval}


def main() -> None:
    baseline = run_short_sim({})
    print(json.dumps(baseline, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
