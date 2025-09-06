import numpy as np

from shel.model.boundary_conditions.registry import get_tracer_bc
from shel.model.solvers.common.stepper import (
    apply_tracer_bc_per_side,
    resolve_bc_sides_from_config,
)


def test_tracer_bc_apply_closed_and_radiative_do_not_crash():
    ny, nx = 16, 24
    C = np.zeros((ny, nx), dtype=float)
    C_old = np.zeros_like(C)
    C[:, nx // 2 :] = 1.0

    # Two configs to exercise both tracer BCs
    cfg_closed = {
        "boundary_conditions": {
            "west": "closed",
            "east": "closed",
            "south": "closed",
            "north": "closed",
        }
    }
    cfg_open = {
        "boundary_conditions": {
            "west": "radiative",
            "east": "radiative",
            "south": "closed",
            "north": "closed",
        }
    }

    for cfg in (cfg_closed, cfg_open):
        bc_sides = resolve_bc_sides_from_config(cfg)
        Cw = C.copy()
        apply_tracer_bc_per_side(
            Cw, bc_sides=bc_sides, C_old=C_old, U=None, V=None, dt=1.0, dx=1.0, dy=1.0
        )
        assert np.isfinite(Cw).all()


def test_tracer_registry_contains_expected():
    assert get_tracer_bc("closed") is not None
    assert get_tracer_bc("radiative") is not None
