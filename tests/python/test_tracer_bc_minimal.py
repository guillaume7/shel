import numpy as np

from shel.model.boundary_conditions.registry import get_tracer_bc, list_tracer_bcs


def test_tracer_registry_contains_closed_and_radiative():
    kinds = list_tracer_bcs()
    assert "closed" in kinds
    assert "radiative" in kinds


ess = ["west", "east", "south", "north"]


def test_tracer_closed_zero_gradient():
    ny, nx = 5, 6
    C = np.arange(ny * nx, dtype=float).reshape(ny, nx)
    cls = get_tracer_bc("closed")
    assert cls is not None
    bc = cls()
    for side in ess:
        C2 = C.copy()
        bc.apply_side_tracer(C2, side)
        if side == "west":
            assert np.allclose(C2[:, 0], C2[:, 1])
        if side == "east":
            assert np.allclose(C2[:, -1], C2[:, -2])
        if side == "south":
            assert np.allclose(C2[0, :], C2[1, :])
        if side == "north":
            assert np.allclose(C2[-1, :], C2[-2, :])


def test_tracer_radiative_one_sided_copy():
    ny, nx = 5, 6
    C = np.random.RandomState(0).rand(ny, nx)
    cls = get_tracer_bc("radiative")
    assert cls is not None
    bc = cls()
    for side in ess:
        C2 = C.copy()
        bc.apply_side_tracer(C2, side)
        if side == "west":
            assert np.allclose(C2[:, 0], C[:, 1])
        if side == "east":
            assert np.allclose(C2[:, -1], C[:, -2])
        if side == "south":
            assert np.allclose(C2[0, :], C[1, :])
        if side == "north":
            assert np.allclose(C2[-1, :], C[-2, :])
