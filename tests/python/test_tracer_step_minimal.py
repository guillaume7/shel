import numpy as np

from shel.model.solvers.tracer.core import tracer_step_minimal


def test_tracer_step_minimal_applies_bcs():
    ny, nx = 10, 14
    C = np.zeros((ny, nx))
    C[:, 1:] = 1.0
    U = np.zeros((ny, nx + 1))
    V = np.zeros((ny + 1, nx))

    cfg = {"boundary_conditions": {"west": "closed", "east": "radiative", "south": "closed", "north": "closed"}}

    Cn = tracer_step_minimal(C, U, V, dt=0.1, dx=1.0, dy=1.0, config=cfg)

    # west boundary copied from interior by closed BC; east boundary copies from interior by radiative
    assert np.allclose(Cn[:, 0], Cn[:, 1])
    assert np.allclose(Cn[:, -1], Cn[:, -2])
