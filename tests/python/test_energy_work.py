import numpy as np

from shel.model.diagnostics.energy import (
    drag_dissipation,
    kinetic_energy,
    potential_energy,
    wind_work,
)
from shel.model.forcings.bottom import linear_drag
from shel.model.forcings.surface import wind_stress


def test_energy_work_rates():
    # Simple domain
    nx, ny = 4, 3
    H = np.ones((ny, nx)) * 10.0
    u = np.ones((ny, nx)) * 2.0
    v = np.ones((ny, nx)) * 0.0
    eta = np.ones((ny, nx)) * 0.5
    # Wind stress: constant wind, no surface current
    tau = wind_stress(10.0, 0.0)
    # Drag: linear drag
    drag = linear_drag(u, 0.01)
    # Energy rates
    KE = kinetic_energy(u, v, H)
    PE = potential_energy(eta, H)
    WW = wind_work(tau, u, v)
    DD = drag_dissipation(drag, u, v)
    # Physical checks
    assert KE > 0
    assert PE > 0
    assert WW > 0  # wind does positive work
    assert DD < 0  # drag dissipates energy
    # Magnitude checks (order of magnitude)
    assert np.isclose(KE, 0.5 * np.sum(H * u**2))
    assert np.isclose(PE, 0.5 * 9.81 * np.sum(eta**2))
    assert np.isclose(WW, np.sum(tau * u))
    assert np.isclose(DD, np.sum(drag * u))
