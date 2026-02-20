import numpy as np

from shel.model.diagnostics.potential_enstrophy import (
    integrated_potential_enstrophy,
    potential_enstrophy_field,
)


def test_potential_enstrophy_field():
    pv = np.array([1.0, -2.0, 0.0])
    penst = potential_enstrophy_field(pv)
    expected = 0.5 * pv**2
    assert np.allclose(penst, expected)


def test_integrated_potential_enstrophy():
    pv = np.ones((3, 3))
    dx = dy = 1.0
    penst = integrated_potential_enstrophy(pv, dx, dy)
    # Should be 0.5 * 9 * 1^2 * 1 * 1 = 4.5
    assert np.isclose(penst, 4.5)
    # With mask
    mask = np.zeros((3, 3), dtype=bool)
    mask[0, 0] = True
    penst_masked = integrated_potential_enstrophy(pv, dx, dy, mask=mask)
    assert np.isclose(penst_masked, 0.5)
