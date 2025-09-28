import numpy as np

from shel.model.diagnostics.enstrophy import enstrophy_field, integrated_enstrophy


def test_enstrophy_field():
    vorticity = np.array([1.0, -2.0, 0.0])
    enst = enstrophy_field(vorticity)
    expected = 0.5 * vorticity**2
    assert np.allclose(enst, expected)


def test_integrated_enstrophy():
    vorticity = np.ones((3, 3))
    dx = dy = 1.0
    enst = integrated_enstrophy(vorticity, dx, dy)
    # Should be 0.5 * 9 * 1^2 * 1 * 1 = 4.5
    assert np.isclose(enst, 4.5)
    # With mask
    mask = np.zeros((3, 3), dtype=bool)
    mask[0, 0] = True
    enst_masked = integrated_enstrophy(vorticity, dx, dy, mask=mask)
    assert np.isclose(enst_masked, 0.5)
