import numpy as np

from shel.model.acceleration.numba import ENABLE_NUMBA, accelerated_sum
from shel.model.diagnostics.profiling import profile_function


def test_profile_function():
    arr = np.arange(100000)
    stats, result = profile_function(np.sum, arr)
    assert result == np.sum(arr)
    assert "ncalls" in stats


def test_accelerated_sum():
    arr = np.arange(100000)
    result = accelerated_sum(arr)
    assert result == np.sum(arr)
    # If Numba is enabled, accelerated_sum should be compiled
    assert isinstance(ENABLE_NUMBA, bool)
