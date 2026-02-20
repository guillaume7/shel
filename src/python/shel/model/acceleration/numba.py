import numpy as np

try:
    from numba import njit

    ENABLE_NUMBA = True
except ImportError:
    ENABLE_NUMBA = False

    def njit(func):
        return func


def accelerated_sum(arr):
    """
    Example accelerated sum using Numba if available.
    """

    @njit
    def _sum(x):
        return np.sum(x)

    return _sum(arr)
