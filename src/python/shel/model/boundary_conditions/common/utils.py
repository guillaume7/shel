from __future__ import annotations

import numpy as np

Array = np.ndarray


def mean_c_along_side(H: Array, g: float, side: str) -> float:
    if side == "west":
        edge = H[:, 0]
    elif side == "east":
        edge = H[:, -1]
    elif side == "south":
        edge = H[0, :]
    elif side == "north":
        edge = H[-1, :]
    else:
        edge = H
    return float(np.sqrt(g * float(np.mean(edge))))


def mean_H_along_side(H: Array, side: str) -> float:
    if side == "west":
        edge = H[:, 0]
    elif side == "east":
        edge = H[:, -1]
    elif side == "south":
        edge = H[0, :]
    elif side == "north":
        edge = H[-1, :]
    else:
        edge = H
    return float(np.maximum(1e-12, np.mean(edge)))
