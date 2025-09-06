"""Bottom forcings placeholder (Phase 1).

Behavior intentionally inert. Future phases (see refactor prompt) will
introduce modules under ``forcings/bottom/`` (e.g. drag.py) and this
module will become a thin backward-compatible wrapper or be removed.
"""

from __future__ import annotations

from typing import Any, Dict

import numpy as np
from numpy.typing import NDArray


def apply_bottom_drag(
    u: NDArray, v: NDArray, H: NDArray, params: Dict[str, Any]
) -> tuple[NDArray, NDArray]:
    """Identity bottom drag (no-op) for Phase 1 skeleton."""
    return u, v


__all__ = ["apply_bottom_drag"]
