"""Output management placeholder (Phase 1).

Future responsibilities:
- Time series accumulation & flush policies
- NetCDF / Zarr / JSON serialization
- Snapshot frequency control
- Optional compression settings

Design principles:
- Stateless functions where possible; minimal class wrappers for policy.
- Separation between diagnostic computation and persistence.
"""
from __future__ import annotations

from typing import Dict, Any


def serialize_state(state_dict: Dict[str, Any]) -> Dict[str, Any]:
    """Identity serialization placeholder.

    Parameters
    ----------
    state_dict : dict
        Output of ModelState.to_dict()

    Returns
    -------
    dict
        Same object (copy in future if mutation safety required).
    """
    return state_dict

__all__ = ["serialize_state"]
