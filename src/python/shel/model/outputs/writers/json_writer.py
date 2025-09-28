"""
Minimal deterministic JSON writer for SHEL outputs.
Ensures stable key ordering and float64 normalization for regression artifacts.
"""

import json
import logging
from typing import Any, Dict

logger = logging.getLogger(__name__)


def write_json(data: Dict[str, Any], file_path: str) -> None:
    """Write a dict to JSON with sorted keys and float64 normalization."""

    def _normalize(obj):
        if isinstance(obj, dict):
            return {k: _normalize(obj[k]) for k in sorted(obj)}
        elif isinstance(obj, list):
            return [_normalize(x) for x in obj]
        elif hasattr(obj, "dtype") and str(getattr(obj, "dtype", "")) == "float64":
            return [float(x) for x in obj.tolist()]
        elif isinstance(obj, float):
            return float(obj)
        return obj

    normed = _normalize(data)
    with open(file_path, "w") as f:
        json.dump(normed, f, sort_keys=True, indent=2, separators=(",", ": "))
    logger.info("Wrote deterministic JSON to %s", file_path)
