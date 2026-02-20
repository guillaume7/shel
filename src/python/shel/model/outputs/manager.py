"""
OutputManager: schedules and manages output writes (snapshots, time-series, diagnostics).
"""

from __future__ import annotations

import os
from typing import Any, Callable, Dict


class OutputManager:
    def __init__(self, output_dir: str, schedule: Dict[str, Any]):
        """
        Args:
            output_dir: directory to write outputs
            schedule: dict with keys (e.g. 'snapshot_interval', 'timeseries_interval')
        """
        self.output_dir = output_dir
        self.schedule = schedule
        self.last_snapshot_step = None
        self.last_timeseries_step = None

    def maybe_write_snapshot(self, state: Dict[str, Any], step: int, writer: Callable):
        interval = self.schedule.get("snapshot_interval", None)
        if interval is not None and (step % interval == 0 or step == 0):
            path = os.path.join(self.output_dir, f"snapshot_step{step}.nc")
            writer(state, path)
            self.last_snapshot_step = step

    def maybe_write_timeseries(
        self, timeseries: Dict[str, Any], step: int, writer: Callable
    ):
        interval = self.schedule.get("timeseries_interval", None)
        if interval is not None and (step % interval == 0 or step == 0):
            path = os.path.join(self.output_dir, "timeseries.parquet")
            writer(timeseries, path)
            self.last_timeseries_step = step

    def maybe_write_json(self, data: Dict[str, Any], step: int, writer: Callable):
        interval = self.schedule.get("json_interval", None)
        if interval is not None and (step % interval == 0 or step == 0):
            path = os.path.join(self.output_dir, f"diagnostics_step{step}.json")
            writer(data, path)


__all__ = ["OutputManager"]
