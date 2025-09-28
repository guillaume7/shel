# SHEL ZeroMQ Pub/Sub Protocol

This document describes the message model and example messages for the SHEL solver-GUI communication using ZeroMQ PUB/SUB.

## Topics
- `state.eta` — Water elevation field
- `state.velocity` — Velocity fields (U, V)
- `diag.global` — Global diagnostics (energy, enstrophy, volume, etc.)
- `diag.field.<name>` — Field diagnostics (e.g., vorticity, PV)
- `event.progress` — Simulation progress/events

## Message Format
All messages are JSON objects, sent as UTF-8 encoded strings.

### Example: Water Elevation Update
**Topic:** `state.eta`
```json
{
  "t": 12.0,
  "shape": [100, 120],
  "dtype": "float64",
  "data": "<base64-encoded array>",
  "units": "meters"
}
```

### Example: Velocity Update
**Topic:** `state.velocity`
```json
{
  "t": 12.0,
  "U": {
    "shape": [100, 121],
    "dtype": "float64",
    "data": "<base64-encoded array>"
  },
  "V": {
    "shape": [101, 120],
    "dtype": "float64",
    "data": "<base64-encoded array>"
  },
  "units": "m/s"
}
```

### Example: Global Diagnostics
**Topic:** `diag.global`
```json
{
  "t": 12.0,
  "energy": 1.234,
  "enstrophy": 0.567,
  "volume": 12345.6
}
```

### Example: Field Diagnostic
**Topic:** `diag.field.vorticity`
```json
{
  "t": 12.0,
  "shape": [100, 120],
  "dtype": "float64",
  "data": "<base64-encoded array>",
  "units": "1/s"
}
```

### Example: Progress Event
**Topic:** `event.progress`
```json
{
  "step": 120,
  "t": 12.0,
  "message": "Simulation running",
  "percent": 60.0
}
```

## Notes
- Data arrays are base64-encoded for compactness and cross-language compatibility.
- Units and dtype are included for clarity.
- Topics allow the GUI to subscribe only to relevant updates.
- Extensible: Add new topics for new diagnostics or events as needed.
