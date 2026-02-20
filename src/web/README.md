# SHEL Web UI

Modern web-based interface for the SHEL shallow water model, replacing the PyQt GUI with HTML5 + TypeScript + SCSS.

## Features

- **Real-time visualization** with Plotly.js
  - Water elevation (η) heatmap
  - Velocity field plots
  - Energy and diagnostics time series
- **WebSocket communication** for low-latency updates
- **REST API** for configuration and control
- **Dark theme** optimized for long debugging sessions
- **Responsive layout** with 2x2 plot grid

## Architecture

### Backend (Python)
- **FastAPI** WebSocket server (`src/python/shel/web/server.py`)
- **REST API** for control (`src/python/shel/web/api.py`)
- Replaces ZeroMQ PUB/SUB with browser-compatible WebSocket

### Frontend (TypeScript)
- **Vite** for fast development and building
- **Plotly.js** for scientific visualization
- **SCSS** for maintainable styling
- **WebSocket client** with automatic reconnection

## Quick Start

### 1. Install Dependencies

```bash
# Python dependencies (if not already installed)
pip install fastapi uvicorn websockets

# Frontend dependencies
cd src/web
npm install
```

### 2. Start the Backend Server

```bash
# From project root
python -m shel.web.run_server

# Or with uvicorn directly
uvicorn shel.web.server:app --reload --host 0.0.0.0 --port 8000
```

### 3. Start the Frontend Dev Server

```bash
# In another terminal
cd src/web
npm run dev
```

### 4. Open in Browser

Navigate to `http://localhost:3000`

## Development

### Backend Development

The backend server provides:
- WebSocket endpoint at `/ws`
- REST API at `/api/*`
- Health check at `/api/health`

To publish updates from the solver:

```python
from shel.web import publish_state_eta, publish_diag_global

# In your solver loop
await publish_state_eta(t, {
    "shape": [ny, nx],
    "dtype": "float64",
    "data": base64_encoded_eta,
    "units": "meters"
})

await publish_diag_global(t, {
    "energy": total_energy,
    "volume": total_volume,
    "enstrophy": enstrophy,
})
```

### Frontend Development

The frontend is organized as:
- `src/main.ts` - Application entry point
- `src/websocket.ts` - WebSocket client
- `src/plots/` - Plotly visualization components
- `src/styles/` - SCSS stylesheets

Hot reload is enabled for rapid development.

## Building for Production

```bash
cd src/web
npm run build
```

The built files will be in `src/web/dist/`.

## Debugging Numerical Instabilities

The web UI is specifically designed to help debug numerical instabilities:

1. **Real-time visualization** shows instabilities as they develop
2. **Browser DevTools** for inspecting WebSocket messages
3. **Time series plots** reveal energy/volume conservation issues
4. **Responsive controls** for quick pause/resume

## Next Steps

- [ ] Add velocity vector field plot
- [ ] Implement parameter panel with live updates
- [ ] Add export functionality (PNG, data)
- [ ] Create configuration save/load
- [ ] Add more diagnostic plots (vorticity, Okubo-Weiss)
