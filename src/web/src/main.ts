import './styles/main.scss';
import { WebSocketClient } from './websocket';
import { EtaPlot } from './plots/EtaPlot';
import { VelocityPlot } from './plots/VelocityPlot';
import { QuiverPlot } from './plots/QuiverPlot';
import { DiagPlot } from './plots/DiagPlot';

interface ConfigState {
    grid_nx: number;
    grid_ny: number;
    dt: number;
    viscosity: number;
}

class ShelApp {
    private wsClient: WebSocketClient;
    private etaPlot: EtaPlot;
    private velocityPlot: VelocityPlot;
    private quiverPlot: QuiverPlot;
    private diagPlot: DiagPlot;
    private config: ConfigState = {
        grid_nx: 100,
        grid_ny: 100,
        dt: 0.5,
        viscosity: 0,
    };

    constructor() {
        // Initialize WebSocket client
        const wsUrl = `ws://${window.location.hostname}:${window.location.port}/ws`;
        this.wsClient = new WebSocketClient(wsUrl);

        // Initialize plots
        this.etaPlot = new EtaPlot('plot-eta');
        this.velocityPlot = new VelocityPlot('plot-velocity');
        this.quiverPlot = new QuiverPlot('plot-quiver');
        this.diagPlot = new DiagPlot('plot-diagnostics');

        // Set up event listeners
        this.setupWebSocketHandlers();
        this.setupControlHandlers();
        this.setupParameterHandlers();

        // Connect to server
        this.wsClient.connect();

        // Load initial config
        this.loadConfig();
    }

    private setupWebSocketHandlers(): void {
        const statusEl = document.getElementById('connection-status');

        this.wsClient.on('open', () => {
            if (statusEl) {
                statusEl.textContent = 'Connected';
                statusEl.className = 'status-online';
            }
        });

        this.wsClient.on('close', () => {
            if (statusEl) {
                statusEl.textContent = 'Disconnected';
                statusEl.className = 'status-offline';
            }
        });

        // Handle eta state updates
        this.wsClient.on('state.eta', (message) => {
            const { t, data, shape } = message;
            if (data && shape) {
                const eta = this.decodeBase64ToArray(data, shape);
                this.etaPlot.update(eta, t);
            }
        });

        // Handle velocity state updates
        this.wsClient.on('state.velocity', (message) => {
            const { t, U, V } = message;
            if (U && V && U.data && V.data && U.shape) {
                const uArr = this.decodeBase64ToArray(U.data, U.shape);
                const vArr = this.decodeBase64ToArray(V.data, U.shape);
                this.velocityPlot.update(uArr, vArr, t);
                this.quiverPlot.update(uArr, vArr, t);
            }
        });

        // Handle global diagnostics
        this.wsClient.on('diag.global', (message) => {
            const { t, total_energy, kinetic_energy, potential_energy, volume, enstrophy } = message;

            // Update diagnostic display
            this.updateDiagnosticValues({
                te: total_energy,
                ke: kinetic_energy,
                pe: potential_energy,
                volume,
                enstrophy
            });

            // Add to time series plot
            this.diagPlot.addPoint(t, {
                te: total_energy,
                ke: kinetic_energy,
                pe: potential_energy,
                volume,
                enstrophy
            });
        });

        // Handle progress updates
        this.wsClient.on('event.progress', (message) => {
            const { step, t, message: msg, percent } = message;
            this.updateProgress(step, t, msg, percent);
        });
    }

    private setupControlHandlers(): void {
        const btnStart = document.getElementById('btn-start') as HTMLButtonElement;
        const btnStop = document.getElementById('btn-stop') as HTMLButtonElement;
        const btnReset = document.getElementById('btn-reset') as HTMLButtonElement;

        btnStart?.addEventListener('click', async () => {
            await this.startSimulation();
            if (btnStart) btnStart.disabled = true;
            if (btnStop) btnStop.disabled = false;
        });

        btnStop?.addEventListener('click', async () => {
            await this.stopSimulation(); // This now calls pause in the backend
            if (btnStart) btnStart.disabled = false;
            if (btnStop) btnStop.disabled = true;
        });

        btnReset?.addEventListener('click', async () => {
            await this.resetSimulation();
            if (btnStart) btnStart.disabled = false;
            if (btnStop) btnStop.disabled = true;
            this.etaPlot.clear();
            this.velocityPlot.clear();
            this.diagPlot.clear();
            this.quiverPlot.clear();
            this.updateDiagnosticValues({
                te: 0, ke: 0, pe: 0, volume: 0, enstrophy: 0
            });
            // Reset headers
            this.updateProgress(0, 0, 'Ready', 0);
        });
    }

    private setupParameterHandlers(): void {
        const inputs = {
            grid_nx: document.getElementById('param-nx') as HTMLInputElement,
            grid_ny: document.getElementById('param-ny') as HTMLInputElement,
            dt: document.getElementById('param-dt') as HTMLInputElement,
            viscosity: document.getElementById('param-nu') as HTMLInputElement,
        };

        const updateParam = async (name: keyof ConfigState, value: number) => {
            this.config[name] = value;
            await this.saveConfig();
        };

        inputs.grid_nx?.addEventListener('change', (e) => updateParam('grid_nx', parseInt((e.target as HTMLInputElement).value)));
        inputs.grid_ny?.addEventListener('change', (e) => updateParam('grid_ny', parseInt((e.target as HTMLInputElement).value)));
        inputs.dt?.addEventListener('change', (e) => updateParam('dt', parseFloat((e.target as HTMLInputElement).value)));
        inputs.viscosity?.addEventListener('change', (e) => updateParam('viscosity', parseFloat((e.target as HTMLInputElement).value)));
    }

    private async loadConfig(): Promise<void> {
        try {
            const response = await fetch('/api/config');
            const result = await response.json();

            // Map backend config to frontend state
            if (result.grid) {
                this.config.grid_nx = result.grid.nx;
                this.config.grid_ny = result.grid.ny;
            }
            if (result.model) {
                this.config.dt = result.model.timestep;
                this.config.viscosity = result.model.viscosity || 0;
            }

            // Update UI inputs
            const inputs = {
                grid_nx: document.getElementById('param-nx') as HTMLInputElement,
                grid_ny: document.getElementById('param-ny') as HTMLInputElement,
                dt: document.getElementById('param-dt') as HTMLInputElement,
                viscosity: document.getElementById('param-nu') as HTMLInputElement,
            };
            if (inputs.grid_nx) inputs.grid_nx.value = (this.config.grid_nx || 100).toString();
            if (inputs.grid_ny) inputs.grid_ny.value = (this.config.grid_ny || 100).toString();
            if (inputs.dt) inputs.dt.value = (this.config.dt || 0.5).toString();
            if (inputs.viscosity) inputs.viscosity.value = (this.config.viscosity || 0).toString();

        } catch (error) {
            console.error('Failed to load config:', error);
        }
    }

    private async saveConfig(): Promise<void> {
        try {
            await fetch('/api/config', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify(this.config),
            });
            console.log('Config saved:', this.config);
        } catch (error) {
            console.error('Failed to save config:', error);
        }
    }

    private async startSimulation(): Promise<void> {
        try {
            const response = await fetch('/api/control/start', { method: 'POST' });
            const result = await response.json();
            console.log('Simulation started:', result);
        } catch (error) {
            console.error('Failed to start simulation:', error);
        }
    }

    private async stopSimulation(): Promise<void> {
        try {
            const response = await fetch('/api/control/stop', { method: 'POST' });
            const result = await response.json();
            console.log('Simulation stopped:', result);
        } catch (error) {
            console.error('Failed to stop simulation:', error);
        }
    }

    private async resetSimulation(): Promise<void> {
        try {
            const response = await fetch('/api/control/reset', { method: 'POST' });
            const result = await response.json();
            console.log('Simulation reset:', result);
        } catch (error) {
            console.error('Failed to reset simulation:', error);
        }
    }

    private decodeBase64ToArray(base64: string, shape: number[]): number[][] {
        // Decode base64 to Float64Array
        const binaryString = atob(base64);
        const bytes = new Uint8Array(binaryString.length);
        for (let i = 0; i < binaryString.length; i++) {
            bytes[i] = binaryString.charCodeAt(i);
        }
        const float64Array = new Float64Array(bytes.buffer);

        // Reshape to 2D array
        const [ny, nx] = shape;
        const array2d: number[][] = [];
        for (let i = 0; i < ny; i++) {
            array2d.push(Array.from(float64Array.slice(i * nx, (i + 1) * nx)));
        }
        return array2d;
    }

    private updateDiagnosticValues(diag: { te?: number; ke?: number; pe?: number; volume?: number; enstrophy?: number }): void {
        const formatSci = (val: number | undefined) =>
            val !== undefined ? val.toExponential(3) : '-';

        const teEl = document.getElementById('diag-energy');
        const keEl = document.getElementById('diag-ke');
        const peEl = document.getElementById('diag-pe');
        const volumeEl = document.getElementById('diag-volume');
        const enstrophyEl = document.getElementById('diag-enstrophy');

        if (teEl) teEl.textContent = formatSci(diag.te);
        if (keEl) keEl.textContent = formatSci(diag.ke);
        if (peEl) peEl.textContent = formatSci(diag.pe);
        if (volumeEl) volumeEl.textContent = formatSci(diag.volume);
        if (enstrophyEl) enstrophyEl.textContent = formatSci(diag.enstrophy);
    }

    private updateProgress(step: number, t: number, message: string, percent: number): void {
        const stepEl = document.getElementById('step-count');
        const timeEl = document.getElementById('sim-time');
        const progressBar = document.getElementById('progress-bar');
        const progressMsg = document.getElementById('progress-message');

        if (stepEl) stepEl.textContent = `Step: ${step}`;
        if (timeEl) timeEl.textContent = `t = ${t.toFixed(2)} s`;
        if (progressBar) progressBar.style.width = `${percent}%`;
        if (progressMsg) progressMsg.textContent = message;
    }
}

// Initialize app when DOM is ready
document.addEventListener('DOMContentLoaded', () => {
    new ShelApp();
});

