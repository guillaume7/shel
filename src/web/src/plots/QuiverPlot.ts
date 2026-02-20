/**
 * Velocity vector (Quiver) plot using Plotly.
 */

import Plotly from 'plotly.js-dist-min';

export class QuiverPlot {
    private plotDiv: HTMLElement;
    private maxArrows = 400; // Limit number of arrows for performance

    constructor(divId: string) {
        const el = document.getElementById(divId);
        if (!el) {
            throw new Error(`Element ${divId} not found`);
        }
        this.plotDiv = el;
        this.initialize();
    }

    private initialize(): void {
        const layout: any = {
            title: 'Velocity Vectors',
            xaxis: { title: 'X', showgrid: false, zeroline: false },
            yaxis: { title: 'Y', showgrid: false, zeroline: false },
            paper_bgcolor: '#1e293b',
            plot_bgcolor: '#0f172a',
            font: { color: '#e2e8f0' },
            margin: { l: 50, r: 50, t: 50, b: 50 },
        };

        const data: any[] = [{
            type: 'scatter',
            mode: 'lines+markers',
            x: [],
            y: [],
            line: { color: '#60a5fa', width: 1.5 },
            marker: { size: 2, color: '#60a5fa' },
            name: 'Flow'
        }];

        Plotly.newPlot(this.plotDiv, data, layout, { responsive: true });
    }

    update(U: number[][], V: number[][], t: number): void {
        const ny = U.length;
        const nx = U[0].length;

        // Downsample for performance
        const skip = Math.max(1, Math.ceil(Math.sqrt((nx * ny) / this.maxArrows)));

        const x: (number | null)[] = [];
        const y: (number | null)[] = [];

        // Huge scale factor for visibility
        // Typical velocity is ~0.005, we want it to look like ~5 pixels
        const scale = 1000.0;

        for (let i = 0; i < ny; i += skip) {
            for (let j = 0; j < nx; j += skip) {
                const u = U[i][j];
                const v = V[i][j];

                // Skip very small velocities
                if (Math.abs(u) < 1e-7 && Math.abs(v) < 1e-7) continue;

                // Start of arrow
                x.push(j);
                y.push(i);

                // End of arrow
                x.push(j + u * scale);
                y.push(i + v * scale);

                // Gap between arrows
                x.push(null);
                y.push(null);
            }
        }

        const data: any[] = [{
            type: 'scatter',
            mode: 'lines+markers',
            x: x,
            y: y,
            line: { color: '#60a5fa', width: 1.5 },
            marker: { size: 2, color: '#60a5fa' },
            name: 'Flow'
        }];

        const layout: any = {
            title: `Velocity Vectors - t = ${t.toFixed(2)} s`,
            xaxis: { range: [0, nx], autorange: false },
            yaxis: { range: [0, ny], autorange: false }
        };

        Plotly.react(this.plotDiv, data, layout);
    }

    clear(): void {
        const data: any[] = [{
            type: 'scatter',
            mode: 'lines+markers',
            x: [],
            y: [],
            line: { color: '#60a5fa', width: 1.5 },
            marker: { size: 2, color: '#60a5fa' },
            name: 'Flow'
        }];
        Plotly.react(this.plotDiv, data, { title: 'Velocity Vectors' });
    }
}
