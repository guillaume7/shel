/**
 * Velocity magnitude plot using Plotly.
 */

import Plotly from 'plotly.js-dist-min';

export class VelocityPlot {
    private plotDiv: HTMLElement;

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
            title: 'Velocity Magnitude',
            xaxis: { title: 'X' },
            yaxis: { title: 'Y' },
            paper_bgcolor: '#1e293b',
            plot_bgcolor: '#0f172a',
            font: { color: '#e2e8f0' },
            margin: { l: 50, r: 50, t: 50, b: 50 },
        };

        const data: any[] = [{
            z: [[0]],
            type: 'heatmap',
            colorscale: 'Viridis',
            colorbar: {
                title: 'Speed (m/s)',
                titlefont: { color: '#e2e8f0' },
                tickfont: { color: '#e2e8f0' },
            },
        }];

        Plotly.newPlot(this.plotDiv, data, layout, { responsive: true });
    }

    update(U: number[][], V: number[][], t: number): void {
        // Calculate magnitude
        const ny = U.length;
        const nx = U[0].length;
        const magnitude: number[][] = [];

        for (let i = 0; i < ny; i++) {
            const row: number[] = [];
            for (let j = 0; j < nx; j++) {
                const u = U[i][j];
                const v = V[i][j];
                row.push(Math.sqrt(u * u + v * v));
            }
            magnitude.push(row);
        }

        const data: any[] = [{
            z: magnitude,
            type: 'heatmap',
            colorscale: 'Viridis',
            zauto: false,
            zmin: 0.0,
            zmax: 0.003,
            colorbar: {
                title: 'Speed (m/s)',
                titlefont: { color: '#e2e8f0' },
                tickfont: { color: '#e2e8f0' },
            },
        }];

        const layout: any = {
            title: `Velocity Magnitude - t = ${t.toFixed(2)} s`,
        };

        Plotly.react(this.plotDiv, data, layout);
    }

    clear(): void {
        Plotly.restyle(this.plotDiv, { z: [[]] }, 0);
        // Reset title to initial state
        Plotly.relayout(this.plotDiv, { title: 'Velocity Magnitude' });
    }
}
