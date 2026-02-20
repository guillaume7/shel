/**
 * Water elevation (eta) plot using Plotly.
 */

import Plotly from 'plotly.js-dist-min';

export class EtaPlot {
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
            title: 'Water Elevation (η)',
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
            colorscale: 'RdBu',
            reversescale: true,
            colorbar: {
                title: 'η (m)',
                titlefont: { color: '#e2e8f0' },
                tickfont: { color: '#e2e8f0' },
            },
        }];

        Plotly.newPlot(this.plotDiv, data, layout, { responsive: true });
    }

    update(eta: number[][], t: number): void {
        const data: any[] = [{
            z: eta,
            type: 'heatmap',
            colorscale: 'RdBu',
            reversescale: true,
            zauto: false,
            zmin: -0.1,
            zmax: 0.1,
            colorbar: {
                title: 'η (m)',
                titlefont: { color: '#e2e8f0' },
                tickfont: { color: '#e2e8f0' },
            },
        }];

        const layout: any = {
            title: `Water Elevation (η) - t = ${t.toFixed(2)} s`,
        };

        Plotly.react(this.plotDiv, data, layout);
    }

    clear(): void {
        Plotly.restyle(this.plotDiv, { z: [[]] }, 0);
        Plotly.relayout(this.plotDiv, { title: 'Water Elevation (η)' });
    }
}
