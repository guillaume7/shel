/**
 * Diagnostics time series plot using Plotly.
 */

import Plotly from 'plotly.js-dist-min';

interface DiagnosticData {
    time: number[];
    te: number[];
    ke: number[];
    pe: number[];
    volume: number[];
    enstrophy: number[];
}

export class DiagPlot {
    private plotDiv: HTMLElement;
    private data: DiagnosticData = {
        time: [],
        te: [],
        ke: [],
        pe: [],
        volume: [],
        enstrophy: [],
    };
    private maxPoints = 1000;

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
            title: 'Simulation Diagnostics',
            xaxis: { title: 'Time (s)', autorange: true },
            yaxis: {
                title: 'Energy (J)',
                type: 'log',
                autorange: true,
                titlefont: { color: '#3b82f6' },
                tickfont: { color: '#3b82f6' },
            },
            yaxis2: {
                title: 'Volume (m³)',
                type: 'linear',
                overlaying: 'y',
                side: 'right',
                rangemode: 'tozero',
                autorange: true,
                titlefont: { color: '#10b981' },
                tickfont: { color: '#10b981' },
                showgrid: false,
            },
            yaxis3: {
                title: 'Enstrophy (s⁻²)',
                type: 'log',
                overlaying: 'y',
                side: 'right',
                anchor: 'free',
                position: 0.1, // Move slightly inside
                autorange: true,
                titlefont: { color: '#f59e0b' },
                tickfont: { color: '#f59e0b' },
                showgrid: false,
            },
            paper_bgcolor: '#1e293b',
            plot_bgcolor: '#0f172a',
            font: { color: '#e2e8f0' },
            margin: { l: 60, r: 120, t: 50, b: 60 },
            showlegend: true,
            legend: {
                font: { color: '#e2e8f0' },
                orientation: 'h',
                y: -0.2,
                x: 0,
            },
        };

        const traces: any[] = [
            {
                x: [], y: [],
                name: 'Total Energy',
                type: 'scatter',
                mode: 'lines',
                line: { color: '#3b82f6', width: 2 },
                yaxis: 'y1',
            },
            {
                x: [], y: [],
                name: 'Kinetic',
                type: 'scatter',
                mode: 'lines',
                line: { color: '#60a5fa', dash: 'dash' },
                yaxis: 'y1',
            },
            {
                x: [], y: [],
                name: 'Potential',
                type: 'scatter',
                mode: 'lines',
                line: { color: '#93c5fd', dash: 'dot' },
                yaxis: 'y1',
            },
            {
                x: [], y: [],
                name: 'Volume',
                type: 'scatter',
                mode: 'lines',
                line: { color: '#10b981' },
                yaxis: 'y2',
            },
            {
                x: [], y: [],
                name: 'Enstrophy',
                type: 'scatter',
                mode: 'lines',
                line: { color: '#f59e0b' },
                yaxis: 'y3',
            },
        ];

        Plotly.newPlot(this.plotDiv, traces, layout, { responsive: true });
    }

    addPoint(t: number, diagnostics: { te?: number; ke?: number; pe?: number; volume?: number; enstrophy?: number }): void {
        this.data.time.push(t);
        this.data.te.push(diagnostics.te ?? 0);
        this.data.ke.push(diagnostics.ke ?? 0);
        this.data.pe.push(diagnostics.pe ?? 0);
        this.data.volume.push(diagnostics.volume ?? 0);
        this.data.enstrophy.push(diagnostics.enstrophy ?? 0);

        // Limit data points
        if (this.data.time.length > this.maxPoints) {
            this.data.time.shift();
            this.data.te.shift();
            this.data.ke.shift();
            this.data.pe.shift();
            this.data.volume.shift();
            this.data.enstrophy.shift();
        }

        this.update();
    }

    private update(): void {
        // Use restyle for data
        Plotly.restyle(this.plotDiv, {
            x: [this.data.time, this.data.time, this.data.time, this.data.time, this.data.time],
            y: [this.data.te, this.data.ke, this.data.pe, this.data.volume, this.data.enstrophy]
        }, [0, 1, 2, 3, 4]);
    }


    clear(): void {
        this.data = {
            time: [],
            te: [],
            ke: [],
            pe: [],
            volume: [],
            enstrophy: [],
        };
        this.update();
    }
}
