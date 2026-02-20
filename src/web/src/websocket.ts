/**
 * WebSocket client for SHEL real-time updates.
 */

export interface StateMessage {
    topic: string;
    t: number;
    [key: string]: any;
}

export type MessageHandler = (message: StateMessage) => void;

export class WebSocketClient {
    private ws: WebSocket | null = null;
    private handlers: Map<string, MessageHandler[]> = new Map();
    private reconnectAttempts = 0;
    private maxReconnectAttempts = 5;
    private reconnectDelay = 1000;

    constructor(private url: string) { }

    connect(): void {
        try {
            this.ws = new WebSocket(this.url);

            this.ws.onopen = () => {
                console.log('WebSocket connected');
                this.reconnectAttempts = 0;
                this.updateConnectionStatus(true);
            };

            this.ws.onmessage = (event) => {
                try {
                    const message: StateMessage = JSON.parse(event.data);
                    this.handleMessage(message);
                } catch (error) {
                    console.error('Failed to parse message:', error);
                }
            };

            this.ws.onerror = (error) => {
                console.error('WebSocket error:', error);
            };

            this.ws.onclose = () => {
                console.log('WebSocket disconnected');
                this.updateConnectionStatus(false);
                this.attemptReconnect();
            };
        } catch (error) {
            console.error('Failed to connect:', error);
            this.attemptReconnect();
        }
    }

    disconnect(): void {
        if (this.ws) {
            this.ws.close();
            this.ws = null;
        }
    }

    on(topic: string, handler: MessageHandler): void {
        if (!this.handlers.has(topic)) {
            this.handlers.set(topic, []);
        }
        this.handlers.get(topic)!.push(handler);
    }

    off(topic: string, handler: MessageHandler): void {
        const handlers = this.handlers.get(topic);
        if (handlers) {
            const index = handlers.indexOf(handler);
            if (index > -1) {
                handlers.splice(index, 1);
            }
        }
    }

    private handleMessage(message: StateMessage): void {
        const handlers = this.handlers.get(message.topic);
        if (handlers) {
            handlers.forEach((handler) => handler(message));
        }
    }

    private attemptReconnect(): void {
        if (this.reconnectAttempts < this.maxReconnectAttempts) {
            this.reconnectAttempts++;
            console.log(`Reconnecting... (attempt ${this.reconnectAttempts})`);
            setTimeout(() => this.connect(), this.reconnectDelay * this.reconnectAttempts);
        } else {
            console.error('Max reconnect attempts reached');
        }
    }

    private updateConnectionStatus(connected: boolean): void {
        const statusEl = document.getElementById('connection-status');
        if (statusEl) {
            statusEl.textContent = connected ? 'Connected' : 'Disconnected';
            if (connected) {
                statusEl.classList.add('connected');
            } else {
                statusEl.classList.remove('connected');
            }
        }
    }
}
