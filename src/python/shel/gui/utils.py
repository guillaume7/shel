"""
Utility functions and classes for the GUI.
"""

import logging
from typing import Any, Dict

import zmq
from PyQt5.QtCore import QThread, pyqtSignal

logger = logging.getLogger(__name__)


class MessageSubscriber(QThread):
    """
    Thread for subscribing to ZeroMQ messages from the model runner.

    This thread subscribes to ZeroMQ messages from the model runner and
    emits a signal when a message is received.

    Attributes:
        message_received: Signal emitted when a message is received
    """

    message_received = pyqtSignal(dict)

    def __init__(self, port=5556, parent=None):
        """
        Initialize the message subscriber.

        Args:
            port: Port to subscribe to
            parent: Parent object
        """
        super().__init__(parent)

        self.port = port
        self.running = False

        # ZeroMQ context and socket
        self.context = None
        self.socket = None

        logger.info("MessageSubscriber initialized with port %s", port)

    def run(self):
        """Run the subscriber thread."""
        self.running = True

        # Initialize ZeroMQ
        self.context = zmq.Context()
        self.socket = self.context.socket(zmq.SUB)
        self.socket.connect(f"tcp://localhost:{self.port}")
        self.socket.setsockopt_string(zmq.SUBSCRIBE, "")

        logger.info("Connected to ZeroMQ publisher on port %s", self.port)

        # Set up polling to allow for thread termination
        poller = zmq.Poller()
        poller.register(self.socket, zmq.POLLIN)

        while self.running:
            # Poll for messages with timeout
            socks = dict(poller.poll(100))  # 100ms timeout

            if self.socket in socks and socks[self.socket] == zmq.POLLIN:
                # Receive and process message
                try:
                    message = self.socket.recv_json()
                    self.message_received.emit(message)
                    logger.debug("Received message: %s", message)
                except Exception as e:
                    logger.error("Error receiving message: %s", e)

        # Clean up
        self.socket.close()
        self.context.term()
        logger.info("MessageSubscriber stopped")

    def stop(self):
        """Stop the subscriber thread."""
        self.running = False
        logger.info("Stopping MessageSubscriber thread...")
