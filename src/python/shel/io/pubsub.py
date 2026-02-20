import base64
import json

import numpy as np
import zmq


class SHELPublisher:
    """
    ZeroMQ PUB socket for SHEL solver state and diagnostics.
    """

    def __init__(self, port=5556):
        self.context = zmq.Context()
        self.socket = self.context.socket(zmq.PUB)
        self.socket.bind(f"tcp://*:{port}")

    def send(self, topic: str, payload: dict):
        msg = json.dumps(payload)
        self.socket.send_string(f"{topic} {msg}")

    @staticmethod
    def encode_array(arr: np.ndarray) -> dict:
        return {
            "shape": list(arr.shape),
            "dtype": str(arr.dtype),
            "data": base64.b64encode(arr.tobytes()).decode("ascii"),
        }

    def send_eta(self, t, eta, units="meters"):
        payload = {"t": t, **self.encode_array(eta), "units": units}
        self.send("state.eta", payload)

    def send_velocity(self, t, U, V, units="m/s"):
        payload = {
            "t": t,
            "U": self.encode_array(U),
            "V": self.encode_array(V),
            "units": units,
        }
        self.send("state.velocity", payload)

    def send_diag_global(self, t, energy, enstrophy, volume):
        payload = {
            "t": t,
            "energy": energy,
            "enstrophy": enstrophy,
            "volume": volume,
        }
        self.send("diag.global", payload)

    def send_diag_field(self, t, name, arr, units):
        payload = {"t": t, **self.encode_array(arr), "units": units}
        self.send(f"diag.field.{name}", payload)

    def send_progress(self, step, t, message, percent):
        payload = {
            "step": step,
            "t": t,
            "message": message,
            "percent": percent,
        }
        self.send("event.progress", payload)


class SHELSubscriber:
    """
    ZeroMQ SUB socket for SHEL GUI to receive state and diagnostics.
    """

    def __init__(self, port=5556, topics=None):
        self.context = zmq.Context()
        self.socket = self.context.socket(zmq.SUB)
        self.socket.connect(f"tcp://localhost:{port}")
        self.topics = topics or [
            "state.eta",
            "state.velocity",
            "diag.global",
            "diag.field.",
            "event.progress",
        ]
        for topic in self.topics:
            self.socket.setsockopt_string(zmq.SUBSCRIBE, topic)

    def recv(self):
        try:
            topic_msg = self.socket.recv_string(flags=zmq.NOBLOCK)
            topic, msg = topic_msg.split(" ", 1)
            payload = json.loads(msg)
            return topic, payload
        except zmq.Again:
            # No message available
            return None, None

    @staticmethod
    def decode_array(arr_dict):
        arr = np.frombuffer(base64.b64decode(arr_dict["data"]), dtype=arr_dict["dtype"])
        return arr.reshape(arr_dict["shape"])
