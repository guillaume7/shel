"""
Integration tests for GUI-solver communication via ZeroMQ.
"""

import logging
import multiprocessing
import os
import tempfile
import time
from multiprocessing import Process
from typing import Any, Dict

import pytest
import yaml

from shel.io.pubsub import SHELSubscriber
from shel.model.model_runner import ModelRunner

# Set up module-level logger
logger = logging.getLogger(__name__)

# Set multiprocessing start method to 'spawn' for better compatibility
multiprocessing.set_start_method("spawn", force=True)


def run_solver_process(config_file: str) -> None:
    """
    Run the solver in a separate process.

    Args:
        config_file: Path to the YAML configuration file
    """

    # Set up basic logging for the solver process
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    try:
        # Load config and run model
        with open(config_file, "r") as f:
            config = yaml.safe_load(f)

        logger.info("Starting solver process")
        print(f"Process starting with config: {config}")

        runner = ModelRunner(config)
        print("ModelRunner created")

        # Give time for publisher to bind
        time.sleep(0.5)
        print("About to initialize ModelRunner...")

        runner.initialize()
        print("ModelRunner initialized")

        print("About to run ModelRunner...")
        runner.run()
        print("ModelRunner completed")

        # Keep process alive a bit longer to ensure messages are sent
        time.sleep(1.0)
        print("Process finishing")
        logger.info("Solver process completed")

    except Exception as e:
        logger.error("Solver process failed: %s", e)
        print(f"ERROR in solver process: {e}")
        import traceback

        traceback.print_exc()
        raise


def create_test_config() -> Dict[str, Any]:
    """Create a minimal test configuration for integration testing."""
    return {
        "grid": {"nx": 20, "ny": 20, "dx": 100.0, "dy": 100.0},
        "model": {
            "timestep": 0.5,
            "num_steps": 10,  # Short run for testing
            "output_interval": 2,
            "solver": "leapfrog",
            "viscosity": 0.0,
            "bottom_drag_coef": 0.0,
            "coriolis_parameter": 0.0,
            "gravity": 9.81,
        },
        "initial_conditions": {"type": "flat", "elevation": 0.1},
        "bathymetry": {"type": "flat", "flat": {"depth": 50.0}},
        "boundary_conditions": {
            "north": "closed",
            "south": "closed",
            "east": "closed",
            "west": "closed",
        },
        "output": {
            "directory": "/tmp/shel_test_output",
            "enabled": False,  # Disable file output for faster testing
        },
        "communication": {
            "zmq_pub_port": 5557,  # Use different port to avoid conflicts
            "enable_zmq": True,
        },
    }


class TestGUIIntegration:
    """Test GUI-solver integration via ZeroMQ."""

    def test_pubsub_communication(self):
        """Test basic pub/sub communication without GUI."""
        config = create_test_config()

        # Create temporary config file
        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            yaml.dump(config, f)
            config_file = f.name

        try:
            # Set up subscriber
            subscriber = SHELSubscriber(port=5557)

            # Start solver process
            solver_process = Process(target=run_solver_process, args=(config_file,))
            solver_process.start()

            # Allow more time for solver to start and bind
            print("TEST: Waiting for solver to start...")
            time.sleep(2.0)

            # Collect messages for a longer time
            messages_received = []
            start_time = time.time()
            timeout = 20.0  # 20 second timeout
            print(f"TEST: Starting message collection with {timeout}s timeout...")

            while time.time() - start_time < timeout:
                try:
                    topic, payload = subscriber.recv()
                    if topic and payload:
                        messages_received.append((topic, payload))
                        logger.info("Received message: topic=%s", topic)
                        print(f"TEST: Received message: {topic}")

                        # Stop collecting after we get a completion message
                        if (
                            topic == "event.progress"
                            and payload.get("percent", 0) >= 100
                        ):
                            print(
                                "TEST: Received completion message, stopping collection"
                            )
                            break
                    else:
                        print("TEST: recv() returned empty topic/payload")

                except Exception as e:
                    logger.debug("No message received: %s", e)
                    print(f"TEST: Exception during recv: {e}")

                time.sleep(0.1)

            # Wait for solver to finish
            print("TEST: Waiting for solver process to finish...")
            solver_process.join(timeout=10.0)
            if solver_process.is_alive():
                print("TEST: Solver process still alive, terminating...")
                solver_process.terminate()
                solver_process.join()
            else:
                print("TEST: Solver process finished normally")
            print(f"TEST: Solver exit code: {solver_process.exitcode}")

            # Verify we received expected message types
            topics_received = {topic for topic, _ in messages_received}
            print(f"TEST: Total messages received: {len(messages_received)}")
            print(f"TEST: Topics received: {topics_received}")

            for i, (topic, payload) in enumerate(messages_received):
                print(f"TEST: Message {i+1}: {topic} -> {type(payload)}")

            assert (
                len(messages_received) > 0
            ), f"No messages received from solver. Process exit code: {solver_process.exitcode}"
            assert "state.eta" in topics_received, "No eta state messages received"
            assert "event.progress" in topics_received, "No progress messages received"

            # Verify message structure
            for topic, payload in messages_received:
                assert isinstance(payload, dict), f"Expected dict payload for {topic}"
                print(f"TEST: Message structure for {topic}: {list(payload.keys())}")

                # Check expected fields based on actual SHEL protocol
                if topic == "state.eta":
                    assert "t" in payload, f"Missing time field in {topic}"
                    assert "data" in payload, f"Missing data field in {topic}"
                    assert "dtype" in payload, f"Missing dtype field in {topic}"
                    assert "shape" in payload, f"Missing shape field in {topic}"
                elif topic == "state.velocity":
                    assert "t" in payload, f"Missing time field in {topic}"
                    assert "U" in payload, f"Missing U field in {topic}"
                    assert "V" in payload, f"Missing V field in {topic}"
                elif topic == "diag.global":
                    assert "t" in payload, f"Missing time field in {topic}"
                    # Global diagnostics have energy, enstrophy, volume fields
                elif topic == "event.progress":
                    assert (
                        "percent" in payload
                    ), "Missing percent field in progress event"
                    assert (
                        "message" in payload
                    ), "Missing message field in progress event"

            logger.info(
                "Pub/sub communication test PASSED: %d messages received",
                len(messages_received),
            )

        finally:
            # Clean up
            os.unlink(config_file)

    @pytest.mark.skip(reason="Requires full Qt application - skip in CI")
    def test_gui_solver_integration(self):
        """Test full GUI-solver integration (interactive test)."""
        # This would test the actual GUI interaction
        # Skipped for automated testing
        pass


if __name__ == "__main__":
    # Allow running this test module directly
    pytest.main([__file__, "-v"])
