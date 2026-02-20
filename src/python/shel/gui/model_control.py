"""
Model control panel for the SHEL GUI.

This module provides the model control panel for the SHEL GUI, allowing users
to start, stop, and pause the model, as well as monitor its progress.
"""

import logging
import os
import sys
import tempfile
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple, Union, cast

import yaml
from PyQt5.QtCore import QProcess, Qt, QTimer
from PyQt5.QtWidgets import (
    QCheckBox,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

if TYPE_CHECKING:  # Avoid runtime import to prevent cyclic import with main_window
    from shel.gui.main_window import MainWindow  # noqa: F401

logger = logging.getLogger(__name__)


class ModelControlPanel(QWidget):
    """
    Panel for controlling the model.

    This class provides a panel for controlling the model, including:
    - Starting, stopping, and pausing the model
    - Monitoring model progress
    - Displaying model status

    Attributes:
        process: QProcess for running the model
        main_window: Main window of the application
    """

    def __init__(self, parent: Optional[QWidget] = None) -> None:
        """Initialize the model control panel."""
        super().__init__(parent)

        # Store reference to main window (typed only for checkers; avoid runtime import)
        self.main_window: Optional["MainWindow"] = cast("MainWindow", parent)

        # Set up the layout
        layout = QVBoxLayout(self)

        # Add status label
        self.status_label = QLabel("Model Status: Ready")
        layout.addWidget(self.status_label)

        # Add progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        layout.addWidget(self.progress_bar)

        # Add control buttons
        button_layout = QHBoxLayout()
        layout.addLayout(button_layout)

        # Run button
        self.run_button = QPushButton("Run Model")
        self.run_button.clicked.connect(self.run_model)
        button_layout.addWidget(self.run_button)

        # Stop button
        self.stop_button = QPushButton("Stop")
        self.stop_button.clicked.connect(self.stop_model)
        self.stop_button.setEnabled(False)
        button_layout.addWidget(self.stop_button)

        # Initialize QProcess
        self.process: Optional[QProcess] = None

        # Initialize attributes
        self.current_step: int = 0
        self.total_steps: int = 0
        self.running: bool = False

    def run_model(self) -> None:
        """Run the model."""

        if not self.main_window:
            logger.error("Main window reference is not set.")
            return

        # Get configuration from the parameter panel
        config: Dict[str, Any] = self.main_window.parameter_panel.get_config()

        # Update total steps
        self.total_steps = config["model"]["num_steps"]

        # Create temporary configuration file
        temp_config_file = tempfile.NamedTemporaryFile(suffix=".yaml", delete=False)
        temp_config_path: str = temp_config_file.name

        try:
            # Write configuration to temporary file
            with open(temp_config_path, "w") as f:
                yaml.dump(config, f, default_flow_style=False)

            # Build command
            cmd = [
                sys.executable,  # Python executable
                os.path.join(
                    os.path.dirname(
                        os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
                    ),
                    "run.py",
                ),
                "--config",
                temp_config_path,
                "--verbose",
            ]

            # Create output directory if it doesn't exist
            output_dir = config["output"]["directory"]
            os.makedirs(output_dir, exist_ok=True)
            cmd.extend(["--output-dir", output_dir])

            # Set ZeroMQ port
            zmq_port: int = config["communication"]["zmq_pub_port"]
            cmd.extend(["--pub-port", str(zmq_port)])

            # Disable output if requested
            if not config["output"]["enabled"]:
                cmd.append("--no-output")

            # Create QProcess if not already created
            if self.process is None:
                self.process = QProcess()
                # Merge stdout and stderr so all logs are captured uniformly
                self.process.setProcessChannelMode(cast(Any, QProcess).MergedChannels)
                self.process.readyReadStandardOutput.connect(self.handle_stdout)
                self.process.readyReadStandardError.connect(self.handle_stderr)
                self.process.finished.connect(self.handle_finished)

            # Start the process
            logger.info("Starting model with command: %s", " ".join(cmd))
            if self.process is not None:
                self.process.start(cmd[0], cmd[1:])

            # Update UI
            self.run_button.setEnabled(False)
            self.stop_button.setEnabled(True)
            self.status_label.setText("Model Status: Running")
            self.running = True

            # Reset progress
            self.current_step = 0
            self.progress_bar.setValue(0)

        except Exception as e:
            logger.error("Failed to start model: %s", e)
            QMessageBox.critical(self, "Error", f"Failed to start model: {str(e)}")

            # Clean up
            os.unlink(temp_config_path)

    def stop_model(self) -> None:
        """Stop the model."""
        if self.process and self.process.state() != cast(Any, QProcess).NotRunning:
            # Confirm with user
            reply = QMessageBox.question(
                self,
                "Stop Model",
                "Are you sure you want to stop the model?",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No,
            )

            if reply == QMessageBox.Yes:
                # Terminate the process
                logger.info("Stopping model")
                self.process.terminate()

                # Give it some time to terminate gracefully
                if not self.process.waitForFinished(3000):
                    logger.warning(
                        "Model did not terminate gracefully, killing process"
                    )
                    self.process.kill()

    def handle_stdout(self) -> None:
        """Handle standard output from the model process."""
        if self.process:
            output: str = self.process.readAllStandardOutput().data().decode()
            # Emit each CLI line at INFO so it shows up by default
            for line in output.splitlines():
                if line.strip():
                    logger.info("CLI: %s", line)

                # Look for progress information
                if "Step" in line and "/" in line:
                    try:
                        parts: List[str] = line.split("/")
                        current_step: int = int(parts[0].split("Step")[-1].strip())
                        total_steps: int = int(parts[1].split()[0])

                        # Update progress
                        self.current_step = current_step
                        self.total_steps = total_steps
                        progress: int = min(100, int(current_step / total_steps * 100))
                        self.progress_bar.setValue(progress)
                    except Exception as e:
                        logger.error("Failed to parse progress: %s", e)

    def handle_stderr(self) -> None:
        """Handle standard error from the model process."""
        if self.process:
            error: str = self.process.readAllStandardError().data().decode()
            logger.error("Model error: %s", error)

    def handle_finished(self, exit_code: int, exit_status: QProcess.ExitStatus) -> None:
        """
        Handle model process completion.

        Args:
            exit_code: Exit code of the process
            exit_status: Exit status of the process
        """
        logger.info("Model process finished with exit code %s", exit_code)

        # Update UI
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.running = False

        if exit_code == 0:
            self.status_label.setText("Model Status: Completed")
            self.progress_bar.setValue(100)
        else:
            self.status_label.setText(f"Model Status: Error (code {exit_code})")
            QMessageBox.warning(
                self, "Model Error", f"The model process exited with code {exit_code}."
            )

    def update_status(self, message: Dict[str, Any]) -> None:
        """
        Update status based on model message.

        Args:
            message: Message from the model
        """
        if not self.running:
            return

        # Update progress based on message
        step: int = message.get("step", 0)
        status: str = message.get("status", "running")

        if status == "running":
            # Update progress
            if self.total_steps > 0:
                progress: int = min(100, int(step / self.total_steps * 100))
                self.progress_bar.setValue(progress)

            # Update status label
            self.status_label.setText(
                f"Model Status: Running (Step {step}/{self.total_steps})"
            )

        elif status == "complete":
            # Update progress to 100%
            self.progress_bar.setValue(100)

            # Update status label
            self.status_label.setText("Model Status: Completed")

    def handle_completion(self) -> None:
        """Handle model completion."""
        # Enable run button
        self.run_button.setEnabled(True)

        # Disable stop button
        self.stop_button.setEnabled(False)

        # Update status
        self.running = False
