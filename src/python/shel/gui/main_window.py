#!/usr/bin/env python
"""
Main GUI for SHEL.

This module provides the main GUI for the SHEL model, allowing users
to interact with the model, visualize results, and set parameters.
"""
import logging
import sys
from functools import partial
from typing import Optional, cast

# Matplotlib widgets are handled within PlotManager; no direct imports needed here
from PyQt5.QtCore import QSettings, QTimer
from PyQt5.QtGui import QCloseEvent
from PyQt5.QtWidgets import (
    QAction,
    QApplication,
    QFileDialog,
    QHBoxLayout,
    QMainWindow,
    QMenu,
    QMenuBar,
    QMessageBox,
    QStatusBar,
    QVBoxLayout,
    QWidget,
)

from shel.gui.model_control import ModelControlPanel
from shel.gui.parameters import ParameterPanel
from shel.gui.visualization import PlotManager
from shel.io.pubsub import SHELSubscriber

logger = logging.getLogger(__name__)


class MainWindow(QMainWindow):
    """
    Main window for the SHEL GUI.

    This class provides the main window for the SHEL GUI, including:
    - Visualization panels
    - Parameter panels
    - Model control
    - Menus

    Attributes:
        plot_manager: Manager for plots
        parameter_panel: Panel for editing parameters
        model_control: Panel for controlling the model
        message_subscriber: Thread for receiving ZMQ messages
    """

    def __init__(self, parent=None):
        """Initialize the main window."""
        super().__init__(parent)

        # Set window properties
        self.setWindowTitle("SHEL - SHallow-water numerical modEL")
        self.setMinimumSize(1200, 800)

        # Set up the central widget
        central_widget = QWidget()
        self.setCentralWidget(central_widget)

        # Set up the main layout
        main_layout = QHBoxLayout(central_widget)

        # Set up the plot area
        plot_layout = QVBoxLayout()
        main_layout.addLayout(plot_layout, 3)  # Plot area takes 3/4 of the space

        # Initialize plot manager
        self.plot_manager = PlotManager(self)
        plot_layout.addWidget(self.plot_manager)

        # Set up parameter and control panels
        control_layout = QVBoxLayout()
        main_layout.addLayout(control_layout, 1)  # Control area takes 1/4 of the space

        # Add parameter panel
        self.parameter_panel = ParameterPanel(self)
        control_layout.addWidget(self.parameter_panel)

        self.model_control = ModelControlPanel(self)
        control_layout.addWidget(self.model_control)

        # Add status bar
        self.status_bar = cast("QStatusBar", self.statusBar())
        self.status_bar.showMessage("Ready")

        # Set up ZeroMQ subscriber for model updates
        self.subscriber = SHELSubscriber(port=5556)
        self._start_subscriber_timer()

        # Create menus
        self.create_menus()

        # Load settings
        self.load_settings()

    def _on_exit_triggered(self):
        """Handle exit action triggered."""
        self.close()

    def create_menus(self):
        """Create menu bar with actions."""
        # Get menu bar with proper typing
        self.menu_bar = cast("QMenuBar", self.menuBar())

        # File menu
        file_menu = cast("QMenu", self.menu_bar.addMenu("&File"))

        # New configuration
        new_action = QAction("&New", self)
        new_action.setShortcut("Ctrl+N")
        new_action.triggered.connect(self.new_configuration)
        file_menu.addAction(new_action)

        # Open configuration
        open_action = QAction("&Open...", self)
        open_action.setShortcut("Ctrl+O")
        open_action.triggered.connect(self.open_configuration)
        file_menu.addAction(open_action)

        # Save configuration
        save_action = QAction("&Save", self)
        save_action.setShortcut("Ctrl+S")
        save_action.triggered.connect(self.save_configuration)
        file_menu.addAction(save_action)

        # Save As configuration
        save_as_action = QAction("Save &As...", self)
        save_as_action.setShortcut("Ctrl+Shift+S")
        save_as_action.triggered.connect(self.save_configuration_as)
        file_menu.addAction(save_as_action)

        file_menu.addSeparator()

        # Export menu
        export_menu = cast("QMenu", file_menu.addMenu("&Export"))

        # Export plot as image
        export_image_action = QAction("&Image...", self)
        export_image_action.triggered.connect(self.export_image)
        export_menu.addAction(export_image_action)

        # Export animation
        export_animation_action = QAction("&Animation...", self)
        export_animation_action.triggered.connect(self.export_animation)
        export_menu.addAction(export_animation_action)

        file_menu.addSeparator()

        # Exit action
        exit_action = QAction("E&xit", self)
        exit_action.setShortcut("Ctrl+Q")
        exit_action.triggered.connect(self._on_exit_triggered)
        file_menu.addAction(exit_action)

        # View menu
        view_menu = cast("QMenu", self.menu_bar.addMenu("&View"))

        # Toggle parameter panel
        toggle_params_action = QAction("&Parameters", self)
        toggle_params_action.setCheckable(True)
        toggle_params_action.setChecked(True)
        toggle_params_action.triggered.connect(self.parameter_panel.setVisible)
        view_menu.addAction(toggle_params_action)

        # Toggle control panel
        toggle_control_action = QAction("&Control Panel", self)
        toggle_control_action.setCheckable(True)
        toggle_control_action.setChecked(True)
        toggle_control_action.triggered.connect(self.model_control.setVisible)
        view_menu.addAction(toggle_control_action)

        # Plot type submenu
        plot_type_menu = cast("QMenu", view_menu.addMenu("&Plot Type"))

        for plot_type in self.plot_manager.available_plots:
            plot_action = QAction(plot_type, self)
            plot_action.setCheckable(True)
            if plot_type == self.plot_manager.current_plot_type:
                plot_action.setChecked(True)
            plot_action.triggered.connect(
                partial(self.plot_manager.set_plot_type, plot_type)
            )
            plot_type_menu.addAction(plot_action)

        # Help menu
        help_menu = cast("QMenu", self.menu_bar.addMenu("&Help"))

        # About action
        about_action = QAction("&About", self)
        about_action.triggered.connect(self.show_about)
        help_menu.addAction(about_action)

    def _start_subscriber_timer(self):
        # Poll for messages every 100ms
        self.timer = QTimer(self)
        self.timer.timeout.connect(self._poll_pubsub)
        self.timer.start(100)

    def _poll_pubsub(self):
        """Poll for pubsub messages."""
        if self.subscriber:
            try:
                topic, payload = self.subscriber.recv()
                if topic and payload:
                    self.handle_model_update(topic, payload)
            except Exception as e:
                logger.error(f"Error receiving pubsub message: {e}")

    def handle_model_update(self, topic, payload):
        """
        Handle model update message from ZeroMQ pubsub.
        Args:
            topic: Topic string
            payload: Decoded message dict
        """
        logger.debug(f"Received topic: {topic}, payload keys: {list(payload.keys())}")

        if topic.startswith("state."):
            # Update plots with state data
            self.plot_manager.update_plots(payload)
        elif topic == "diag.global":
            # Update diagnostics display (TODO: implement diagnostics panel)
            logger.info(f"Global diagnostics: {payload}")
        elif topic.startswith("diag.field."):
            # Update field diagnostics (TODO: implement field diagnostics)
            logger.debug(f"Field diagnostics: {topic} - {payload}")
        elif topic == "event.progress":
            # Update progress and status
            self.model_control.update_status(payload)
        else:
            logger.warning(f"Unknown topic: {topic}")

        # If model is complete, handle completion
        if payload.get("status") == "complete":
            self.handle_model_completion()

    def handle_model_completion(self):
        """Handle model completion."""
        # Update UI state
        self.model_control.handle_completion()

        # Show completion message
        QMessageBox.information(
            self, "Model Complete", "The model run has completed successfully."
        )

    def new_configuration(self):
        """Create a new configuration."""
        # Check if current configuration should be saved
        if not self.check_save_current():
            return

        # Reset parameters to defaults
        self.parameter_panel.reset_to_defaults()

        # Clear current file path
        self.current_config_path = None

        # Update status
        self.status_bar.showMessage("New configuration created")

    def open_configuration(self):
        """Open a configuration file."""
        # Check if current configuration should be saved
        if not self.check_save_current():
            return

        # Show file dialog
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Open Configuration", "", "YAML Files (*.yaml);;All Files (*)"
        )

        if file_path:
            try:
                # Load configuration
                self.parameter_panel.load_from_file(file_path)

                # Update current path
                self.current_config_path = file_path

                # Update status
                self.status_bar.showMessage(f"Loaded configuration from {file_path}")
            except Exception as e:
                QMessageBox.critical(
                    self, "Error", f"Failed to load configuration: {str(e)}"
                )

    def save_configuration(self):
        """Save the current configuration."""
        if self.current_config_path:
            self.parameter_panel.save_to_file(self.current_config_path)
            self.status_bar.showMessage(
                f"Saved configuration to {self.current_config_path}"
            )
        else:
            self.save_configuration_as()

    def save_configuration_as(self):
        """Save the current configuration to a new file."""
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save Configuration", "", "YAML Files (*.yaml);;All Files (*)"
        )

        if file_path:
            self.parameter_panel.save_to_file(file_path)
            self.current_config_path = file_path
            self.status_bar.showMessage(f"Saved configuration to {file_path}")

    def export_image(self):
        """Export the current plot as an image."""
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Export Image",
            "",
            "PNG Files (*.png);;JPEG Files (*.jpg);;PDF Files (*.pdf);;EPS Files (*.eps);;All Files (*)",
        )

        if file_path:
            self.plot_manager.export_current_plot(file_path)
            self.status_bar.showMessage(f"Exported image to {file_path}")

    def export_animation(self):
        """Export an animation of the model run."""
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Export Animation",
            "",
            "MP4 Files (*.mp4);;GIF Files (*.gif);;All Files (*)",
        )

        if file_path:
            try:
                # Export animation with default settings (10 fps)
                self.plot_manager.export_animation(file_path, fps=10)
                self.status_bar.showMessage(f"Exported animation to {file_path}")
            except Exception as e:
                QMessageBox.warning(
                    self,
                    "Export Failed",
                    f"Failed to export animation: {str(e)}\n\n"
                    "Make sure you have ffmpeg installed for MP4 export.",
                )

    def check_save_current(self):
        """
        Check if the current configuration should be saved.

        Returns:
            True if it's okay to proceed, False if the operation should be cancelled
        """
        # Check if there are unsaved changes
        if self.parameter_panel.has_unsaved_changes():
            reply = QMessageBox.question(
                self,
                "Unsaved Changes",
                "The current configuration has unsaved changes. Do you want to save them?",
                QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel,
                QMessageBox.Save,
            )

            if reply == QMessageBox.Save:
                self.save_configuration()
                return True
            if reply == QMessageBox.Cancel:
                return False

        return True

    def show_about(self):
        """Show the about dialog."""
        QMessageBox.about(
            self,
            "About SHEL",
            """
            <h1>SHEL - SHallow-water numerical modEL</h1>
            <p>Version 1.0.0</p>
            <p>A finite volume, free-surface, variable bottom, shallow-waters equations numerical solver.</p>
            <p><b>Author:</b> Guillaume Riflet</p>
            <p><b>AI Assistant:</b> GitHub Copilot</p>
            <p>Copyright © 2025 SHEL Developers</p>
            <p>This program is free software: you can redistribute it and/or modify
            it under the terms of the GNU General Public License as published by
            the Free Software Foundation, either version 3 of the License, or
            (at your option) any later version.</p>
            """,
        )

    def load_settings(self):
        """Load application settings."""
        settings = QSettings("SHEL", "SHELApp")

        # Load window geometry
        geometry = settings.value("geometry")
        if geometry:
            self.restoreGeometry(geometry)

        # Load window state
        state = settings.value("windowState")
        if state:
            self.restoreState(state)

        # Load current plot type
        plot_type = settings.value("plotType")
        if plot_type and plot_type in self.plot_manager.available_plots:
            self.plot_manager.set_plot_type(plot_type)

    def save_settings(self):
        """Save application settings."""
        settings = QSettings("SHEL", "SHELApp")

        # Save window geometry
        settings.setValue("geometry", self.saveGeometry())

        # Save window state
        settings.setValue("windowState", self.saveState())

        # Save current plot type
        settings.setValue("plotType", self.plot_manager.current_plot_type)

    def closeEvent(self, a0: Optional["QCloseEvent"]) -> None:
        """Handle window close event."""
        if a0 is None:
            return

        # Check if there are unsaved changes
        if not self.check_save_current():
            a0.ignore()
            return

        # Save settings
        self.save_settings()

        # Stop polling timer
        if hasattr(self, "timer"):
            self.timer.stop()

        # Accept the event
        a0.accept()


def main():
    """
    Main entry point for the GUI.
    Returns:
        Exit code
    """
    app = QApplication(sys.argv)
    # Set application properties
    app.setApplicationName("SHEL")
    app.setOrganizationName("SHEL")

    # Create and show the main window
    window = MainWindow()
    window.show()

    return app.exec_()


if __name__ == "__main__":
    sys.exit(main())
