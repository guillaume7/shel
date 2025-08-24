"""
Tests for the GUI components.
"""

import pytest
from PyQt5.QtWidgets import QApplication
from PyQt5.QtTest import QTest
from PyQt5.QtCore import Qt

from shel.gui.main_window import MainWindow
from shel.gui.parameters import ParameterPanel
from shel.gui.visualization import PlotManager


@pytest.fixture
def app():
    """Fixture for Qt application."""
    return QApplication([])


@pytest.fixture
def main_window(app):
    """Fixture for main window."""
    return MainWindow()


def test_main_window_creation(main_window):
    """Test that the main window can be created."""
    assert main_window is not None
    assert main_window.windowTitle() == "SHEL - SHallow-water numerical modEL"


def test_parameter_panel(app):
    """Test that the parameter panel can be created and initialized."""
    panel = ParameterPanel()

    # Check that default values are set
    assert panel.nx_spinbox.value() == 100
    assert panel.ny_spinbox.value() == 100
    assert panel.dt_spinbox.value() == 1.0

    # Check that config is initialized
    config = panel.get_config()
    assert config is not None
    assert "grid" in config
    assert "model" in config
    assert "initial_conditions" in config
    assert "boundary_conditions" in config


def test_plot_manager(app):
    """Test that the plot manager can be created and initialized."""
    plot_manager = PlotManager()

    # Check available plot types
    assert "Water Elevation" in plot_manager.available_plots
    assert "Velocity Field" in plot_manager.available_plots
    assert "Vorticity" in plot_manager.available_plots

    # Check that the current plot type is set
    assert plot_manager.current_plot_type == "Water Elevation"


def test_parameter_changes(app):
    """Test that parameter changes are tracked."""
    panel = ParameterPanel()

    # Check initial state
    assert panel.has_unsaved_changes() is False

    # Change a parameter
    panel.nx_spinbox.setValue(200)

    # Check that changes are tracked
    assert panel.has_unsaved_changes() is True

    # Get updated config
    config = panel.get_config()
    assert config["grid"]["nx"] == 200
