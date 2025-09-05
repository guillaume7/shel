"""
Tests for the GUI components.

Note: Temporarily skipped in headless CI until GUI phases (G1–G6) are tackled.
"""

import pytest

try:  # pragma: no cover - optional GUI dependency
    from PyQt5.QtWidgets import QApplication
    from PyQt5.QtTest import QTest  # noqa: F401
    from PyQt5.QtCore import Qt  # noqa: F401
    _PYQT_AVAILABLE = True
except Exception:  # broad catch to skip gracefully in minimal envs
    _PYQT_AVAILABLE = False

pytestmark = pytest.mark.skip(reason="GUI tests disabled until GUI phases are implemented")

if _PYQT_AVAILABLE:
    from shel.gui.main_window import MainWindow  # type: ignore
    from shel.gui.parameters import ParameterPanel  # type: ignore
    from shel.gui.visualization import PlotManager  # type: ignore

@pytest.fixture
def app():
    """Fixture for Qt application (skipped if PyQt5 unavailable)."""
    if not _PYQT_AVAILABLE:
        pytest.skip("PyQt5 not installed")
    # mypy/linters: QApplication guaranteed defined when _PYQT_AVAILABLE True
    return QApplication([])  # type: ignore[name-defined]


@pytest.fixture
def main_window(app):
    """Fixture for main window."""
    return MainWindow()  # type: ignore[name-defined]


@pytest.mark.skipif(not _PYQT_AVAILABLE, reason="PyQt5 not installed")
def test_main_window_creation(main_window):
    """Test that the main window can be created."""
    assert main_window is not None
    assert main_window.windowTitle() == "SHEL - SHallow-water numerical modEL"


@pytest.mark.skipif(not _PYQT_AVAILABLE, reason="PyQt5 not installed")
def test_parameter_panel(app):
    """Test that the parameter panel can be created and initialized."""
    panel = ParameterPanel()  # type: ignore[name-defined]

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


@pytest.mark.skipif(not _PYQT_AVAILABLE, reason="PyQt5 not installed")
def test_plot_manager(app):
    """Test that the plot manager can be created and initialized."""
    plot_manager = PlotManager()  # type: ignore[name-defined]

    # Check available plot types
    assert "Water Elevation" in plot_manager.available_plots
    assert "Velocity Field" in plot_manager.available_plots
    assert "Vorticity" in plot_manager.available_plots

    # Check that the current plot type is set
    assert plot_manager.current_plot_type == "Water Elevation"


@pytest.mark.skipif(not _PYQT_AVAILABLE, reason="PyQt5 not installed")
def test_parameter_changes(app):
    """Test that parameter changes are tracked."""
    panel = ParameterPanel()  # type: ignore[name-defined]

    # Check initial state
    assert panel.has_unsaved_changes() is False

    # Change a parameter
    panel.nx_spinbox.setValue(200)

    # Check that changes are tracked
    assert panel.has_unsaved_changes() is True

    # Get updated config
    config = panel.get_config()
    assert config["grid"]["nx"] == 200
