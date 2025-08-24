"""
Parameter panel for the SHEL GUI.

This module provides the parameter panel for the SHEL GUI, allowing users
to set model parameters, initial conditions, and boundary conditions.
"""

import logging
import yaml
from typing import Dict, Any, Optional, List

from PyQt5.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QComboBox,
    QTabWidget,
    QSpinBox,
    QDoubleSpinBox,
    QGroupBox,
    QFormLayout,
    QLineEdit,
    QFileDialog,
    QMessageBox,
    QCheckBox,
)
from PyQt5.QtCore import Qt

logger = logging.getLogger(__name__)


class ParameterPanel(QWidget):
    """
    Panel for editing model parameters.

    This class provides a panel for editing model parameters, including:
    - Grid parameters
    - Time stepping parameters
    - Initial conditions
    - Boundary conditions
    - Physical parameters

    Attributes:
        tabs: Tab widget for organizing parameter groups
        modified: Flag indicating if parameters have been modified
        config: Current configuration dictionary
    """

    def __init__(self, parent=None):
        """Initialize the parameter panel."""
        super().__init__(parent)

        # Set up the layout
        layout = QVBoxLayout(self)

        # Create tab widget for parameters
        self.tabs = QTabWidget()
        layout.addWidget(self.tabs)

        # Create tabs for different parameter groups
        self.create_grid_tab()
        self.create_time_tab()
        self.create_initial_conditions_tab()
        self.create_boundary_conditions_tab()
        self.create_physical_params_tab()
        self.create_output_tab()

        # Add buttons for loading/saving
        button_layout = QHBoxLayout()
        layout.addLayout(button_layout)

        self.reset_button = QPushButton("Reset to Defaults")
        self.reset_button.clicked.connect(self.reset_to_defaults)
        button_layout.addWidget(self.reset_button)

        button_layout.addStretch(1)

        # Initialize state
        self.modified = False
        self.config = self.get_default_config()

        # Apply initial configuration
        self.apply_config_to_ui()

    def create_grid_tab(self):
        """Create the grid parameters tab."""
        tab = QWidget()
        layout = QFormLayout(tab)

        # Grid dimensions
        self.nx_spinbox = QSpinBox()
        self.nx_spinbox.setRange(10, 1000)
        self.nx_spinbox.setValue(100)
        self.nx_spinbox.valueChanged.connect(self.mark_as_modified)
        layout.addRow("Grid Points X:", self.nx_spinbox)

        self.ny_spinbox = QSpinBox()
        self.ny_spinbox.setRange(10, 1000)
        self.ny_spinbox.setValue(100)
        self.ny_spinbox.valueChanged.connect(self.mark_as_modified)
        layout.addRow("Grid Points Y:", self.ny_spinbox)

        # Grid spacing
        self.dx_spinbox = QDoubleSpinBox()
        self.dx_spinbox.setRange(0.1, 10000.0)
        self.dx_spinbox.setValue(100.0)
        self.dx_spinbox.setSuffix(" m")
        self.dx_spinbox.valueChanged.connect(self.mark_as_modified)
        layout.addRow("Grid Spacing X:", self.dx_spinbox)

        self.dy_spinbox = QDoubleSpinBox()
        self.dy_spinbox.setRange(0.1, 10000.0)
        self.dy_spinbox.setValue(100.0)
        self.dy_spinbox.setSuffix(" m")
        self.dy_spinbox.valueChanged.connect(self.mark_as_modified)
        layout.addRow("Grid Spacing Y:", self.dy_spinbox)

        self.tabs.addTab(tab, "Grid")

    def create_time_tab(self):
        """Create the time stepping parameters tab."""
        tab = QWidget()
        layout = QFormLayout(tab)

        # Time step
        self.dt_spinbox = QDoubleSpinBox()
        self.dt_spinbox.setRange(0.01, 100.0)
        self.dt_spinbox.setValue(1.0)
        self.dt_spinbox.setSuffix(" s")
        self.dt_spinbox.setDecimals(3)
        self.dt_spinbox.valueChanged.connect(self.mark_as_modified)
        layout.addRow("Time Step:", self.dt_spinbox)

        # Number of steps
        self.num_steps_spinbox = QSpinBox()
        self.num_steps_spinbox.setRange(1, 1000000)
        self.num_steps_spinbox.setValue(1000)
        self.num_steps_spinbox.valueChanged.connect(self.mark_as_modified)
        layout.addRow("Number of Steps:", self.num_steps_spinbox)

        # Output interval
        self.output_interval_spinbox = QSpinBox()
        self.output_interval_spinbox.setRange(1, 10000)
        self.output_interval_spinbox.setValue(10)
        self.output_interval_spinbox.valueChanged.connect(self.mark_as_modified)
        layout.addRow("Output Interval:", self.output_interval_spinbox)

        # Solver type
        self.solver_combo = QComboBox()
        self.solver_combo.addItems(["leapfrog"])
        self.solver_combo.currentTextChanged.connect(self.mark_as_modified)
        layout.addRow("Solver:", self.solver_combo)

        self.tabs.addTab(tab, "Time")

    def create_initial_conditions_tab(self):
        """Create the initial conditions tab."""
        tab = QWidget()
        layout = QVBoxLayout(tab)

        # Water elevation group
        elevation_group = QGroupBox("Water Elevation")
        elevation_layout = QFormLayout(elevation_group)
        layout.addWidget(elevation_group)

        # Elevation type
        self.elevation_type_combo = QComboBox()
        self.elevation_type_combo.addItems(["flat", "gaussian", "custom"])
        self.elevation_type_combo.currentTextChanged.connect(
            self.on_elevation_type_changed
        )
        elevation_layout.addRow("Type:", self.elevation_type_combo)

        # Gaussian parameters
        self.gaussian_amplitude_spinbox = QDoubleSpinBox()
        self.gaussian_amplitude_spinbox.setRange(-100.0, 100.0)
        self.gaussian_amplitude_spinbox.setValue(1.0)
        self.gaussian_amplitude_spinbox.setSuffix(" m")
        self.gaussian_amplitude_spinbox.valueChanged.connect(self.mark_as_modified)
        elevation_layout.addRow("Amplitude:", self.gaussian_amplitude_spinbox)

        self.gaussian_x0_spinbox = QDoubleSpinBox()
        self.gaussian_x0_spinbox.setRange(0.0, 10000.0)
        self.gaussian_x0_spinbox.setValue(5000.0)
        self.gaussian_x0_spinbox.setSuffix(" m")
        self.gaussian_x0_spinbox.valueChanged.connect(self.mark_as_modified)
        elevation_layout.addRow("X Center:", self.gaussian_x0_spinbox)

        self.gaussian_y0_spinbox = QDoubleSpinBox()
        self.gaussian_y0_spinbox.setRange(0.0, 10000.0)
        self.gaussian_y0_spinbox.setValue(5000.0)
        self.gaussian_y0_spinbox.setSuffix(" m")
        self.gaussian_y0_spinbox.valueChanged.connect(self.mark_as_modified)
        elevation_layout.addRow("Y Center:", self.gaussian_y0_spinbox)

        self.gaussian_sigma_spinbox = QDoubleSpinBox()
        self.gaussian_sigma_spinbox.setRange(1.0, 5000.0)
        self.gaussian_sigma_spinbox.setValue(500.0)
        self.gaussian_sigma_spinbox.setSuffix(" m")
        self.gaussian_sigma_spinbox.valueChanged.connect(self.mark_as_modified)
        elevation_layout.addRow("Sigma:", self.gaussian_sigma_spinbox)

        # Custom file input (for custom type)
        self.custom_elevation_file_layout = QHBoxLayout()
        self.custom_elevation_file_edit = QLineEdit()
        self.custom_elevation_file_edit.setReadOnly(True)
        self.custom_elevation_file_button = QPushButton("Browse...")
        self.custom_elevation_file_button.clicked.connect(self.browse_elevation_file)

        self.custom_elevation_file_layout.addWidget(self.custom_elevation_file_edit)
        self.custom_elevation_file_layout.addWidget(self.custom_elevation_file_button)
        elevation_layout.addRow("File:", self.custom_elevation_file_layout)

        # Bathymetry group
        bathymetry_group = QGroupBox("Bathymetry")
        bathymetry_layout = QFormLayout(bathymetry_group)
        layout.addWidget(bathymetry_group)

        # Bathymetry type
        self.bathymetry_type_combo = QComboBox()
        self.bathymetry_type_combo.addItems(["flat", "linear_slope", "bump", "custom"])
        self.bathymetry_type_combo.currentTextChanged.connect(
            self.on_bathymetry_type_changed
        )
        bathymetry_layout.addRow("Type:", self.bathymetry_type_combo)

        # Flat depth
        self.flat_depth_spinbox = QDoubleSpinBox()
        self.flat_depth_spinbox.setRange(0.1, 10000.0)
        self.flat_depth_spinbox.setValue(100.0)
        self.flat_depth_spinbox.setSuffix(" m")
        self.flat_depth_spinbox.valueChanged.connect(self.mark_as_modified)
        bathymetry_layout.addRow("Depth:", self.flat_depth_spinbox)

        # Slope parameters
        self.slope_min_depth_spinbox = QDoubleSpinBox()
        self.slope_min_depth_spinbox.setRange(0.1, 10000.0)
        self.slope_min_depth_spinbox.setValue(10.0)
        self.slope_min_depth_spinbox.setSuffix(" m")
        self.slope_min_depth_spinbox.valueChanged.connect(self.mark_as_modified)
        bathymetry_layout.addRow("Min Depth:", self.slope_min_depth_spinbox)

        self.slope_max_depth_spinbox = QDoubleSpinBox()
        self.slope_max_depth_spinbox.setRange(0.1, 10000.0)
        self.slope_max_depth_spinbox.setValue(100.0)
        self.slope_max_depth_spinbox.setSuffix(" m")
        self.slope_max_depth_spinbox.valueChanged.connect(self.mark_as_modified)
        bathymetry_layout.addRow("Max Depth:", self.slope_max_depth_spinbox)

        # Bump parameters
        self.bump_amplitude_spinbox = QDoubleSpinBox()
        self.bump_amplitude_spinbox.setRange(-1000.0, 1000.0)
        self.bump_amplitude_spinbox.setValue(-50.0)
        self.bump_amplitude_spinbox.setSuffix(" m")
        self.bump_amplitude_spinbox.valueChanged.connect(self.mark_as_modified)
        bathymetry_layout.addRow("Bump Height:", self.bump_amplitude_spinbox)

        self.bump_x0_spinbox = QDoubleSpinBox()
        self.bump_x0_spinbox.setRange(0.0, 10000.0)
        self.bump_x0_spinbox.setValue(5000.0)
        self.bump_x0_spinbox.setSuffix(" m")
        self.bump_x0_spinbox.valueChanged.connect(self.mark_as_modified)
        bathymetry_layout.addRow("Bump X Center:", self.bump_x0_spinbox)

        self.bump_y0_spinbox = QDoubleSpinBox()
        self.bump_y0_spinbox.setRange(0.0, 10000.0)
        self.bump_y0_spinbox.setValue(5000.0)
        self.bump_y0_spinbox.setSuffix(" m")
        self.bump_y0_spinbox.valueChanged.connect(self.mark_as_modified)
        bathymetry_layout.addRow("Bump Y Center:", self.bump_y0_spinbox)

        self.bump_sigma_spinbox = QDoubleSpinBox()
        self.bump_sigma_spinbox.setRange(1.0, 5000.0)
        self.bump_sigma_spinbox.setValue(500.0)
        self.bump_sigma_spinbox.setSuffix(" m")
        self.bump_sigma_spinbox.valueChanged.connect(self.mark_as_modified)
        bathymetry_layout.addRow("Bump Sigma:", self.bump_sigma_spinbox)

        # Custom file input (for custom type)
        self.custom_bathymetry_file_layout = QHBoxLayout()
        self.custom_bathymetry_file_edit = QLineEdit()
        self.custom_bathymetry_file_edit.setReadOnly(True)
        self.custom_bathymetry_file_button = QPushButton("Browse...")
        self.custom_bathymetry_file_button.clicked.connect(self.browse_bathymetry_file)

        self.custom_bathymetry_file_layout.addWidget(self.custom_bathymetry_file_edit)
        self.custom_bathymetry_file_layout.addWidget(self.custom_bathymetry_file_button)
        bathymetry_layout.addRow("File:", self.custom_bathymetry_file_layout)

        # Initialize visibility based on current selections
        self.on_elevation_type_changed(self.elevation_type_combo.currentText())
        self.on_bathymetry_type_changed(self.bathymetry_type_combo.currentText())

        self.tabs.addTab(tab, "Initial Conditions")

    def create_boundary_conditions_tab(self):
        """Create the boundary conditions tab."""
        tab = QWidget()
        layout = QFormLayout(tab)

        # Boundary condition types
        bc_types = ["closed", "free_slip", "radiative"]

        # North boundary
        self.north_bc_combo = QComboBox()
        self.north_bc_combo.addItems(bc_types)
        self.north_bc_combo.currentTextChanged.connect(self.mark_as_modified)
        layout.addRow("North:", self.north_bc_combo)

        # South boundary
        self.south_bc_combo = QComboBox()
        self.south_bc_combo.addItems(bc_types)
        self.south_bc_combo.currentTextChanged.connect(self.mark_as_modified)
        layout.addRow("South:", self.south_bc_combo)

        # East boundary
        self.east_bc_combo = QComboBox()
        self.east_bc_combo.addItems(bc_types)
        self.east_bc_combo.currentTextChanged.connect(self.mark_as_modified)
        layout.addRow("East:", self.east_bc_combo)

        # West boundary
        self.west_bc_combo = QComboBox()
        self.west_bc_combo.addItems(bc_types)
        self.west_bc_combo.currentTextChanged.connect(self.mark_as_modified)
        layout.addRow("West:", self.west_bc_combo)

        self.tabs.addTab(tab, "Boundary Conditions")

    def create_physical_params_tab(self):
        """Create the physical parameters tab."""
        tab = QWidget()
        layout = QFormLayout(tab)

        # Viscosity
        self.viscosity_spinbox = QDoubleSpinBox()
        self.viscosity_spinbox.setRange(0.0, 1000.0)
        self.viscosity_spinbox.setValue(1.0)
        self.viscosity_spinbox.setSuffix(" m²/s")
        self.viscosity_spinbox.setDecimals(3)
        self.viscosity_spinbox.valueChanged.connect(self.mark_as_modified)
        layout.addRow("Viscosity:", self.viscosity_spinbox)

        # Bottom friction
        self.bottom_friction_spinbox = QDoubleSpinBox()
        self.bottom_friction_spinbox.setRange(0.0, 1.0)
        self.bottom_friction_spinbox.setValue(0.002)
        self.bottom_friction_spinbox.setDecimals(5)
        self.bottom_friction_spinbox.setSingleStep(0.001)
        self.bottom_friction_spinbox.valueChanged.connect(self.mark_as_modified)
        layout.addRow("Bottom Friction:", self.bottom_friction_spinbox)

        # Coriolis parameter
        self.coriolis_spinbox = QDoubleSpinBox()
        self.coriolis_spinbox.setRange(-1e-3, 1e-3)
        self.coriolis_spinbox.setValue(0.0)
        self.coriolis_spinbox.setSuffix(" s⁻¹")
        self.coriolis_spinbox.setDecimals(8)
        self.coriolis_spinbox.setSingleStep(1e-5)
        self.coriolis_spinbox.valueChanged.connect(self.mark_as_modified)
        layout.addRow("Coriolis Parameter:", self.coriolis_spinbox)

        # Gravity
        self.gravity_spinbox = QDoubleSpinBox()
        self.gravity_spinbox.setRange(1.0, 20.0)
        self.gravity_spinbox.setValue(9.81)
        self.gravity_spinbox.setSuffix(" m/s²")
        self.gravity_spinbox.valueChanged.connect(self.mark_as_modified)
        layout.addRow("Gravity:", self.gravity_spinbox)

        self.tabs.addTab(tab, "Physical Parameters")

    def create_output_tab(self):
        """Create the output parameters tab."""
        tab = QWidget()
        layout = QFormLayout(tab)

        # Output directory
        self.output_dir_layout = QHBoxLayout()
        self.output_dir_edit = QLineEdit("./output")
        self.output_dir_button = QPushButton("Browse...")
        self.output_dir_button.clicked.connect(self.browse_output_dir)

        self.output_dir_layout.addWidget(self.output_dir_edit)
        self.output_dir_layout.addWidget(self.output_dir_button)
        layout.addRow("Output Directory:", self.output_dir_layout)

        # Enable output
        self.enable_output_checkbox = QCheckBox("Enable Output Files")
        self.enable_output_checkbox.setChecked(True)
        self.enable_output_checkbox.stateChanged.connect(self.mark_as_modified)
        layout.addRow("", self.enable_output_checkbox)

        # ZeroMQ publishing port
        self.zmq_port_spinbox = QSpinBox()
        self.zmq_port_spinbox.setRange(1024, 65535)
        self.zmq_port_spinbox.setValue(5556)
        self.zmq_port_spinbox.valueChanged.connect(self.mark_as_modified)
        layout.addRow("ZeroMQ Port:", self.zmq_port_spinbox)

        # Enable ZeroMQ
        self.enable_zmq_checkbox = QCheckBox("Enable ZeroMQ Communication")
        self.enable_zmq_checkbox.setChecked(True)
        self.enable_zmq_checkbox.stateChanged.connect(self.mark_as_modified)
        layout.addRow("", self.enable_zmq_checkbox)

        self.tabs.addTab(tab, "Output")

    def get_default_config(self) -> Dict[str, Any]:
        """
        Get the default configuration.

        Returns:
            Default configuration dictionary
        """
        return {
            "grid": {"nx": 100, "ny": 100, "dx": 100.0, "dy": 100.0},
            "model": {
                "dt": 1.0,
                "num_steps": 1000,
                "output_interval": 10,
                "solver": "leapfrog",
            },
            "initial_conditions": {
                "type": "gaussian",
                "gaussian": {
                    "amplitude": 1.0,
                    "x0": 5000.0,
                    "y0": 5000.0,
                    "sigma": 500.0,
                },
                "custom": {"file": ""},
            },
            "bathymetry": {
                "type": "flat",
                "flat": {"depth": 100.0},
                "linear_slope": {"min_depth": 10.0, "max_depth": 100.0},
                "bump": {
                    "amplitude": -50.0,
                    "x0": 5000.0,
                    "y0": 5000.0,
                    "sigma": 500.0,
                    "base_depth": 100.0,
                },
                "custom": {"file": ""},
            },
            "boundary_conditions": {
                "north": "radiative",
                "south": "radiative",
                "east": "radiative",
                "west": "radiative",
            },
            "physical_parameters": {
                "viscosity": 1.0,
                "bottom_friction": 0.002,
                "coriolis": 0.0,
                "gravity": 9.81,
            },
            "output": {"directory": "./output", "enabled": True},
            "communication": {"zmq_pub_port": 5556, "enable_zmq": True},
        }

    def reset_to_defaults(self):
        """Reset all parameters to defaults."""
        # Confirm with user
        reply = QMessageBox.question(
            self,
            "Reset Parameters",
            "Are you sure you want to reset all parameters to their default values?",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No,
        )

        if reply == QMessageBox.Yes:
            # Reset to default configuration
            self.config = self.get_default_config()

            # Apply to UI
            self.apply_config_to_ui()

            # Mark as modified
            self.modified = True

    def apply_config_to_ui(self):
        """Apply the current configuration to the UI controls."""
        # Grid tab
        self.nx_spinbox.setValue(self.config["grid"]["nx"])
        self.ny_spinbox.setValue(self.config["grid"]["ny"])
        self.dx_spinbox.setValue(self.config["grid"]["dx"])
        self.dy_spinbox.setValue(self.config["grid"]["dy"])

        # Time tab
        self.dt_spinbox.setValue(self.config["model"]["dt"])
        self.num_steps_spinbox.setValue(self.config["model"]["num_steps"])
        self.output_interval_spinbox.setValue(self.config["model"]["output_interval"])
        self.solver_combo.setCurrentText(self.config["model"]["solver"])

        # Initial conditions tab
        ic_type = self.config["initial_conditions"]["type"]
        self.elevation_type_combo.setCurrentText(ic_type)

        if ic_type == "gaussian":
            self.gaussian_amplitude_spinbox.setValue(
                self.config["initial_conditions"]["gaussian"]["amplitude"]
            )
            self.gaussian_x0_spinbox.setValue(
                self.config["initial_conditions"]["gaussian"]["x0"]
            )
            self.gaussian_y0_spinbox.setValue(
                self.config["initial_conditions"]["gaussian"]["y0"]
            )
            self.gaussian_sigma_spinbox.setValue(
                self.config["initial_conditions"]["gaussian"]["sigma"]
            )
        elif ic_type == "custom":
            self.custom_elevation_file_edit.setText(
                self.config["initial_conditions"]["custom"]["file"]
            )

        bathy_type = self.config["bathymetry"]["type"]
        self.bathymetry_type_combo.setCurrentText(bathy_type)

        if bathy_type == "flat":
            self.flat_depth_spinbox.setValue(self.config["bathymetry"]["flat"]["depth"])
        elif bathy_type == "linear_slope":
            self.slope_min_depth_spinbox.setValue(
                self.config["bathymetry"]["linear_slope"]["min_depth"]
            )
            self.slope_max_depth_spinbox.setValue(
                self.config["bathymetry"]["linear_slope"]["max_depth"]
            )
        elif bathy_type == "bump":
            self.bump_amplitude_spinbox.setValue(
                self.config["bathymetry"]["bump"]["amplitude"]
            )
            self.bump_x0_spinbox.setValue(self.config["bathymetry"]["bump"]["x0"])
            self.bump_y0_spinbox.setValue(self.config["bathymetry"]["bump"]["y0"])
            self.bump_sigma_spinbox.setValue(self.config["bathymetry"]["bump"]["sigma"])
        elif bathy_type == "custom":
            self.custom_bathymetry_file_edit.setText(
                self.config["bathymetry"]["custom"]["file"]
            )

        # Boundary conditions tab
        self.north_bc_combo.setCurrentText(self.config["boundary_conditions"]["north"])
        self.south_bc_combo.setCurrentText(self.config["boundary_conditions"]["south"])
        self.east_bc_combo.setCurrentText(self.config["boundary_conditions"]["east"])
        self.west_bc_combo.setCurrentText(self.config["boundary_conditions"]["west"])

        # Physical parameters tab
        self.viscosity_spinbox.setValue(self.config["physical_parameters"]["viscosity"])
        self.bottom_friction_spinbox.setValue(
            self.config["physical_parameters"]["bottom_friction"]
        )
        self.coriolis_spinbox.setValue(self.config["physical_parameters"]["coriolis"])
        self.gravity_spinbox.setValue(self.config["physical_parameters"]["gravity"])

        # Output tab
        self.output_dir_edit.setText(self.config["output"]["directory"])
        self.enable_output_checkbox.setChecked(self.config["output"]["enabled"])
        self.zmq_port_spinbox.setValue(self.config["communication"]["zmq_pub_port"])
        self.enable_zmq_checkbox.setChecked(self.config["communication"]["enable_zmq"])

    def update_config_from_ui(self):
        """Update the configuration from UI controls."""
        # Grid tab
        self.config["grid"]["nx"] = self.nx_spinbox.value()
        self.config["grid"]["ny"] = self.ny_spinbox.value()
        self.config["grid"]["dx"] = self.dx_spinbox.value()
        self.config["grid"]["dy"] = self.dy_spinbox.value()

        # Time tab
        self.config["model"]["dt"] = self.dt_spinbox.value()
        self.config["model"]["num_steps"] = self.num_steps_spinbox.value()
        self.config["model"]["output_interval"] = self.output_interval_spinbox.value()
        self.config["model"]["solver"] = self.solver_combo.currentText()

        # Initial conditions tab
        ic_type = self.elevation_type_combo.currentText()
        self.config["initial_conditions"]["type"] = ic_type

        if ic_type == "gaussian":
            self.config["initial_conditions"]["gaussian"] = {
                "amplitude": self.gaussian_amplitude_spinbox.value(),
                "x0": self.gaussian_x0_spinbox.value(),
                "y0": self.gaussian_y0_spinbox.value(),
                "sigma": self.gaussian_sigma_spinbox.value(),
            }
        elif ic_type == "custom":
            self.config["initial_conditions"]["custom"] = {
                "file": self.custom_elevation_file_edit.text()
            }

        bathy_type = self.bathymetry_type_combo.currentText()
        self.config["bathymetry"]["type"] = bathy_type

        if bathy_type == "flat":
            self.config["bathymetry"]["flat"] = {
                "depth": self.flat_depth_spinbox.value()
            }
        elif bathy_type == "linear_slope":
            self.config["bathymetry"]["linear_slope"] = {
                "min_depth": self.slope_min_depth_spinbox.value(),
                "max_depth": self.slope_max_depth_spinbox.value(),
            }
        elif bathy_type == "bump":
            self.config["bathymetry"]["bump"] = {
                "amplitude": self.bump_amplitude_spinbox.value(),
                "x0": self.bump_x0_spinbox.value(),
                "y0": self.bump_y0_spinbox.value(),
                "sigma": self.bump_sigma_spinbox.value(),
                "base_depth": 100.0,  # Default base depth
            }
        elif bathy_type == "custom":
            self.config["bathymetry"]["custom"] = {
                "file": self.custom_bathymetry_file_edit.text()
            }

        # Boundary conditions tab
        self.config["boundary_conditions"] = {
            "north": self.north_bc_combo.currentText(),
            "south": self.south_bc_combo.currentText(),
            "east": self.east_bc_combo.currentText(),
            "west": self.west_bc_combo.currentText(),
        }

        # Physical parameters tab
        self.config["physical_parameters"] = {
            "viscosity": self.viscosity_spinbox.value(),
            "bottom_friction": self.bottom_friction_spinbox.value(),
            "coriolis": self.coriolis_spinbox.value(),
            "gravity": self.gravity_spinbox.value(),
        }

        # Output tab
        self.config["output"] = {
            "directory": self.output_dir_edit.text(),
            "enabled": self.enable_output_checkbox.isChecked(),
        }

        self.config["communication"] = {
            "zmq_pub_port": self.zmq_port_spinbox.value(),
            "enable_zmq": self.enable_zmq_checkbox.isChecked(),
        }

    def mark_as_modified(self):
        """Mark the configuration as modified."""
        self.modified = True

    def has_unsaved_changes(self) -> bool:
        """
        Check if there are unsaved changes.

        Returns:
            True if there are unsaved changes, False otherwise
        """
        return self.modified

    def on_elevation_type_changed(self, elevation_type):
        """
        Handle changes to the elevation type.

        Args:
            elevation_type: New elevation type
        """
        # Show/hide parameters based on the selected type
        show_gaussian = elevation_type == "gaussian"
        show_custom = elevation_type == "custom"

        self.gaussian_amplitude_spinbox.setVisible(show_gaussian)
        self.gaussian_x0_spinbox.setVisible(show_gaussian)
        self.gaussian_y0_spinbox.setVisible(show_gaussian)
        self.gaussian_sigma_spinbox.setVisible(show_gaussian)

        self.custom_elevation_file_edit.setVisible(show_custom)
        self.custom_elevation_file_button.setVisible(show_custom)

        # Mark as modified
        self.mark_as_modified()

    def on_bathymetry_type_changed(self, bathymetry_type):
        """
        Handle changes to the bathymetry type.

        Args:
            bathymetry_type: New bathymetry type
        """
        # Show/hide parameters based on the selected type
        show_flat = bathymetry_type == "flat"
        show_slope = bathymetry_type == "linear_slope"
        show_bump = bathymetry_type == "bump"
        show_custom = bathymetry_type == "custom"

        self.flat_depth_spinbox.setVisible(show_flat)

        self.slope_min_depth_spinbox.setVisible(show_slope)
        self.slope_max_depth_spinbox.setVisible(show_slope)

        self.bump_amplitude_spinbox.setVisible(show_bump)
        self.bump_x0_spinbox.setVisible(show_bump)
        self.bump_y0_spinbox.setVisible(show_bump)
        self.bump_sigma_spinbox.setVisible(show_bump)

        self.custom_bathymetry_file_edit.setVisible(show_custom)
        self.custom_bathymetry_file_button.setVisible(show_custom)

        # Mark as modified
        self.mark_as_modified()

    def browse_elevation_file(self):
        """Browse for an elevation file."""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Open Elevation File", "", "NetCDF Files (*.nc);;All Files (*)"
        )

        if file_path:
            self.custom_elevation_file_edit.setText(file_path)
            self.mark_as_modified()

    def browse_bathymetry_file(self):
        """Browse for a bathymetry file."""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Open Bathymetry File", "", "NetCDF Files (*.nc);;All Files (*)"
        )

        if file_path:
            self.custom_bathymetry_file_edit.setText(file_path)
            self.mark_as_modified()

    def browse_output_dir(self):
        """Browse for an output directory."""
        dir_path = QFileDialog.getExistingDirectory(
            self, "Select Output Directory", self.output_dir_edit.text()
        )

        if dir_path:
            self.output_dir_edit.setText(dir_path)
            self.mark_as_modified()

    def load_from_file(self, file_path: str):
        """
        Load configuration from a file.

        Args:
            file_path: Path to the configuration file
        """
        # Load YAML file
        with open(file_path, "r") as f:
            self.config = yaml.safe_load(f)

        # Apply to UI
        self.apply_config_to_ui()

        # Reset modified flag
        self.modified = False

    def save_to_file(self, file_path: str):
        """
        Save configuration to a file.

        Args:
            file_path: Path to save the configuration file
        """
        # Update config from UI
        self.update_config_from_ui()

        # Save to YAML file
        with open(file_path, "w") as f:
            yaml.dump(self.config, f, default_flow_style=False)

        # Reset modified flag
        self.modified = False

    def get_config(self) -> Dict[str, Any]:
        """
        Get the current configuration.

        Returns:
            Current configuration dictionary
        """
        # Update config from UI
        self.update_config_from_ui()

        return self.config
