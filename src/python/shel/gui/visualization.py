"""
Visualization components for the SHEL GUI.

This module provides the plotting and visualization components for the SHEL GUI,
including 2D plots of model variables and time series plots.
"""

import logging

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_qt import NavigationToolbar2QT as NavigationToolbar
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PyQt5.QtWidgets import QComboBox, QHBoxLayout, QLabel, QVBoxLayout, QWidget

logger = logging.getLogger(__name__)


class PlotManager(QWidget):
    """
    Manager for plots in the GUI.

    This class manages the plotting area in the GUI, including:
    - 2D plots of model variables
    - Time series plots of energy, volume, etc.

    Attributes:
        figure: Matplotlib figure
        canvas: Matplotlib canvas
        toolbar: Matplotlib toolbar
        available_plots: List of available plot types
        current_plot_type: Current plot type
    """

    def __init__(self, parent=None):
        """
        Initialize the plot manager.

        Args:
            parent: Parent widget
        """
        super().__init__(parent)

        # Set up the layout
        layout = QVBoxLayout(self)

        # Create toolbar area at the top
        toolbar_layout = QHBoxLayout()
        layout.addLayout(toolbar_layout)

        # Add plot type selector
        self.plot_type_label = QLabel("Plot Type:")
        toolbar_layout.addWidget(self.plot_type_label)

        self.plot_type_selector = QComboBox()
        self.available_plots = [
            "Water Elevation",
            "Velocity Field",
            "Vorticity",
            "Energy Time Series",
            "Volume Time Series",
        ]
        self.plot_type_selector.addItems(self.available_plots)
        self.plot_type_selector.currentTextChanged.connect(self.set_plot_type)
        toolbar_layout.addWidget(self.plot_type_selector)

        toolbar_layout.addStretch(1)  # Add stretch to push the selector to the left

        # Create matplotlib figure and canvas
        self.figure = Figure(figsize=(8, 6), dpi=100)
        self.canvas = FigureCanvas(self.figure)
        layout.addWidget(self.canvas)

        # Add matplotlib toolbar
        self.toolbar = NavigationToolbar(self.canvas, self)
        layout.addWidget(self.toolbar)

        # Initialize attributes
        self.current_plot_type = "Water Elevation"
        self.last_message = None

        # Predefine plot object handles and data containers to satisfy linters
        self.ax = None  # type: ignore[assignment]
        self.eta_image = None  # set by setup_plot/update
        self.vort_image = None  # set by setup_plot/update
        self.quiver_obj = None  # set by update on velocity field
        self.kinetic_line = None  # type: ignore[assignment]
        self.potential_line = None  # type: ignore[assignment]
        self.total_line = None  # type: ignore[assignment]
        self.volume_line = None  # type: ignore[assignment]
        self.time_data = []
        self.kinetic_data = []
        self.potential_data = []
        self.total_data = []
        self.volume_time_data = []
        self.volume_data = []
        self._animation_frames = []

        # Set up the initial empty plot
        self.setup_plot()

    def setup_plot(self):
        """Set up the plot based on the current plot type."""
        # Clear the figure
        self.figure.clear()

        if self.current_plot_type in ["Water Elevation", "Velocity Field", "Vorticity"]:
            # 2D plots
            self.ax = self.figure.add_subplot(111)
            self.ax.set_aspect("equal")

            if self.current_plot_type == "Water Elevation":
                # Prepare image handle for elevation
                self.eta_image = self.ax.imshow(
                    np.zeros((10, 10)),
                    cmap="coolwarm",
                    interpolation="bilinear",
                    origin="lower",
                )
                self.figure.colorbar(self.eta_image, ax=self.ax, label="Elevation (m)")
                self.ax.set_title("Water Elevation")

            elif self.current_plot_type == "Velocity Field":
                self.ax.set_title("Velocity Field")
                # Quiver will be created on first update
                self.quiver_obj = None  # type: ignore[attr-defined]

            elif self.current_plot_type == "Vorticity":
                # Prepare image handle for vorticity
                self.vort_image = self.ax.imshow(
                    np.zeros((10, 10)),
                    cmap="RdBu_r",
                    interpolation="bilinear",
                    origin="lower",
                )
                self.figure.colorbar(
                    self.vort_image, ax=self.ax, label="Vorticity (s⁻¹)"
                )
                self.ax.set_title("Vorticity")

            self.ax.set_xlabel("X (m)")
            self.ax.set_ylabel("Y (m)")

        elif self.current_plot_type in ["Energy Time Series", "Volume Time Series"]:
            # Time series plots
            self.ax = self.figure.add_subplot(111)

            if self.current_plot_type == "Energy Time Series":
                self.ax.set_title("Energy vs. Time")
                self.ax.set_xlabel("Time (s)")
                self.ax.set_ylabel("Energy (J)")

                # Initialize empty lines for different energy components
                (self.kinetic_line,) = self.ax.plot([], [], "b-", label="Kinetic")
                (self.potential_line,) = self.ax.plot([], [], "g-", label="Potential")
                (self.total_line,) = self.ax.plot([], [], "r-", label="Total")

                self.ax.legend()

            elif self.current_plot_type == "Volume Time Series":
                self.ax.set_title("Water Volume vs. Time")
                self.ax.set_xlabel("Time (s)")
                self.ax.set_ylabel("Volume (m³)")

                # Initialize empty line for volume
                (self.volume_line,) = self.ax.plot([], [], "b-")

        # Update the canvas
        self.canvas.draw()

    def set_plot_type(self, plot_type):
        """
        Set the current plot type.

        Args:
            plot_type: Plot type to set
        """
        self.current_plot_type = plot_type
        self.setup_plot()

        # Update the plot if we have data
        if self.last_message:
            self.update_plots(self.last_message)

    def update_plots(self, message):
        """
        Update the plots with new data.

        Args:
            message: Message containing model state
        """
        # Store the last message
        self.last_message = message

        # Check if we have field data
        fields = message.get("fields", {})
        if not fields:
            return

        # Update the current plot
        ax = self.ax
        if ax is None:
            return
        if self.current_plot_type == "Water Elevation":
            # Update water elevation plot
            eta = np.array(fields.get("eta", []))
            if eta.size > 0:
                # Initialize image if needed (in case plot type switched after init)
                if not hasattr(self, "eta_image") or self.eta_image is None:
                    self.eta_image = ax.imshow(
                        eta,
                        cmap="coolwarm",
                        interpolation="bilinear",
                        origin="lower",
                    )
                else:
                    self.eta_image.set_data(eta)
                    self.eta_image.set_clim(float(np.min(eta)), float(np.max(eta)))

                # Update axes limits and ticks
                ny, nx = eta.shape
                stride_y = fields.get("stride_y", 1)
                stride_x = fields.get("stride_x", 1)

                # Set extent based on grid size
                extent = (0.0, float(nx * stride_x), 0.0, float(ny * stride_y))
                self.eta_image.set_extent(extent)

                # Update title with time
                ax.set_title(f"Water Elevation (t = {message.get('time', 0):.2f} s)")

        elif self.current_plot_type == "Velocity Field":
            # Update velocity field plot
            u = np.array(fields.get("u", []))
            v = np.array(fields.get("v", []))

            if u.size > 0 and v.size > 0:
                # Clear previous quiver plot
                ax.clear()

                # Get grid dimensions
                ny, nx = u.shape
                stride_y = fields.get("stride_y", 1)
                stride_x = fields.get("stride_x", 1)

                # Create grid for quiver plot
                X, Y = np.meshgrid(
                    np.arange(0, nx * stride_x, stride_x),
                    np.arange(0, ny * stride_y, stride_y),
                )

                # Create quiver plot
                self.quiver_obj = ax.quiver(X, Y, u, v, scale=0.2)  # type: ignore[attr-defined]

                # Add colorbar for velocity magnitude
                vel_mag = np.sqrt(u**2 + v**2)
                contour = ax.contourf(X, Y, vel_mag, cmap="viridis", alpha=0.3)
                self.figure.colorbar(contour, ax=ax, label="Velocity (m/s)")

                # Set up axes
                ax.set_aspect("equal")
                ax.set_xlabel("X (m)")
                ax.set_ylabel("Y (m)")
                ax.set_title(f"Velocity Field (t = {message.get('time', 0):.2f} s)")

        elif self.current_plot_type == "Vorticity":
            # Update vorticity plot
            # We need to compute vorticity from u and v
            u = np.array(fields.get("u", []))
            v = np.array(fields.get("v", []))

            if u.size > 0 and v.size > 0:
                # Get grid dimensions
                ny, nx = u.shape
                stride_y = fields.get("stride_y", 1)
                stride_x = fields.get("stride_x", 1)

                # Compute vorticity (simple central difference)
                vort = np.zeros((ny, nx))
                for j in range(1, ny - 1):
                    for i in range(1, nx - 1):
                        dvdx = (v[j, i + 1] - v[j, i - 1]) / (2 * stride_x)
                        dudy = (u[j + 1, i] - u[j - 1, i]) / (2 * stride_y)
                        vort[j, i] = dvdx - dudy

                # Initialize image if needed
                if not hasattr(self, "vort_image") or self.vort_image is None:
                    self.vort_image = ax.imshow(
                        vort,
                        cmap="RdBu_r",
                        interpolation="bilinear",
                        origin="lower",
                    )
                else:
                    # Update plot
                    self.vort_image.set_data(vort)
                vmax = np.max(np.abs(vort))
                self.vort_image.set_clim(float(-vmax), float(vmax))

                # Set extent based on grid size
                extent = (0.0, float(nx * stride_x), 0.0, float(ny * stride_y))
                self.vort_image.set_extent(extent)

                # Update title with time
                ax.set_title(f"Vorticity (t = {message.get('time', 0):.2f} s)")

        elif self.current_plot_type == "Energy Time Series":
            # Update energy time series plot
            # Get the stored time series data
            time_data = getattr(self, "time_data", [])
            kinetic_data = getattr(self, "kinetic_data", [])
            potential_data = getattr(self, "potential_data", [])
            total_data = getattr(self, "total_data", [])

            # Add new data point
            time_data.append(message.get("time", 0))
            kinetic_data.append(message.get("kinetic_energy", 0))
            potential_data.append(message.get("potential_energy", 0))
            total_data.append(message.get("total_energy", 0))

            # Store data for future updates
            self.time_data = time_data
            self.kinetic_data = kinetic_data
            self.potential_data = potential_data
            self.total_data = total_data

            # Update the lines
            if (
                self.kinetic_line is not None
                and self.potential_line is not None
                and self.total_line is not None
            ):
                self.kinetic_line.set_data(time_data, kinetic_data)
                self.potential_line.set_data(time_data, potential_data)
                self.total_line.set_data(time_data, total_data)

            # Adjust axes limits
            if time_data:
                ax.set_xlim(0, max(time_data) * 1.1)

                # Determine y-axis limits
                all_data = kinetic_data + potential_data + total_data
                if all_data:
                    ymin = min(all_data) * 0.9
                    ymax = max(all_data) * 1.1
                    ax.set_ylim(ymin, ymax)

        elif self.current_plot_type == "Volume Time Series":
            # Update volume time series plot
            # Get the stored time series data
            time_data = getattr(self, "volume_time_data", [])
            volume_data = getattr(self, "volume_data", [])

            # Add new data point
            time_data.append(message.get("time", 0))
            volume_data.append(message.get("volume", 0))

            # Store data for future updates
            self.volume_time_data = time_data
            self.volume_data = volume_data

            # Update the line
            if self.volume_line is not None:
                self.volume_line.set_data(time_data, volume_data)

            # Adjust axes limits
            if time_data:
                ax.set_xlim(0, max(time_data) * 1.1)

                if volume_data:
                    # Calculate initial volume
                    initial_volume = volume_data[0] if volume_data else 0

                    # Calculate relative change from initial volume
                    rel_change = [
                        (v - initial_volume) / initial_volume for v in volume_data
                    ]

                    # Set y-axis limits based on relative change
                    max_rel_change = (
                        max(abs(min(rel_change)), abs(max(rel_change)))
                        if rel_change
                        else 0
                    )

                    if max_rel_change < 1e-10:
                        # If volume is nearly constant, show a narrow range around the initial value
                        ax.set_ylim(initial_volume * 0.9999, initial_volume * 1.0001)
                    else:
                        # Otherwise show the full range of variation
                        ymin = min(volume_data) * 0.99
                        ymax = max(volume_data) * 1.01
                        ax.set_ylim(ymin, ymax)

        # Store frame for potential animation export
        self.store_frame_for_animation(message)

        # Update the canvas
        self.canvas.draw()

    def update_timeseries_from_diag(self, diag_msg: dict):
        """Update time-series lines from a diag.global message.

        Expected keys: t (time), energy (total), enstrophy, volume.
        """
        ax = self.ax
        if ax is None:
            return
        t = diag_msg.get("t", 0.0)
        total = diag_msg.get("energy", None)
        volume = diag_msg.get("volume", None)

        # Energy plot: we only have total energy here; add it if the plot is active
        if self.current_plot_type == "Energy Time Series" and total is not None:
            time_data = getattr(self, "time_data", [])
            kinetic_data = getattr(self, "kinetic_data", [])
            potential_data = getattr(self, "potential_data", [])
            total_data = getattr(self, "total_data", [])

            time_data.append(t)
            # Keep kinetic/potential as previous values or 0 when unknown
            kinetic_data.append(kinetic_data[-1] if kinetic_data else 0.0)
            potential_data.append(potential_data[-1] if potential_data else 0.0)
            total_data.append(total)

            self.time_data = time_data
            self.kinetic_data = kinetic_data
            self.potential_data = potential_data
            self.total_data = total_data

            if self.kinetic_line is not None:
                self.kinetic_line.set_data(time_data, kinetic_data)
            if self.potential_line is not None:
                self.potential_line.set_data(time_data, potential_data)
            if self.total_line is not None:
                self.total_line.set_data(time_data, total_data)

            if time_data:
                ax.set_xlim(0, max(time_data) * 1.1)
                all_data = kinetic_data + potential_data + total_data
                if all_data:
                    ymin = min(all_data) * 0.9
                    ymax = max(all_data) * 1.1
                    ax.set_ylim(ymin, ymax)
            self.canvas.draw()

        # Volume plot
        if self.current_plot_type == "Volume Time Series" and volume is not None:
            time_data = getattr(self, "volume_time_data", [])
            volume_data = getattr(self, "volume_data", [])
            time_data.append(t)
            volume_data.append(volume)
            self.volume_time_data = time_data
            self.volume_data = volume_data
            if self.volume_line is not None:
                self.volume_line.set_data(time_data, volume_data)
            if time_data:
                ax.set_xlim(0, max(time_data) * 1.1)
                if volume_data:
                    initial_volume = volume_data[0]
                    rel_change = [
                        (v - initial_volume) / initial_volume
                        for v in volume_data
                        if initial_volume
                    ]
                    if (
                        not rel_change
                        or max(abs(min(rel_change)), abs(max(rel_change))) < 1e-10
                    ):
                        ax.set_ylim(initial_volume * 0.9999, initial_volume * 1.0001)
                    else:
                        ymin = min(volume_data) * 0.99
                        ymax = max(volume_data) * 1.01
                        ax.set_ylim(ymin, ymax)
            self.canvas.draw()

    def export_current_plot(self, file_path):
        """
        Export the current plot to a file with scientific publication quality.

        Args:
            file_path: Path to save the plot to
        """

        # Determine format from file extension
        file_ext = file_path.lower().split(".")[-1]

        # Set format-specific options for publication quality
        if file_ext == "png":
            # High-resolution raster for digital publications
            dpi = 300
            format_options = {
                "facecolor": "white",
                "edgecolor": "none",
                "bbox_inches": "tight",
                "pad_inches": 0.1,
                "metadata": {
                    "Title": f"SHEL {self.current_plot_type}",
                    "Author": "SHEL Model",
                    "Description": f"Shallow water model {self.current_plot_type.lower()} visualization",
                },
            }
        elif file_ext == "eps":
            # Vector format for print publications
            dpi = 300  # EPS can handle high DPI
            format_options = {
                "facecolor": "white",
                "edgecolor": "none",
                "bbox_inches": "tight",
                "pad_inches": 0.1,
                "format": "eps",
            }
            # Ensure fonts are embedded for vector output
            plt.rcParams["ps.useafm"] = True
            plt.rcParams["pdf.use14corefonts"] = True
        elif file_ext == "pdf":
            # Vector format alternative
            dpi = 300
            format_options = {
                "facecolor": "white",
                "edgecolor": "none",
                "bbox_inches": "tight",
                "pad_inches": 0.1,
                "metadata": {
                    "Title": f"SHEL {self.current_plot_type}",
                    "Author": "SHEL Model",
                    "Subject": f"Shallow water model {self.current_plot_type.lower()} visualization",
                },
            }
        else:
            # Default high-quality settings
            dpi = 300
            format_options = {
                "facecolor": "white",
                "edgecolor": "none",
                "bbox_inches": "tight",
                "pad_inches": 0.1,
            }

            # Save with publication-quality settings
        self.figure.savefig(file_path, dpi=dpi, **format_options)

        # Reset matplotlib rcParams if we changed them
        if file_ext == "eps":
            plt.rcParams["ps.useafm"] = False
            plt.rcParams["pdf.use14corefonts"] = False

    def export_animation(self, file_path, frame_data=None, fps=10):
        """
        Export an animation of the time series data to MP4.

        Args:
            file_path: Path to save the animation to (should end with .mp4)
            frame_data: List of message dictionaries for each frame, or None to use stored data
            fps: Frames per second for the animation
        """
        import matplotlib.animation as animation

        if not frame_data and not hasattr(self, "_animation_frames"):
            logger.warning("No frame data available for animation export")
            return

        # Use provided frame data or stored frames
        frames = frame_data or getattr(self, "_animation_frames", [])

        if not frames:
            logger.warning("No frames to animate")
            return

        # Create a temporary figure for animation
        temp_fig = Figure(figsize=(8, 6), dpi=150)

        def animate_frame(frame_idx):
            """Animation function for each frame."""
            temp_fig.clear()
            ax = temp_fig.add_subplot(111)

            message = frames[frame_idx]
            fields = message.get("fields", {})

            if self.current_plot_type == "Water Elevation":
                eta = np.array(fields.get("eta", []))
                if eta.size > 0:
                    im = ax.imshow(
                        eta,
                        cmap="coolwarm",
                        interpolation="bilinear",
                        origin="lower",
                        aspect="equal",
                    )
                    ax.set_title(
                        f"Water Elevation (t = {message.get('time', 0):.2f} s)"
                    )
                    ax.set_xlabel("X")
                    ax.set_ylabel("Y")
                    plt.colorbar(im, ax=ax, label="Elevation (m)")

            elif self.current_plot_type == "Velocity Field":
                u = np.array(fields.get("u", []))
                v = np.array(fields.get("v", []))
                if u.size > 0 and v.size > 0:
                    ny, nx = u.shape
                    X, Y = np.meshgrid(np.arange(nx), np.arange(ny))
                    ax.quiver(X, Y, u, v, scale=0.2)
                    ax.set_aspect("equal")
                    ax.set_title(f"Velocity Field (t = {message.get('time', 0):.2f} s)")
                    ax.set_xlabel("X")
                    ax.set_ylabel("Y")

            elif self.current_plot_type == "Vorticity":
                # Compute vorticity from u and v
                u = np.array(fields.get("u", []))
                v = np.array(fields.get("v", []))
                if u.size > 0 and v.size > 0:
                    # Simple central difference vorticity
                    dv_dx = np.gradient(v, axis=1)
                    du_dy = np.gradient(u, axis=0)
                    vorticity = dv_dx - du_dy

                    im = ax.imshow(
                        vorticity,
                        cmap="RdBu_r",
                        interpolation="bilinear",
                        origin="lower",
                        aspect="equal",
                    )
                    ax.set_title(f"Vorticity (t = {message.get('time', 0):.2f} s)")
                    ax.set_xlabel("X")
                    ax.set_ylabel("Y")
                    plt.colorbar(im, ax=ax, label="Vorticity (s⁻¹)")

            return (ax,)

        # Create animation
        anim = animation.FuncAnimation(
            temp_fig, animate_frame, frames=len(frames), interval=1000 / fps, blit=False
        )

        # Save animation with high quality settings
        writer = animation.FFMpegWriter(
            fps=fps,
            metadata={
                "title": f"SHEL {self.current_plot_type} Animation",
                "artist": "SHEL Model",
                "comment": f"Shallow water model {self.current_plot_type.lower()} time series",
            },
            bitrate=1800,  # High bitrate for quality
            extra_args=["-vcodec", "libx264", "-preset", "slow", "-crf", "22"],
        )

        try:
            anim.save(file_path, writer=writer, dpi=150)
            logger.info("Animation saved to %s", file_path)
        except Exception as e:
            logger.error("Failed to save animation: %s", e)
            # Try alternative writer if FFMpeg fails
            try:
                writer_alt = animation.FFMpegWriter(fps=fps, bitrate=1800)
                anim.save(file_path, writer=writer_alt, dpi=150)
            except Exception as e2:
                logger.error("Alternative animation save also failed: %s", e2)
                raise

    def store_frame_for_animation(self, message):
        """
        Store a frame for potential animation export.

        Args:
            message: Message containing model state
        """
        if not hasattr(self, "_animation_frames"):
            self._animation_frames = []

        # Limit stored frames to prevent memory issues (keep last 1000 frames)
        if len(self._animation_frames) >= 1000:
            self._animation_frames.pop(0)

        self._animation_frames.append(message.copy())
