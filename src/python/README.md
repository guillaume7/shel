# SHEL Python Port

SHEL (SHallow-water numerical modEL) is a finite volume, free-surface, variable bottom, shallow-waters equations numerical solver. This is the Python port of the original MATLAB implementation.

## Features

- Solves the shallow water equations using an Arakawa C-grid staggered mesh
- Implements leapfrog time-stepping scheme
- Supports various boundary conditions (closed, free-slip, radiative)
- Includes several initial condition types for water elevation and bathymetry
- Provides visualization and analysis tools for model results
- Decouples numerical solver from GUI for better performance and flexibility
- Uses modern data formats (NetCDF, Parquet) for efficient I/O
- Implements ZeroMQ publish/subscribe pattern for real-time communication between solver and GUI

## Installation

### Prerequisites

- Python 3.8 or higher
- NumPy, SciPy
- PyQt5 for the GUI
- NetCDF4 for gridded data I/O
- Pyarrow for Parquet support
- PyYAML for configuration
- ZeroMQ for real-time communication

### Setup

1. Clone the repository:
   ```
   git clone https://github.com/your-username/shel.git
   cd shel
   ```

2. Run the setup script to create a virtual environment and install dependencies:
   ```
   bash setup_env.sh
   ```

3. Activate the virtual environment:
   ```
   source shel-env/bin/activate
   ```

## Usage

### Command-Line Interface

To run the model from the command line:

```bash
python src/python/run.py --config path/to/config.yaml --output-dir ./output
```

Command-line options:
- `--config`: Path to the YAML configuration file (required)
- `--output-dir`: Directory for output files (default: ./output)
- `--verbose`: Increase verbosity (can be specified multiple times)
- `--pub-port`: Port for ZeroMQ publisher (default: 5556)
- `--no-output`: Do not write output files (useful for benchmarking)

### Graphical User Interface

To run the model with the GUI:

```bash
python src/python/gui.py
```

The GUI allows you to:
- Configure model parameters
- Set initial and boundary conditions
- Run, pause, and stop the model
- Visualize results in real-time
- Export plots and animations

## Configuration

The model is configured using YAML files. Here's an example configuration:

```yaml
grid:
  nx: 100
  ny: 100
  dx: 100.0
  dy: 100.0

model:
  dt: 1.0
  num_steps: 1000
  output_interval: 10
  solver: leapfrog

initial_conditions:
  type: gaussian
  gaussian:
    amplitude: 1.0
    x0: 5000.0
    y0: 5000.0
    sigma: 500.0

bathymetry:
  type: flat
  flat:
    depth: 100.0

boundary_conditions:
  north: radiative
  south: radiative
  east: radiative
  west: radiative

physical_parameters:
  viscosity: 1.0
  bottom_friction: 0.002
  coriolis: 0.0
  gravity: 9.81

output:
  directory: ./output
  enabled: true

communication:
  zmq_pub_port: 5556
  enable_zmq: true
```

## Project Structure

```
src/python/
├── run.py                  # CLI entry point
├── gui.py                  # GUI entry point
└── shel/
    ├── __init__.py
    ├── cli.py              # Command line interface
    ├── config/             # Configuration handling
    │   └── config_loader.py
    ├── gui/                # GUI components
    │   ├── __init__.py
    │   ├── main_window.py  # Main GUI window
    │   ├── model_control.py # Model control panel
    │   ├── parameters.py   # Parameter editing panel
    │   ├── utils.py        # Utility functions
    │   └── visualization.py # Plotting components
    ├── io/                 # I/O components
    │   ├── netcdf_reader.py
    │   └── parquet_reader.py
    └── model/              # Core model components
        ├── __init__.py
        ├── grid.py         # Arakawa C-grid implementation
        ├── model_runner.py # Model execution coordinator
        ├── state.py        # Model state management
        ├── boundary_conditions/ # Boundary condition implementations
        │   └── boundary.py
        ├── initial_conditions/ # Initial condition implementations
        │   ├── bottom.py
        │   └── waterlevel.py
        └── solvers/        # Numerical solver implementations
            ├── base.py
            ├── factory.py
            └── leapfrog.py
```

## Development

### Running Tests

To run the test suite:

```bash
pytest tests/python/
```

### Adding New Components

#### Adding a New Solver

1. Create a new file in `src/python/shel/model/solvers/` that extends the `BaseSolver` class
2. Implement the required methods
3. Add the solver to the factory in `solvers/factory.py`

#### Adding a New Initial Condition Type

1. Extend the appropriate factory method in `initial_conditions/waterlevel.py` or `initial_conditions/bottom.py`
2. Add UI components to the parameter panel in `gui/parameters.py`

## License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
