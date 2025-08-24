# SHEL Python Port Guidelines

## Overview
This document provides guidelines for porting the SHEL (SHallow-water numerical modEL) from MATLAB to Python 3. The port should include both the GUI interface and the numerical solver components, with a focus on the restructured parts of the codebase.

## Core Requirements

### Code Structure and Design
- Follow clean code principles throughout the implementation
- Apply the DRY (Don't Repeat Yourself) principle rigorously
- Implement factory patterns to maintain and enhance extensibility for custom numerical solvers
- Maintain the mathematical accuracy and validation properties of the original model
- Preserve all functionality from the MATLAB version
- Structure code in a modular, maintainable way with clear separation of concerns
- Decouple the GUI from the numerical solver to allow headless CLI operation
- Design a clean, well-defined API between the solver and GUI components

### Performance Considerations
- Leverage NumPy/SciPy for vectorized operations to maximize performance
- Optimize computational bottlenecks, especially in the main solution loops
- Profile code to identify and address performance issues
- Consider using Numba for performance-critical sections if needed

### GUI Implementation
- Implement the GUI using a modern Python framework (options include PyQt, Tkinter, or other suitable libraries)
- Maintain all visualization capabilities from the original MATLAB GUI
- Ensure responsive real-time visualization of simulation results
- Support all parameter configuration options from the original interface
- Enable export of results in various formats (PNG, EPS, MP4/AVI)

### Numerical Solver Implementation
- Port the Arakawa C-grid staggered mesh system with proper indexing
- Implement the leapfrog and central differences schemes
- Support all boundary conditions (Dirichelet, Neumann, Sommerfeld)
- Correctly handle variable bottom bathymetry
- Maintain implementation of viscosity and drag coefficients
- Preserve conservation properties (volume, momentum, energy)

## Development Standards

### Code Quality
- Use consistent code formatting (PEP 8 compliance)
- Add comprehensive docstrings (NumPy or Google style) for all functions, classes, and modules
- Keep functions small and focused on a single responsibility
- Use type hints throughout the codebase
- Apply proper naming conventions for variables, functions, and classes
- Use static typing strictly. Strong typing. Every argument and variable MUST have its type specified.
- Use pylint after each new file and fix it until it achieves a good score
- Trim whitespaces. Apply black formatter with line length of 88 characters.
- Avoid using Any type in type hints; use Union types instead when needed.

### Project Configuration and Environment
- Use a Python virtual environment named 'venv' for development
- Configure mypy.ini for strict type checking with appropriate settings
- Set up .pylintrc for consistent code quality enforcement
- Create .editorconfig to enforce consistent coding styles across editors
- Configure .vscode settings for proper code navigation and debugging
- Define pyproject.toml with project metadata and tool configurations
- Maintain an up-to-date setup.py for package installation
- Create .gitignore to exclude appropriate files from version control
- Use .gitattributes to enforce LF line endings and proper file handling
- Place all support and automation scripts in the devops/scripts/ folder
- Set up proper pytest configuration in pyproject.toml or conftest.py

### DevOps and CI/CD
- Include comprehensive CI/CD configurations
- Set up GitHub Actions for automated testing and linting
- Configure pre-commit hooks for quality enforcement
- Ensure all tests run automatically on push/PR
- Validate type checking as part of the CI process
- Create automation scripts for common tasks (place in devops/scripts/)

### Error Handling and Logging
- Implement robust exception handling throughout the codebase
- Use appropriate custom exceptions where needed
- Set up a clean, configurable logging system
- Include appropriate debug, info, warning, and error messages
- Ensure errors are handled gracefully with informative user feedback

### Testing
- Develop comprehensive unit tests for all components
- Include integration tests for the overall system
- Create validation tests to verify mathematical correctness
- Add regression tests for the test cases from the original model
- Ensure tests verify conservation properties

## Implementation Approach
1. Start by porting the core numerical solver (`model_handles.m`)
2. Refactor the code following the modular structure in `run_modular.m`
3. Implement and test the boundary conditions
4. Develop the visualization components
5. Build the GUI interface
6. Integrate all components
7. Validate against test cases

## Data Formats and Communication

### Data Storage Strategy
- Use NetCDF for structured gridded data (bathymetry, initial conditions, simulation results)
- Use YAML for configuration parameters and model settings
- Use Parquet for timeseries data (metrics, point measurements over time)
- Ensure all formats support appropriate metadata

### CLI-GUI Communication
- Implement a publish-subscribe pattern using ZeroMQ for real-time data streaming
- CLI solver publishes state updates during simulation
- GUI subscribes to updates for real-time visualization
- Design protocol to handle different update frequencies and data types
- Support multiple GUI clients connecting to a single solver instance

### Project Structure
```
shel/
├── README.md
├── pyproject.toml
├── setup.py
├── mypy.ini
├── .pylintrc
├── .gitignore
├── .gitattributes
├── .editorconfig
├── requirements.txt
├── .vscode/
│   ├── settings.json
│   └── launch.json
├── devops/
│   └── scripts/
│       ├── setup_env.sh
│       ├── create_init_files.sh
│       └── ... (other utility scripts)
├── docs/
│   ├── latex/
│   └── markdown/
├── src/
│   ├── matlab/      # Original MATLAB code
│   └── python/
│       ├── shel/
│       │   ├── __init__.py
│       │   ├── cli.py        # Command-line interface
│       │   ├── config/       # Configuration handling
│       │   │   └── __init__.py
│       │   ├── io/           # File I/O operations
│       │   │   └── __init__.py
│       │   ├── model/        # Core numerical model
│       │   │   ├── __init__.py
│       │   │   ├── grid.py
│       │   │   ├── state.py
│       │   │   ├── model_runner.py
│       │   │   ├── solvers/
│       │   │   │   └── __init__.py
│       │   │   ├── boundary_conditions/
│       │   │   │   └── __init__.py
│       │   │   └── initial_conditions/
│       │   │       └── __init__.py
│       │   ├── utils/        # Utility functions
│       │   │   └── __init__.py
│       │   └── gui/          # GUI implementation
│       │       ├── __init__.py
│       │       ├── main_window.py
│       │       ├── model_control.py
│       │       ├── parameters.py
│       │       ├── visualization.py
│       │       ├── utils.py
│       │       ├── widgets/
│       │       │   └── __init__.py
│       │       └── visualizers/
│       │           └── __init__.py
│       ├── run.py            # CLI entry point
│       └── gui.py            # GUI entry point
└── tests/
    └── python/
        ├── __init__.py
        ├── conftest.py
        ├── test_grid.py
        ├── test_state.py
        ├── test_gui.py
        └── ... (other test files)
```

## Dependencies
- NumPy/SciPy for numerical calculations
- Matplotlib for plotting and visualization
- PyQt5 for GUI implementation
- pytest for testing
- pylint for code quality
- mypy for static type checking
- black for code formatting
- isort for import sorting
- logging for structured logs
- netCDF4/xarray for gridded data storage and manipulation
- PyYAML for configuration file handling
- pyarrow for Parquet file support
- ZeroMQ (pyzmq) for publish-subscribe communication

## Configuration Files

### pyproject.toml
The pyproject.toml file should include:
- Build system requirements
- Project metadata
- Dependencies with version constraints
- Development dependencies
- Tool configurations (black, isort, pytest)
- Mypy configuration options

### mypy.ini
Configure with strict type checking:
- python_version = 3.8+
- disallow_untyped_defs = true
- disallow_incomplete_defs = true
- check_untyped_defs = true
- disallow_untyped_decorators = true
- no_implicit_optional = true
- strict_optional = true
- warn_redundant_casts = true
- warn_unused_ignores = true
- warn_no_return = true
- warn_unreachable = true
- Appropriate overrides for third-party libraries

### .pylintrc
Include comprehensive linting rules:
- Enforce type checking
- Configure naming conventions
- Set appropriate max-line-length (88 characters)
- Enable relevant error checks
- Disable false positives for specific use cases
- Configure proper Python path handling

### .editorconfig
Define consistent coding styles:
- Set indentation to 4 spaces for Python
- Use LF line endings
- Trim trailing whitespace
- Insert final newline
- Set character encoding to UTF-8
- Configure specific rules for different file types

### .vscode/
Include VS Code configuration files:
- settings.json: Configure Python path, linters, formatters
- launch.json: Set up debugging configurations for different components
- extensions.json: Recommend useful extensions for the project

### Virtual Environment
Use a standard virtual environment structure:
- Place in project root as 'venv/'
- Use Python 3.8+ as the baseline
- Install all dependencies in the venv
- Document activation procedures in README.md
- Include scripts for setting up the environment in devops/scripts/

### Scripts
All support scripts should be placed in the devops/scripts/ directory:
- Environment setup scripts
- Utility scripts for code generation
- Test runners
- Linting and formatting automation
- Documentation generation scripts

## Documentation
- Include comprehensive README with installation and usage instructions
- Provide API documentation for all major components
- Add examples demonstrating key features
- Include developer documentation explaining the code structure
- Document mathematical background and implementation details
- Update documentation when code changes

Remember that the goal is to create a clean, maintainable, and extensible Python implementation that preserves all the capabilities of the original MATLAB version while following modern Python development practices.
