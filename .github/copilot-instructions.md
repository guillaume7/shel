# GitHub Copilot Instructions for SHEL (SHallow-water numerical modEL)

## Project Overview
SHEL (SHallow-water numerical modEL) is a finite volume, free-surface, variable bottom, shallow-waters equations numerical solver implemented in MATLAB. It includes a built-in graphical interface for:
- Loading, editing, and saving simulation parameters and forcings
- Running simulations
- Visualizing results
- Exporting images (eps, png) and movies (avi)

The model is designed to be compact, efficient, and extensible, allowing developers to replace core solver files with custom numerical schemes or contribute to the stack of available numerical schemes.

## Core Mathematical Model
The model solves the shallow water equations (SWE) which describe 2D barotropic motion of water masses. These equations include:
- Momentum equations (x and y components)
- Continuity equation
- Various forcings (Coriolis, wind stress, bottom stress)

### Key Components
- Uses Arakawa C-grid staggered mesh for numerical stability
- Implements leapfrog and central differences schemes (second-order accurate in time and space)
- Supports various boundary conditions (Dirichelet, Neumann, Sommerfeld)
- Includes variable bottom bathymetry
- Incorporates viscosity and drag coefficients

## Code Organization
The codebase is organized as follows:

### Main Entry Points
- `src/matlab/run.m` - Main entry point for the model
- `src/matlab/gui/GUI_ControlPanel2D.m` - Main GUI implementation

### Model Implementation
- `src/matlab/model/model_handles.m` - Core model functionality (primary working implementation)
- `src/matlab/model/run_modular.m` - Modular implementation of the model (refactoring in progress, not fully tested)
- `src/matlab/model/initial_conditions/` - Initial condition implementations (part of refactoring)
- `src/matlab/model/solvers/` - Numerical solver implementations (part of refactoring)
- `src/matlab/model/boundary_conditions/` - Boundary condition implementations (part of refactoring)
- `src/matlab/model/grid/` - Grid-related functionality (part of refactoring)
- `src/matlab/model/state/` - Model state management (part of refactoring)
- `src/matlab/model/outputs/` - Output and visualization functions (part of refactoring)

> **Important Note**: Currently, `model_handles.m` is the de-facto working implementation of the numerical model. The other files represent an ongoing refactoring effort aimed at better structuring the codebase, but they have not been fully tested. The ultimate goal is to port both the GUI and the restructured model to Python to enable running the model without requiring a MATLAB license.

## Key Features
1. **Multiple Numerical Schemes**: Implements the shallow water equations using an Arakawa C-grid with a leapfrog and central differences scheme, and provides extensibility for other schemes.

2. **Boundary Conditions**: Supports various types including:
   - Null-flux (land mask)
   - No-slip
   - Radiative (Sommerfeld and Flather methods)

3. **Validation Test Cases**: Includes test cases like Gaussian bell-shaped initial conditions to validate:
   - Conservation of volume
   - Conservation of momentum
   - Energy conservation and dissipation
   - Vorticity analysis

4. **Visualization**: Built-in GUI for real-time visualization of:
   - Water elevation
   - Velocity fields
   - Energy metrics
   - Vorticity

## GUI Functionality
The GUI allows users to:
- Set simulation parameters (grid size, time step, etc.)
- Choose initial conditions (Gaussian elevation, flat, etc.)
- Select boundary conditions
- Adjust physical parameters (viscosity, Coriolis parameter, etc.)
- Run/pause/stop simulations
- Visualize results in real-time
- Export results as images or videos

## Common Development Tasks
When working with this codebase, you might need to:

1. **Implement new numerical schemes**: Add new files to the solvers directory following the existing patterns.

2. **Add new initial conditions**: Extend the initial condition factories with new types.

3. **Modify boundary conditions**: Update or add new boundary condition implementations.

4. **Add visualization options**: Extend the GUI plotting capabilities for new metrics or views.

5. **Optimize performance**: Look for computational bottlenecks in the main solution loops.

## Key Mathematical Concepts
Important concepts to understand when working with this code:

1. **Shallow Water Equations**: The mathematical foundation that assumes horizontal length scales are much larger than vertical ones.

2. **Arakawa C-grid**: A staggered grid system where variables are evaluated at different points to improve numerical stability.

3. **Leapfrog Scheme**: A time-stepping method that uses information from two previous time steps.

4. **Conservation Properties**: The model tracks conservation of volume, momentum, and monitors energy dissipation.

5. **Okubo-Weiss Scalar**: Used for analyzing vorticity in the flow.

## Implementation Notes
- The code uses MATLAB's object-oriented features and factory patterns
- Modular design allows for easy replacement of components
- GUI is implemented using MATLAB's GUI development environment (GUIDE)
- The model implements various numerical filters for stability
- Performance optimizations include vectorized operations where possible

## Future Development Direction
- The primary development goal is to port the model to Python
- This includes porting both the numerical model and the GUI components
- The refactored structure (in progress) will serve as the basis for the Python implementation
- Maintaining the same mathematical accuracy and validation properties is critical
- The Python version should maintain feature parity with the MATLAB version
- The port should leverage appropriate Python libraries for numerical computation (NumPy, SciPy) and visualization (Matplotlib, possibly PyQt or other GUI frameworks)

## Scientific Background
The model is based on well-established oceanographic and hydrodynamic principles, suitable for:
- Educational purposes in oceanography and fluid dynamics
- Research on numerical schemes for shallow water equations
- Simple simulations of coastal and ocean dynamics

This implementation follows scientific approaches documented in:
- Kantha and Clayson (2000)
- Arakawa (1966)
- Gill (1982)
- Other standard references in physical oceanography
