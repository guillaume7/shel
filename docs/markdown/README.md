# SHEL, a SHallow water equations modEL: Scientific Documentation

**Author: Guillaume Riflet**

## Table of Contents

### [Abstract and Keywords](swe-abstract.md)
- [Abstract](swe-abstract.md#abstract)
- [Keywords](swe-abstract.md#keywords)

### [Part 1: Model Description](swe-part1.md)
- [Introduction](swe-part1.md#introduction)
- [Mathematical Model of the Shallow-Water Equations](swe-part1.md#mathematical-model-of-the-shallow-water-equations)
  - [The Mathematical Model](swe-part1.md#the-mathematical-model)
  - [The Mesh](swe-part1.md#the-mesh)
  - [Boundary Conditions](swe-part1.md#boundary-conditions)
    - [Null-flux](swe-part1.md#null-flux)
    - [No-slip](swe-part1.md#no-slip)
    - [Radiative Boundary Conditions](swe-part1.md#radiative-boundary-conditions)
  - [The Numerical Scheme](swe-part1.md#the-numerical-scheme)

### [Part 2: Model Validation](swe-part2.md)
- [Validation](swe-part2.md#validation)
  - [Gaussian Bell-Shaped Geometry](swe-part2.md#gaussian-bell-shaped-geometry)
  - [Energy](swe-part2.md#energy)
  - [Geometric and Similitude Considerations](swe-part2.md#geometric-and-similitude-considerations)
  - [Basic Results](swe-part2.md#basic-results)

### [Part 3: Advanced Analysis](swe-part3.md)
- [Energy Decay Study](swe-part3.md#energy-decay-study)
- [Radiation Boundary Condition](swe-part3.md#radiation-boundary-condition)
- [Geostrophic Equilibrium](swe-part3.md#geostrophic-equilibrium)
- [Applying the Okubo-Weiss Scalar to Assess the Open-Boundary Condition](swe-part3.md#applying-the-okubo-weiss-scalar-to-assess-the-open-boundary-condition)

### [Conclusions and Bibliography](swe-references.md)
- [Conclusions](swe-references.md#conclusions)
- [References](swe-references.md#references)

## Document Overview

This documentation is split into multiple parts for better readability:

1. **Abstract and Keywords** provides a brief overview of the model and its purpose.
2. **Part 1** covers the mathematical and numerical fundamentals of the model, including the equations, mesh design, and boundary conditions.
3. **Part 2** focuses on model validation using a Gaussian bump test case, analyzing conservation properties and energy behavior.
4. **Part 3** explores advanced topics like energy decay, radiation boundary conditions, geostrophic equilibrium, and the application of the Okubo-Weiss scalar.
5. **Conclusions and Bibliography** contains the overall conclusions of the study and all references cited throughout the documentation.

All equations are rendered using standard Markdown equation syntax compatible with both GitHub and VS Code viewers.

# SHEL Python Port

SHEL (SHallow-water numerical modEL) is a finite volume, free-surface, variable bottom, shallow-water equations solver. This Python port preserves the MATLAB model's numerical fidelity and extends it with modular architecture and modern tooling.

## Features
- Arakawa C-grid staggered mesh
- Leapfrog and central difference schemes
- Multiple boundary condition strategies
- Conservation diagnostics (volume, energy, enstrophy, potential enstrophy)
- Modular solver and IO subsystems
- Built-in regression and conservation tests
- Extensible for new schemes and diagnostics

## Solver Documentation
- Operators: pressure, advection, diffusion, friction, continuity
- Time integration: leapfrog, Asselin filter (see `src/python/shel/model/solvers/time/`)
- Boundary conditions: strategy/registry pattern
- Diagnostics: field and integrated metrics

## IO Documentation
- NetCDF: grid, state, diagnostics snapshots (`model/outputs/writers/netcdf.py`)
- Parquet: time-series diagnostics (`model/outputs/writers/parquet.py`)
- JSON: deterministic serialization for regression (`model/outputs/writers/json_writer.py`)
- OutputManager: schedules output hooks
- YAML: configuration for reproducible runs

## Developer Guide
See `DEV_GUIDE.md` for best practices, workflow, and contribution checklist.

## Getting Started
1. Clone the repository
2. Install dependencies: `pip install -r requirements.txt`
3. Run tests: `pytest`
4. Explore examples in `examples/python/`

## Documentation
- Full API docs and guides in `docs/markdown/`
- Status table in `.github/prompts/MASTER_PORT_PROMPT.md`

---
For scientific background and implementation notes, see `.github/copilot-instructions.md`.
