# SHEL - SHallow-waters numerical modEL

<div align="center">
  <img src="docs/markdown/figs/bump-4.svg" alt="Water level simulation" width="600"/>
  <p><em>Simulation of water elevation propagation in a shallow water domain</em></p>
</div>

## What is SHEL?

The SHEL (SHallow-waters numerical modEL) is a finite volume, free-surface, variable bottom, shallow-waters equations numerical solver.

SHEL has an original MATLAB implementation (with a GUI) and an ongoing Python port engineered for modularity, testability, and open execution (no MATLAB license required).

<div align="center">
  <img src="docs/markdown/figs/arakawaCgrid.svg" alt="Arakawa C-grid" width="400"/>
  <p><em>SHEL uses the Arakawa C-grid staggered mesh system for numerical stability</em></p>
</div>

The code is compact, efficient and extensible. In the Python port, numerical kernels are decomposed into pure functions with strategy/registry patterns for boundary conditions and interchangeable algorithms.

SHEL uses an Arakawa C grid type over a land-mask. Time integration includes explicit Euler (for ministeps) and leapfrog + Robert–Asselin filter. Boundary conditions include Closed, Free-slip, Radiative (Sommerfeld), Flather, and Dirichlet (eta). Tracer BCs include Closed and Radiative.

## Repository Structure

The repository is organized as follows:

- `src/matlab/`: MATLAB source code (original, feature-complete)
  - `run.m`: The main entry point for running the model
  - `data/`: Input data and simulation results
  - `gui/`: Graphical user interface components
  - `model/`: Core implementation of the numerical model
- `src/python/`: Python port (active development)
  - `shel/model/solvers/time/`: Time drivers (`asselin_filter`, `leapfrog_stepper`, `leapfrog_step_with_config`)
  - `shel/model/solvers/common/`: Tendencies, ministep, and config-aware step helpers
  - `shel/model/boundary_conditions/`: Strategy-based BCs (momentum, waterlevel, tracer) + registry
  - `shel/model/initial_conditions/`, `grid/`, `diagnostics/`, `state/`, `outputs/` (modular domains)
- `examples/python/`: Minimal runnable Python examples
- `tests/python/`: Unit/integration tests for the Python port
- `docs/`: Documentation
  - `latex/`: Original LaTeX documentation and figures
  - `markdown/`: Converted markdown documentation
- `COPYING`: License information

## Python Quickstart

Requirements: Python 3.10+.

Install (editable):

```bash
pip install -e .
```

Run tests:

```bash
pytest -q
```

### Time API (Python)

- Explicit ministep (Euler): `shel.model.solvers.common.ministep.explicit_step(eta, H, U, V, ...)`.
- Config-aware explicit step: `shel.model.solvers.common.stepper.explicit_step_with_config(eta, H, U, V, config=...)`.
- Leapfrog + Asselin:
  - `from shel.model.solvers.time import leapfrog_stepper, leapfrog_step_with_config`
  - `leapfrog_stepper` advances one step given (n-1, n) states; `leapfrog_step_with_config` also applies per-side BCs and sponge according to a config mapping.

Array staggering (C-grid): `eta: (ny,nx)`, `U: (ny,nx+1)`, `V: (ny+1,nx)`, `H: (ny,nx)`.

### Boundary Conditions (config)

Set per-side types and options under a config dict. Example (Flather east + cosine sponge):

```python
cfg = {
  "boundary_conditions": {"west": "closed", "east": "flather", "south": "closed", "north": "closed"},
  "boundary_eta_ext": {"east": eta_ext},     # external eta along open side
  "boundary_eta_relax": 0.5,                  # 0..1 blend on eta
  "boundary_momentum_relax": 0.3,             # 0..1 blend on Flather correction
  "eta_bc_stage": "post",                     # apply eta BCs post-step (default)
  "sponge": {"enabled": True, "width": 4, "alpha": 0.2, "taper": "cosine", "apply_to": "both"},
}
```

Supported names: `closed`, `freeslip` (momentum); `radiative`, `flather`, `dirichlet` (eta); `closed`, `radiative` (tracer).

### Examples (Python)

- Flather + sponge conservative config: `examples/python/flather_sponge_config_example.py`
- Leapfrog with Flather (east), Dirichlet (north), sponge: `examples/python/leapfrog_flather_dirichlet_sponge_example.py`

Run an example:

```bash
python examples/python/leapfrog_flather_dirichlet_sponge_example.py
```

## Documentation

Comprehensive documentation is available in the [docs/markdown](docs/markdown) directory, including:
- [Abstract and Keywords](docs/markdown/swe-abstract.md)
- [Part 1: Model Fundamentals](docs/markdown/swe-part1.md)
- [Part 2: Model Implementation](docs/markdown/swe-part2.md)
- [Part 3: Validation and Results](docs/markdown/swe-part3.md)
- [Conclusions and References](docs/markdown/swe-references.md)

<div align="center">
  <img src="docs/markdown/figs/radiate-coriolis-energy.svg" alt="Energy conservation" width="600"/>
  <p><em>SHEL tracks energy conservation during simulations, showing kinetic, potential, and total energy</em></p>
</div>

## How to Use (MATLAB)

1. Open MATLAB
2. Set the workfolder to the `src/matlab` directory of the SHEL repository
3. Type `run` and press enter

## How to Cite

If you use SHEL in your work, please cite the scientific documentation as follows:

```
Riflet, G., 2010. SHEL, a Shallow-Water Numerical 
Model: Scientific Documentation. Instituto Superior Técnico, 
Universidade Técnica de Lisboa.
```

## Notes & Contact

- The legacy class-based Python solver `shel.model.solvers.leapfrog.LeapfrogSolver` is deprecated; use the time API in `shel.model.solvers.time` instead.
- Issues and contributions welcome via GitHub.

- Email: guillaume.riflet at gmail.com
- Last Updated: 2025-09-06

<div align="center">
  <img src="docs/markdown/figs/radiate-coriolis-velocity-modulus-sam2p.svg" alt="Velocity field with Coriolis effect" width="600"/>
  <p><em>Visualization of velocity field modulus depicting the evolution of motion from a gaussian elevation as the initial waterlevel</em></p>
</div>
