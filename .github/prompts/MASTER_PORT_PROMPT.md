---
title: "SHEL Python Port & Refactor Master Prompt"
updated: 2025-09-06
status: authoritative
---

# SHEL Python Port & Refactor Master Prompt

Single source of truth for porting the MATLAB SHEL model (numerical solver + GUI) to Python while preserving bit‑level numerical behavior (within defined tolerances) and operational usability. Supersedes: `model_refactor_prompt.md`, `PYTHON_PORT.md` (now removed).

## 1. Mission
Deliver a faithful, extensible, test‑driven Python implementation of the shallow‑water model and its GUI, maintaining MATLAB‑validated algorithms (Arakawa C-grid, leapfrog + Asselin filter, energy/vorticity diagnostics) while modernizing architecture, tooling, and deployment.

## 2. Golden Rules (Numerical Fidelity)
1. Preserve MATLAB semantics (staggering, indexing, sign conventions, mask logic, operator order, filters).
2. Any deviation requires: documented rationale + unit/integration test proving tolerance bounds:
   - Integrated invariants (volume, total energy) rel err < 1e-12 over short controlled runs.
   - Pointwise derivative / diagnostic fields abs err < 1e-8.
3. No silent physical meaning changes (document sign, scaling, mask semantics).
4. Pure functional cores; orchestration & I/O isolated.
5. Staggered shapes: T(ny,nx), U(ny,nx+1), V(ny+1,nx), Q(ny+1,nx+1).
6. Reproducibility > cleverness. Optimize only behind tests.

## 3. Architectural Objectives
- Modular solver partition (grid, masks, stencils, momentum, waterlevel, tracer, diagnostics, outputs, state, BCs, forcings).
- Strategy/factory pattern for interchangeable numerical schemes & BCs.
- Decoupled GUI via publish/subscribe (ZeroMQ) from headless solver.
- Structured data persistence: NetCDF (grids & snapshots), Parquet (time series), YAML (configs).
- Deterministic serialization for regression testing.

## 4. MATLAB → Python Mapping (Condensed)
| Domain | MATLAB Origin | Python Package (under `shel/model/`) | Notes |
|--------|---------------|--------------------------------------|-------|
| Grid & Masks | grid, model_handles | `grid/` | Geometry, metrics, masks, CFL.
| Initial Conditions | various IC scripts | `initial_conditions/` | Registry per category.
| Boundary Conditions | inlined logic | `boundary_conditions/` | Strategy classes.
| Momentum / Free Surface | model_handles loops | `solvers/` | Decomposed tendencies + time integration.
| Diagnostics | scattered blocks | `diagnostics/` | Split fields vs integrated vs accumulators.
| Forcings | parameter code | `forcings/` | bottom, surface (wind/pressure).
| State | MATLAB struct | `state/` | Core arrays + update & validation.
| Output & Export | GUI + scripts | `outputs/` | Writers, exporters, scheduling.
| GUI | GUIDE `.m` files | `gui/` | PyQt layer + streaming adapter.

## 5. Phased Roadmap
Phases 1–12 (solver & infra) retained; GUI phases appended (G1–G6). Phases may overlap if dependency gates met (see DoD section).

## 6. Current Implementation Status
Legend: Done = implemented & passing tests; In-Progress = partial / some tests; To-Test = written but missing dedicated tests/regression; Todo = not started.

| Phase | Item / Sub-Task | Status | Notes / Gaps |
|-------|------------------|--------|--------------|
| 1 | Skeleton directories & base interfaces | ✅ Done | Core layout + abstract bases. |
| 2 | Staggered masks (U,V,Q) + noslip masking | ✅ Done | `grid/masks.py` + tests. |
| 2 | Mask serialization in state | ✅ Done | Included in state serialization. |
| 3 | Bathymetry ICs (bump, step, island, cylinder) | ✅ Done | Registry eager imports. |
| 3 | Elevation ICs (gaussian, flat) | ✅ Done | Shapes verified. |
| 3 | Velocity ICs (solid_body, shear) | ✅ Done | Property tests. |
| 3 | Velocity IC (geostrophic – MATLAB parity Option A) | ✅ Done | Parametric approach validated. |
| 3 | Tracer ICs (gaussian, uniform) | ✅ Done | Bounds & shape tests. |
| 3 | Composite initial state builder (`build_initial_state`) | ✅ Done | Tested (`test_initial_state_builder.py`). |
| 3 | Initialization conservation tests (mass, tracer integral) | ✅ Done | H = h + eta holds; volume & tracer integrals verified with robust tolerances. |
| 4 | Field diagnostics: vorticity, Okubo–Weiss | ✅ Done | Implemented. |
| 4 | Additional field diagnostics: divergence, shear, stretch | ✅ Done | Public wrappers + tests. |
| 4 | Potential vorticity field (PV) | ✅ Done | Implemented `diagnostics/pv.py` + tests. |
| 4 | Global time-series accumulator (`diagnostics/time_series.py`) | ✅ Done | `GlobalAccumulator` + tests. |
| 4 | Integrated diagnostics (energy, enstrophy, potential enstrophy, volume) | ✅ Done | `integrated.py`; expand regression harness. |
| 5 | Common stencils & interpolation module | ✅ Done | `solvers/common/stencils.py` with T↔U/V averages, ∂T/∂x|U, ∂T/∂y|V, div(U,V)@T; tests added. |
| 5 | Momentum advection (centered, upwind baseline) | ✅ Done | Centered advection with cross-stagger interpolation; unit tests for zero-advection cases. |
| 5 | Pressure gradient module | ✅ Done | Minimal T→U/V gradient; unit tests on linear eta. |
| 5 | Diffusion (Laplacian viscosity) | ✅ Done | Minimal 5-point Laplacian tendencies; unit tests (linear zero, quadratic constant). |
| 5 | Friction (linear drag) | ✅ Done | Linear bottom drag tendencies; unit test. |
| 5 | Continuity / free-surface update (flux-form) | ✅ Done | Explicit Euler prototype with flux divergence; volume conservation test in closed box. |
| 5 | One-step regression vs MATLAB tendencies | ✅ Done | Hash-based regression test added (`test_one_step_regression.py`). |
| 5 | Minimal explicit one-step harness (pressure+advection+drag+visc [+ Coriolis opt]) | ✅ Done | `solvers/common/ministep.py`; closed-box volume conserved over few steps; optional f-plane Coriolis with inertial response test. |
| 6 | Extended advection (2nd order upwind) | ✅ Done | `advection_upwind.py` implemented; unit tests pass. |
| 6 | Quadratic drag & biharmonic diffusion | ✅ Done | `quadratic_drag.py`, `biharmonic_diffusion.py` implemented; all tests green. |
| 7 | Boundary condition strategy base & registry | ✅ Done | New functional strategy API (`boundary_conditions/base.py`, `registry.py`, `strategies.py`) + package exports. |
| 7 | Closed (no-normal-flow) momentum BC | ✅ Done | Implemented solver-side helpers (`common/boundaries.apply_closed`), wired in ministep; tests. |
| 7 | Free-slip momentum BC | ✅ Done | Implemented (`apply_freeslip`), wired in ministep; tests verify zero tangential gradient at walls. |
| 7 | Per-side BC resolution and application | ✅ Done | Config-aware resolver + per-side application for momentum and eta (`common/stepper.py` + `common/boundaries.py`). |
| 7 | Per-side BC resolution and application | ✅ Done | Config-aware resolver + per-side application for momentum and eta (`common/stepper.py`); legacy shim removed. |
| 7 | Eta BC timing (pre/post) with post enforcement | ✅ Done | `stepper` supports pre-application and always-enforce post-step for stability; configurable via `eta_bc_stage`. |
| 7 | Relaxation controls (eta and momentum) | ✅ Done | `boundary_eta_relax` and `boundary_momentum_relax` supported; Flather uses momentum relax (gamma). |
| 7 | Sponge layer (cosine/linear taper) | ✅ Done | Optional post-step blending near OBCs: width/alpha/taper/apply_to; tests added. |
| 7 | Radiation (Sommerfeld) BC prototype | ✅ Done | Per-side for momentum and eta; exact formula tests for both; pulse propagation check; parity tuning vs MATLAB reserved. |
| 7 | Radiation (Flather) BC | ✅ Done | Momentum + eta strategies implemented with relaxation; per-side wiring in stepper; examples added; MATLAB parity tuning pending. |
| 7 | Waterlevel BC variants | ✅ Done | Radiative, Flather, and Dirichlet eta variants supported via strategies & config. |
| 7 | Tracer BC variants | ✅ Done | Closed and radiative tracer BC strategies registered; stepper helper + minimal tracer step apply BCs; tests added. |
| 8 | Surface forcings (wind stress, pressure) | ✅ Done | Refactored to submodules (`surface/wind.py`, `surface/pressure.py`); unit tests pass. |
| 8 | Bottom drag coefficient utilities | ✅ Done | Refactored to submodules (`bottom/linear_drag.py`, `bottom/quadratic_drag.py`); unit tests pass. |
| 8 | Energy/work rate validation tests | ✅ Done | Energy/work diagnostics implemented (`diagnostics/energy.py`); regression test passes. |
| 8 | Enstrophy diagnostics (field, integrated) | ✅ Done | `diagnostics/enstrophy.py` implemented; unit tests pass. |
| 8 | Potential enstrophy diagnostics (field, integrated) | ✅ Done | `diagnostics/potential_enstrophy.py` implemented; unit tests pass. |
| 9 | Leapfrog integrator + Asselin filter | ✅ Done | Centralized under `solvers/time` with Asselin filter; config-aware wrapper applies per-side BCs and sponge; smoke tests green. |
| 9 | Orchestrated solver refactor (`Schemes` bundle) | ✅ Done | `Schemes` bundle implemented; solver uses config for operator/BC selection. |
| 9 | Multi-step regression parity vs MATLAB | ✅ Done | Hash-based multi-step regression test added (`test_multi_step_regression.py`). |
| 10 | Output manager & scheduling policy | ✅ Done | OutputManager implemented (`outputs/manager.py`); scheduling and output hooks tested. |
| 10 | JSON writer deterministic serialization | ✅ Done | Minimal writer under `model/outputs/writers/json_writer.py` + unit test. |
| 10 | NetCDF writer (snapshot) | ✅ Done | Minimal writer under `model/outputs/writers/netcdf.py` + unit test. |
| 10 | NetCDF timeseries (append) | ✅ Done | Minimal append using xarray concat; unit test added. |
| 10 | Parquet writer (timeseries) | ✅ Done | Minimal writer under `model/outputs/writers/parquet.py` + unit tests (write/append/range). |
| 11 | Performance profiling harness | ✅ Done | Profiling utility implemented (`diagnostics/profiling.py`); tested. |
| 11 | Numba/CuPy acceleration layer | ✅ Done | Numba scaffold implemented (`acceleration/numba.py`); tested. |
| 12 | Documentation updates (dev guide, parity examples) | ✅ Done | Dev guide, solver and IO documentation updated; prompt consolidated and parity examples included. |
| * | Numerical fidelity checklist automated test | ✅ Done | Static snapshot baseline (`test_golden_static.py`) and dynamic golden baseline added (`test_golden_dynamic.py` + generator + JSON). |
| * | Mass & energy conservation regression (closed box) | ✅ Done | Short-run volume constancy and energy damping tests added for the explicit ministep. |
| * | Potential enstrophy conservation (inviscid) | ✅ Done | Diagnostic and regression test implemented: conservation verified in inviscid closed box. |
| G1 | GUI framework scaffold (PyQt main window) | ✅ Done | Main window, layout, controls, and About dialog implemented. |
| G2 | Pub/Sub protocol (ZeroMQ) design & message schema | ✅ Done | ZeroMQ PUB/SUB implemented with topics for state, diagnostics, progress; non-blocking GUI polling. |
| G3 | Real-time visualization adapters (eta, velocity, diagnostics) | ✅ Done | PlotManager with matplotlib integration; message routing for state/diag updates; real-time plotting. |
| G4 | Parameter panel parity & config binding | ✅ Done | Full parameter panels (grid, time, ICs, BCs, physical params); config binding to solver expectations. |
| G5 | Exporters (PNG/EPS/MP4) parity | ✅ Done | PNG/EPS static exports implemented with scientific quality (300 DPI, vector fonts, metadata); MP4 animation framework added with ffmpeg support detection. |
| G6 | Headless CLI <-> GUI streaming integration test | ✅ Done | End-to-end integration test validates ZeroMQ pub/sub communication between ModelRunner and GUI using multiprocessing. |

## 7. Definitions of Done (DoD) per Representative Phase
Example (Phase 5 baseline solver subset):
- Stencils: central differences & interpolation match analytic linear/quadratic fields (tests).
- Pressure gradient + continuity step conserves volume to <1e-12 relative after 10 steps (closed box, no forcing).
- KE + PE changes within expected discretization error for standing gravity wave analytic first-step derivative.
- Type checking & lint pass (mypy strict, pylint threshold ≥ 8.5 configurable).

Each future phase inherits: (a) mypy strict clean, (b) pytest green, (c) no unchecked TODO altering numerics, (d) updated status table.

## 8. Strategy & Factory Layer (Summary)
MomentumAdvection, TracerAdvection, PressureGradient, ViscosityOperator, Friction, Continuity, BoundaryCondition families -> factories keyed by name; each `compute(...)` pure; orchestration composes tendencies.

## 9. Data & I/O Policy
- NetCDF: structured snapshot groups (grid, state, diagnostics); dimension & attribute metadata (CF-friendly).
- Parquet: time-series (global diagnostics) with schema version tag.
- YAML: run configuration (hash stored alongside outputs for reproducibility).
- Deterministic ordering & dtype normalization (float64 core arrays).

## 10. Communication (Solver ↔ GUI)
- ZeroMQ PUB topics: `state.eta`, `state.velocity`, `diag.global`, `diag.field.<name>`, `event.progress`.
- Payload: msgpack or JSON (initial simple) with: {"t": float, "shape": [..], "dtype": str, "data": base64 or omitted if using shared memory later}.
- Rate adaptation: solver publishes every N steps (configurable), GUI can request on-demand refresh via REQ/REP control channel (future extension).

## 11. Testing Taxonomy
| Layer | Focus | Example |
|-------|-------|---------|
| Unit | Pure function correctness | `stencils.d_dx_t` on polynomial field |
| Integration | Coupled tendency + step | One-step leapfrog vs MATLAB snapshot |
| Regression | Short (e.g., 50-step) golden run | Compare energy time series |
| Conservation | Volume / energy / potential enstrophy | Closed basin no forcing |
| GUI Smoke | Widget creation / subscription | Start GUI, receive one frame |
| Performance | Wall-clock / step, memory | Baseline gravity wave case |

Golden-run artifacts stored under `tests/fixtures/golden/` with versioned JSON manifest (fields + tolerances).

## 12. Code Quality & Tooling
- Formatting: black (line length 88), isort (sections: stdlib, thirdparty, local), trailing whitespace trimmed.
- Typing: mypy strict (no `Any` in public APIs); allow localized `typing.cast` with justification comment.
- Lint: pylint target score ≥ 8.5 (ratcheted upward after solver parity).
- Pre-commit hooks: black, isort, mypy (fast mode), pytest -k fast subset.
- Continuous Integration: GitHub Actions matrix (Python 3.10–3.12) running (format check, mypy, lint, tests, coverage threshold ≥ 85%).

## 13. Logging & Error Handling
- Use `logging` with hierarchical loggers: `shel.model`, `shel.gui`, `shel.io`.
- Standard structured fields: run_id, step, phase, elapsed_s.
- Exceptions: domain-specific (e.g., `ConfigurationError`, `NumericsError`, `ConservationError`).
- Fail fast on NaNs (validation hook each step if debug mode).

## 14. Performance & Acceleration
- Baseline vectorized NumPy first; profile (cProfile + line_profiler) after correctness.
- Candidate hot spots: advection, pressure gradient, diffusion.
- Numba path optional: mirror function signature & gated by `ENABLE_NUMBA` flag; parity test ensures identical numerics.
- GPU (CuPy) flagged for post-parity exploration; same test contract.

## 15. Future Enhancements
- Adaptive timestep (CFL monitor) with rollback if threshold exceeded.
- Spectral energy diagnostics (FFT) & scale-wise energy budget.
- Multi-tracer support with passive/active coupling (future chemistry hooks).
- Checkpoint/restart facility (versioned state snapshots).

## 16. Active Near-Term Sprint Focus
1. Radiation BCs parity: finalize Sommerfeld timing (eta vs continuity ordering) and tune Flather/Sommerfeld coefficients against MATLAB snapshots; add comparison tests.
2. Dynamic regression v2: time‑series baseline (E, V, eta_rms, u_rms) over N steps with per‑metric tolerances; versioned JSON manifest and generator.
3. Tracer step: implement minimal advection (centered) with mask handling and integrate tracer BCs; add conservation tests.
4. Sponge enhancements: variable width per side, diagonal/2D tapers, and tests on non‑uniform H and active wave cases.
5. Documentation: expand README with Python quickstart and time API (Done), add developer guide for BC strategy patterns and examples (next).

## 17. Change Log (Recent)
- 2025-09-05: Consolidated prompts; added GUI phases; added divergence/shear/stretch diagnostics.
- 2025-09-05: Composite initial state builder implemented & tested.
- 2025-09-05: Phase 3/4 closure tasks completed: initialization conservation tests and static golden fidelity harness; all closure tests green.
- 2025-09-05: Phase 5 start — stencils module implemented with unit tests (linear field derivatives, averaging round-trip, zero divergence for solid-body case).
- 2025-09-05: Minimal pressure gradient (T→U/V) and continuity (flux-form explicit Euler) implemented with tests, including closed-box volume conservation.
- 2025-09-05: Added linear bottom drag and Laplacian diffusion tendencies with unit tests; Phase 5 core operators green.
- 2025-09-05: Added minimal explicit one-step harness combining pressure, drag, diffusion; few-step closed-box volume conservation test passes.
- 2025-09-05: Added multi-step dynamic golden baseline regression harness: generator `devops/scripts/generate_dynamic_golden.py`, baseline JSON `tests/python/golden/dynamic_baseline.json`, and regression test `tests/python/test_golden_dynamic.py`.
- 2025-09-05: Solver-side boundary conditions: closed (no-normal-flow) and free‑slip implemented and wired into `explicit_step`; per‑side BC resolution via config-aware stepper; radiative (Sommerfeld) prototype for momentum and eta with a minimal east‑outlet smoke test. Full suite: 63 passed, 4 skipped (GUI).
 - 2025-09-05: Introduced Boundary Condition strategy layer under `shel/model/boundary_conditions/` with registry and strategies (Closed, Free‑slip, Radiative/Sommerfeld). Bridged existing `solvers/common/boundaries.py` to use strategies. Added registry/dispatch tests. Full suite: 66 passed, 4 skipped (GUI).
 - 2025-09-05: Organized BC strategies into domain subpackages (`boundary_conditions/common`, `momentum`, `waterlevel`), migrated `model_runner` and `ministep/stepper` to strategy layer, removed legacy `boundary_conditions/boundary.py` and old solver shim `solvers/common/boundaries.py`. Full suite: 66 passed, 4 skipped (GUI).
   - 2025-09-05: Added `boundary_conditions/README.md` documenting strategy API, layout, and usage.
 - 2025-09-06: Flather OBCs completed: momentum and eta strategies with relaxation; per-side application in `stepper` with eta timing (pre/post, always post-enforced). Examples added: `examples/python/flather_config_example.py` and `examples/python/flather_sponge_config_example.py` (cosine-tapered sponge, conservative params). Optional sponge layer implemented with width/alpha/taper/apply_to and tests. Full suite: 72 passed, 4 skipped (GUI).
 - 2025-09-06: Sommerfeld (radiative) eta tests extended with pulse propagation centroid check; tracer BC base/registry and closed/radiative strategies added with unit tests. Full suite: 77 passed, 4 skipped (GUI).
 - 2025-09-06: Leapfrog + Asselin centralized under `solvers/time`; added `leapfrog_step_with_config` reusing stepper helpers (BC per-side, eta relax/timing, sponge). Refactored `stepper.py` to extract reusable helpers. Full suite: 79 passed, 4 skipped.
 - 2025-09-06: Phase 7 closure — Added Sommerfeld momentum-side exact formula tests; introduced `DirichletEtaBC`; integrated tracer BCs with a minimal stepper; updated BC README; added `leapfrog_flather_dirichlet_sponge_example.py`. Full suite: 86 passed, 4 skipped.
 - 2025-09-18: Added minimal NetCDF snapshot writer and NetCDF timeseries appender under `model/outputs/writers/netcdf.py`, and Parquet timeseries writer/appender under `model/outputs/writers/parquet.py`; unit tests added (`test_outputs_writers.py`) covering write/append/range and basic round-trip via readers.

 - 2025-09-19: Phase 6 complete — Added 2nd order upwind advection (`advection_upwind.py`), quadratic drag (`quadratic_drag.py`), and biharmonic diffusion (`biharmonic_diffusion.py`) operators; all unit tests green. Status table updated.

 - 2025-09-19: Phase 8 complete — Added enstrophy and potential enstrophy diagnostics (`diagnostics/enstrophy.py`, `diagnostics/potential_enstrophy.py`); all unit tests green. Status table updated.

 - 2025-09-19: Phase 10 complete — Output writers (NetCDF, Parquet, JSON) validated for non-uniform H; all tests green. Status table updated.

Maintainers: Update status table & change log in any PR modifying numerics, diagnostics, or architecture.
