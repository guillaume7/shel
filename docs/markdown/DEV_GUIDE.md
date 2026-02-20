# SHEL Python Development Guide

This guide provides best practices, architectural notes, and implementation details for contributing to the SHEL Python port.

## Project Structure
- `src/python/shel/`: Main model codebase
- `tests/python/`: Unit, integration, and regression tests
- `examples/python/`: Example configurations and usage
- `docs/markdown/`: Documentation and guides

## Key Concepts
- **Arakawa C-grid**: Staggered grid for stability
- **Leapfrog Scheme**: Time-stepping for second-order accuracy
- **Strategy/Factory Pattern**: Interchangeable operators and boundary conditions
- **Conservation Diagnostics**: Volume, energy, enstrophy, potential enstrophy

## Development Workflow
1. Fork and clone the repository
2. Create a feature branch
3. Write code and tests (unit, integration, regression)
4. Run `pytest` and ensure all tests pass
5. Update documentation and status table
6. Submit a pull request

## Coding Standards
- Use type hints and docstrings for all public functions
- Follow PEP8 and use `black` for formatting
- Use lazy logging formatting: `logger.info("Message: %s", value)`
- Write pure functions for all core operators
- Isolate I/O and orchestration logic

## Testing
- Unit tests for pure functions and operators
- Integration tests for coupled solver steps
- Regression tests for numerical fidelity (golden runs)
- Conservation tests for volume, energy, potential enstrophy

## Documentation
- Update `README.md` for new features
- Add API docs and usage examples in `docs/markdown/`
- Document new strategies, operators, and diagnostics

## Solver Implementation
- Operators: pressure, advection, diffusion, friction, continuity
- Time integration: leapfrog, Asselin filter (in `solvers/time/`)
- Boundary conditions: strategy/registry pattern
- Diagnostics: field and integrated metrics

## IO Implementation
- NetCDF: grid, state, diagnostics snapshots
- Parquet: time-series diagnostics
- JSON: deterministic serialization for regression
- OutputManager: schedules output hooks
- YAML: configuration for reproducible runs

## GUI Implementation (Future)
- PyQt main window scaffold
- ZeroMQ pub/sub protocol for solver-GUI communication
- Real-time visualization adapters

## Contribution Checklist
- [ ] Code passes all tests
- [ ] Documentation updated
- [ ] Status table updated
- [ ] Follows coding standards
- [ ] Pull request submitted

---
For more details, see the full documentation in `docs/markdown/` and the status table in `.github/prompts/MASTER_PORT_PROMPT.md`.
