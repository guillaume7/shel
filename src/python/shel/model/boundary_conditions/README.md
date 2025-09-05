# Boundary Conditions (strategy layer)

Purpose: domain‑oriented, functional BC strategies used by the solver/ministep path. Legacy OO BCs and factory have been removed.

## Layout
- base.py: abstract strategy interfaces
  - MomentumBC.apply_uniform(U,V)
  - MomentumBC.apply_side(U,V, side, ...)
  - EtaBC.apply_side_eta(eta_next, side, ...)
- registry.py: register/get/list strategies by name
- strategies.py: registrar wiring domain strategies under common names
- common/utils.py: shared helpers (e.g., mean_c_along_side)
- momentum/strategies.py: ClosedBC, FreeSlipBC, RadiativeSommerfeldBC, FlatherBC
- waterlevel/strategies.py: RadiativeSommerfeldEtaBC
- tracer/: placeholder for future tracer BCs

## Supported names
- closed: momentum only (no normal flow)
- freeslip: momentum only (no normal flow + zero tangential gradient)
- radiative: momentum + eta (Sommerfeld prototype)
- flather: momentum only (requires external free-surface eta_ext along the open boundary)

## Usage (solver side)
- Uniform BC in a step: get_bc(name)[0]().apply_uniform(U, V)
- Per‑side momentum: m = get_bc(name)[0](); m.apply_side(U,V, side, U_old=U0, V_old=V0, H=H, g=g, dt=dt, dx=dx, dy=dy)
- Per‑side eta: e = get_bc(name)[1](); e.apply_side_eta(eta1, side, eta_old=eta0, H=H, g=g, dt=dt, dx=dx, dy=dy)

For Flather, pass eta_ext and eta_old to apply_side:

- m = get_bc("flather")[0](); m.apply_side(U,V, side, U_old=U0, V_old=V0, H=H, g=g, dt=dt, dx=dx, dy=dy, eta_old=eta0, eta_ext=eta_boundary)

Notes
- Strategies are stateless; instantiate and reuse as needed.
- Radiative timing for eta vs continuity is approximate in the ministep path and may be refined.
- Add new BCs by implementing the appropriate interface in the domain folder, then register in strategies.py via register_bc.
