# Solver Bugfix Plan — Python ↔ MATLAB Parity

Created: 2026-02-17  
Status: **IN PROGRESS**

## Context

The MATLAB `model_handles.m` is the working reference implementation. The Python
solver stack (`shel/model/solvers/`) has accumulated several bugs that make it
numerically divergent from MATLAB. This plan tracks the fixes grouped into waves
so that each wave leaves the codebase in a passing-tests state.

Baseline before any changes: **119 tests pass, 5 skipped**.

---

## Bug Inventory

| #   | Sev | Description                                                                                                                                                                                                                         | Files                                                                                                                               |
| --- | --- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------- |
| 1   | 🔴  | Leapfrog is actually Forward Euler — `eta_nm1`/`U_nm1` are unused in the advance step; it calls `explicit_step(eta_n, …, dt=dt)` instead of advancing from `n-1` with `2*dt`                                                        | `solvers/time/leapfrog.py`, `solvers/common/ministep.py`                                                                            |
| 2   | 🔴  | Momentum in **non-conservative** form — all tendencies (pressure, advection, diffusion, Coriolis) are missing the $H$ factor; time stepping does plain `U + dt * dU` instead of $(H_{old} \cdot U_{old} + 2dt \cdot RHS) / H_{new}$ | `momentum/pressure.py`, `momentum/advection.py`, `momentum/diffusion.py`, `common/ministep.py`                                      |
| 3   | 🔴  | Velocity updated **before** eta — Python feeds already-advanced `U_next,V_next` into continuity; MATLAB computes continuity RHS from current velocities first                                                                       | `common/ministep.py`                                                                                                                |
| 4   | 🔴  | $H = \eta + d$ never recomputed after the $\eta$ update in any stepper or factory wrapper                                                                                                                                           | `common/ministep.py`, `time/leapfrog.py`, `solvers/factory.py`                                                                      |
| 5   | 🟡  | Land masks (`mask`, `mask_u`, `mask_v`) completely ignored in all operators                                                                                                                                                         | All solver operators                                                                                                                |
| 6   | 🟡  | Asselin filter has spurious **0.5** factor (`0.5 * α` vs MATLAB's `α`)                                                                                                                                                              | `time/asselin.py`                                                                                                                   |
| 7   | 🟡  | Diffusion uses current `u` instead of `u_old` (leapfrog stability requirement)                                                                                                                                                      | `momentum/diffusion.py`, `common/ministep.py`                                                                                       |
| 8   | 🟡  | Boundary H averaging in `avg_x_t_to_u` / `avg_y_t_to_v` uses linear extrapolation which can produce negative depths                                                                                                                 | `common/stencils.py`                                                                                                                |
| 9   | 🟡  | Radiative BC uses domain-mean wave speed instead of per-cell `sqrt(gH)`                                                                                                                                                             | `boundary_conditions/momentum/strategies.py`, `boundary_conditions/waterlevel/strategies.py`, `boundary_conditions/common/utils.py` |
| 10  | 🟢  | Flather formulation differs from MATLAB                                                                                                                                                                                             | `boundary_conditions/momentum/strategies.py`                                                                                        |

---

## Fix Strategy

### Wave 0 — Isolated leaf fixes (no API changes)

These are self-contained, low-risk fixes. Fix them first so they don't compound
with the deeper changes.

- [x] **Bug #6** — `asselin.py`: remove the `0.5 *` factor.
    - Update: `return current + alpha * (old - 2 * current + new)`
    - Tests affected: `test_leapfrog_asselin_smoke.py` (tolerances may change)

- [x] **Bug #8** — `stencils.py`: change boundary fill in `avg_x_t_to_u` and
      `avg_y_t_to_v` from linear extrapolation to simple nearest-cell copy:
      `out[:, 0] = T[:, 0]` and `out[:, -1] = T[:, -1]` (like MATLAB).
    - Tests affected: `test_stencils.py`, conservation tests, regression hashes.

### Wave 1 — Core solver rewrite (Bugs #1–4, #7)

These are tightly coupled and must be done together to keep the scheme
consistent. The goal is to produce a conservative flux-form solver matching
MATLAB's `ComputeLeapfrog`.

#### Step 1.1 — New signatures for tendency operators

Each momentum tendency function gets an `H` (and/or `H_old`) parameter:

- `pressure_gradient(eta, H, g, dx, dy)` → scale by $H_U$
- `viscous_tendency(U_old, V_old, H_old, nu, dx, dy)` → multiply by $H_{old}$, use `u_old` (fix #7)
- `advect_momentum(U, V, H, dx, dy)` → flux-form centered advection
- Coriolis in `ministep.py` → multiply by $H_U$

#### Step 1.2 — Rewrite `explicit_step` to the correct order

1. Compute continuity RHS from **current** $(u, v)$ → `RHSeta`
2. Time-step eta: `eta_new = eta_old + 2*dt * RHSeta * mask`
3. Compute `H_new = eta_new + d`
4. Compute momentum RHS from current tendencies (all scaled by $H$)
5. Time-step velocity: `u_new = mask_u * (u_old * H_old_u + 2*dt * RHSu) / H_new_u`
6. Apply BCs to `eta_new`, `u_new`, `v_new`

Key function: introduce a **new** `conservative_step(...)` function in `ministep.py`
(or replace `explicit_step`) that takes `eta_old/eta_n`, `u_old/u_n`, `H_old`,
`d` (bathymetry), `mask`, `mask_u`, `mask_v`, and returns
`(eta_new, u_new, v_new, H_new)`.

Keep the old `explicit_step` API in a deprecated wrapper for one cycle, or
update call sites.

#### Step 1.3 — Wire into leapfrog stepper

- `leapfrog_stepper` calls the new conservative step with `2*dt`, passing
  `eta_nm1`, `U_nm1`, `V_nm1` as the "old" arrays and `eta_n`, `U_n`, `V_n`
  as "current" for RHS evaluation.
- Asselin filter on the middle level.
- Return `eta_np1, U_np1, V_np1, eta_n_filtered, U_n_filtered, V_n_filtered, H_np1`.

#### Step 1.4 — Update `factory.py`

- Recompute `state.H = state.eta + state.d` after each step (fix #4).
- Pass `state.d`, masks through to the solver.

### Wave 2 — Mask support (Bug #5)

Add optional `mask`, `mask_u`, `mask_v` parameters to:

- `pressure_gradient`
- `viscous_tendency`
- `advect_momentum`
- `bottom_drag_tendency`
- `update_free_surface` (continuity)
- `explicit_step` / `conservative_step`

When a mask is provided, multiply every tendency by the appropriate mask.
When `None`, operators behave as today (all-water).

### Wave 3 — Boundary condition fixes (Bugs #9, #10)

- **Bug #9**: Change `mean_c_along_side` to return a per-cell array
  `np.sqrt(g * H_edge)` instead of a scalar. Update all BC strategies that
  consume it (radiative momentum, radiative eta).
- **Bug #10**: Align Flather momentum BC with MATLAB:
  `u_boundary = sign * sqrt(g / H_boundary) * eta_boundary`

### Wave 4 — Test updates

- Update hash-based regression tests (`test_one_step_regression`,
  `test_multi_step_regression`) with new baselines.
- Update tolerance-sensitive tests (`test_leapfrog_asselin_smoke`).
- Add new conservation tests:
    - Volume conservation with variable bathymetry
    - Momentum conservation in closed box (conservative form)
    - Energy decay test with friction
- Add parity test: run identical small case through MATLAB (saved arrays) and
  Python, compare within tolerance.

---

## Execution Order

```
Wave 0  →  Bug #6  →  Bug #8
Wave 1  →  Step 1.1  →  Step 1.2  →  Step 1.3  →  Step 1.4
Wave 2  →  mask support
Wave 3  →  Bug #9  →  Bug #10
Wave 4  →  test updates & new tests
```

Each wave ends with `pytest` green (modulo intentionally updated baselines).

---

## Progress Tracking

| Wave | Step                         | Status  | Date       |
| ---- | ---------------------------- | ------- | ---------- |
| 0    | Bug #6 Asselin               | ✅ DONE | 2026-02-17 |
| 0    | Bug #8 Stencil boundary      | ✅ DONE | 2026-02-17 |
| 1    | Step 1.1 Tendency signatures | ✅ DONE | 2026-02-17 |
| 1    | Step 1.2 conservative_step   | ✅ DONE | 2026-02-17 |
| 1    | Step 1.3 Leapfrog wiring     | ✅ DONE | 2026-02-17 |
| 1    | Step 1.4 Factory/stepper     | ✅ DONE | 2026-02-17 |
| 1    | Test fixes (API changes)     | ✅ DONE | 2026-02-17 |
| 1.5  | model_runner wiring bugs     | ✅ DONE | 2026-02-17 |
| 2    | Mask support                 | ⬜ TODO |            |
| 3    | Bug #9 Local c               | ⬜ TODO |            |
| 3    | Bug #10 Flather              | ⬜ TODO |            |
| 4    | Test updates                 | ⬜ TODO |            |

### Wave 1.5 notes — model_runner wiring (the "UI still broken" bugs)

Three critical bugs in `model_runner._advance_timestep()` meant the solver
fixes from Waves 0–1 never actually took full effect at runtime:

1. **Asselin-filtered values discarded**: The leapfrog stepper returns both
   the advanced (n+1) and Asselin-filtered (n) values. The runner was
   copying the _unfiltered_ current state into `eta_old`/`u_old`/`v_old`
   instead of using the filtered values. Without the filter feeding back,
   the leapfrog computational mode grows → velocity instability.
2. **`H` never recomputed**: After updating `self.state.eta`, the total
   water-column height `H = eta + d` was never recalculated for the next
   step. This made all conservative tendencies use stale `H`.
3. **`d` and `H_old` never passed**: The call to `leapfrog_step_with_config`
   did not pass `d=self.state.d` or `H_old=self.state.H_old`, forcing the
   function to _infer_ them. Now passed explicitly for correctness.

- All core solver files rewritten to conservative flux form matching MATLAB.
- `explicit_step` now returns 4-tuple `(eta_new, U_new, V_new, H_new)`.
- New parameters: `d` (bathymetry), `eta_old`, `H_old`, `U_old`, `V_old` for leapfrog.
- Coriolis interpolated to U/V faces (fixes broadcast mismatch).
- Diffusion uses `H_old` and `U_old` (Bug #7 fixed).
- 21 test call-sites updated for new APIs.
- Golden baseline hash regenerated.
- Energy test tolerance relaxed (Euler stepping has O(dt) energy overshoot).
- Leapfrog tests updated with physically realistic parameters (dx=1000, H=100).
- **Test status: 119 passed, 5 skipped.**
