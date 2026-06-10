# `engines.fluid` — 2-D rotating shallow-water solver — CONTRACT

> **Status: LIVING CONTRACT — versioned; the test suite *is* the contract.** Built
> in Planet Phase 3, extended for the GCM climb (the advected tracer, rung 1 — see
> below + the Changelog). Guarded by its passing validation suite
> (`engines/fluid/tests/`, run via `./run_tests.ps1`); this one page is the unit of
> context downstream code loads — Planet's circulation coupler and the documented
> GCM climb depend on *this*, never on `planet/` internals. **Extend it directly
> when a consumer needs it:** keep the suite green, add tests covering the new
> surface, and record the change in the Changelog. Editing a shared engine triggers
> the **full-repo gate** + the import-drift guard (ADR 0003) — that is the guardrail.
> Engines are no longer frozen artifacts (ADR 0005, which supersedes the freeze
> clauses of ADR 0001 / 0003).
>
> **The program's second shared engine** — and deliberately a *different solver
> class* from `engines.diffusion`: hyperbolic & **explicit** (CFL-limited), where the
> diffusion spine is parabolic & implicit (unconditionally stable). They share no
> machinery; that contrast is the point (ARCHITECTURE.md §5).

## What it solves

The **rotating-frame shallow-water equations** on a doubly-periodic β-plane, in
**vector-invariant** form:

```
∂u/∂t = (f + ζ) v − ∂_x B
∂v/∂t = −(f + ζ) u − ∂_y B
∂h/∂t = −∇·(h u)
       ζ = ∂_x v − ∂_y u            (relative vorticity)
       B = g(h + h_b) + ½(u² + v²)  (Bernoulli function)
       f = f₀ + β·(y − y_ref)       (β-plane Coriolis)
       q = (f + ζ)/h                (potential vorticity — the rotation-aware invariant)
```

`h` is layer thickness, `(u, v)` the depth-independent horizontal velocity, `g`
gravity (or a *reduced* gravity `g′` for a baroclinic mode), `h_b` optional bottom
topography. The engine is **physics-agnostic** in the diffusion-engine sense: the
consumer supplies `g, H, f₀, β` (Planet's `circulation.py` pins the planetary
values). Derived scales: gravity-wave speed `c = √(gH)`, deformation radius
`L_R = √(gH)/|f₀|`.

## Discretization (fixed)

- **Arakawa C-grid (staggered).** `h`, `h_b` at cell *centers*; `u` on E–W *faces*;
  `v` on N–S *faces*; `ζ`, `f`, `q` at cell *corners*. The staggering is what lets a
  centered scheme carry geostrophic balance and `L_R` without a checkerboard mode — and
  is where a wrong Coriolis averaging hides (hence the PV seal is the discriminating leg).
- **Vector-invariant momentum, symmetric corner-PV Coriolis flux** (Sadourny-1975
  lineage). This realization conserves **energy** semi-discretely (not potential
  enstrophy — a single such scheme conserves one or the other; Arakawa–Lamb conserves
  both, at a large step up in complexity, *not* built here).
- **SSP-RK3 (Shu–Osher) explicit time stepping.** Stable for the non-dissipative
  centered operator up to a CFL on the fastest signal `√(gH) + |u|`. Explicit Coriolis
  is fine (`f·dt ≪ 1` at the gravity-wave CFL — no semi-implicit treatment needed).

## API

```python
from engines.fluid import ShallowWater, SWState, Grid2D, uniform_grid

grid = uniform_grid(Lx, Ly, nx, ny)             # uniform, doubly-periodic
sw   = ShallowWater(grid, g, mean_depth,         # g, H  → c=√(gH), L_R
                    f0=0.0, beta=0.0,            # β-plane: f = f0 + β(y − y_ref)
                    y_ref=None, bottom=None)     # y_ref defaults to domain centre

state = SWState(h, u, v, tracer=None)            # stacked plain (ny,nx) arrays
state = sw.step(state, dt)                        # one SSP-RK3 step; returns new state
state = sw.solve(state, t_end, dt=None)           # march to t_end (dt=None → CFL step)
dt    = sw.max_dt(state, safety=0.3)              # recommended CFL step

# diagnostics (all consume a state, return scalars / corner|center arrays)
sw.mass(state)                  # ∫h dA            — machine-precision conserved
sw.energy(state)                # ∫[½h(u²+v²)+½gh²]dA — bounded, dt³-convergent drift
sw.potential_vorticity(state)   # q=(f+ζ)/h at corners
sw.potential_enstrophy(state)   # ∫½q²h dA         — bounded (spatial-limited) drift
sw.relative_vorticity(state)    # ζ at corners     (∫ζ dA = 0 to machine precision)
sw.gravity_wave_speed           # √(gH);   sw.rossby_radius  √(gH)/|f0|
```

- **`SWState`** is the stable (ADR-0001) data boundary: three plain `(ny, nx)` `ndarray` fields
  `(h, u, v)` — *stacked fields and nothing more* — plus an **optional `tracer`** slot.
- **`step`** does not mutate its input; it **raises** on non-positive `dt`, on `dt`
  above the CFL stability limit (the explicit analogue of the diffusion engine's
  stability promise — here *conditional*, so enforced), and on a set `tracer` (below).

### The stable data boundary (ADR 0001)

`state` is an `SWState` of plain 2-D arrays — and only it crosses the per-step boundary:
`step`/`solve` consume and return it, the diagnostics consume it. No live object crosses.
`Grid2D`, `g`, `H`, `f₀`, `β`, and `bottom` are **construction-time configuration** that
reduce to numbers/arrays during RHS assembly; a compiled reimplementation parameterizes
them natively and exposes the same `SWState`. The viz layer (ADR 0002) consumes the same
arrays — never a live solver object.

### The GCM-climb seam (extension-ready, ARCHITECTURE.md §8 / Planet §5)

The stacked-field shape is the growth axis: **N vertical layers** is a leading-axis
extension (rung 3 — baroclinic), and an **advected tracer** (moisture/temperature for
rungs 1–2) is one more field. The `tracer` slot is **declared on `SWState` but not
advected in v1** — `step` raises `NotImplementedError` naming rung 1, the same "build the
seam, not the machinery" idiom as the map's unpainted `vector_overlay`. Adding either is a
contract *extension* that does not change the per-step array boundary.

## Guaranteed invariants (what the test suite enforces — = the contract)

1. **Gravity-wave speed `√(gH)`** and **Poincaré dispersion `ω² = f₀² + gH·k²`**
   reproduced to ~1e-3 (`test_waves.py`). The Poincaré leg is the *rotation* check — the
   frequency rises with `f₀` exactly as predicted; a wrong Coriolis fails it loudly.
2. **Rossby-wave dispersion `ω = −βk/(k²+l²+1/L_R²)`** — **westward**, **dispersive**
   (longer waves faster), reproduced to a **loose band that converges to analytic as the
   grid refines** (`test_waves.py`, slow). The slow balanced mode carries a few-percent
   numerical-dispersion error a gravity wave does not — asserted loose, named.
3. **Geostrophic balance is steady** — a balanced zonal jet sits still (`test_geostrophy.py`),
   and **geostrophic adjustment** of a height bump settles to the analytic Helmholtz state
   `(1 − L_R²∇²)η_adj = η_init` over `L_R`, to ~few % (`test_geostrophy.py`, slow) — Rossby's
   benchmark, with `L_R` (hence `f₀`) the discriminating scale.
4. **Mass `∫h` conserved to machine precision** for any state/step (`test_conservation.py`);
   `∫ζ dA = 0` to machine precision (a discrete curl).
5. **Energy (KE+PE) conserved to a bounded, dt→0 convergent drift** — the RK3
   time-truncation of the semi-discrete energy invariant, shrinking as `dt³`
   (`test_conservation.py`, the convergence test is slow).
6. **Potential vorticity / enstrophy bounded at FINITE amplitude** (`test_conservation.py`)
   — the discriminating Coriolis seal, run on a balanced vortex at Rossby number ~0.5 (where
   advection genuinely moves PV around): enstrophy holds to a small spatial-discretization
   bound and PV grows no spurious extrema (PV is materially conserved). At small amplitude
   this leg would be near-vacuous; the finite-amplitude requirement is what gives it teeth.

## Validation boundary (what Phase 3 does *not* claim)

- **Energy *or* enstrophy, not both** (semi-discretely). This realization conserves
  energy; potential enstrophy drifts at a small, *spatial*-discretization-limited bound
  (it shrinks with Δx, not dt). Machine-precision conservation is claimed only for **mass**.
- **No exact material PV at the grid scale.** Material PV conservation is asserted as
  bounded (no spurious extrema) over a few inertial periods at moderate Rossby number — not
  to machine precision. Long, fully-turbulent, high-Rossby integrations can cascade enstrophy
  to the grid (the energy-conserving scheme's known behaviour); add hyperviscosity or an
  enstrophy-conserving operator for that regime (named, not built).
- The *physical constants* (`g, H, f₀, β`, planetary values) are the **consumer's** to supply
  and validate (Planet `circulation.py`, against geostrophic balance + the jet benchmark).

## Units & scope

- **SI throughout** (`h, H` m; `u, v` m/s; `g` m/s²; `f₀` 1/s; `β` 1/(m·s); time s).
- **v1 scope (named rungs on the §5 staircase, not built):** single layer (dry — no
  thermodynamic variable, *the* fact making Phase-4 coupling one-way); **doubly-periodic
  β-plane** (rigid channel walls in y are a named BC extension); **explicit** (CFL-limited);
  no vertical structure (multi-layer = rung 3), no moisture/tracer advection (rung 1), no
  sphere (rung 5 — the pole problem). The `SWState` array boundary is the seam where each is
  slotted without touching consumers.

## Changelog

- **2026-06-09 — Planet Phase 3.** Initial build & validation (the guaranteed invariants above).
- **2026-06-10 — doctrine.** Status FROZEN → living, versioned contract (ADR 0005); no API change.
