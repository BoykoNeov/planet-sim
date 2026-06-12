# `engines.fluid` — 2-D rotating shallow-water solver — CONTRACT

> **Status: LIVING CONTRACT — versioned; the test suite *is* the contract.** Built
> in Planet Phase 3, extended for the GCM climb (the advected tracer, rung 1; **N vertical
> layers → baroclinic instability, rung 3** — see below + the Changelog). Guarded by its passing
> validation suite
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

state = SWState(h, u, v, tracer=None)            # stacked plain (ny,nx) arrays; tracer optional
state = sw.step(state, dt)                        # one SSP-RK3 step; advects the tracer if set
state = sw.solve(state, t_end, dt=None)           # march to t_end (dt=None → CFL step)
dt    = sw.max_dt(state, safety=0.3)              # recommended CFL step

# diagnostics (all consume a state, return scalars / corner|center arrays)
sw.mass(state)                  # ∫h dA            — machine-precision conserved
sw.energy(state)                # ∫[½h(u²+v²)+½gh²]dA — bounded, dt³-convergent drift
sw.potential_vorticity(state)   # q=(f+ζ)/h at corners
sw.potential_enstrophy(state)   # ∫½q²h dA         — bounded (spatial-limited) drift
sw.relative_vorticity(state)    # ζ at corners     (∫ζ dA = 0 to machine precision)
sw.tracer_mass(state)           # ∫hθ dA           — machine-precision conserved (passive tracer)
sw.tracer_variance(state)       # ∫½hθ² dA         — bounded drift (NOT monotone — no limiter)
sw.gravity_wave_speed           # √(gH);   sw.rossby_radius  √(gH)/|f0|
```

- **`SWState`** is the stable (ADR-0001) data boundary: three plain `(ny, nx)` `ndarray` fields
  `(h, u, v)` — *stacked fields and nothing more* — plus an **optional `tracer`** slot.
- **`step`** does not mutate its input; it **raises** on non-positive `dt` and on `dt`
  above the CFL stability limit (the explicit analogue of the diffusion engine's
  stability promise — here *conditional*, so enforced). A set **`tracer`** is advected as a
  **passive** scalar (below); `tracer=None` runs the dry dynamics bit-for-bit.

### N-layer baroclinic API (rung 3 — the leading-axis seam, built)

A **sibling** engine for the vertical structure that single-layer SW categorically cannot carry
(no APE, no shear ⇒ no baroclinic instability). The single-layer `ShallowWater`/`SWState` are
**untouched** — existing consumers are unaffected.

```python
from engines.fluid import (LayeredShallowWater, LayeredState,
                           ThermalWindBackground, TwoLayerStability)

lay = LayeredShallowWater(grid, g, mean_depth, reduced_gravities,   # mean_depth: per-layer H_k (len nl)
                          f0=0.0, beta=0.0, y_ref=None, bottom=None, #  reduced_gravities: nl-1 internal g'_k
                          background=None)                           #  default-off basic state
bg    = lay.thermal_wind([U0, U1, ...])         # per-layer mean zonal flow → ThermalWindBackground(U, G)
state = LayeredState(h, u, v)                    # stacked (nl, ny, nx) arrays, layer 0 = top
state = lay.step(state, dt)                       # one SSP-RK3 step; layers couple ONLY via Montgomery pressure
dt    = lay.max_dt(state, safety=0.3)             # CFL on the fast EXTERNAL wave √(g·H_tot)

lay.layer_mass(state)            # ∫h_k dA per layer — machine-precision conserved (background=None)
lay.perturbation_energy(state)   # Σ_k ∫(u²+v²+h'²) dA — the linear-growth diagnostic (grows as e^{2σt})
lay.stability()                  # the matching TwoLayerStability analytic operator (nl==2)
TwoLayerStability(f0,g,gp,H1,H2).growth_rate(k, l, Us)   # σ = max Im ω from the 6×6 dispersion matrix
```

- **Leading layer axis.** `LayeredState` stacks the same C-grid fields with a leading axis
  `(nl, ny, nx)` — the GCM-climb seam, *not* a new per-step boundary. `h` is the **total**
  thickness `H_k + h'`; with a `background`, `(u, v)` are the **perturbation** velocities (the
  mean `U_k` rides in the engine's Doppler coefficient, not the field).
- **Coupling is pressure-only.** Per-layer stencils are *identical* to the single-layer engine
  (vectorized over the layer axis); the only inter-layer term is the **Montgomery pressure stack**
  `M₀ = g·η_top`, `M_k = M_{k-1} + g'_k·η_{k-1}`.
- **The background is default-off** and injects an unstable thermal-wind basic state as **constant
  coefficients** (a doubly-periodic domain cannot carry a meridional gradient *as a field*): a
  Doppler `−U_k ∂_x` and a baroclinic `−G_k v'` in continuity, with the prognostic fields the
  perturbations. `None` → the plain nonlinear dynamics (and the bit-for-bit single-layer reduction).

### The stable data boundary (ADR 0001)

`state` is an `SWState` of plain 2-D arrays — and only it crosses the per-step boundary:
`step`/`solve` consume and return it, the diagnostics consume it. No live object crosses.
`Grid2D`, `g`, `H`, `f₀`, `β`, and `bottom` are **construction-time configuration** that
reduce to numbers/arrays during RHS assembly; a compiled reimplementation parameterizes
them natively and exposes the same `SWState`. The viz layer (ADR 0002) consumes the same
arrays — never a live solver object.

### The GCM-climb seam (extension-ready, ARCHITECTURE.md §8 / Planet §5)

The stacked-field shape is the growth axis: **N vertical layers** is a leading-axis
extension (rung 3 — baroclinic), and an **advected tracer** (a temperature/moisture proxy
for rungs 1–2) is one more field. The **tracer is built (rung 1):** when `state.tracer` is
set, `step` advects the tracer mass `h·θ` in flux form (conserving `∫hθ` to machine
precision), strictly **passively** — no back-reaction on `(h, u, v)` (an *active* buoyancy
tracer feeding `h` would be a different reduced-gravity / two-layer model). **N vertical layers
are built (rung 3, Phase A):** the sibling `LayeredShallowWater` stacks the fields on a leading
layer axis, coupled *only* through the Montgomery pressure, with the **interface displacement now
the active dynamical buoyancy** (the leap past the *passive* tracer — exactly the "active buoyancy
tracer feeding `h` = a two-layer model" this contract named). The single-layer engine is left
untouched (it is a sibling), so the `nl=1` reduction is byte-identical and existing consumers are
unaffected. Adding either is a contract *extension* that does not change the per-step array
boundary (ADR 0005 — engines are living contracts).

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
7. **Passive tracer: `∫hθ` machine-precision; `∫½hθ²` bounded** (`test_tracer.py`).
   Flux-form advection telescopes like mass (the clean anchor); the tracer is strictly passive,
   so the dry `(h, u, v)` trajectory is **byte-identical** to a `tracer=None` run (the re-seal),
   and a uniform tracer stays uniform (consistency — no spurious source). Tracer **variance** is
   a near-invariant held to a small, **bounded** drift (the enstrophy honesty class — round-off for
   smooth fields, cascading toward the grid only under strong filamentation) — **not** monotone
   (no flux limiter → over/undershoot on sharp gradients).
8. **N-layer baroclinic (rung 3, Phase A): the linear growth rate + the byte-identical reduction**
   (`test_layered.py`, `test_stability.py`). The single-layer reduction (`nl=1`, no background) is
   **`np.array_equal`** to `ShallowWater` over many steps (the by-construction rung-0 reduction; the
   single-layer engine is untouched, so this is a *meaningful cross-engine* check). Per-layer mass
   `∫h_k` is **machine-precision** conserved on the `background=None` path (flux-form telescoping).
   The **tight anchor**: with a supercritical thermal-wind background, a small `l=0` perturbation at
   the analytic most-unstable wavenumber grows at the **two-layer SW linear rate** — `σ = max Im ω`
   from the 6×6 dispersion matrix (`TwoLayerStability`, itself validated by zero-shear neutrality to
   machine precision, two-layer Poincaré recovery, the Eady coefficient `≈0.31`, and the f-plane
   no-critical-shear law) — within **a few %, converging with resolution** (slow). The two-layer
   **Poincaré dispersions** (external `√(gH_tot)`, internal `√(g'H_e)`) are reproduced too (slow) —
   the tight check on the Montgomery coupling.

## Validation boundary (what Phase 3 does *not* claim)

- **Energy *or* enstrophy, not both** (semi-discretely). This realization conserves
  energy; potential enstrophy drifts at a small, *spatial*-discretization-limited bound
  (it shrinks with Δx, not dt). Machine-precision conservation is claimed only for **mass**.
- **No exact material PV at the grid scale.** Material PV conservation is asserted as
  bounded (no spurious extrema) over a few inertial periods at moderate Rossby number — not
  to machine precision. Long, fully-turbulent, high-Rossby integrations can cascade enstrophy
  to the grid (the energy-conserving scheme's known behaviour); add hyperviscosity or an
  enstrophy-conserving operator for that regime (named, not built).
- **The tracer is passive and not monotone.** Advecting `θ` does not feed back on the dynamics
  (no buoyancy/active-tracer model — that is a reduced-gravity / multi-layer change). The
  centered scheme has **no flux limiter**, so a sharp tracer front **over/undershoots** (Gibbs);
  `∫hθ` (machine-exact) and `∫½hθ²` (bounded) are the claims, **not** boundedness/monotonicity of
  `θ`. A TVD/WENO limiter or hyperdiffusion is the named, unbuilt upgrade (it would break the
  energy-conserving symmetry — a deliberate non-goal).
- The *physical constants* (`g, H, f₀, β`, planetary values) are the **consumer's** to supply
  and validate (Planet `circulation.py`, against geostrophic balance + the jet benchmark).
- **N-layer Phase A is the *linear* growth rate, not the saturated flux.** `LayeredShallowWater`
  validates that baroclinic instability **exists and grows at the analytic rate** at small
  amplitude. The **saturated, irreversible eddy heat flux** (Phase B — the headline payoff: the
  rung-1 reduction-to-EBM made non-vacuous) is **not** claimed here: it needs the named-but-unbuilt
  **hyperviscosity** (a turbulent layered run cascades enstrophy to the grid — the inherited
  no-limiter behaviour, now load-bearing) and long post-saturation integrations. The background
  basic state is **f-plane** (`β=0`, matching the f-plane stability operator); a *finite* critical
  shear needs a β-capable PV-gradient treatment (named, not built). The layered engine carries **no
  passive tracer** (the interface displacement *is* the buoyancy at this rung) and is **not**
  energy-conserving *with a background* (it extracts APE — that growth is the signal). The free
  surface's fast **external** mode `√(g·H_tot)` sets the explicit step (the named compute cost; the
  rigid-lid elliptic-solve fork is the named within-rung upgrade).

## Units & scope

- **SI throughout** (`h, H` m; `u, v` m/s; `g` m/s²; `f₀` 1/s; `β` 1/(m·s); time s).
- **Scope (named rungs on the §5 staircase):** the single-layer engine is dry (no thermodynamic
  variable — *the* fact making Phase-4 coupling one-way); a **passive** scalar tracer *is* advected
  (rung 1, built); **vertical structure is built** as the sibling `LayeredShallowWater` (rung 3,
  Phase A — *active* buoyancy via the interface displacement), with the **saturated eddy flux**
  (Phase B) and **hyperviscosity** named-but-unbuilt. Still unbuilt: **rigid channel walls in y**
  (the classic baroclinic-lifecycle geometry — both engines are doubly-periodic β-plane; the layered
  growth test uses the `l=0` mode to sidestep them), the **rigid-lid / external-mode** fork, and the
  **sphere** (rung 5 — the pole problem). The stacked-field array boundary is the seam where each is
  slotted without touching consumers.

## Changelog

- **2026-06-09 — Planet Phase 3.** Initial build & validation (the guaranteed invariants above).
- **2026-06-10 — doctrine.** Status FROZEN → living, versioned contract (ADR 0005); no API change.
- **2026-06-10 — rung 1.** Added **passive flux-form tracer advection** (`SWState.tracer`): `step`
  advects `h·θ` through the same SSP-RK3, with `tracer_mass` / `tracer_variance` diagnostics and
  `test_tracer.py`. Additive — dry dynamics bit-for-bit unchanged; `step` no longer raises on a set
  tracer.
- **2026-06-12 — rung 3, Phase A (baroclinic).** Added the **N-layer** sibling engine
  `LayeredShallowWater` + `LayeredState` (`layered.py`) and the analytic anchor `TwoLayerStability`
  (`stability.py`), with `test_layered.py` / `test_stability.py`. Layers couple **only** through the
  Montgomery pressure; an optional default-off thermal-wind `background` injects the unstable basic
  state as Doppler `−U_k∂_x` + baroclinic `−G_k v'` coefficients. The single-layer `ShallowWater` /
  `SWState` are **untouched** (sibling design), so the `nl=1` reduction is **byte-identical**. Phase A
  validates the **linear growth rate** against the 6×6 dispersion matrix; the saturated Phase-B eddy
  flux + hyperviscosity are named-but-unbuilt (above).
