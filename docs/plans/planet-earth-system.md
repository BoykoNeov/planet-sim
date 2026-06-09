# Earth-System / Planet Simulator — Project Plan

> Per-project plan **#3** of the educational-simulator program — the **capstone**.
> Built to the **Section 10 template** of `ARCHITECTURE.md`; inherits Sections
> 2–9 as fixed invariants (compliance check in §8 below). This is the **third and
> last** project in build order (Steel → Microchip → Planet), and the **first
> project to build a *second* shared engine** (the fluid/PDE solver) and the first
> whose per-project-gate `uses` entry is genuinely multi-engine (ARCHITECTURE.md
> §4, §11).

---

## 1. One-line vision & the dramatic early win

**Vision.** *Planetary knobs in, climate & habitability out:* set a world's solar
constant, CO₂, obliquity, and land/ocean character, and watch a **climate state**
— the temperature field, the ice line, the wind systems, and ultimately a **map of
biomes** — emerge from energy balance and fluid dynamics, not a lookup table. This
is the planet analogue of Steel's *cooling curve → microstructure* and Chip's
*process → device*: here it is **forcing → climate → habitability**.

**The anchor demo (Phase 1's banked artifact).** *Snowball Earth — the hysteresis
of a frozen planet.* A 1-D latitudinal **energy-balance model** (EBM) with
diffusive heat transport and an **ice-albedo feedback** is driven by a slowly
**ramped solar constant**. As the sun dims, the ice line creeps equatorward — then,
past a threshold, the feedback runs away and the whole planet **freezes over in a
catastrophic jump** (Snowball Earth). Brighten the sun back up and the planet
**stays frozen far past where it froze** — a wide **hysteresis loop**, two stable
climates for one sun. One planet, one knob, a dramatic counter-intuitive
bifurcation — the cheapest, most teachable payoff in the project, and
simultaneously the integration test for every Phase-1 module. Its 0-D limit has an
**exact analytic equilibrium temperature**, and global energy balance is a **free
conservation check**.

**Why this is the right early win for the capstone.** It reuses the **frozen
diffusion spine verbatim** (latitudinal heat transport *is* a diffusion equation),
banks a complete planet-scale artifact before any new engine exists, and
establishes the validation-triad habit on ground with an exact analytic anchor —
exactly the build-order logic of ARCHITECTURE.md §4 ("the phasing-and-validation
discipline is internalized before the coupler is attempted").

---

## 2. Shared engines consumed

| Engine | Status here | Contract pointer |
|---|---|---|
| **Diffusion/heat (Fick / erfc)** — the program spine | **`[reuse frozen ✓ — Steel Phase 1a]`** | `engines/diffusion/CONTRACT.md`. Loaded as the **one-page contract**, never Steel's internals (ARCHITECTURE.md §6/§11). Planet instantiates **heat mode on the sphere**: `u = T`, the latitudinal transport `D·∂/∂x[(1−x²)∂T/∂x]` maps onto the engine's spatially-varying array diffusivity `D_eng(x) = D·(1−x²)` (x = sin φ); insulated **Neumann(0)** at both poles (no flux through a pole — the symmetry/conservation BC). The radiative source/sink is **not** carried by the engine (see the §3 Phase-1 splitting note). |
| **Fluid / PDE (shallow-water on a rotating frame)** | **`[to build & freeze HERE — Phase 3]`** | **`engines/fluid/CONTRACT.md` — the project's central new contribution.** A rotating-frame **shallow-water** solver: hyperbolic, **explicit** (CFL-limited), staggered **C-grid** — sharing *no* machinery with the parabolic-implicit diffusion engine (that *is* the point). Built standalone, validated against geostrophic balance / wave speeds / PV conservation, then **frozen** behind its `CONTRACT.md` before Phase 4 couples to it. Designed **extension-ready** (§3 Phase-3, §5): `state` as *stacked fields* so 1→N layers is a contract *extension*, and *tracer-ready* so moisture/temperature advection slots in for the deep-end rungs. |

**The one new shared engine in the whole portfolio-trio is built here.** Steel built
and froze the diffusion spine; Chip added only a *chip-local* module (Fourier
optics) and reused the spine; **Planet is where the second `engines/` member is
born** (ARCHITECTURE.md §5). It is promoted to `engines/` (not kept project-local
like Chip's `litho.py`) because it is a genuine general solver the wider portfolio
draws on (seismic, traffic, glacier, acoustics — §5 table) and Planet itself uses
it across Phases 3–4 + the documented climb.

> **No third shared engine.** The §5 toolkit table lists ODE integrators against
> Planet, but the time-stepping Planet needs is **folded into `engines/fluid`** (the
> shallow-water explicit integrator) and a small EBM relaxation loop that stays
> **planet-local** (`projects/planet/ebm.py`). Per invariant 5 / rule-of-three, a
> standalone ODE engine is **not** built here (1 use ≠ 3) — saying so explicitly so a
> future session does not read the §5 table as a mandate for a third engine build.

**Language & performance.** Python + NumPy/SciPy, per ADR 0001. Rung-0 compute is a
laptop-second job (1-D EBM relaxation; a coarse 2-D β-plane shallow-water run).
The array-out contract is kept strictly data-oriented (ADR 0001) because — uniquely
for Planet — that seam is *also* the runway for the GCM climb (§5): the heavier
rungs (3-D, many timesteps) are exactly where the ADR-0001 escalation
(vectorize → Numba/JAX/GPU → compiled core → separate process) gets exercised, each
localized behind one engine's contract.

---

## 3. Phases — each a complete, demonstrable artifact

Every phase names its **validation triad** concretely — an *analytical limit*, a
*conservation law*, and a *published benchmark* (invariant 3 / ARCHITECTURE.md §7) —
plus, per the program discipline (the way Steel/Chip recorded it), the
**non-circularity split** (what is *validated* tight vs *calibrated* and flagged) and
the **scope edge** (the regime where the model is honestly wrong, named not papered).

> **Phase ordering — the payoff is banked early (user decision, 2026-06-09).** The
> dramatic end-to-end win is the **biome map**, and it needs only the Phase-1 EBM
> temperature + a diagnostic precipitation field — **not** the fluid engine or the
> coupler. So biomes are **Phase 2**, before the riskier new-engine work, and the
> project banks a *complete, demonstrable planet* (climate + habitability map) before
> the shallow-water engine is attempted. This is invariant 2 ("fail gracefully into
> something real") used as deliberate risk management: if the new engine proves hard,
> a real artifact is already banked. The portfolio arc lists biomes *after* the
> coupler, but that is a *complexity* gradient, not a *dependency* one (§8 lists
> biomes as a peer module).

### Phase 1 — The EBM & the Snowball bifurcation (the foundation, spine reuse)

The 0-D global energy balance, then the **1-D latitudinal EBM**

```
C ∂T/∂t = D ∂/∂x[(1−x²) ∂T/∂x] + S(x)(1−α(x,T)) − (A + B·T)      x = sin φ
```

with **ice-albedo feedback** (`α` jumps high where `T < T_freeze`), driven by a
ramped solar constant to produce the **Snowball hysteresis**.

> **The reuse is real but not glib — it reuses the *Jominy splitting pattern* too
> (advisor, the crux).** The frozen engine solves `∂u/∂t = ∂/∂x(D ∂u/∂x) + S(x,t)`.
> The transport term maps cleanly onto the engine's array diffusivity `D·(1−x²)`. But
> the radiative terms **cannot be the engine's source**: `−(A+B·T)` is **linear in the
> state** and the albedo source `S(x)(1−α(T))` is a **nonlinear pointwise function of
> `T`**, while the engine's `S` is only `S(x,t)`, *not* `S(u)`. So the radiation is
> composed *around* the frozen engine by **Strang operator splitting** — **exactly the
> Phase-2a Jominy precedent**, which split a state-dependent lateral sink `−h(T−T_air)`
> around the same engine for this identical reason. The linear `−B·T` relaxation is an
> exact analytic exponential half-step (as Jominy's lateral sink was); the **albedo
> threshold makes the local step a pointwise *nonlinear* relaxation — which is what
> *creates* the multiple equilibria.** Phase 1 therefore reuses **both** the frozen
> engine *and* the frozen splitting idiom — a *stronger* reuse story than "the EBM is
> the engine," and the plan states it so the builder does not hit a wall on day one.

**Validation triad — Phase 1**
- *Analytical limit.* (a) The **0-D** global-mean equilibrium `T̄ = [(S₀/4)(1−ᾱ) − A]/B`
  — exact. (b) The **North (1975) two-mode** analytic EBM solution: with constant
  albedo and linear OLR, `T(x) = T₀ + T₂·P₂(x)` solves a 2×2 linear system exactly —
  the headline check that the **transport + linear radiation** are assembled
  correctly (validates the splitting in the *no-feedback* limit, where it is exact).
- *Conservation.* At equilibrium **net top-of-atmosphere flux integrates to zero**
  over the sphere: `∫(S(1−α) − (A+BT)) dx = 0` (absorbed solar = OLR globally). The
  **diffusive transport conserves structurally** — it only redistributes — guaranteed
  by the frozen engine's exact no-flux invariant (poles are Neumann(0)), re-confirmed
  for this BC pair.
- *Benchmark.* **climlab's EBM** (`EBM`/`EBM_seasonal`) — the present-day **ice-line
  latitude** (~70°), the **Snowball threshold** (a solar dimming of ~2–10% triggers
  runaway glaciation), and the **hysteresis width** (the small-/large-ice-cap
  instabilities of North 1975 / Budyko–Sellers). climlab consumed as a **reference
  tool** (the pycalphad pattern, §7), never copied.

**Non-circularity split.** What is *validated* (asserted tight): the **structural
reuse** (frozen-engine transport + the analytic two-mode solution it must reproduce)
and the **existence and qualitative structure of the hysteresis** (two stable branches,
a catastrophic jump — emergent, nothing but the feedback produces it). What is
*calibrated* (cited, flagged): the radiation/albedo constants `A, B, D, α_ice, α_open,
T_freeze` (Budyko/North/climlab values — `[[ebm-radiation-source]]`, pinned at build),
so the *exact threshold numbers* are calibration-dependent and asserted only in **loose
bands**, the way Steel's 1045 knee and Chip's contrast curve were.

**Scope edge, named.** The **linear OLR `A+B·T`** is a *parameterization* of
radiation, accurate only near the present climate; far from it (deep snowball, hot
states) it is increasingly wrong — this is the **rung-4** deferral (real radiative
transfer) on the §5 staircase, named not modeled. v1 is **annual-mean** (no
seasonal cycle / obliquity-driven insolation) unless the obliquity knob is added as a
seasonal-insolation option; the planet is **zonal-mean** (no land/ocean contrast,
no orography).

**Banked artifact:** the **Snowball hysteresis loop** — global-mean temperature (and
ice-line latitude) vs solar constant, traced **up and back down** via a
**parameter-continuation sweep**, the two branches and the catastrophic jumps marked;
beside the equilibrium `T(φ)` profile with its ice line. *Knob in, frozen-or-temperate
planet out.*

### Phase 2 — Biomes & habitability (the consequence payoff, banked early)

The end-to-end win: map the climate state to a **biome map**. Two inputs — the
**Phase-1 EBM temperature** `T(φ)` and a **diagnostic precipitation field** `P(φ)` —
feed a **Whittaker biome classifier** (`(T, P) → biome`), producing latitude bands of
biome that **migrate as the knobs turn**.

> **The precipitation is a *parameterization*, not a simulated water cycle (advisor,
> the honesty flag).** Neither the EBM (energy only) nor single-layer shallow-water
> (dry dynamics, no thermodynamic variable) produces precipitation; a real
> hydrological cycle is moist thermodynamics = the GCM tar pit. So `P(φ)` is a
> **prescribed kinematic diagnostic**: physically-motivated latitude bands — wet
> equatorial **ITCZ**, dry **subtropics ~25–30°** (the world's deserts), wet
> **midlatitude** storm tracks, dry **poles** — **scaled by the EBM temperature via
> Clausius–Clapeyron** (~7%/K), so warming shifts and intensifies the bands. It is
> **named as a parameterization that encodes the *known* observed precip-by-latitude
> structure, not as a derived field** — or the phase over-claims. (The deep-end hook:
> once the Phase-3 circulation exists, the pattern can be made **circulation-informed**
> — rain where the flow converges — a rung-1/2 enhancement at the array seam, still
> without real moisture physics.)

**Validation triad — Phase 2**
- *Analytical limit.* A **uniform-climate planet** (constant `T, P`) maps to a single
  biome; and the classifier's decision boundaries reproduce the **published Whittaker
  diagram** *exactly* at chosen `(T, P)` probe points — the classifier is an exact,
  testable partition of the `(T, P)` plane.
- *Conservation* (honestly weaker — a classifier has no energy/mass law; named as
  such). The **partition/consistency** check: every `(T, P)` maps to exactly one biome,
  the map **tiles the planet with no unclassified gaps**, and biome **area fractions sum
  to 1**. Plus the precip param's **global-water budget** (the C–C-scaled global mean
  moves as designed). This leg is a *consistency* check, not a conservation law — stated
  plainly, not dressed up as one.
- *Benchmark.* **Present-day Earth** (Phase-1 `T(φ)` + the diagnostic `P(φ)` at present
  insolation) reproduces the **observed major biome bands**: tropical rainforest at the
  equator, savanna/desert at ~15–30°, temperate forest/grassland in midlatitudes,
  boreal forest then tundra toward the poles. The **Whittaker biome diagram** is the
  reference fact.

**Non-circularity split.** *Validated* tight: the classifier as an exact `(T, P)`
partition and the present-Earth band ordering. *Calibrated/flagged:* the precip
band shape and the C–C scaling exponent (cited — `[[whittaker-biome-source]]`,
`[[precip-parameterization-source]]`); the **absolute biome latitudes** depend on the
calibrated precip param, asserted in loose bands.

**Scope edge, named.** **Whittaker** (annual `T, P`) is chosen over **Köppen**
(which needs seasonal/monthly precip the v1 annual-mean model does not produce) — named.
No continentality / orography (zonal-mean planet); precip prescribed (above); no
dynamic vegetation / carbon feedback.

**Banked artifact (the showcase):** the **biome-band map of the planet**, with a
**knob panel** — raise CO₂ (lower `A`) → biomes shift poleward, ice retreats; dim the
sun → the Snowball **white planet**; change obliquity → bands redistribute. This is the
project's dramatic end-to-end demo and the centerpiece of the **interactive map** (§9).

### Phase 3 — The shallow-water engine (build & freeze the new shared engine)

The project's **central new contribution** and its **risk phase**: a rotating-frame
**shallow-water** solver on a **β-plane** channel —

```
∂u/∂t + (u·∇)u + f×u = −g∇h            (+ optional forcing/drag)
∂h/∂t + ∇·(h u) = 0                     f = f₀ + βy   (β-plane)
```

— hyperbolic, **explicit** time stepping (CFL), staggered **C-grid**. Built standalone
in `engines/fluid/`, validated, then **frozen** behind its own `CONTRACT.md` before
Phase 4 depends on it (freeze-before-reuse, invariant 5).

**Validation triad — Phase 3**
- *Analytical limit.* **Geostrophic balance** (steady state: `f×u = −g∇h`) reproduced;
  **gravity-wave speed `√(gH)`**; **Rossby-wave dispersion** `ω = −βk/(k²+l²+1/L_R²)`;
  and **geostrophic adjustment** to the **Rossby radius** `L_R = √(gH)/f`.
- *Conservation.* **Mass `∫h`**, **energy** (KE + PE), and — the discriminating leg —
  **potential vorticity** `q = (ζ+f)/h` materially conserved.
- *Benchmark.* Classic shallow-water test cases: **Rossby's geostrophic-adjustment**
  problem (an initial height step relaxes to a geostrophic jet over `L_R`) and the
  **westward Rossby-wave phase speed** — vs the analytic dispersion relation and a
  standard SW reference (`[[shallow-water-source]]`, pinned at build).

> **The discriminating legs catch a wrong Coriolis (advisor).** Gravity-wave speed
> `√(gH)` alone does **not** exercise rotation; **PV conservation** and the **Rossby
> radius / geostrophic adjustment** are the legs that fail loudly if `f`, `β`, or the
> C-grid Coriolis averaging is wrong — so those are the seal, not `√(gH)`.

**Extension-ready contract (the GCM-climb foresight, §5).** The frozen
`engines/fluid/CONTRACT.md` declares `state` as a **stack of plain 2-D fields**
(`h`, `u`, `v`, and an *optional tracer slot*) so that (a) **1→N layers** is a
contract **extension**, not a rewrite (rung 3 — baroclinic), and (b) an **advected
tracer** (moisture/temperature) slots in (rungs 1–2). The β-plane geometry is the
rung-0 choice; the **sphere** is the rung-5 swap (it brings the pole problem — a real
numerical tar pit — named, not built). Like the diffusion contract, **only plain
arrays cross the per-step boundary** (ADR 0001), which is what makes both the layer
extension and a future compiled core slot in without touching consumers.

**Scope edge, named.** Single-layer, **dry** (no thermodynamic variable → no heat
transport on its own — *the* fact that makes Phase 4's coupling one-way and the
two-way version a stretch), **β-plane** (not sphere), **explicit** (CFL-limited). All
are named rungs on the §5 staircase.

**Banked artifact:** a **geostrophic-adjustment animation** (height step → balanced jet
over `L_R`) and a **Rossby wave** propagating westward — beside the conservation
diagnostics (mass / energy / PV) holding flat.

### Phase 4 — The one-way coupler (climate forces circulation; jets emerge)

Couple the two engines: the **EBM's meridional temperature gradient sets a
pressure/height field that forces the shallow-water flow → midlatitude jets emerge.**
**One-way** in v1 (climate → circulation); the two-way feedback is **rung 1** of the
§5 climb, seamed not built.

> **Loose coupling at cadence, never one time loop (§8 mandate).** The EBM (slow,
> radiative-thermal) and the shallow-water flow (fast, dynamical) **exchange boundary
> conditions through a coupler at appropriate cadences** — the EBM relaxes to
> quasi-equilibrium, hands the flow a height field, the flow spins up its jets; they do
> **not** share one time step. This respects the timescale separation ARCHITECTURE.md
> §8 makes non-negotiable, and is the structural seam the GCM climb extends.

**Validation triad — Phase 4**
- *Analytical limit.* **Geostrophic balance of the emergent jet.** The EBM's meridional
  temperature gradient maps to a height field `h(φ)`, and the steady forced flow must
  settle into `f×u = −g∇h` — a jet at the latitude the climate gradient puts it. This is
  *distinct* from Phase 3's seal (which tests geostrophic *adjustment in isolation*):
  Phase 4 tests that the **coupled** system produces a balanced jet in the right place.
- *Conservation.* The **shallow-water invariants** (mass `∫h`, PV, energy) preserved
  under the steady EBM forcing — the frozen engine's guarantees re-confirmed in the
  forced/coupled configuration.
- *Benchmark.* **Jet latitude and strength** vs the observed midlatitude jet
  (~30–45°, tens of m/s — loose).

> **Why one-way banks cleanly, and what it honestly *cannot* claim (advisor).** One-way
> has a clean **geostrophic-balance anchor** (a coupled jet that must balance the imposed
> height field) — checkable, so the phase banks a *validated* artifact. Crucially, three
> legs a richer coupler *would* add are **not available here and are not claimed**:
> **poleward heat transport** and the **reduction-to-diffusive-EBM limit** need an
> advected thermodynamic variable the **dry** layer lacks (the very fact that makes the
> coupling one-way — there is no genuine `v·T` flux, and any diagnosed value would scale
> with the arbitrary forcing amplitude = tuning, not validation), and **thermal-wind
> balance** needs vertical shear a **single layer** does not have (it is a multi-layer
> concept). Those are **rung-1 / rung-3 anchors** (§5), not v1 ones — moving them up
> *strengthens* the climb rather than padding Phase 4 with legs it cannot support. The
> full two-way coupler ("derive the EBM's diffusivity from resolved eddies") needs that
> same advected tracer and has no clean exact anchor — the "doesn't bank cleanly" failure
> mode §7 exists to prevent. So two-way is **rung 1**, designed-for at the tracer seam,
> not built in v1.

**Scope edge, named.** One-way only; the geometry bridge (1-D zonal-mean EBM ↔ 2-D
β-plane channel) is handled by forcing the channel's meridional structure from the
zonal-mean gradient — a deliberate reduced coupling, not a full 2-D climate.

**Banked artifact:** the **EBM temperature field forcing the shallow-water → an
emergent, geostrophically-balanced jet**, drawn on the interactive map (circulation
streamlines over the temperature field), with the mass/PV/energy invariants holding
flat and the jet-latitude diagnostic.

---

## 4. Module map & contracts

Small files, so any single task loads with its neighbours' *contracts*, not their
internals (ARCHITECTURE.md §6). Mirrors `projects/steel/` and `projects/chip/`.

```
BigSim/
  engines/
    diffusion/CONTRACT.md     # the FROZEN spine Planet reuses (load this, not steel/)
    fluid/                    # the NEW shared engine, built & frozen in Phase 3
      shallowwater.py         #   rotating β-plane shallow water (C-grid, explicit)
      CONTRACT.md             #   the FROZEN one-page API (extension-ready: stacked fields, tracer slot)
      tests/                  #   geostrophic balance, wave speeds, PV/mass/energy — the seal
  projects/planet/
    ebm.py                    # 0-D + 1-D latitudinal EBM; frozen-diffusion transport + Strang-split radiation  (Phase 1)
    albedo.py                 # ice-albedo feedback + the Snowball hysteresis continuation sweep                 (Phase 1)
    precip.py                 # diagnostic precipitation parameterization (named kinematic; C–C-scaled)          (Phase 2)
    biomes.py                 # Whittaker classifier: (T, P) → biome map                                          (Phase 2)
    circulation.py            # planet's instantiation of engines/fluid (β-plane channel, planetary params)      (Phase 3)
    coupler.py                # one-way EBM → shallow-water forcing (cadence-based, timescale-separated)          (Phase 4)
    plots.py                  # planet-local static figures (→ promote to viz/ by rule-of-three)
    planetmap.py              # the deep-end INTERACTIVE map (Plotly/web) — §9 chosen surface
    planet.ipynb              # single teaching notebook (sliders → climate → biome map)
    demo_snowball.py / demo_biomes.py / demo_shallowwater.py / demo_coupler.py    # banked artifacts
    climate_reference.py      # frozen climlab reference table (keeps the triad green without the [climate] extra)
    README.md                 # per-module map + per-session load pointer
    tests/                    # the validation triads (the seal)
  pyproject.toml              # testpaths += projects/planet, engines/fluid; [climate] + [webviz] extras
```

**Contracts kept short.** Each module's docstring is its contract (the steel/chip
convention). `ebm.py` loads `engines/diffusion/CONTRACT.md`; `circulation.py` loads
`engines/fluid/CONTRACT.md`; the rest exchange **plain arrays** — a temperature
profile `T(φ)`, a precip field `P(φ)`, a biome map, a height/velocity field — the
loose-coupling currency (§5 / ADR 0001).

---

## 5. Scope ceiling — the GCM staircase (consequence now, climb documented)

**The named tar pit (ARCHITECTURE.md §8):** a monolithic **GCM-grade** 3-D
primitive-equation Earth-system model — full moist physics, radiative transfer,
sub-grid convection/cloud parameterizations, an ocean GCM, sea ice, carbon cycle. It
is squarely **infeasible-tier**: a compute wall (clusters, simulated decades), a
**context-coherence** wall (a GCM codebase vastly exceeds an agent's context window),
and a **validation** wall (the exact analytic anchors fade; you rely on inter-model
comparison + reanalysis).

**What v1 targets instead — the consequence (rung 0):** *a loosely-coupled reduced
climate* — 1-D EBM + 2-D dry shallow-water + a one-way coupler + diagnostic precip +
a biome map. A learner sees *"ice line at 70°, a midlatitude jet, tropical rainforest
to polar tundra,"* and *why* (which knob moved it) — not a meshed primitive-equation
field. v1 is decided to be **rung 0** (user, 2026-06-09).

**The staircase — the scope ceiling is a documented *growth axis*, not a wall**
(ARCHITECTURE.md §8: "deferrals, not permanent foreclosures... the same module seam is
where a deferred heavy regime is later slotted"). Each rung slots at the ADR-0001
array seam; each banks an artifact with a (progressively weakening) anchor:

| Rung | Adds | Unlocks | Anchor | Wall |
|---|---|---|---|---|
| **0 (v1)** | EBM + dry SW + one-way coupler + diagnostic precip + biomes | hysteresis, jets, biome map | 0-D temp; geostrophic balance; climlab | laptop |
| **1** | two-way coupler (advected T tracer; close budget to EBM; circulation-informed precip) | emergent transport | poleward ~5–6 PW; reduction-to-EBM | laptop; validation loosens |
| **2** | moist dynamics (moisture var; evap/condense; latent heat) | **emergent precipitation** | global `E=P`; C–C scaling | first sub-grid closure |
| **3** | vertical structure (stack layers → multi-level) | **baroclinic instability** = real storms | Eady/Charney growth rates; Held–Suarez | compute climbs (3-D); compiled-core escalation |
| **4** | real radiation (gray/band over the column) | computed CO₂/H₂O forcing + feedbacks | gray-atmosphere (Schwarzschild); RCE | a radiation engine |
| **5** | idealized GCM (sphere core; convection/PBL/cloud params; slab ocean + sea ice) | aquaplanet GCM | inter-model (Held–Suarez); reanalysis | leaves the laptop |
| **6** | full GCM / ESM (topography, full moist physics, clouds, ocean GCM, carbon cycle) | operational climate model | CMIP + observations | the infeasible end |

**The seam carries the climb.** `engines/fluid`'s **stacked-field, tracer-ready**
state (built in Phase 3) turns rungs 1–3 from rewrites into contract extensions; the
ADR-0001 escalation is the compute runway for rungs 3+; the §6 contract/test
discipline is the *only* thing that keeps the upper rungs coherent for an
agent — which is exactly why the climb is **one validated rung at a time**, never a
leap. The one-way→two-way coupler progression of Phase 4 → rung 1 is the first step
of this climb. Concretely, the validation legs Phase 4 **cannot honestly support**
under one-way + dry single-layer live on the staircase: **poleward heat transport**
and the **reduction-to-diffusive-EBM limit** arrive at **rung 1** (with the advected
tracer), and **thermal wind** at **rung 3** (with vertical structure) — so the climb is
where those anchors belong, not a v1 phase.

---

## 6. Terms-of-use status

**Clean per ARCHITECTURE.md §9 — no export-control dimension** (climate/planetary
science is published fundamental science, like steel; the §9 carve-out is not even
needed). Copyright a non-issue (equations and physical facts; original code/prose; no
verbatim listings/figures).

**The dataset diligence item — the recurring §9 one.** Reference facts (EBM
radiation constants; Whittaker biome boundaries; shallow-water test-case parameters)
are **cited and pinned at build** (`[[…-source]]` notes), used for comparison, not
redistributed. The one thing to **license-check before redistributing** is any
**observed dataset** used as a benchmark — a **biome / land-cover map**, **topography**,
or **climate reanalysis** (§9 names topography/reanalysis explicitly). *Action:* use
only openly-licensed datasets (e.g. Köppen-Geiger/Whittaker reference figures as
*facts*, not copied; any reanalysis as validation-only), and **never commit** a
redistribution-restricted dataset — the analogue of Steel's CALPHAD-TDB rule.
**climlab** (the reference tool) is for **validation**, never copied; it is consumed
behind the optional `[climate]` extra with a **frozen reference table** so the triad
stays green without it (§7).

---

## 7. Test runner — the per-project gate (and engine #2)

The **per-project gate** (ADR 0003 → *Successor*, `tools/gate.py`, built at the end of
Microchip):

```powershell
python -m tools.gate planet -m "not slow"   # routine commit gate: planet's own tests
                                             #   + the tests of every engine/module it uses
python -m tools.gate planet                  # full gate for planet — EXCEPTIONAL
./run_tests.ps1                              # whole-repo full gate — release / CI
```

**Planet is the manifest's most informative entry yet.** The `tools/gate.py` `GATES`
manifest gets a `planet` row whose **`uses` is genuinely multi-engine** —
`{engines/diffusion, engines/fluid}` — the **first** entry that is not a second *row*
of the single-engine pattern (Steel and Chip both used only `engines/diffusion`). This
is the case the per-project gate was *designed to validate* (ADR 0003); a `planet`
commit runs planet's tests + **both** engines' seals, while a Steel/Chip commit still
need not run `engines/fluid`'s.

**Two gate-infrastructure consequences of building engine #2:**
- **`engines/fluid` joins the manifest as a shared engine**, and editing it is the
  **cross-cutting case that triggers the full gate** (its seal is the contract Phase 4
  + the GCM climb rely on), exactly as editing `engines/diffusion` does.
- **The import-drift guard, deferred to "engine #2," is now built** (memory
  `[[test-execution-policy]]`): with two engines, the guard that checks the
  hand-declared `uses` manifest against the actual imports a project makes finally has a
  second engine to discriminate (a `planet` that imports `engines.fluid` but forgot to
  declare it is the failure the guard catches). This lands with Phase 3 (when
  `engines/fluid` first exists and `circulation.py` imports it).

`pyproject.toml` gains `projects/planet` and `engines/fluid` in `testpaths`; the
existing `pythonpath = ["."]` lets planet import both engines with no install step. New
tests that drive a live external solver / kernel / subprocess / the `[climate]` or
`[webviz]` stack get the **`slow`** marker (the climlab live cross-check; the notebook
smoke-test; any long shallow-water integration).

---

## 8. Invariant-compliance check (against ARCHITECTURE.md §2–9 — not re-litigated)

| Program invariant | How this plan honors it |
|---|---|
| 1 — build toolkit once, solver-heavy first | Reuses the frozen diffusion spine (Phases 1–2, 4) **and builds the program's one remaining shared engine** (`engines/fluid`, Phase 3) — the capstone's structural job. |
| 2 — phase so each stage banks a working artifact | Four phases, each an explicit banked artifact (hysteresis loop, biome map, geostrophic-adjustment animation, emergent jet). **Payoff banked early** (biomes = Phase 2) so the project degrades gracefully if the new engine proves hard. |
| 3 — validation triad from day one | Instantiated *concretely per phase* in §3 (analytic + conservation + benchmark), each with its non-circularity split + scope edge. The Phase-2 "conservation" leg is honestly flagged as a consistency/partition check, not a fabricated law. |
| 4 — target consequence where mechanism is a wall | §5: a reduced coupled climate + biome map instead of a GCM — and the ceiling is written as a **documented staircase**, the consequence-now / climb-later form of §8. |
| 5 — reuse only frozen modules | Reuses `engines/diffusion/CONTRACT.md` (sealed in Steel 1a); **freezes `engines/fluid` behind its own seal before Phase 4 couples to it.** |
| 6 — updating docs is part of every change | This plan + per-module READMEs + `engines/fluid/CONTRACT.md` + `docs/decisions/` + the ARCHITECTURE.md §11 pointer are maintained per change. |
| Terms of use (§9) | §6: clean (no export dimension); the lone diligence item is an **observed biome/topography/reanalysis dataset** license-check — flagged, the CALPHAD-DB analogue. |

---

## 9. Visualization & UX

Per ARCHITECTURE.md §12 / ADR 0002: compute stays headless; views consume the engines'
plain arrays; a figure is never in the correctness path.

- **Floor (universal):** the §3 banked figures — the Snowball hysteresis loop, the
  `T(φ)` profile + ice line, the biome-band map, the geostrophic-adjustment / Rossby
  animation, the emergent jet — as **static matplotlib figures** (the `[viz]` extra),
  testable against numeric output.
- **Mechanism views** (the "teach *why*" target, ADR 0002 §5): the **hysteresis traced
  by the continuation sweep** (why a dimming sun *jumps* to a frozen planet, and why it
  stays frozen on the way back — the bistability *seen*, not stated), the **biome bands
  migrating** as a knob turns (why deserts sit at 30°), and the **jet emerging from the
  EBM-imposed height field** (geostrophic balance made visible — thermal wind itself is
  the rung-3 upgrade).
- **Experimentation — the deep-end interactive map (user decision, 2026-06-09).**
  Planet's payoff is inherently a **map**, exactly the case ADR 0002 §4 reserves the
  *selective deep-end* (Plotly / web) for ("planet maps"). So Planet ships, **beyond**
  the floor + a single teaching notebook (`planet.ipynb`, the chip-style thin skin):
  an **interactive planet map** (`planetmap.py`) — knob sliders (solar constant, CO₂,
  obliquity) drive a live **temperature / ice-line / biome map**, with the Phase-4
  **circulation overlay** when coupled. This is the one project in the trio that earns
  the deep-end viz; it is **opt-in behind a `[webviz]` extra**, consumes only the
  validated array outputs (never a live solver object — the ADR-0001/0002 data
  boundary), and its test is an **execution smoke-test**, not a physics check (the
  triads in §3 are the validation).
- **Toolkit:** plot primitives start planet-local in `plots.py`; the **2-D field /
  heatmap** and **time-animation** primitives are exactly the ADR-0002 §3 candidates a
  third reuse (after steel/chip) would **promote to the shared `viz/`** by
  rule-of-three (ARCHITECTURE.md §6).

Responsiveness is free at rung 0: compute is laptop-seconds, so knob → re-run →
re-map needs no special engineering (ADR 0001 scope). (Rungs 3+ are where the
ADR-0001 escalation is exercised — §5.)

---

## 10. Immediate next step

**Plan banked (this document).** The build order (ARCHITECTURE.md §4) advances from a
100%-complete Steel + Microchip to the **Planet capstone**, with all four scope
decisions locked (2026-06-09): **biomes banked early** (Phase 2, before the new
engine), **one-way coupler in v1** (two-way = rung 1), **v1 = rung 0** with the **GCM
staircase documented** as the growth axis (§5), and the **deep-end interactive map** as
the chosen viz surface (§9).

**Phase 1 — the EBM & the Snowball bifurcation.** Build `ebm.py` (0-D global balance →
1-D latitudinal EBM, the **frozen-diffusion transport + Strang-split radiation** — the
Jominy-2a idiom reused) and `albedo.py` (ice-albedo feedback + the **continuation-sweep
hysteresis**). Validation triad: the **0-D + North two-mode** analytic anchors, **net-TOA
global energy balance**, and the **climlab** ice-line / Snowball-threshold / hysteresis
benchmark (behind the `[climate]` extra + a frozen reference table). This is **low-risk
spine reuse** that banks a complete planet-scale artifact and re-proves the program
thesis — the diffusion spine reuses a *third* time, now as a sphere's heat transport.

**Phase-1 reference sources — to pin at build (the `[[…-source]]` discipline, not
carried from memory):** the **EBM radiation/albedo constants** (`A, B, D, α, T_freeze` —
Budyko 1969 / North 1975 / climlab defaults → `[[ebm-radiation-source]]`); later phases
pin **`[[whittaker-biome-source]]`** + **`[[precip-parameterization-source]]`** (Phase 2)
and **`[[shallow-water-source]]`** (Phase 3 test-case parameters). Gather and pin each
when its phase is built, as Steel/Chip did.
