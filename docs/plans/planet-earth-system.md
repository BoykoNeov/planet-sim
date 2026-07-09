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

**Why this is the right early win for the capstone.** It reuses the **diffusion
spine verbatim** (latitudinal heat transport *is* a diffusion equation),
banks a complete planet-scale artifact before any new engine exists, and
establishes the validation-triad habit on ground with an exact analytic anchor —
exactly the build-order logic of ARCHITECTURE.md §4 ("the phasing-and-validation
discipline is internalized before the coupler is attempted").

---

## 2. Shared engines consumed

| Engine | Status here | Contract pointer |
|---|---|---|
| **Diffusion/heat (Fick / erfc)** — the program spine | **`[reuse unchanged ✓ — Steel Phase 1a]`** | `engines/diffusion/CONTRACT.md`. Loaded as the **one-page contract**, never Steel's internals (ARCHITECTURE.md §6/§11). Planet instantiates **heat mode on the sphere**: `u = T`, the latitudinal transport `D·∂/∂x[(1−x²)∂T/∂x]` maps onto the engine's spatially-varying array diffusivity `D_eng(x) = D·(1−x²)` (x = sin φ); insulated **Neumann(0)** at both poles (no flux through a pole — the symmetry/conservation BC). The radiative source/sink is **not** carried by the engine (see the §3 Phase-1 splitting note). |
| **Fluid / PDE (shallow-water on a rotating frame)** | **`[to build HERE — Phase 3]`** | **`engines/fluid/CONTRACT.md` — the project's central new contribution.** A rotating-frame **shallow-water** solver: hyperbolic, **explicit** (CFL-limited), staggered **C-grid** — sharing *no* machinery with the parabolic-implicit diffusion engine (that *is* the point). Built standalone, validated against geostrophic balance / wave speeds / PV conservation, then **sealed** behind its `CONTRACT.md` before Phase 4 couples to it. Designed **extension-ready** (§3 Phase-3, §5): `state` as *stacked fields* so 1→N layers is a contract *extension*, and *tracer-ready* so moisture/temperature advection slots in for the deep-end rungs. |

**The one new shared engine in the whole portfolio-trio is built here.** Steel built
and sealed the diffusion spine; Chip added only a *chip-local* module (Fourier
optics) and reused the spine; **Planet is where the second `engines/` member is
born** (ARCHITECTURE.md §5). It is promoted to `engines/` (not kept project-local
like Chip's `litho.py`) because it is a genuine general solver the wider portfolio
draws on (seismic, traffic, glacier, acoustics — §5 table) and Planet itself uses
it across Phases 3–4 + the documented climb.

> **No third shared engine.** The §5 toolkit table lists ODE integrators against
> Planet, but the time-stepping Planet needs is **folded into `engines/fluid`** (the
> shallow-water explicit integrator) and a small EBM relaxation loop that stays
> **planet-local** (`planet/ebm.py`). Per invariant 5 / rule-of-three, a
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
> (advisor, the crux).** The shared engine solves `∂u/∂t = ∂/∂x(D ∂u/∂x) + S(x,t)`.
> The transport term maps cleanly onto the engine's array diffusivity `D·(1−x²)`. But
> the radiative terms **cannot be the engine's source**: `−(A+B·T)` is **linear in the
> state** and the albedo source `S(x)(1−α(T))` is a **nonlinear pointwise function of
> `T`**, while the engine's `S` is only `S(x,t)`, *not* `S(u)`. So the radiation is
> composed *around* the engine by **Strang operator splitting** — **exactly the
> Phase-2a Jominy precedent**, which split a state-dependent lateral sink `−h(T−T_air)`
> around the same engine for this identical reason. The linear `−B·T` relaxation is an
> exact analytic exponential half-step (as Jominy's lateral sink was); the **albedo
> threshold makes the local step a pointwise *nonlinear* relaxation — which is what
> *creates* the multiple equilibria.** Phase 1 therefore reuses **both** the
> engine *and* the splitting idiom — a *stronger* reuse story than "the EBM is
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
  by the engine's exact no-flux invariant (poles are Neumann(0)), re-confirmed
  for this BC pair.
- *Benchmark.* **climlab's EBM** (`EBM`/`EBM_seasonal`) — the present-day **ice-line
  latitude** (~70°), the **Snowball threshold** (a solar dimming of ~2–10% triggers
  runaway glaciation), and the **hysteresis width** (the small-/large-ice-cap
  instabilities of North 1975 / Budyko–Sellers). climlab consumed as a **reference
  tool** (the pycalphad pattern, §7), never copied.

**Non-circularity split.** What is *validated* (asserted tight): the **structural
reuse** (engine transport + the analytic two-mode solution it must reproduce)
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

### Phase 3 — The shallow-water engine (build the new shared engine)

The project's **central new contribution** and its **risk phase**: a rotating-frame
**shallow-water** solver on a **β-plane** channel —

```
∂u/∂t + (u·∇)u + f×u = −g∇h            (+ optional forcing/drag)
∂h/∂t + ∇·(h u) = 0                     f = f₀ + βy   (β-plane)
```

— hyperbolic, **explicit** time stepping (CFL), staggered **C-grid**. Built standalone
in `engines/fluid/`, validated, then sealed behind its own `CONTRACT.md` (its passing suite)
before Phase 4 depends on it (build-and-validate before reuse; engines are *living* contracts,
ADR 0005).

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

**Extension-ready contract (the GCM-climb foresight, §5).** The
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
- *Conservation* (**reframed at build, 2026-06-09 — advisor-blessed honesty call**). The
  original wording ("mass, PV, energy preserved under the steady forcing") is internally
  inconsistent with this phase's own *forced* prose: a **forced–dissipative** system does **not**
  conserve energy or PV — the forcing (thermal relaxation) injects available potential energy and the
  drag removes kinetic energy, and *that balance is precisely what selects the steady jet*. So the leg
  is reframed to what is actually true (the Phase-3 "energy *or* enstrophy, as measured" honesty
  class): **(a)** mass `∫h` is **machine-exact under forcing** (a discretely zero-mean height target +
  the engine's exact mass invariant), and **(b)** a **release test** — switch the forcing & drag *off*
  and run the bare engine: mass / energy / potential enstrophy are conserved (the engine's
  Phase-3 guarantees re-confirmed in the coupled configuration) **and the jet persists**, proving it is
  a genuine balanced state, not a forcing-propped artifact. *That* is "the engine's guarantees
  re-confirmed", read honestly.
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
    diffusion/CONTRACT.md     # the diffusion spine Planet reuses (load this, not steel/)
    fluid/                    # the NEW shared engine, built & validated in Phase 3
      shallowwater.py         #   rotating β-plane shallow water (C-grid, explicit)
      CONTRACT.md             #   the one-page API (living/versioned; extension-ready: stacked fields, tracer slot)
      tests/                  #   geostrophic balance, wave speeds, PV/mass/energy — the seal
  planet/
    ebm.py                    # 0-D + 1-D latitudinal EBM; diffusion-spine transport + Strang-split radiation  (Phase 1)
    albedo.py                 # ice-albedo feedback + the Snowball hysteresis continuation sweep                 (Phase 1)
    precip.py                 # diagnostic precipitation parameterization (named kinematic; C–C-scaled)          (Phase 2)
    biomes.py                 # Whittaker classifier: (T, P) → biome map                                          (Phase 2)
    circulation.py            # planet's instantiation of engines/fluid (β-plane channel, planetary params)      (Phase 3)
    coupler.py                # one-way EBM → shallow-water forcing (cadence-based, timescale-separated)          (Phase 4)
    plots.py                  # planet-local static figures (a future planet/viz/ once a 3rd geometry consumer appears — §9.4)
    planetmap.py              # the deep-end INTERACTIVE map: a LAYER REGISTRY painted by Plotly+ipywidgets (ADR 0004 #1, §9.1)
    planet_spec.py            # the planet-spec interchange schema: export/import the layer stack; round-trip-identity tested (ADR 0004 #3-4, §9.3)
    flow_serialize.py         # a vector flow field (u,v)+coverage+provenance THROUGH the planet-spec schema — the producer-agnostic viz/output seam (R1, §11.2)
    planet.ipynb              # single teaching notebook (sliders → climate → biome map)
    demo_snowball.py / demo_biomes.py / demo_shallowwater.py / demo_coupler.py    # banked artifacts
    climate_reference.py      # frozen climlab reference table (keeps the triad green without the [climate] extra)
    README.md                 # per-module map + per-session load pointer
    tests/                    # the validation triads (the seal)
  pyproject.toml              # testpaths += planet, engines/fluid; [climate] + [webviz] extras
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

**Two staircase consequences for the interactive map (ADR 0004 / §9).** (1) **The
interaction model changes with the rung, the renderer does not.** At rung 0 a slider
drives an instant recompute-and-remap; as compute climbs (rungs 3+: 3-D, many
timesteps), that live loop becomes impossible and the model shifts to **set parameters →
launch a run → view the result** — but the map renderer, a pure array consumer, is
unchanged (named so a future session does not engineer to preserve the live loop past
where compute allows). (2) **"Editable planet" is itself a climb on this staircase, not a
v1 feature.** The cheap tier — elevation as a **lapse-rate** map diagnostic, land/ocean
as an **albedo difference**, land/ocean fraction per band as continentality-lite — rides
**rungs 0–1** (no engine change); **ocean heat capacity → seasonality** needs the
seasonal cycle the annual-mean v1 lacks (the §3 scope edge); and **true longitudinal
geography → regional climate, orographic precip, rain shadows** (the user's north star)
is the **rung-5** exit from the zonal-mean planet — new transport that leaves the
1-D engine. Until those rungs, an imported/edited geography is **inert** (§9.3).

**A sibling growth axis — off Earth (scoped, not built).** The shallow-water and two-layer-QG engines
(Phase 3 / rung 3) are the *same idealized model family* the literature uses for **gas-giant
atmospheres** (banded jets, the Great Red Spot). A feasibility sketch — three tiers (β-plane mechanism =
~one rung on `baroclinic_qg.py`; sphere-correct globe = a new geometry engine; deep convective interior
= a *steeper reach, not out of scope* — the Busse-annulus QG model and rotating Rayleigh–Bénard are
reduced laptop-scale entries, only the realistic anelastic-deep-shell + dynamo regime is the wall), the
load-bearing "isotropic condensate ≠ zonal jets" correction, and off-the-shelf options
(pyqg / Dedalus / EPIC / MagIC) — is recorded in
[`docs/explorations/gas-giant-atmosphere.md`](../explorations/gas-giant-atmosphere.md).

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

`pyproject.toml` gains `planet` and `engines/fluid` in `testpaths`; the
existing `pythonpath = ["."]` lets planet import both engines with no install step. New
tests that drive a live external solver / kernel / subprocess / the `[climate]` or
`[webviz]` stack get the **`slow`** marker (the climlab live cross-check; the notebook
smoke-test; any long shallow-water integration).

---

## 8. Invariant-compliance check (against ARCHITECTURE.md §2–9 — not re-litigated)

| Program invariant | How this plan honors it |
|---|---|
| 1 — build toolkit once, solver-heavy first | Reuses the shared diffusion spine (Phases 1–2, 4) **and builds the program's one remaining shared engine** (`engines/fluid`, Phase 3) — the capstone's structural job. |
| 2 — phase so each stage banks a working artifact | Four phases, each an explicit banked artifact (hysteresis loop, biome map, geostrophic-adjustment animation, emergent jet). **Payoff banked early** (biomes = Phase 2) so the project degrades gracefully if the new engine proves hard. |
| 3 — validation triad from day one | Instantiated *concretely per phase* in §3 (analytic + conservation + benchmark), each with its non-circularity split + scope edge. The Phase-2 "conservation" leg is honestly flagged as a consistency/partition check, not a fabricated law. |
| 4 — target consequence where mechanism is a wall | §5: a reduced coupled climate + biome map instead of a GCM — and the ceiling is written as a **documented staircase**, the consequence-now / climb-later form of §8. |
| 5 — reuse validated modules behind their `CONTRACT.md` | Reuses `engines/diffusion/CONTRACT.md` (validated in Steel 1a); **builds & validates `engines/fluid` behind its own suite before Phase 4 couples to it** (engines are *living* contracts — ADR 0005). |
| 6 — updating docs is part of every change | This plan + per-module READMEs + `engines/fluid/CONTRACT.md` + `docs/decisions/` + the ARCHITECTURE.md §11 pointer are maintained per change. |
| Terms of use (§9) | §6: clean (no export dimension); the lone diligence item is an **observed biome/topography/reanalysis dataset** license-check — flagged, the CALPHAD-DB analogue. |

---

## 9. Visualization & UX

Per ARCHITECTURE.md §12 / ADR 0002: compute stays headless; views consume the engines'
plain arrays; a figure is never in the correctness path. Planet is the **first project
to reach the deep end** of ADR 0002 §4, so its interactive map is governed by the
doctrine `docs/decisions/0004-interactive-maps-and-state-interchange.md` (ADR 0004 —
the layer registry, the tier-dependent interaction model, and the state-interchange
schema, all first-instanced here). The map's design was converged with the user
(2026-06-09); the decisions and their rationale are recorded inline below.

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

### 9.1 The deep-end interactive map (the chosen surface)

Planet's payoff is inherently a **map**, exactly the case ADR 0002 §4 reserves the
*selective deep-end* for. Beyond the floor + a single teaching notebook (`planet.ipynb`,
the chip-style thin skin), Planet ships an **interactive planet map** (`planetmap.py`),
**opt-in behind a `[webviz]` extra**, consuming only validated array outputs (never a
live solver object — the ADR 0001/0002 boundary). Its render test is an **execution
smoke-test**, not a physics check (the §3 triads are the validation) — *except* the
state-interchange round-trip, which is a real invariant and is tested as one (below).

- **Tech (D5):** **Plotly + ipywidgets** — a genuinely interactive globe (rotate / zoom /
  hover) driven by knob-sliders, runnable in a notebook *and* as a thin standalone page.
  Matplotlib is too weak for the deep-end map; a full custom web app is the *editable*
  future (§9.3), not v1. All behind the `[webviz]` extra (the `[viz]`/`[climate]`
  pattern) so the core sim and the test suite never depend on it.
- **Layers, not a monolith (D4 / ADR 0004 #1).** The renderer is a generic painter over
  a stack of `(name, kind, data-array, style, z-order)` layers, `kind ∈ {scalar field,
  vector/line overlay, annotation}`. Each phase **registers** layers — temperature & ice
  mask (Phase 1) → biome & precipitation fields (Phase 2) → circulation streamlines /
  jet axis (Phase 4) → elevation/bathymetry & coastlines (the geography seam, §9.3) —
  and **never edits the renderer**. This *is* the user's "show more features as the
  phases progress" requirement made structural. The registry stays planet-local; it would
  be the **third consumer** that eventually promotes the 2-D-field / animation primitives
  to a planet-internal `planet/viz/` by rule-of-three — *not* a cross-repo `viz/` (that
  pre-split shared package stayed in the monorepo archive; §9.4 post-split note). It does
  not pre-empt that (ADR 0004): the geometry helpers are still two-consumer (§9.4).
- **Knobs — the "knob in, climate out" panel, and the exoplanet sandbox.** v1 sliders:
  **solar constant `S₀`** (the *amount* of radiation — already the Snowball lever),
  **CO₂** (→ lower `A`), **obliquity**, and transport `D`. Two further "exoplanet" knobs
  are *named here, built when their phase/source is pinned* (the `[[…-source]]`
  discipline): **stellar spectrum / type** enters as an **albedo modifier** (a redder /
  cooler star → near-IR ice albedo drops → the ice-albedo feedback weakens → harder to
  snowball — a real, modest, citeable effect; full spectral radiative transfer is the
  rung-4 deferral); and **planet size**, which in a rung-0 EBM **leaves the 0-D
  global-mean `T` untouched** (that is `S₀, α, A, B` only) and enters *only the
  transport* — a larger planet transports heat less effectively per unit area, so it
  **sharpens the equator-to-pole gradient**. Size's richer effects route through
  **rotation** (the real circulation lever — Rossby radius / β-plane), which lives in the
  Phase-3 fluid engine, not the EBM. So both deepen up the staircase rather than being
  faked at rung 0.

### 9.2 The interaction model is tier-dependent; the renderer is not (ADR 0004 #2)

Because the map only consumes arrays, it is **invariant up the whole §5 staircase** —
only the *trigger* changes: **live slider → instant remap** (rung 0, compute is
laptop-seconds, no special engineering — ADR 0001 scope) **→ set parameters → launch a
run → view the result** (the heavy upper rungs, where a live loop is impossible). The
slider is a *driver of compute*, the map a *consumer of arrays*; only the former is
tier-dependent. This is recorded as an explicit staircase consequence (§5) so a future
session does not try to keep the live loop alive past where compute allows.

### 9.3 Editable planet & state interchange — preplanned, not built (D1/D2 + ADR 0004 #3–4)

The user's forward requirements (2026-06-09): the planet should be **editable** —
elevation/depth, water/surface areas — and the map should **export/import** layer-by-layer
in a simple format, so a *future* custom map-editing app can author a world the current
app imports and runs models on. Both are **preplanned in the architecture, not built in
v1** — satisfied by a written contract + a non-foreclosing seam, the same discipline as
the rest of this plan. Two seams, kept distinct:

- **Renderer-input seam (built 2-D-ready now).** The map consumes a **2-D lat×lon field**;
  v1, being zonal-mean (§3 scope edge), simply **broadcasts `T(φ)` across longitude** — so
  the globe paints **latitude bands**, honestly. This is the natural shape for a globe,
  not premature 2-D.
- **Geography-physics seam (design-only now).** A documented **"geography spec"** —
  elevation grid + bathymetry + land/ocean mask, as plain *input* arrays — that the
  compute layer *will* consume. v1 writes the contract and names the consuming rungs; it
  **builds no machinery**. **Honesty flag (the user is aware):** an imported/edited
  geography is **inert** at v1 — carried, displayed, round-tripped, but it **does not
  change the climate** until the climb below. The round-trip guarantee is therefore
  *array identity*, not a changed climate.

**The "editable planet" is a climb on the §5 staircase, not a v1 feature** — with a real
cost gradient (mirrored into the §5 table):

| Tier | What it edits | What responds | Cost |
|---|---|---|---|
| **cheap — stays 1-D, annual-mean** | elevation `z`; land/ocean fraction per band | elevation → a **lapse-rate** map diagnostic `T_surf = T_band − Γ·z`; land/ocean → an **albedo difference** into the existing per-band albedo; fraction-per-band = continentality-lite | no engine change |
| **needs seasonality** | ocean depth / **heat capacity** | thermal lag, continentality | **v1 cannot** — annual-mean drops `C` at equilibrium; named, not built |
| **needs 2-D** (the north star) | true longitudinal geography | **regional climate, orographic precip, rain shadows** | new transport — **leaves the 1-D engine** |

**State interchange (ADR 0004 #3–4): pin a schema, not a format.** A documented
**planet-spec** schema — grid geometry, **explicit units** (self-describing), the
**layer list** (the §9.1 registry *is* the export manifest — one structure, not two), the
knob values, and a **`schema_version`**. Encoding per consumer, behind the schema:
v1 lean = a **JSON manifest + `.npz`** arrays; **editable-geography heightmaps**
(elevation/bathymetry/mask) interchange as **16-bit grayscale PNG** (the native currency
of paint tools / web canvases — 8-bit is too coarse; this is the round-trip a future
editor needs); **NetCDF** is a documented future encoding for climate-tool interop
(deliberately *not* the v1 choice — NetCDF is browser-hostile, and the editor is the
consumer). **Round-trip identity (`import(export(s)) == s`) is a genuine correctness
property and gets a genuine test**, unlike the map's smoke-tests.

**Inert-honesty has two forms — geometry and disclosure.** Through v1 the discipline is **inertness**:
an imported geography is *carried but not consumed*, the eddy globe is *one honest band, not a wrap*. The
showcase flow renderer (§9.5, Rung C, decided 2026-06-12) introduces the second form, **honest-by-disclosure**:
for the *showcase layer only*, an illustration may depart from what the model literally computes —
a global-looking flow, "currents carrying heat" — **so long as a visible on-screen disclaimer documents the
departure**, and a test machine-checks that the disclaimer is present (the documentation is verified even
when the physics is not). The carve-out is narrow: it never touches the science layer or the
honest-by-construction A/B rungs. Its data hook is this same manifest — a **vector-field-on-a-globe layer
type** (grid + `(u,v)` + scalar + frames + **coverage-extent** + **provenance/honesty label**) joins the
§9.1 layer registry, the coverage-extent carrying the band-vs-globe truth into the data itself.

### 9.4 Toolkit promotion

Plot primitives start planet-local (`plots.py` static floor; the `planetmap.py` layer
registry for the deep end). The **2-D field / heatmap** and **time-animation** primitives
are the ADR 0002 §3 candidates whose third reuse (after steel/chip) would **promote to
the shared `viz/`** by rule-of-three (ARCHITECTURE.md §6); the layer registry is that
third consumer-in-waiting. Promotion is **not** done pre-emptively (the existing three
`plots.py` share conventions, not copy-pasted code — the thin-extraction finding,
2026-06-09). **Building visualization rung A (§9.5) is what finally trips this trigger** —
the eddy life-cycle animation is the time-animation primitive's first real consumer.

> **Post-split status (R2, 2026-06-14).** Two things this paragraph assumed are no longer
> live, and one promotion has now landed:
> - **No cross-repo `viz/`.** When the monorepo split into standalone repos (2026-06-10),
>   `ARCHITECTURE.md` and the program-shared `viz/` stayed in the archive; *steel* / *chip*
>   are separate repos now. So "promote to the **shared `viz/`** across steel/chip/planet by
>   rule-of-three (ARCHITECTURE.md §6)" is a **dangling pre-split reference** — the relevant
>   scope is now planet-sim-**internal** (a future `planet/viz/`, if ever), and the third
>   consumer must be an in-repo one.
> - **The serialization/interchange machinery did meet rule-of-three — and is now documented
>   (R1, 2026-06-14).** `planet_spec` now serves a **third** consumer-class: the biome-map
>   export, the two-world diff, **and** the vector-field interchange (`flow_serialize`, R1 /
>   §11.2 — a producer-agnostic `VECTOR_OVERLAY` layer carrying coverage-extent + provenance,
>   round-trip-pinned on two producers). That contract — *the* §9.3 schema + the §9.1 layer
>   registry — **is** the "shared, documented contract" R2 names; it is already a clean shared
>   module, so the promotion is **documentation, not extraction** (this note + the §11.2 R1/R2
>   record).
> - **The globe-geometry viz *helpers* are still two-consumer — deliberately not promoted.**
>   `_sphere_xyz` (planetmap → eddy_globe) and `_band_geometry` / `_earth_radius` (eddy_globe
>   → flow_globe) each have **exactly two** consumers (R1's `flow_serialize` renders through
>   `planetmap.render`, so it added a serialization consumer, **not** a geometry-helper one).
>   They are clean single-source imports, not copy-paste, so extracting them to a `planet/viz/`
>   now would be the **pre-emptive promotion this very paragraph forbids** — named here, held,
>   to re-trip when a genuine third geometry consumer appears.

### 9.5 Animated flow — the visualization rungs (decided 2026-06-11; **rungs A+B+C BUILT** — A 2026-06-11, B 2026-06-12, C 2026-06-13)

The emergent eddy life cycle (`eddy_flux.eddy_life_cycle`) is the only genuinely
time-varying, longitudinally-structured 2-D flow the project produces — the instability
grows, saturates, and stirs the tracer. The user's forward requirement (2026-06-11,
referencing NASA's *Perpetual Ocean* and Ventusky as *broad* visual references — weather,
not climate) is to **animate** it. Decision (user, 2026-06-11): **build all three
visualization rungs, A → B → C** — rising cost, *falling* pedagogical return, *rising*
overclaim risk, so each is a deliberate step, not an automatic one. (These are
**visualization** rungs A/B/C; do **not** conflate them with the §5 GCM staircase rungs
0–6.)

**The shared prerequisite (all three): bank the flow frames.** `eddy_life_cycle` currently
discards the 2-D state each step, keeping only the scalar diagnostics + the κ life-cycle
integrals. The one data change every rung needs is an **opt-in `n_frames`** that snapshots
`(h, u, v, θ)` into an `EddyFrames` side-channel — `h` included (the verification anchors
need it), on **even *time* thresholds over the full release `[0, t_end]`** (the step `dt`
is adaptive, and the full span — not the κ window — is where growth→saturation, *the
mechanism*, lives). It is **diagnostic-pure**: `n_frames=0` leaves the κ result bit-for-bit
unchanged — the inert-seam discipline (§9.3) applied to motion.

| Viz rung | Renderer | Tech | What it is |
|---|---|---|---|
| **A** | matplotlib `FuncAnimation` | in-repo, `[viz]` | the **mechanism artifact** + the repo's first time-animation primitive (`plots.eddy_life_animation` + `demo_eddy_life.py`) |
| **B** ✅ | Plotly globe animation | existing `[webviz]` stack | the **globe view** (frames → play/slider on `planetmap`'s sphere), no new tech — **BUILT 2026-06-12** |
| **C** ✅ | three.js / WebGL perspective globe | self-contained inline HTML (three.js vendored, no CDN) | the **showcase** — a general flow-on-a-globe **particle-streaming** renderer; *honest-by-disclosure* — **BUILT 2026-06-13** (decided 2026-06-12, amended same day to a true 3-D sphere, below) |

- **Rung A is meaningful on its own** (not merely a step to C): it validates the banked
  frames, builds the §9.4 primitive, and **teaches the mechanism honestly**. It is a
  **two-panel** artifact — the θ field stirring *beside* a running cumulative `∫F̄dt` vs the
  `|F̄|` throughput — so that `eddy_flux`'s headline finding, that the instantaneous flux is
  **~90 % reversible** (`irreversible_fraction ~0.1`), is made *visible*: the swirls rage
  while the net stays near zero, then settles to the small down-gradient residual. Without
  the second panel a stirring movie *contradicts* the module's own finding (the overclaim
  this repo polices). "Visualize the mechanism, not the output" (ADR 0002 §5). GIF via
  Pillow is the CI-safe default; MP4 via ffmpeg is optional, guarded like the `[viz]`
  figures. Frame-fidelity tests: `∫hθ` machine-exact across banked frames, `eddy_ke`
  recomputed-from-a-frame matches the series, `n_frames=0`-vs-`N` bit-for-bit.
- **Rung B** reuses the existing Plotly globe — a globe view for a fraction of C's cost. **(BUILT — see below.)**
- **Rung C** is the **showcase** — build approach **decided (2026-06-12), amended the same day to a true
  3-D sphere (below):** a **three.js / WebGL perspective globe** with particles streaming on a *real*,
  rotatable sphere, architected as a general flow-on-a-globe renderer and governed by the
  *honest-by-disclosure* carve-out. **B and C are now both 3-D globes, so their roles — not their
  dimensionality — separate them:** rung B paints the one true band as a faithful **scalar field**
  (honest-by-construction); rung C is the **immersive particle-streaming** view (illustrative, global-
  looking, honest-by-disclosure). It is **reach / delivery, not new teaching**. (Two earlier renderer
  ideas are superseded: first, vendoring **mapbox/webgl-wind** + **cambecc/earth**; then an *original
  Canvas2D orthographic* globe — the latter now demotes to a lighter considered-alternative / fallback,
  see the Tech bullet.)

**Named scope edges — carried through all three, *hardest to preserve at C*.** (1) The flow
is a **doubly-periodic midlatitude β-plane band patch, not a global field** (the same edge as
`circulation_layer`): on a globe it is one honest latitude band, not a planet-wide
circulation. (2) The flux is **~90 % reversible**: particles stream, but the streaming
mostly *sloshes* — genuine net transport is only the small κ residual. The prettier the
rung, the louder the medium itself whispers "global ocean currents carrying heat" — exactly
the two things the model lacks — so **B keeps a flux indicator and a band label**
(the §9.3 inert-honesty discipline, applied to motion — *honest-by-construction*), and **C, the
showcase, keeps the truth in a documented on-screen disclaimer instead** (*honest-by-disclosure*,
decided 2026-06-12 — see the Rung C subsection below). This is *why* the recommended order
is A first (honest by construction), then judge B-vs-C *after seeing real frames move*.

**Verification answers the user's own screenshot concern.** Judging animation from
screenshots of a *rotating* globe is unreliable; rung A is a **fixed-camera flat field**
(frames directly comparable), and the real proof is **numerical** — `∫hθ` machine-exact, θ
bounded, `eddy_ke` grows→saturates at `saturation_period`, and the cumulative-flux trace
lands on the diagnosed `kappa_bulk` — not the eye. The banked `(h,u,v,θ)` frames are exactly
what rung C consumes, so the cheap rung A de-risks the data *before* the WebGL investment.

**Built — rung A + the shared frame prerequisite (2026-06-11).** The `n_frames` side-channel is built
on `eddy_flux.eddy_life_cycle`: an opt-in `EddyFrames` snapshotting `(h, u, v, θ)` (plus the cumulative
transport traces) at even *time* thresholds over the full release `[0, t_end]`. It is **diagnostic-pure
— `n_frames=0` is bit-for-bit** (a test asserts `==` on `kappa_bulk`/`F_int`/`G_int`/`jet_speed`/the
`eddy_ke` series), the inert-seam discipline (§9.3) applied to motion. The two-panel animation
(`plots.eddy_life_animation` → `demo_eddy_life.py`, banked `docs/figures/planet-eddy-life.gif`, GIF via
Pillow) is the program's **first time-animation primitive**: left = the θ field stirred by the released
eddies with the eddy-velocity `(u−ū, v−v̄)` overlaid (fixed colour range + fixed quiver scale so the eye
reads the *growth*, not autoscale); right = the cumulative meridional transport — the **throughput**
`Σ∫|F̄|dt` raging upward while the **net** `Σ|∫F̄dt|` stays a small fraction, so `eddy_flux`'s
~90 %-reversible finding is made *visible* (a bare stirring movie would silently overclaim "ocean
currents carrying heat" — the two things the channel lacks). **The advisor caught one design error:** the
net trace must integrate the flux **per latitude first** (then `Σ|·|`), or *spatial* cancellation across
the band leaks into a curve whose only job is *temporal* reversibility — making the flux look *more*
reversible than it is; integrating-per-latitude-first makes the endpoint ratio land on the diagnosed
`irreversible_fraction` by construction, and a marked `window_start` line reconciles the full-release
curve with the banked windowed number. Frame-fidelity tests: `∫hθ` machine-exact across frames; `eddy_ke`
recomputed *from a banked frame* reproduces the series exactly; the render is an `importorskip`-gated
execution smoke-test. **C remains pending** — judge it after seeing the globe (rung B) land.

**Built — rung B, the Plotly-globe animation (2026-06-12).** The *same* banked frames lifted onto the
existing `planetmap` sphere — *no new stack* (`planet/eddy_globe.py`: `eddy_globe_figure` +
`save_eddy_globe_html`; `demo_eddy_globe.py` banks `docs/figures/planet-eddy-globe.html`, ~3.8 MB on par
with the comparison-triptych globes). Plotly-free at import (Plotly lazy, like `planetmap`; a headless
subprocess guard asserts no plotly/matplotlib/ipywidgets is pulled), reusing `planetmap._sphere_xyz`.
**Both honesty edges are carried geometrically, not just captioned:** (1) the doubly-periodic channel is
laid as a **bounded ~55° longitude sector at its true physical width** — `Δlon = Lx/(a·cosφ_c)`, with `a`
recovered from the frames' own linear β-plane `(y, φ)` metric (no new constant/schema field) — a *single*
NH band (measured 18.9°–61.1° lat × 55.1° lon) on an otherwise-bare base sphere, **not** a 360° wrap and
**not** mirrored to the SH (the `circulation_layer` two-hemisphere broadcast is valid only for a
*zonal-mean* jet; the eddy field is the project's only longitudinally-structured field, so it does not
transfer). (2) The Rung-A flux-budget panel rides beside the globe (throughput-rages / net-stays-small),
so the ~90 %-reversible finding stays on screen; a moving cursor + `κ diagnosed` line tie the two panels.
**The advisor's gate** (before push): a figure that *builds* is not a figure that *animates* — the frames
update only the band's `surfacecolor` (lean HTML; `x/y/z` static-merged from the base trace), which is a
known Plotly-3D soft spot, so (a) `cmin/cmax` are re-stated per frame to forbid colour-autoscale, (b)
`redraw=True` on play + every slider step forces the gl3d repaint, (c) a structural test pins each
frame's `surfacecolor` to its θ snapshot + the `traces` indices + that the band changes across frames.
The one thing not self-verifiable headlessly (no kaleido/browser here) is the actual play-through, handed
to the user to eyeball. Geometry pin (always-green): the band is a bounded midlat sector, never
pole-to-pole / 360°.

**Rung C — the showcase: decided 2026-06-12; amended the same day to a true 3-D sphere (three.js / WebGL
perspective; *honest-by-disclosure*); BUILT 2026-06-13 (build record at the end of this subsection).**
The build approach was locked, then built as locked. Rung C is what §9.5 already calls
it — *reach / delivery, not new teaching* — and the user's forward framing (2026-06-12) widens it into a
**general-purpose flow-on-a-globe renderer**: its design aim is to *one day* animate a full **GCM / ESM**
wind-or-current field; until (and after) then it renders whatever lesser model we have — today, the one eddy
band. The renderer is therefore architected around the *data*, not the eddy. **The amendment touches only
the renderer.** The renderer-agnostic data contract (below) was built to *"commit to nothing about
projection or particle representation"* — so swapping the originally-planned Canvas2D-orthographic globe for
a three.js perspective sphere changes nothing in the contract, the honest-by-disclosure carve-out, §9.3, or
the ADR 0002 note. That the swap is *local* is the proof the boundary was drawn in the right place.

- **Tech — a three.js / WebGL perspective globe (three.js vendored inline).** Particles stream on a
  *real* 3-D sphere — three.js `PerspectiveCamera` + an orbit control (rotate / zoom / tilt) with correct
  back-face occlusion — fed by the banked `(u, v)` frames: the immersive *Perpetual-Ocean* look the user
  is after, now on an actual planet rather than an orthographic disk. The **particle technique stays
  original** (our own advection updating a three.js `BufferGeometry`, CPU-side to start); three.js itself
  is the only vendored piece — the WebGL scene/camera framework we build *on*, not a particle library we
  copy.
   - **§6 reverses from *sidestepped* to *owed* — the real cost of this choice.** The prior
     original-Canvas2D plan existed precisely to avoid vendoring; three.js is a **library** (you vendor it,
     you do not reimplement it), and its MIT licence **requires attribution** — so a `NOTICE` file
     (currently absent in the repo) becomes a **named build deliverable**, not an afterthought. Forgetting
     attribution in a published artifact is exactly the kind of thing that bites later.
   - **Vendor three.js *inline*, not via CDN — decided.** Inlining the minified three.js into the emitted
     HTML preserves the `interactive.py` property the repo relies on: **self-contained, works straight off
     `file://`, deterministic, golden-able, shareable offline**. A CDN `<script src>` is the **rejected**
     alternative — lighter to write, but it breaks offline self-containment (the file would silently need
     the network), the one property that pattern exists to guarantee. The artifact grows by three.js'
     bundle weight; that is the accepted price of the real 3-D sphere.
   - **The original Canvas2D orthographic globe demotes to a lighter considered-alternative / fallback** —
     a no-dependency, §6-free renderer kept in reserve (e.g. if a future build wants a zero-vendor floor),
     no longer the v1.

- **The renderer-agnostic data contract — the "one day a GCM" hook.** The renderer consumes a generic
  **vector-field-on-a-globe** layer, *not* the eddy band specifically. The contract carries only: a
  lat×lon **grid**, per-cell `(u, v)` components, an optional **scalar field** for colour, optional **time
  frames**, a **coverage-extent**, and a **provenance / honesty label** (the disclaimer text). It commits
  to **nothing** about projection or particle representation — those stay renderer-side, which is exactly
  why *this very amendment* (original-Canvas2D → three.js / WebGL sphere) touches it not at all, and why a
  later GPU-advection swap (below) will reuse it unchanged. The **coverage-extent carries the band-vs-globe truth into
  the data itself**: the same renderer *labels the band* today (coverage = the ~55° NH sector measured in
  rung B) and *illustrates the globe with disclosure* tomorrow (coverage = global, fed by GCM
  winds/currents). This is a new **layer type** for the §9.1 registry / planet-spec schema (§9.3) — a
  vector field joining the existing scalar layers — so the showcase renderer and the interactive map share
  one manifest.

- **Honest-by-disclosure — the scoped honesty carve-out (the one real policy change).** Rungs A and B are
  **honest-by-construction**: the geometry cannot lie — one band, the flux-budget panel on screen. Rung C,
  the *showcase*, is **honest-by-disclosure**: it **may** render a global-looking, continuously-streaming
  field the model does not literally produce — a band extended toward a globe; particles that *imply* net
  currents though the instantaneous flux is **~90 % reversible** — **provided the departure is documented
  visibly in the artifact itself.** The user's condition (2026-06-12) is explicit — illustrate freely *"if
  it['s] documented … currents carry heat, when they do not."* Two consequences, deliberately asymmetric:
   - **Physics-fidelity verification relaxes.** Approximate is fine for the showcase; **no byte-golden, no
     numerical transport proof.** The model/science layer keeps *all* of that discipline (the validation
     triad, the scope edges, the inert seams) — the figure was never in the correctness path (ADR 0002 #2:
     *a figure is never evidence of validity*), which is exactly what *licenses* the illustration.
   - **Documentation verification tightens.** Because the on-screen disclaimer **is the entire license**,
     it is the one thing **machine-checked**: a structural test asserts the artifact carries the honesty
     caption, **on-screen and legible to a casual viewer** — not a code comment, not buried in this plan.
     **Principle: the documentation is machine-checked even when the physics is not.**
   The carve-out is **narrow — the Rung-C showcase renderer only.** A/B's honest-by-construction guarantees
   and the two named scope edges (as they bind A/B) are untouched; the science layer's overclaim-policing
   is untouched.

- **What the on-screen disclaimer must say (minimum).** That the flow is **illustrative**; **which**
  model/coverage it depicts (today: one midlatitude β-plane band; when fed a richer model: that model,
  named); that the streaming **mostly sloshes** (genuine net transport is only the small κ residual); and
  that *"currents carrying heat"* is an artistic reading the numbers do **not** validate. This is the §9.3
  inert-honesty discipline re-expressed as **disclosure** (a caption) rather than **geometry** (a literal
  band).

- **GPU-advection upgrade seam + a concrete trigger (the seam shifted down one level).** WebGL is no
  longer the *future* upgrade — three.js makes it the v1 renderer — so the named seam moves *inside* it:
  **CPU-side particle advection** (updating the `BufferGeometry` each frame) is v1; a **GPU ping-pong**
  integrator (the *webgl-wind* technique, run in a fragment shader) is the documented upgrade **behind the
  same data contract** — swap the integrator, keep the contract and the disclaimer test. **Trigger (named,
  not "someday"):** when grid resolution × particle count pushes the interactive frame-rate below ~30 fps —
  the GCM-resolution regime — move advection onto the GPU. Until then, CPU advection on the three.js sphere
  stays responsive, self-contained (inline three.js, no CDN), and golden-friendly. **Built 2026-06-13
  (follow-up, user-requested ahead of the trigger): the GPU ping-pong integrator now ships as the default
  with the CPU loop demoted to a runtime fallback — see the "GPU ping-pong advection" build paragraph in §9.5.**

- **Deliverables for the build session (not executed this session).** A generic renderer module
  (`planet/flow_globe.py` — named for the role, the eddy being its first consumer) emitting a three.js
  scene with **three.js vendored inline**; a demo (`planet/demo_eddy_particles.py`) banking a
  self-contained `docs/figures/planet-eddy-particles.html`; a `catalog.py` `DEMOS` entry + the
  drift-guarded `python -m planet site` regenerate; a **`NOTICE` / attribution file** carrying three.js'
  MIT licence (the §6 deliverable this renderer now *owes* — see the Tech bullet); and **structural +
  disclaimer-presence tests** (mirroring `test_eddy_globe.py` *minus* the byte-golden, *plus* the caption
  assertion). The disclaimer is now a **DOM overlay element** over the WebGL canvas — *easier* to
  machine-check in the HTML source than canvas-drawn text would have been — so the disclaimer-presence
  test asserts that overlay's honesty text is present. Browser play-through is the one thing handed to the
  user to eyeball.

**Built — rung C, the three.js particle flow-globe (2026-06-13).** Built exactly as locked. `planet/flow_globe.py`
is the **generic** renderer: a renderer-agnostic `FlowField` contract (lat×lon grid + per-cell `(u,v)` +
optional `scalar` + `Coverage` extent + `honesty` disclaimer string) + `flow_field_from_eddy` (its first
consumer) + `flow_globe_html`/`save_flow_globe_html`. It is **NumPy-only at import** (builds an HTML *string*;
the `eddy_globe` import is local to the builder), mirroring the headless discipline. `demo_eddy_particles.py`
banks `docs/figures/planet-eddy-particles.html` (~758 KB — three.js inlined) reusing `demo_eddy_life.compute`
(one shared life cycle, three views); a `catalog.py` entry (`extras=()` — generation needs only core, the
artifact is self-contained) auto-adds it to the menu/CLI/landing page, the golden site test enforcing
regeneration. **Decisions, advisor-blessed at the done-check:** (1) **three.js r137 UMD vendored inline**
(`planet/vendor/three.min.js`, global `THREE`, plain `<script>`) — *not* ESM, because ES-module imports are
blocked over `file://` (CORS), which would break the self-contained-off-disk property; the orbit camera +
particle advection are **hand-rolled** (original), so three.js core is the only vendored piece. The §6
deliverable it owes — a repo `NOTICE` with three.js' **full MIT body** (not just the SPDX tag) + its inlined
`@license` banner — ships, so attribution travels with both repo and artifact. (2) **Band-confined coverage**
(seed particles only within the true ~55° NH sector recovered from Rung-B's shared `_band_geometry`/
`_earth_radius` — extracted so B and C can't drift) — fabricating a global `(u,v)` from a 55° patch would be
*inventing* data, not illustrating a richer model. Band-confinement does **not** downgrade C to
honest-by-construction: streaming particles still imply persistent currents the ~90 %-reversible flux does not
produce, so the disclaimer is mandatory and carries the *mostly-sloshes / net-is-the-small-κ-residual* clause.
(3) **Steady stream from the saturated frame** (nearest `saturation_period`) — the *Perpetual-Ocean* look, a
lean one-field payload. **Verification, per the carve-out** (`planet/tests/test_flow_globe.py`, 6 fast + 1
slow): physics-fidelity **relaxes** (no byte-golden, no transport proof — a figure is never in the correctness
path); **documentation verification tightens** — the disclaimer is a *visible* DOM element (`<div
class="disclaimer">`), and the machine-checked test asserts it carries **both** honesty clauses and is never
hidden (`display:none`/`visibility:hidden`/`opacity:0`). Plus the always-green guards: headless-import,
coverage-is-a-bounded-midlat-band-not-global, self-contained-with-three.js-inlined (no external `src=`). The
GPU ping-pong advection seam stayed **named-not-built** at this initial build (CPU `BufferGeometry` v1 is
responsive at this resolution; trigger = <~30 fps at GCM scale) — **now built 2026-06-13, see the "GPU
ping-pong advection" paragraph below.** **The one thing not headlessly self-verifiable** (no WebGL here)
is the actual browser play-through — handed to the user to eyeball (acceptance: particles stream along one
tilted band, rest of globe bare, disclaimer legible, drag/zoom responsive); same hand-off Rung B took. Gate:
303 fast-lane tests pass.

**Polish pass (2026-06-13, after the user's browser play-through).** Three cheap appearance wins, advisor-blessed:
(1) **Dark-spawn/death fix.** The fade lived in *RGB* (colour multiplied toward black over the first/last fraction
of a particle's life), so fresh and dying particles read as distracting **dark dots** against the dark globe. Moved
the fade into a **4th alpha channel** (vertex colour → RGBA, `USE_COLOR_ALPHA` confirmed in the vendored r137):
RGB now stays full-brightness `cmap(t)` *always*, alpha = `fadeIn × fadeOut` (the `+0.15` death-floor dropped too,
so particles fade fully to transparent at both ends). (2) **Round soft particles** via a generated radial-gradient
`CanvasTexture` as `material.map` — square GL points were the amateur tell; the white sprite preserves the
per-vertex temperature colour and rounds the dot (the single biggest *showcase* upgrade, answering the "shape"
ask as a better *default*, not a knob). (3) **Three live sliders** (particle size + opacity + edge sharpness) —
the first two are cheap because `material.size`/`.opacity` are live-mutable; **edge sharpness** (added the same day
after the user found the round sprite *too* defocused) rebuilds the cheap 64² `CanvasTexture` live and disposes the
old one — it sets the sprite's opaque-core radius (`0` = soft bloom, `1` = near-hard disc, capped at 0.97 for a 1-px
anti-aliased rim). All three initial positions flow from `particle_size`/`particle_opacity`/`particle_sharpness`
Python kwargs (one source feeding both the material init via `FLOW_DATA` and the slider `value=`, no drift), so a
notebook can ship different defaults while a viewer fine-tunes in-browser without regenerating. **Deferred as a
named seam:** the open-ended remainder — colour ramps, particle density presets, trail length, shape menus — is
speculative with no second consumer, so it stays named-not-built. *(2026-07-06: the second consumer is now
scoped — the §9.6 ocean producer; O3 names trail length + density as the first knobs to build.)* No test changes (the carve-out keeps the
disclaimer the only machine-checked thing; the 6 structural tests pass untouched); artifact re-banked, gate still
303 fast-lane pass.

**GPU ping-pong advection (2026-06-13, user-requested ahead of the trigger; advisor green-lit + done-checked).**
The named seam closes: particle advection now runs **entirely on the GPU** by default, with the original CPU
`step()` loop demoted to a **runtime fallback**. State lives in an **RGBA32F float texture**, one texel per
particle = `(lon, lat, age, life)`; each frame an off-screen fragment shader (`UPDATE_FS`) reads the current
state, advects by the *same* `dλ/dt = u/(a cosφ)`, `dφ/dt = v/a` metric and `accel` the CPU path used, and
writes the next state into a second target — the two **ping-pong**. A `Points` cloud then draws them, its vertex
shader (`DRAW_VS`) reading each particle's position straight from the state texture (sphere transform, point-size
attenuation, RdBu_r colour, alpha fade all moved into GLSL); sliders became uniforms. The velocity (+θ) field
rides along as a **half-float `DataTexture`** (`(u, v, θ, 0)`, linear-filtered — core in WebGL2; the state
texture is `Nearest`, so no float-linear extension is needed). **Hand-rolled** with core three.js
`WebGLRenderTarget`×2 + `RawShaderMaterial` — `GPUComputationRenderer` was *rejected* (an ESM addon, CORS-blocked
over `file://`, and a new vendored dep owing its own NOTICE). **No new vendored library → NOTICE untouched.**
The binding constraint is that **WebGL cannot run in CI**, so the design optimises for *a blind hand-off staying
recoverable* (the advisor's framing): (1) the CPU loop ships as a fallback so a GPU failure degrades to a working
globe, never a blank one; (2) the path is chosen at runtime by **feature-detect** (`isWebGL2` +
`EXT_color_buffer_float`) + **raw-compile-validating the GLSL against the live context** (three logs but does not
*throw* on a link failure) + a `try/catch` around init; (3) **console diagnostics** name the active path and the
fallback reason, and a **diagnostic read-back** logs particle 0's round-tripped state (the one residual gap the
gate can't catch — a driver that advertises the extension yet renders an *incomplete* float target → frozen
particles while the console still says "GPU active"). r137 landmines front-loaded: state texture `Nearest`,
`gl_PointSize` attenuation replicated by hand, `depthTest:true` against the opaque base (far-side occlusion),
`RGBAFormat` not `RGBFormat`, GLSL1 on `RawShaderMaterial`, `frustumCulled=false` on both the update quad and the
Points (texture-resident bounds). A 7th structural test pins that the artifact carries **both** the GPU shader
source and the CPU fallback, so a future edit can't silently gut the safety net. The data contract (`FlowField`),
`_build_data`, the disclaimer, and the carve-out are **all unchanged** — the swap touched only the renderer,
exactly as §9.3 predicted. Gate **304 fast-lane pass**; the browser play-through is again the one thing handed to
the user (acceptance unchanged; **"frozen particles" now has a second GPU cause** — an incomplete float target,
visible in the read-back log — and GPU vs CPU default point size may differ by up to the pixel-ratio, ≤2×, which
is slider-correctable, not a bug).

### 9.6 Beautiful ocean currents — the real-data showcase rungs O1–O5 (scoped 2026-07-06; **O1 + O2 + O3 + O4 + O5 all built + banked 2026-07-06 (O4 browser-verified; O5 = the QG emergent producer, the rule-of-three re-affirmed HOLD)**)

**The ask (user, 2026-07-06): "visualize beautiful ocean currents."** The project already owns both halves
of the machinery: the Rung-C particle globe (`flow_globe.py` — GPU-advected by default, honest-by-disclosure)
and the producer-agnostic serialization seam (R1 `flow_serialize.py`, whose own docstring names an
ECCO-shaped field as the target grid). What is missing is an **ocean producer** — every consumer today is
the ~55° emergent eddy band. These rungs feed a real global surface-current field through the pipe that was
built for exactly this, then make the render earn the word *beautiful* (the NASA *Perpetual-Ocean* look the
§9.5 requirement always referenced).

**The scope amendment this section makes (§11.4's living-staircase rule, exercised).** §11.2 declared
*"planet-sim stays atmosphere-only — it never ships an ocean visual,"* and §11.4 settled ECCO-ingest into
the spin-out as S1 **while explicitly keeping the alternative named** ("ECCO as planet-sim's last rung, to
de-risk against real data *before* the split"). **The named alternative is now taken — for the viz half
only** (user-ordered, 2026-07-06): planet-sim ships a real-data ocean **visual** (ingest → viz seam →
globe), because the renderer and the seam live *here* and the deliverable is wanted *now*. The **engine**
half of the boundary stands untouched: no ocean physics in this repo, ClimaOcean stays the spin-out, and
the forcing/input seam stays undesigned until the API is seen (§11.1's two-seams rule unbroken). Net
roadmap effect: **spin-out S1 narrows** — its data pipeline and Perpetual-Ocean visual are inherited from
O2 instead of rebuilt, leaving S1 as "ECCO as the *validation anchor* for S3" (§11.3 status note).

**Rung ladder — execution order O1 → O2 → O3 → O4; O5 independent/opportunistic.** Per the living-staircase
rule each rung is provisional until its predecessor lands; the `Retarget-when-done` notes say what moves.

- **O1 — the mask increment (the first genuine contract growth past R1).** `Coverage` is a lat/lon *box*;
  an ocean field needs a **per-cell validity mask** (land). `FlowField` gains `mask: Optional[(ny, nx)
  bool]` (`None` = all-valid = **bit-for-bit today** — the default-off discipline), `Coverage` stays as the
  bounding box. Particles are **seeded only in masked cells** — the honesty style of the R1 band-zeros
  embedding: *where-the-data-is carried in the data, not a caption*; land = no flow, honest-by-construction.
  Concrete renderer hook: the velocity `DataTexture` is already `(u, v, θ, 0)` — the free 4th channel
  carries the mask, so the GPU respawn/advection logic rejects land texels with **no new texture** (CPU
  fallback mirrors it in `step()`). Serialization: the mask rides as one more array on the `VECTOR_OVERLAY`
  layer's `.npz`; the round-trip `==` extends to cover it, and the no-mask spec stays loadable (schema is
  additive). *Retarget-when-done:* judge the mask's shape against the real ocean product at O2 (partial-ice
  cells, NaN-vs-fill conventions, staggered grids) before freezing it. **Retargeted 2026-07-06 by the O2
  spike (concrete now, not hypothetical):** OSCAR's land cells are `NaN` in `u`/`v`, and the existing
  `_bilinear` in `flow_serialize.py` runs plain `np.interp` — a `NaN` neighbor poisons interpolated ocean
  cells within roughly one source-cell of every coastline. So the mask can't just ride along passively:
  `flow_field_from_ocean` must (a) fill `NaN`→0 in `u`/`v` **before** resampling, (b) resample the boolean
  mask **separately with nearest-neighbor** (not bilinear — keeps it boolean, and nearest-neighbor doesn't
  smear the land/ocean boundary), then (c) re-apply the mask after resampling so filled-zero land cells
  read as "no data" rather than "zero current." Pole handling is the other retarget: OSCAR's grid is
  cell-centered (`-89.75…89.75`, 719 rows) and does **not** reach the true poles, while the interchange's
  `_globe_grid` is pole-inclusive (`-90…90` exactly) — for an `is_global=True` field the poles must come
  out masked (`mask=False`), not edge-clamp-extrapolated open ocean, or the render paints current at a
  point the source never measured. **BUILT 2026-07-06, all three spike retargets in.** Contract:
  `FlowField.mask` (`True` = valid cell; `None` = all-valid — every pre-O1 producer/spec byte-unchanged,
  the default-off discipline). Renderer (`flow_globe.py`): the mask rides the velocity `DataTexture`'s
  formerly-free 4th channel `(u, v, θ, mask)` — no new texture; the GPU update shader recycles a particle
  that drifts onto (or respawns on) a masked texel via the *invisible-retry* idiom (age set past life ⇒
  zero fade, re-rolled next frame — a fragment shader cannot loop a rejection sample; converges in ~2
  frames at OSCAR's 44% land), while the CPU fallback rejection-samples spawns (40 tries) and recycles in
  `step()`; both paths share the bilinear-sample-vs-0.5-threshold rule so they agree at the coastline to
  half a source cell. Serialization (`flow_serialize.py`): `_nearest_mask` resamples the mask
  **nearest-neighbour** (boolean-preserving, coastline-crisp), destination latitudes beyond the source's
  range come out `False` (the pole honesty; longitude edge-clamps instead — periodic axis), and the mask
  is **applied** to the embedded arrays (zero `(u, v)`, NaN scalar on invalid cells) so filled-zero land
  can never read back as measured zero current; it rides as one more *additive* categorical 0/1 `mask`
  layer in the same `.npz` (bonus: `render(active="mask")` paints an honest data-coverage globe through
  the unchanged renderer). Round-trip `==` extended to a masked, OSCAR-shaped probe (cell-centred lat
  stopping short of ±90, a continent strip inside a global box); the NaN→0 fill-before-resample stays
  **O2's producer job** (`flow_field_from_ocean` — an O1 probe's velocities are already finite). Tests:
  6 new in `test_flow_serialize.py` (§6 mask block) + 2 in `test_flow_globe.py`; full fast gate green.
- **O2 — the real-ocean producer (the deliverable, and the S1 de-risk).** Ingest one **OSCAR** global
  surface-current snapshot (0.25°, NASA PO.DAAC; ECCO acceptable if friction is lower) →
  `flow_field_from_ocean(...)`: full-globe grid (**`is_global=True` — the contract's first true-global
  consumer**), ocean mask, scalar = speed (or SST), provenance in `style`. Serialize through `planet_spec`
  (**round-trip identity on real data** — the R1 proof, now non-synthetic) and render through the
  **unchanged** Rung-C globe → banked `docs/figures/planet-ocean-currents.html`. **Honesty flips class:**
  the field is *real* but *not this model's output* — a new **provenance clause** under the
  honest-by-disclosure carve-out ("real reanalysis-class ocean currents [product, version]; **not** computed
  by this project's models"), machine-checked as a visible DOM element exactly like Rung C's reversibility
  clause. Data discipline: the raw netCDF is **never committed** (size); a small subsampled `.npz` fixture
  pins the tests; the netCDF reader is a demo-side **optional dep** (`flow_globe`/`flow_serialize` stay
  NumPy-only at import). **The one external unknown = data acquisition/auth** (PO.DAAC needs an Earthdata
  login) — **spike-first: download one file by hand before building anything.** **Spike DONE 2026-07-06**
  (`OSCAR_L4_OC_FINAL_V2.0`, granule `oscar_currents_final_20200601.nc`, one day/global/0.25°, 33 MB;
  fetched to the local temp workspace, never committed): **auth settled** — an EDL bearer token
  (`Authorization: Bearer <token>` header) against
  `archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-protected/...` works directly, no `.netrc`/URS
  redirect handshake needed, so the demo downloader is a one-header `curl`/`requests` call, token supplied
  via env var (never committed, never logged). **Format settled:** netCDF4/HDF5 (`h5netcdf`+`h5py` read it;
  `cftime` additionally needed to decode the **`julian`** calendar — `pandas` can't, relevant for O4's time
  axis); both `lat` (`-89.75…89.75`, 719 rows) and `lon` (`0…359.75`, 1440 cols, **0–360 convention, not
  ±180**) are already ascending — no north-to-south flip to handle, but the lon **rewrap**
  (`((lon+180)%360)-180`) makes the array non-monotonic and needs a re-sort/roll before
  `flow_field_from_ocean` hands it to `_bilinear`. **Dim order is `(time, longitude, latitude)`** — lon
  before lat, opposite of `FlowField`'s `(n_lat, n_lon)` — a transpose is required, not optional. Variables:
  `u`/`v` (total = geostrophic+Ekman, the fuller "beautiful" signal) and `ug`/`vg` (geostrophic-only,
  named as a future knob); units already m/s, range ±3 — no conversion needed. `NaN` = land/missing (44%
  of the global grid) — this is what makes O1's fill-before-resample/nearest-neighbor-mask/pole-masking
  retarget (above) concrete rather than speculative. **Forward flag for O3/render, not yet acted on:** the
  interchange's round-trip proof grid is coarse (2°, `GLOBE_N_LAT=91`/`GLOBE_N_LON=181`) but OSCAR is
  0.25° — confirm at O2/O3 build time that the *rendered* texture consumes native resolution (or a
  deliberately-chosen finer grid) and the 2° spec stays only the serialization round-trip proof, or
  "beautiful" dies at 2°; a native-res global texture also costs far more than the eddy showcase's 758 KB
  and that tradeoff needs a conscious call, not a default. *Retarget-when-done:* the
  real field's dims/mask/units retarget O1's mask and R1's schema (the §11.2 R1 note foresaw exactly this
  revisit), and the O3/O4 payloads are re-judged against what the product actually carries.
  **BUILT 2026-07-06.** Producer = `planet/ocean_currents.py` (NumPy-only at import): `load_oscar`
  (lazy `h5netcdf` — the new `[ocean]` extra; owns the file-layout knowledge: the `(time, lon, lat)`
  transpose, `_FillValue`→NaN, a `stride` subsampler, the `geostrophic_only` `ug`/`vg` knob) →
  `OceanSnapshot` (loader-normalized but convention-raw: 0–360 lon, NaN land, provenance strings read
  from the granule's own attrs) → `flow_field_from_ocean` (pure arrays: the ±180 rewrap with the
  `argsort` column re-sort that restores a monotone axis, mask-from-finiteness *then* NaN→0 fill — the
  O1 rule that the mask, never a filled zero, carries "no data" — scalar = speed, `is_global=True`
  with the coverage box honestly cell-centred at ±89.75°). **Both flagged conscious calls made:**
  (1) *resolution* — the banked render consumes **0.5°** (stride 2 off native 0.25°; the 2° interchange
  grid stays proof-only): boundary currents stay sharp, the self-contained HTML lands at **5.2 MB**
  (native would be ~20 MB for detail a 20 000-particle render can't resolve); (2) *pace* — the
  auto-accel's "fastest particle crosses the span in 6 s" was band tuning, so `flow_globe_html` grew an
  **additive default-off `crossing_seconds` knob** (eddy artifact byte-unchanged; the ocean demo passes
  45 s for 360°). The provenance clause ships in `field.honesty` ("REAL data … OSCAR L4 v2.0, PO.DAAC,
  DOI 10.5067/OSCAR-25F20 … **NOT computed by planet-sim's models**") and is machine-checked as a
  visible DOM element — including on the **committed artifact itself** (a regeneration can't silently
  drop it). Demo = `planet/demo_ocean_currents.py` (catalogued): bearer-token downloader
  (`EARTHDATA_TOKEN` env var, one `Authorization: Bearer` header, token never logged; raw granule lives
  in gitignored `outputs/`, never committed), asserts the **R1 round-trip identity on the real field**
  (the proof, now non-synthetic) before banking `docs/figures/planet-ocean-currents.html`. Tests =
  `test_ocean_currents.py` (17): real-data checks run off a committed **14 KB 5° fixture**
  (`planet/tests/fixtures/oscar_subsample.npz` — raw conventions kept on purpose so rewrap/mask/pole
  handling is exercised; OSCAR is freely-distributable, provenance in `fixtures/README.md`; NB
  `.gitignore`'s `data/` rule is why the dir is named `fixtures/`), the loader against a synthetic
  granule with OSCAR's exact layout (`importorskip`-gated). Renderer + contract byte-level untouched
  apart from the pacing knob — the §9.3 win recurring. Full fast gate 531 green. *O3 forward flags
  observed while here (renderer-side, not acted on):* particles recycle (not wrap) at the ±180° seam;
  uniform-in-lat spawning over-densifies high latitudes on a global field (cos-weighted/equal-area
  seeding belongs in the O3 beauty pass, alongside speed colouring — the shipped RdBu_r speed scalar
  reads, but a dedicated ramp + speed-weighted seeding is where "beautiful" gets earned).
- **O3 — the beauty pass (renderer-only; `FlowField` untouched — the recurring §9.3 win). BUILT 2026-07-06**
  (all three sub-parts + both §9.5 knobs; `flow_globe.py` only, three commits). The contract held: not one
  line of `FlowField` moved — the entire beauty pass is three.js-side, the §9.3 win recurring a *fourth*
  time. Advisor-shaped on three load-bearing calls: **(a)** the base layer is **honest-by-construction**,
  not a `CanvasTexture` — a base fragment shader inverts each surface point to `(lat, lon)` with the SAME
  `sph()` mapping the particles use (`atan(n.z, n.x)`) and samples the SAME mask on the SAME 0.5 coastline
  rule as `validAt()`, so the coast under the particles can never drift from the coast under the base; no
  mask (the eddy band) or a compile miss degrades to the plain solid sphere. **(b)** trails default-OFF
  behind a kwarg (no WebGL CI + blind hand-off ⇒ the ocean globe is the first to exercise them, the shipped
  eddy artifact can't silently regress); the **depth-only occluder prepass is load-bearing** — it kills
  back-hemisphere particles *before* they enter the accumulation buffer, so nothing bleeds through the
  planet at composite time; **additive** accumulation (One+One), not alpha-`over`, sidesteps
  premultiplied-alpha fringing and *is* the ocean glow; and the **rotation-smear fix** = `decay=0 while
  dragging` (a screen-space buffer smears when the projection moves, so history pauses during rotation and
  resumes when still — trails degrade to a clean fade mid-drag). New fullscreen shaders gated through
  `compileOK`, RT creation `try/catch`-wrapped, screen-sized targets realloc on resize — any miss drops to
  the plain single-pass render, never the CPU fallback (which stays fade-only). **(c)** the colour path was
  *already* producer-driven (eddy=θ, ocean=speed) so there was nothing to recolour; the genuinely new part
  is **speed-weighted seeding in the RESPAWN path** (weighting only the initial seed relaxes back to uniform
  as particles age out), composing with the mask reject via the same invisible-retry idiom, floored so calm
  water keeps an ambient fill — plus a **sequential** speed ramp as an opt-in default (RdBu_r is diverging
  and bleaches a 0→max field; RdBu_r stays default for signed θ, the ocean demo opts in). Both §9.5 knobs
  shipped: **density** (GPU rank-cut uniform / CPU tail-hide) and **trail length** (the decay). The ocean
  demo (`demo_ocean_currents.py`) opts into `colormap="speed"` + `trails=True`; the eddy artifact re-banked
  on the same renderer (trails off, no mask, θ ramp — the visible change is jet-core concentration from
  speed seeding). Verified end-to-end on the committed **5° OSCAR fixture** (masked pipeline, both artifacts
  node-`--check` clean); full fast gate **538 green**. *Honesty edge, named:* trails pause-to-fade during a
  drag is a deliberate property, not a bug — the smear fix, worth saying in the browser hand-off. The 0.5°
  banked ocean artifact re-bank needs an `EARTHDATA_TOKEN` (user hand-off; the code + fixture prove it).
  *Original scoping (kept for the record):* Three upgrades, all on the three.js side, each independently
  landable:
  - **(a) Land/ocean base layer.** A lat/lon two-tone (dark land / deep ocean) `CanvasTexture` generated
    from the O1 mask and draped on the sphere — currents are unreadable on a bare globe, and coastlines are
    what make the Gulf Stream *look* like the Gulf Stream. Reuse-not-invent: the mask is already in the
    payload; no new data source.
  - **(b) Motion trails.** The signature *Perpetual-Ocean* element: an accumulate-and-fade offscreen
    render-target pair (previous accumulation × ~0.96 decay + this frame's particles on top, then composited
    to screen) — a **third ping-pong** alongside the state textures, so it fits the architecture already in
    the file. The CPU fallback simply keeps today's fade-only look — trails degrade gracefully to a working
    globe, the same blind-hand-off discipline as the GPU advection build.
  - **(c) Speed styling.** Colour particles by `|u, v|` (the `speed` layer is already universal to every
    producer) + **speed-weighted seeding** so western-boundary currents (Gulf Stream, Kuroshio, Agulhas)
    visually dominate the way they physically do.
  - The §9.5 **control-surface seam unlocks**: the ocean producer is the named "second consumer," so
    **trail length + particle density** become the first two knobs (colour ramps / shape menus stay named).
- **O4 — frames: the time axis (seasonal currents). BUILT + VERIFIED + BANKED 2026-07-06.** The
  R1-deferred schema increment, all four pieces built and fast-gate-green (552 pass/1 skip), **default-off so
  every pre-O4 producer/consumer is bit-for-bit** and the single-snapshot shaders are byte-untouched.
  **Contract:** `FlowFrames` side-channel on `FlowField` (`u`/`v` `(nt, ny, nx)` + `labels`; optional per-frame
  `scalar` left `None`), `frames=None` = the exact pre-O4 path — the same additive-default-off discipline as
  O1's `mask`/O3's `trails`. **Producer** `flow_field_from_ocean_series` (`ocean_currents.py`): the same
  rewrap→mask→fill pipeline over N snapshots, stacked; **static mask = finite-in-EVERY-frame** (advisor:
  conservative — a cell measured only part of the year is left bare, never blinking; sea-ice seasonality
  folded into the "no data in any frame → bare" rule, not animated) applied across the whole stack;
  **acquisition-agnostic** (advisor: the producer never claims "climatology" — the caller passes the honest
  `period` phrase, which rides verbatim into the provenance clause). **Renderer** (`flow_globe.py`,
  GPU-path-only like trails): **separate two-texture crossfade shaders** `UPDATE_FS_F`/`DRAW_VS_F` (a `velAt()`
  = `mix(velA, velB, uMix)` substitution of the working shaders — the single path's `UPDATE_FS`/`DRAW_VS`
  untouched), N frame textures built once, a `stepSeason` tick advancing wall-clock `yearTime` through the
  frames with a **cyclic Dec→Jan wrap** `(k+1)%NT` (no hard cut), a `seconds_per_year` pace knob, and a
  **live month-label time badge** (the showpiece read-out). **Advisor payload calls, both taken:** colour is
  the **in-shader mixed speed** so the **per-frame scalar is dropped** (−⅓ payload), and the animation grid is
  **coarser than the O2 still** (demo `STRIDE=6` = 1.5°; motion hides resolution) → a 12-frame HTML ~5 MB,
  comparable to O2's single-frame 5.2 MB (12 frames of `(u,v)` at 0.5° would be ~30 MB). **Serializer**
  (`flow_serialize.py`): the R1 frames deferral now acted on — the `(nt, 2, n_lat, n_lon)` stack rides as ONE
  additive `FRAMES_LAYER` `VECTOR_OVERLAY` (each frame embedded on the globe grid exactly like the primary
  snapshot, labels in `style`), round-trip `==` extends for free; the interactive Plotly map **skips the 4-D
  stack** (a one-line `planetmap._overlay_traces` `ndim>=4` guard — it paints the primary snapshot; the
  flow-globe animates the stack). **Demo** `demo_ocean_seasonal.py` (bearer-token, 12 mid-month-day granules
  of 2020, gitignored `outputs/`, round-trip proof on the framed field → `docs/figures/…-seasonal.html`).
  **Frame-data choice (user, 2026-07-06): twelve monthly snapshots — ONE day per month of 2020, NOT means and
  NOT a climatology** (the Somali reversal is a large monsoon signal that reads clearly in a day-per-month
  series; the honesty label says exactly that). **The honesty gap (advisor, load-bearing):** `node --check`
  validates the JS template syntax but **not GLSL**, and there is no WebGL CI — so the *entire frames GPU path*
  (the crossfade shaders, `stepSeason`, the badge, trails+frames composition) **has compiled/run nowhere**; a
  frames-shader compile error degrades **silently** to the CPU fallback (static frame-0, which reads as "just
  not animating"). That gap is now **closed by the owed browser play-through (user, 2026-07-06, PASS):** the
  gold month badge cycles Jan→Dec and the particles stream — the GPU crossfade path compiled and ran, no CPU
  fallback. The 12-frame artifact is **banked** (`docs/figures/planet-ocean-currents-seasonal.html`, **4.23 MB**,
  ~5 MB as predicted) off a real 2020 OSCAR series (52% valid-in-every-frame ocean, |current| max 2.83 m/s,
  round-trip identity OK), the demo is now **catalogued** (`ocean_seasonal`, "Interactive globes") and on the
  landing page. *Retarget-when-done:* whether frames also back-port to the eddy band (Rung B/C animation parity)
  is still decided now-that-it's-banked, not promised — deferred, not owed.
- **O5 — the QG producer (independent bonus, does not gate O1–O4). BUILT + BANKED 2026-07-06.**
  `flow_field_from_qg` (`planet/flow_globe.py`, beside `flow_field_from_eddy`): the rung-3 two-layer QG
  condensate — coherent vortices + rolled-up PV filaments streaming in a box — as the **second emergent**
  producer. Cheap as scoped: `(u, v) = (−∂ψ/∂y, ∂ψ/∂x)` recovered by the model's own `invert`→`velocities`,
  and — advisor-confirmed — the axes already match the contract (`u` eastward / `v` northward; row=`y`→lat,
  col=`x`→lon), so **no transpose** (unlike OSCAR). Colour = the **upper-layer PV anomaly** `q₁` (the
  vortex-filament field the demo headlines), a *signed* scalar → the diverging RdBu_r ramp (like θ; the
  eddy's twin) — **nondimensionalised by f₀** (advisor caught this *pre-commit*): the raw QG PV anomaly is
  `O(1e-4 /s)`, which the renderer's 3-dp payload rounding (`_build_data`'s `flat(scalar, 3)`) collapses to
  a constant 0 ⟹ every particle one flat colour, *erasing* the vortex structure the producer exists to
  show; scaling by the positive constant `f₀` gives an `O(1)` Rossby-like field that survives the rounding
  (monotone, so RdBu_r still centres on 0). The fix lives in the **producer**, not the shared `_build_data`
  (touching its rounding would ripple to eddy/ocean). A payload-fidelity regression test pins that a
  realistic-magnitude state yields a non-flat rendered scalar (the class the earlier tests — which checked
  only the full-precision Python `FlowField` + grepped the disclaimer — missed). No mask (box coverage, no land), no frames (single snapshot) = the plain pre-O1 contract
  shape. **Two advisor content calls taken:** (1) the display latitude is **explicit `center_lat_deg=45°`,
  NOT derived from f₀/β** — the idealized `(f₀, β)` are independent numerical knobs, not a consistent
  `(sinφ, cosφ)` pair (the demo's f₀ implies ~43°, its β ~44°), so `atan(f₀/βa)` would *manufacture* a
  latitude never put in; the box maps by `Δlon=Δx/(a cosφ_c)`, `Δlat=Δy/a`, centred, its honest ~box-width
  sector, never 360°-wrapped. (2) the honesty string is **fresh, NOT the eddy's reversibility clause** —
  the whole rung-3 win is that this saturated flux is *irreversible* (a persistent turbulent field), so
  "~90 %-reversible / mostly-sloshes" is *false* here; it carries the QG edges instead (idealized model /
  not real data / box-not-planet-wide / inverse-cascade condensate) and names the colour (PV). A
  QG-specific honesty test pins the fresh clauses + asserts the eddy's are absent (the generic
  disclaimer machine-check runs only on the eddy field, so it does not misfire). **The architectural
  side-effect landed as documentation, re-affirming HOLD:** O5 is the third *geometry* consumer, and it
  confirms the two-consumer hold — `flow_field_from_qg` **cannot** call `_band_geometry` (that takes a
  frames object with `.phi/.x/.y` the metre-space box lacks) and never touches `_sphere_xyz`
  (renderer-side), so the one-line sector formula is *inlined*; extracting a shared helper would force the
  banked eddy path to recompute its latitude from `y` (ULP-risk on a banked artifact) — the pre-emptive
  promotion §9.4/R2 forbids. **The R1 round-trip is a named deliverable and passes** (the QG box is the
  plain no-mask/no-frames shape, so `vector_spec_from_flow_field`→`save`/`load ==` covers it for free; a
  `test_flow_serialize` case pins it as the third producer). Demo `demo_qg_particles.py` (catalogued
  `qg_particles`, "runs a short sim"; reuses `demo_baroclinic_qg`'s idealized params/box-sizing/dissipation
  but keeps the **model** so the renderer has the grid metric) runs the shear to a saturated `v'≫U_s`
  condensate (nx=96) and banks `docs/figures/planet-qg-particles.html`. Producer tests fast-lane (cheap
  `random_state`, no spin-up); the end-to-end bank is a slow smoke-test. Verification is **lighter than
  O4** (advisor): O5 reuses the already-browser-verified single-snapshot GPU path — no new shaders, so the
  "GLSL never runs in CI" hand-off gap O4 carried does not reapply.

**What stays out (the boundary, restated so this section can't scope-creep).** No ocean *engine* (ClimaOcean
= the spin-out); no forcing/input seam (S2 — designed at the boundary against the *seen* API); no real-time
what-if on real data (a reanalysis snapshot has no knobs — interactivity stays the model side's job). The
honesty ceiling carries: Earth reanalysis is the anchor; nothing here validates a custom world's ocean.

---

## 10. Immediate next step

**Where Planet stands (2026-06-09).** **All four phases are built — the capstone is complete** (the
build records for Phases 3–4 + the interactive map + the teaching notebook are below; the Phase-1/2
detail that follows is retained as the foundation narrative). *Phase 1* — the latitudinal EBM &
the Snowball bifurcation, the diffusion spine's third reuse with Strang-split radiation (banked:
present-day ~73° ice line, freeze at ~8 % dimming). *Phase 2* — the **climate→biome map** (the
payoff, banked early): `precip.py` (a diagnostic precipitation parameterization — a circulation-set
Gaussian latitude pattern in cm/yr × a Clausius–Clapeyron global-moisture amplitude) + `biomes.py`
(an **original, total** Whittaker `(T,P)→biome` partition — the Irvin precedent: the copyrighted
diagram reproduced by an independent calibrated computation, *not* digitized — with cold biomes
temperature-limited and warm biomes moisture-limited on diagonal thresholds). Banked
(`docs/figures/planet-biomes.png`): the biome-band map + the Whittaker (T,P) plane with the planet's
climate trajectory + the poleward migration under a CO₂ warming knob. No new engine; the 20-test
triad is green. The **interactive-map design is converged** and written up (ADR 0004 + §9 here): a
Plotly+ipywidgets **layer registry**, a tier-dependent interaction model, the **planet-spec**
export/import schema, and the editable-geography seam. **The interactive map is now BUILT** (2026-06-09):
`planetmap.py` (the layer registry — `LayerKind`/`Layer`/`Grid`/`PlanetView` + the biome-map builders +
the Plotly globe renderer + the ipywidgets live loop) and `planet_spec.py` (the pin-the-schema
JSON+`.npz` interchange with the round-trip-identity test). Its first version *is* the biome map; the
banked artifact is `docs/figures/planet-map.html`.

**Built — the interactive map (user decision, D3-B): `planetmap.py` v1 + `planet_spec.py`** (2026-06-09).
The map's first version *is* the **biome map** (the §9 centerpiece): a Plotly 3-D globe painted from
the **layer registry** (temperature / precipitation / biome scalar fields + the ice-line annotation +
an inert elevation seam), with knob-sliders **S₀ / CO₂→A / transport `D`** driving an instant
recompute-and-remap (rung-0 live loop, a pure consumer of `demo_biomes.compute` — no new physics). The
**planet-spec** export/import schema ships the v1 lean encoding (JSON manifest + `.npz`) with the
**round-trip-identity test** (the deep end's one real correctness property), and the geography-spec
seam is **written but inert** (the elevation layer is carried/round-tripped, not consumed) per §9.3.
Two non-obvious calls, advisor-blessed: **obliquity is a named, deferred slider** (wiring it needs the
`s₂(obliquity)` annual-mean-insolation relation pinned to a source — the same `[[…-source]]` discipline
as the exoplanet knobs) — **SINCE BUILT (2026-06-10): `obliquity.py` computes `s₂(ε)` from the pinned
daily-mean-insolation geometry (projected onto P₂, ratio-anchored so Earth recovers `−0.48` exactly;
validated by the exact `s₂(0)=−5/8` limit + the `≈−0.48` climlab cross-check, `[[obliquity-insolation-source]]`),
the slider is now live, banked `docs/figures/planet-obliquity.png`; obliquity is no longer deferred** —
and **`vector_overlay` is a declared-but-unpainted `LayerKind`** (the renderer
raises `NotImplementedError` naming Phase 4 — build the seam, not the machinery; the extensibility proof
is the inert elevation scalar the existing renderer paints for free). No new engine, **no gate-manifest
change**; opt-in behind a new `[webviz]` extra (Plotly + ipywidgets). 31-test pair added (round-trip +
registry always-green; render smoke-tests `importorskip` on Plotly — fast, not `slow`).

**Built — the teaching notebook (`planet.ipynb`, 2026-06-09)** then **Phase 3 (the shallow-water
engine, 2026-06-09).** *(Notebook extended 2026-06-10 to the full arc: §4 the winds
(`demo_shallowwater`) + §5 the emergent jet (`demo_coupler`) added as static banked renders, the §6
deep-end globe de-staled — obliquity + `vector_overlay` now-built, "what's next" repointed up the §5
staircase. Reach-not-evidence, so `nx=48` and the honest ~2 %-converges-to-~1 % Helmholtz framing;
no module/engine change. Pedagogy tiered 2026-06-10: §1–§5 each gained an expert `<details>` "Going
deeper" collapsible — the chip/steel convention planet alone lacked — making the three reading depths
explicit (narrative/sliders/collapsible = novice/intermediate/expert); markdown-only, banked outputs
byte-identical.)*

**Built — Phase 3: `engines/fluid` (the program's SECOND shared engine) + `circulation.py`** (2026-06-09).
A rotating shallow-water solver on a doubly-periodic β-plane — **Arakawa C-grid, vector-invariant
form, explicit SSP-RK3** — built standalone in `engines/fluid/` and **sealed behind its passing
suite in `engines/fluid/CONTRACT.md`** before any coupling (build-and-validate before reuse). It deliberately shares *no* machinery
with the parabolic-implicit `engines/diffusion` (a hyperbolic, CFL-limited wave solver). **Validation
triad:** gravity-wave `√(gH)` + **Poincaré dispersion `ω²=f₀²+gHk²`** to ~1e-3 (the rotation check),
**Rossby waves** westward & dispersive (loose band, *converging to analytic as the grid refines* — a
named numerical-dispersion edge), **geostrophic balance steady** + **geostrophic adjustment** to the
analytic Helmholtz state over `L_R` (~1%, the published Rossby benchmark), **mass to machine
precision**, **energy bounded & dt³-convergent**, and — the discriminating Coriolis seal — **potential
vorticity / enstrophy bounded at FINITE amplitude** (a vortex at Rossby ~0.5, where advection
genuinely moves PV around; the design review's key point — at small amplitude this leg is vacuous).
**Two honest design calls baked in:** (1) the symmetric vector-invariant scheme conserves **energy**
semi-discretely, *not* potential enstrophy (one Sadourny-class scheme conserves one or the other;
Arakawa–Lamb conserves both, not built) — so claims are stated *as measured*, not aspirational; (2) the
finite-amplitude PV seal came out crisp, so the engine stayed **nonlinear** (the fallback was a linear
core; not needed). The `tracer` slot is declared on `SWState` but **not advected** (rung 1 — "seam,
not machinery", `step` raises). `circulation.py` pins the planetary numbers (`f₀=2Ω sinφ`, `β=2Ω
cosφ/a`, equivalent depth → `L_R ≈ 960 km` at 45°, [[shallow-water-source]]) and banks the artifact
(`docs/figures/planet-shallowwater.png`: geostrophic adjustment + a westward Rossby wave +
conservation diagnostics). **Gate/infra consequences realized:** `planet`'s `uses` is now
`{engines/diffusion, engines/fluid}` — the manifest's **first genuinely multi-engine row** — and the
**import-drift guard** (deferred to engine #2) is **built and live** in `tools/tests/test_gate.py`.

**Built — Phase 4: the one-way EBM→circulation coupler (2026-06-09 — the capstone complete).**
`planet/coupler.py` (+ `demo_coupler.py`, `plots.coupler_figure`, the 9-test `test_coupler`
triad + 2-test `test_demo_coupler`) **couples the two shared engines**: the diffusion spine's
EBM equilibrium hands its meridional temperature gradient to the shallow-water engine, and a
**geostrophically-balanced midlatitude westerly jet emerges** (banked: ~16.5 m/s @ ~42°, core
geostrophic residual ~0.6%, `docs/figures/planet-coupler.png`). The forcing is composed *around* the
bare engine by **operator splitting — the third reuse of the EBM/Jominy idiom**: exact-exponential
**thermal relaxation** of `h` toward the EBM-derived target + weak **Rayleigh drag** (`τ_drag ≫ 1/f`,
near-geostrophic), *half-forcing / full bare-engine-step / half-forcing*. Three calls baked in:
**(1) y-periodicity** — the engine is doubly-periodic (walls are its named, unbuilt BC extension), so a
monotonic warm→cold target would jump at the seam and force a spurious boundary jet; the target is a
**windowed (Tukey), discretely zero-mean** height anomaly → C¹-periodic + mass-neutral, and the
periodic channel exacts a real **flanking easterly return** (∮u≈0). Its east–west–east *sign* banding
qualitatively resembles the general circulation, but the single-layer periodic channel does **not**
reproduce the observed westerly-dominant magnitudes — here the poleward easterly is actually the
strongest band, and its concentration is window-shaped, not observed — so this is a **named scope
edge**, not a faithful trade-wind/polar-easterly reconstruction. **(2) emergence, not
placement** — the channel brackets the *smooth* midlatitude baroclinic zone (excluding the ice-line
albedo cliff) and the jet lands at the **EBM gradient maximum** (~45°), *poleward of the channel
centre* (40°); a synthetic off-centre gradient makes the jet **follow it** (the decisive emergence
test). **(3) the conservation reframe** above (mass forced-exact + the release test). Non-circularity:
the jet *latitude* and geostrophic balance are amplitude-independent (validated); the height-per-Kelvin
`α` → jet *speed* is calibrated (tuning). Scope edge: **one-way, dry single layer** — no poleward heat
transport / reduction-to-EBM (rung 1, needs the tracer) and no thermal wind (rung 3, needs vertical
shear). The **interactive map now paints the `vector_overlay`** seam (`planetmap.circulation_layer` +
`_vector_overlay_trace` Plotly cones; `build_view(jet=…)`) — the deferred Phase-4 machinery built, the
jet registered as a `circulation` layer over the temperature field (`docs/figures/planet-coupler-map.html`),
round-tripping through `planet_spec`; it is **computed-then-viewed**, not in the live-slider loop (the
first compute too heavy for the rung-0 instant remap, §9.2). No engine modified; planet's `uses`
unchanged (`{diffusion, fluid}`); full planet gate **140 passed, 1 skipped**. Two-way coupling is rung 1
(seamed at the engine's `tracer` slot, not built).

**Started — rung 1 (the two-way coupler), step 1: the passive tracer-advection engine extension
(2026-06-10).** With engines now **living contracts** (ADR 0005 — the freeze ceremony dropped), the
long-declared `engines/fluid` `tracer` slot is **built**: `SWState.tracer` is advected in flux form
(`∂(hθ)/∂t = −∇·(hθ u)`) through the same SSP-RK3, strictly **passive** (no back-reaction on `h,u,v`),
with `tracer_mass` / `tracer_variance` diagnostics and `engines/fluid/tests/test_tracer.py`. Triad:
`∫hθ` machine-exact (telescopes like mass — the anchor), dry dynamics **bit-for-bit** unchanged (the
re-seal), uniform-flow translation (analytic), uniform-tracer consistency; variance is **bounded**
(a build-time honesty correction to the plan-table "dt-convergent" — measured round-off/spatially
limited, the enstrophy class, *not* dt-truncation limited) and the scheme is **not monotone** (no flux
limiter → over/undershoot on sharp fronts — the named scope edge). **Step-0 de-risking
(`outputs/rung1_stability_probe.py`, gitignored):** the existing Phase-4 jet is **barotropically
unstable** (Rayleigh–Kuo met; a v-perturbation grows ~200× then saturates), so a passive tracer on the
meandering flow gets an **emergent** `⟨v'θ'⟩` flux — no imposed stationary wave needed. **Anchor
reframe (for step 2):** the rung-1 anchor is **reduction-to-EBM** (resolved flux → `D_eff`; down-gradient
limit recovers the rung-0 diffusive EBM), **not** the plan-table "~5–6 PW" (an eddy/baroclinic = rung-3
number; here magnitude is window/forcing-tuned). Editing the shared engine triggered the
**full-repo gate** (bare `pytest`, ADR 0003) — green (168 fast + the slow engine seals + slow planet;
the live-climlab and notebook smoke-tests skip, as in CI). (The import-drift guard is a monorepo-only
mechanism — see ADR 0003's standalone-repo extraction note.)

**Rung 1, step 2 — Phase A BUILT (the two-way feedback *machinery*; 2026-06-10).** `planet/transport.py`
closes the loop given *any* meridional eddy heat flux `⟨v'θ'⟩(φ)`: the **κ→D bridge** `D = C_atm·κ/a²`
(physical/citable — `C_atm = c_p·p_s/g`; rung-0 `D=0.555` ⟺ `κ≈2.2×10⁶ m²/s`, the observed midlatitude
eddy-diffusivity order) maps the **band-bulk** down-gradient diffusivity (least-squares over the
window-flat interior) to an EBM transport coefficient, which **re-equilibrates** the EBM (uniform
`D_eff` is the headline; `EnergyBalanceModel` now also accepts a callable `D(x)` for the band-limited
diagnostic). The *design* anchor is **reduction-to-EBM:** the closure `⟨v'θ'⟩=−D_eff·∂θ̄/∂y` has the
*same form* as the EBM transport term, so the two-way model with a constant flow-diagnosed `D_eff` *is*
a rung-0 diffusive EBM with that `D`. **What Phase A actually validates (advisor-corrected honesty):**
(i) the **bridge pinned absolutely** (`C_atm≈1.037e7`, `κ₀≈2.17e6` — *not* just round-tripped, which
would let a wrong `a²`/`C_atm` cancel); (ii) the **right-signed response** (stronger flux ⇒ flatter
contrast — the EBM's genuine physical response, non-tautological); (iii) the **plumbing** (κ recovered
+ `D_eff` routed into re-equilibration, sign rejected if up-gradient). The reduction itself **reduces
to rung-0 by construction** in Phase A (the re-equilibration re-runs the scalar-`D` EBM) — plumbing,
not an independent test; the **genuinely tight reduction** (an independent two-way budget whose
flux-divergence matches the EBM operator, *plus* the **Cartesian-channel ↔ spherical-EBM geometry
correspondence** — the bridge is derived for uniform κ on the sphere) needs the emergent flux and
arrives in **Phase B**. **De-risked the A/B split** (advisor): Phase A drives the machinery with a
**synthetic** down-gradient flux (the Phase-4 synthetic-gradient playbook), so it lands *independent*
of the (tuned) eddy sim; the `flux_fn` argument is the seam Phase B plugs into. Planet gate
**151 passed, 1 skip** (no shared-engine edit → planet gate, not full-repo). **Next (step 2, Phase B):**
the *emergent* eddy flux — advect θ (relaxed to the EBM target) on the barotropically-unstable jet,
diagnose `⟨v'θ'⟩` **post-saturation** via a life-cycle integral (release mode), magnitude named tuned,
+ a **`D_eff`-tracks-climate** non-circularity test (warm/flatter ⇒ smaller `D_eff`, else the loop is
cosmetic); the loop-closes claim scoped to **one feedback pass** (a converged fixed-point is a `slow`
demo if it converges cleanly). Then **step 3** circulation-informed precip.

**Rung 1, step 2 — Phase B BUILT (the EMERGENT eddy flux; 2026-06-11).** `planet/eddy_flux.py` fills
Phase A's `flux_fn` seam with the real thing: the meridional eddy heat flux `⟨v'θ'⟩(φ)` diagnosed from
a passive temperature tracer advected on the **released** barotropically-unstable Phase-4 jet
(`eddy_life_cycle`: spin up the jet **dry** → init `θ` = windowed-EBM profile → deterministic `cos(kx)`
v-perturbation → forcing-**OFF** release → life-cycle integral `κ_eff = −∫F̄ dt / ∫θ̄_y dt` over the
window-flat interior). **Genuinely emergent** — resolved barotropic instability (Rayleigh–Kuo met, jet
~20 m/s), no imposed stationary wave, no down-gradient closure assumed. **What is banked, and its
honesty class** (the spike-driven build — see below): **(headline, DIRECTION) the eddy diffusivity is
STATE-DEPENDENT** — across two climates with the forcing amplitude `α` held **fixed**, a flatter EBM
gradient (high-obliquity-like, `s₂=−0.32`, jet ~14 m/s) gives `κ_eff ≈ 0.5–0.6×` the steep climate's
(`s₂=−0.48`, jet ~20 m/s): the loop is a real, right-signed feedback (flatter mean → weaker jet →
weaker eddies → smaller `κ`), **not** a fixed-`D` re-labelling. (Mechanism, advisor: the gradient
*cancels* in `−∫F̄/∫θ̄_y`, so `κ_eff ≈ v'·ℓ` tracks climate **only** through the jet — which is exactly
why `α` must be held fixed; renormalizing to fix the jet speed would make the test cosmetic by
construction. Warming via `CO₂`/`S₀` was found to barely move the *channel* gradient — the linear OLR
shifts `T₀` uniformly and the ice retreat is poleward of the channel — so `s₂`/obliquity, which
genuinely flattens it, is the non-circularity knob.) **(MAGNITUDE — named, NOT banked)** `κ_eff ~ 10³
m²/s`, ~1000× below rung-0's `2.2×10⁶`: the *instantaneous* `⟨v'θ'⟩` is largely **reversible**
(oscillates sign with the meander; **irreversible fraction ~0.1**), and the value is **resolution-
converged** (`nx=80`≈`96`) but suppressed by configuration choices that can't be cleanly separated from
intrinsic physics (a single coherent seeded wavenumber mixes more reversibly than broadband turbulence;
the band-bulk estimator smears the jet-centred peak) — so it is **window/forcing/configuration-tuned**,
the honesty class of Phase 4's jet *speed*; the **sign** and the **climate ordering** carry validation,
not the number (chasing a bigger number is unbankable regardless). **(the tight reduction — a FINDING,
not a manufactured match)** the barotropic flux does **not** *tightly* reduce to the EBM operator at
rung 1: `reduction_to_ebm_operator` **tests** (not assumes) the resolved flux-divergence's *shape*
against smooth down-gradient diffusion built from the band-bulk **scalar** `κ` (not circular — the
pointwise `κ` is not reused; only the normalised shape) → a **partial** correlation (`~0.6`), and the
comparison is itself **near-vacuous** (a uniform-`κ` diffusion of the near-linear midlatitude gradient
produces ~0 transport-divergence); the tight reduction becomes **non-vacuous only at rung 3** (a strong
baroclinic flux). **(the geometry correspondence — DELIVERED, the genuinely tight part)**
`planet/transport.py` gains the **spherical transport operator** `(1/cosφ)∂/∂y[κ·cosφ·∂θ/∂y]` (the EBM
operator `D·∂/∂x[(1−x²)∂T/∂x]` written in β-plane coordinates), anchored on the **P₂ eigenvalue** (both
forms → `−6·(κ/a²)·P₂`, the analytic Legendre check — **not** a self-comparison of two finite-difference
operators), with the flat-vs-spherical `cos φ` metric gap shown **order-unity (~0.6)** over the wide
channel (`φ≈19°–61°`, `cos φ` varying ~2×): so the bridge's "uniform κ on the sphere" derivation is made
rigorous — the geometry is **not inherited for free**, and the β-plane tangent across a ~42° band is a
named scope edge. **(seam) `close_loop`** routes the emergent `D_eff` through the Phase-A-validated
bridge + re-equilibration: it **converges** with the right sign (weaker transport ⇒ equator-to-pole
contrast 56→122 °C, steeper), but the degenerate near-radiative-equilibrium climate is **not banked**
(magnitude). **Build discipline (advisor-gated):** the physics was **de-risked in throwaway spikes
*before* the module** (`outputs/`, gitignored) — the spikes caught that a *single-transient* diagnosis
is noisy and that *tracer relaxation during release* **breaks** the climate ordering (a τ-artifact), so
the settled **pure-release** life-cycle integral over a fixed post-saturation window is what landed
(flat<steep across every window tried). Tests split by cost: the geometry legs (P₂ eigenvalue +
order-unity metric) are **fast/always-green** in `test_transport.py`; the eddy-sim legs (unstable +
down-gradient + mostly-reversible, state-dependence, non-tight reduction, `close_loop` sign) are
**`slow`** in `test_eddy_flux.py`. No engine edit (planet remains a *consumer* of `engines/fluid` +
`engines/diffusion`); `uses` unchanged.

**Rung 1, step 3 — BUILT (circulation-informed precip; 2026-06-11). RUNG 1 COMPLETE.**
`planet/circ_precip.py` wires the precip **storm-track band centre** to the **emergent jet latitude**
instead of the prescribed constant (the §3 deep-end hook: "rain where the flow puts the storm track,
still without moisture physics"). `precip.precip_pattern` gains a `midlat_center_deg` (default = the
cited 50° → the rung-0 field **bit-for-bit by construction**, the `two_way_pass`-style plumbing
reduction); `circ_precip.circulation_informed_precip(state, jet)` feeds it `jet_lat`. **De-risked in
two throwaway spikes first** (`outputs/rung1_circprecip*`, gitignored — this project's discipline), and
the spikes set the headline + killed an anchor. **What is banked:** (1) the **seam** (centre ← emergent
circulation latitude) + (2) the **reduction** (jet at 50° ⇒ rung-0 exactly) + (3) the **migration
mechanism** — the band tracks a *dynamically-selected* latitude, shown via the coupler's
synthetic-gradient playbook, **anchored to `jet_lat` not the EBM `gradient_peak_lat`** (the two nearly
coincide at present-day, within ~2–3°; the gap opens only off-centre — so it is a *flow* response, not an
EBM-gradient diagnostic). **The rung-1 FINDING (named, NOT an accuracy gain):** the dry circulation
**cannot refine** the rain location — it is a **trade, not a ranking**: the model's own jet sits at ~42°,
**~8° equatorward** of
Earth's observation-calibrated 50° (the Phase-4 channel's known equatorward bias — it excludes the ice
cliff), and across *realistic* obliquity/CO₂/S₀ the gradient peak is pinned ~43–46° (spike #1), so
realistic **migration is mechanism-only** (decisive only under a synthetic gradient). The literal "rain
where the flow **converges**" anchor (centre on the eddy heat-flux convergence `−∂F̄/∂y`) was **tested
in spike #2 and rejected**: the resolved convergence is **near-vacuous in the channel interior + a
window-taper edge dipole**, not a physical storm-track convergence — the **same rung-3 boundary**
`eddy_flux` already found (non-vacuous only under a strong baroclinic flux, where the *shape*-resolving
precip refinement belongs, with the geometry already delivered in `transport.py`). So rung 1 wires the
**position** seam; the **shape/amplitude** refinement (wet-get-wetter) stays deferred, the **ITCZ/
subtropics** stay prescribed (Hadley is out of the midlatitude channel), and **rung-0 `precip.py`
remains the default** in the biome map/demos (circ-informed is **opt-in** — the relocation trades
observational calibration for internal consistency, and would regress the calibrated map). Scope edge:
a large *equatorward* displacement shallows the subtropical trough toward merging with the ITCZ (the
band-tracking is asserted only where the structure survives, centre ≳ 36°). Tests: fast
reduction/migration/structure + **one** `slow` composition (band follows the emergent jet on a synthetic
off-centre gradient — the coupler's own jet-tracks-gradient proof is not re-tested), `test_circ_precip.py`.
No engine edit; `uses` unchanged.

**Rung 1 within-rung — wet-get-wetter, dry-get-drier BUILT (2026-06-14; the §12.2 amplitude slice).**
`planet/moist.py` gains `wet_get_wetter_precip_field` (+ `_amplify_contrast`, `WetGetWetter`/`wet_get_wetter`):
the thermodynamic **contrast sharpening** the rung-0 uniform `CC(T̄)` omits. **Home corrected by the advisor:
`moist.py`, NOT the §12-stated `circ_precip.py`** — the §12 pointer predated the rung-2 build, and this is
literally the **generalization** of the one-line `energy_constrained_precip_field` (which scales mean+anomaly
*together* at one rate; this **splits** them — `P = ⟨P⟩·M(T̄) + (pattern−⟨P⟩)·W(T̄)`, the mean `M` at the
energy-constrained ~2.5 %/K, the anomaly `W` at C–C ~7 %/K). Held & Soden 2006 "rich-get-richer" on the
precip pattern: under warming the ITCZ/storm-track bands **intensify** while the **deserts dry** (where the
rung-0 uniform `CC` wrongly wettens them too). **Honesty class (advisor) = a better *prescribed*
parameterization, NOT derived:** the split *direction* and the two rates are **calibrated/cited** (the
energy slope is the named sub-grid wall, C–C the moisture-capacity rate); what is **structurally exact** is
the **mean-zero anomaly split** (`⟨pattern−⟨pattern⟩⟩=0` to machine precision on the equal-area grid ⟹ ⟨P⟩
scales at the energy rate — **plumbing, by-construction, not a finding**) and the **reduction to BOTH**
existing fields when the rates coincide (`M=W=CC`→rung-0; `M=W=energy`→energy field). **Opt-in/default-off**
(rung-0 `precip.py` unchanged), and **deliberately not fused** with the storm-track *position* seam (the
moisture-budget non-composition rule). Named edges: **global-`T̄`** anomaly factor (local-`q_sat(T(φ))` =
richer named upgrade); **`P≥0` floor** (deep warming → the dry minima to total aridity, the linearization
breaks); **thermodynamic only** (the *dynamic* circulation-driven amplification = the moisture-convergence
path, rung 3+). Demo `planet/demo_wet_get_wetter.py` → `docs/figures/planet-wet-get-wetter.png` (the field
panel + the ΔP warming response: rung-0 positive everywhere vs wet-get-wetter negative in the subtropics).
Tests `test_moist.py` (7 fast + 1 `slow` demo guard); no engine edit, `uses` unchanged.
[[planet-rung1-two-way-coupler]] [[precip-parameterization-source]]; §12.2.

**Rung 2 — SCOPED (design + fork settled; NOT built), 2026-06-11.** The next staircase rung (moist
dynamics → emergent precipitation) is **scoped and its central fork settled empirically**, ahead of any
code (the project's spike-first discipline; advisor-pressure-tested twice). **The fork — *where moisture
lives* — is decided: a *column moist budget*, not fluid-channel moisture transport.** A throwaway spike
(`outputs/rung2_moisture_convergence_spike.py`, gitignored) advected a steep, C–C-shaped **moisture**
tracer on the *same* released barotropically-unstable Phase-4 jet that :mod:`planet.eddy_flux` uses, and
measured its meridional flux **convergence** ``−∂⟨v'q'⟩/∂y`` over the window-flat interior: it is
**near-vacuous in the interior + dominated by the window-taper edge artifact** (interior/edge RMS ratio
**0.32**, *worse* than the temperature tracer's 0.50 — the steeper moisture gradient does **not** rescue
it), and the **mean** balanced-jet flow is near-nondivergent (eddy/mean convergence 11–23×, since
``∇·(hu)=−∂h/∂t≈0``). This **corroborates** rung-1 step-3's rigorous `reduction_to_ebm_operator` finding:
the single layer has no genuine ascent→condensation, so a *resolved storm-track precip pattern* is **rung
3** (the vertical), **not rung 2** — confirmed, not just predicted (one realization suffices given the
convergent evidence; no further spikes — "don't gold-plate the spike", advisor).

**So the bankable rung-2 core is a column moist budget — a moist-EBM *diagnostic* that reuses the
diffusion spine a 4th time and does NOT perturb the validated Phase-1 climate.** Design: a diagnostic
moisture field ``q(φ)=RH·q_sat(T)`` (fixed relative humidity over the rung-0 temperature; ``q_sat`` from
Clausius–Clapeyron) with a **down-gradient latent transport** whose diffusivity is **tied to rung-1's
eddy ``κ`` through the existing :mod:`planet.transport` κ→D bridge** (the *same* eddies stir heat and
moisture — reuse + a consistency check, **not** a new free ``D_q``), and a precipitation ``P``
**diagnosed** from the atmospheric water budget ``P = E − ∇·(moisture transport)``. **The headline unlock
(the banked physics) is the RATE, not the pattern:** it replaces :mod:`planet.precip`'s prescribed global
Clausius–Clapeyron **7 %/K** amplitude with the **energy-constrained ~2–3 %/K** global precipitation rate
— closing the gap `precip.py` *already names in its own scope-edge #3* (moisture-capacity 7 %/K vs the
slower energy-constrained global rate). **The first sub-grid closure (the staircase's named rung-2 wall)
is the atmospheric-energy closure that sets that rate** — an atmospheric radiative-cooling term
``R_atm(T)`` with ``L⟨P⟩ ≈ R_atm − SH`` — **not** the benign fixed-RH assumption; the rate claim's honesty
lives there, named explicitly. **``P`` is a pure diagnostic — it does *not* enter the temperature equation
— so the Phase-1 climate and its triad stay green** ("energy-limited evaporation" would *overstate* it:
this is an atmospheric-cooling *constraint* on ``⟨P⟩``, **not** a closed surface energy budget, which would
change ``T``).

**The honest emergent-pattern finding (named, NOT a win):** down-gradient moisture diffusion captures the
**extratropical** budget (midlat/polar ``P>E`` convergence, subtropical ``E>P`` evaporative source) but
gets the **deep-tropical ITCZ backwards** — diffusion *exports* moisture from the moist equator, whereas
the real ITCZ is *up-gradient* Hadley convergence (a mean-circulation feature, deferred). So Phase A is
again a **trade**: emergent + rate-correct in the extratropics, but the prescribed `precip.py` ITCZ is
*better* in the deep tropics — Phase A banks the moist **energetics + extratropical budget + the global
rate**, not a wholesale-better precip map (rung-0 `precip.py` stays the default, as with circ-informed
precip). **Triad (re-classed for honesty):** *tight* — ``q_sat`` = the **exact Clausius–Clapeyron
function** (the Whittaker-partition precedent: an exact testable function, not a fit); *the real-but-loose
physics result (the unlock)* — the energy-constrained ~2–3 %/K rate; *consistency / plumbing (named as
such)* — global ``∫E=∫P`` (falls out of moisture-mass conservation in steady state) and the reduction to
rung-0 as ``L·q→0`` (a vanishing moisture layer — **by-construction** plumbing, the rung-1 Phase-A
`two_way_pass` honesty class, **not** an independent test); *benchmark (loose)* — the observed
extratropical ``E−P`` belts. **Named scope edges:** ITCZ/Hadley deferred (mean circulation); the spatial
storm-track precip pattern is rung 3 (spike-confirmed above); the ``R_atm`` atmospheric-cooling closure is
the sub-grid wall. **A fuller moist EBM — diffuse *moist static energy* ``m=c_pT+L·q`` so ``T`` itself
responds (emergent polar amplification via moisture) — is a named rung-2.5 extension, NOT Phase A:** it
re-opens the Phase-1 ``(A,B,D)`` calibration (rung-0's ``D=0.555`` is an *effective* diffusivity already
absorbing latent transport — explicit MSE diffusion double-counts the latent heat implicit in the linear
OLR ``A+B·T``), so it is deferred behind the clean diagnostic. **Sources to pin at build** (the
``[[…-source]]`` discipline — named now, pinned when Phase A is built, **not** carried from memory): the
energy-constrained ~2–3 %/K-vs-moisture-capacity-7 %/K rate → **Held & Soden 2006 / Allen & Ingram 2002**
(extending `[[precip-parameterization-source]]`, which already cites the gap); the diffusive-moist-EBM
formulation (MSE / latent diffusion, fixed RH) → **Flannery 1984 / Hwang & Frierson 2010 /
Siler–Roe–Armour 2018**.

**Rung 2, Phase A — BUILT (the column moist-EBM diagnostic; `planet/moist.py`, 2026-06-11).** Built
spike-first (`outputs/rung2_phaseA_moistebm_spike.py`, gitignored) and advisor-pressure-tested, which
**reshaped the deliverable and overturned one scoped claim** (see below). Sources pinned at build (the
`[[…-source]]` discipline, **not** carried from memory): C–C `q_sat` (integrated form + `e₀, L_v, R_v,
ε`) → **Hartmann *GPC* / Bohren & Albrecht**; the energy-constrained rate → **Held & Soden 2006 /
Allen & Ingram 2002**; the diffusive moist EBM (fixed RH, latent diffusion) → **Flannery 1984 / Hwang &
Frierson 2010 / Siler–Roe–Armour 2018**. **The advisor split the deliverable in two** (do not fuse — a
full emergent `P` field forces an unphysical evaporation pattern, see below):
1. **The RATE — the headline unlock, robust, opt-in.** `energy_constrained_factor(T̄)` replaces
   `precip.py`'s C–C **7 %/K** amplitude with the **energy-constrained ~2.5 %/K** (`L⟨P⟩ = R_atm − SH`,
   normalised by `⟨P⟩₀≈100 cm/yr`) — closing `precip.py`'s own scope-edge #3. The energy budget is
   **linear** in `T̄` (so the factor is **linear**, *not* a smaller C–C exponent — the honest functional
   difference), `= 1` at the present reference (present map unchanged, **bit-for-bit**), floored ≥ 0.
   **Opt-in** (`energy_constrained_precip_field`); rung-0 `precip.py` stays the default. `precip.py` is
   **untouched** — its warmed-climate 7 %/K tests stay green.
2. **The emergent budget `P − E` — the trade, a diagnostic, NOT the default.** `P − E = (D/c_p)·∂/∂x[(1−
   x²)∂q/∂x]` with `q=RH·q_sat(T)` over the rung-0 `T`. **The latent heat `L` CANCELS** (the same eddies
   stir heat and moisture → the latent transport is the EBM operator on `L·q` with the **same `D`** via
   the rung-1 κ→D bridge; converting the latent-energy convergence back to a moisture mass flux divides
   out `L`) → **no new free `D_q`**, rung-1's eddy diffusivity reused. Built in **conservative face-flux
   form** on the EBM `x`-grid (the diffusion spine's structure, a **4th** reuse) → `∫(P−E)=0` machine-
   exact (the `∫E=∫P` plumbing leg, under the area-mean rule — *not* `np.trapezoid`, which breaks the
   telescoping). A **pure diagnostic** — does NOT enter the T-equation, so Phase-1 stays green
   (asserted). **Why P−E not full-P (advisor catch):** equatorial export is ~−2.4 m/yr but `⟨P⟩~1` m/yr,
   so *no* honest zonal `E` keeps `P≥0` (uniform `E`→`P(eq)<0`; `E∝q_sat`→absurd ~6 m/yr equatorial
   evaporation) → report `P−E`, skip full-P. **THE OVERTURNED CLAIM (the spike corrected the scope):**
   the budget is **extratropical-ONLY**, not "extratropics good + subtropics good". The deep equator is
   **backwards** (diffusion *exports* from the moist equator; real ITCZ is up-gradient Hadley, deferred)
   **AND the subtropical evaporative belt is NOT reproduced** — the steep equator–pole contrast
   hyper-peaks C–C `q` at the equator, pushing the moisture-flux maximum equatorward (~20° zero-crossing)
   so the subtropics come out **`P>E`** (production: eq export ≈ −267 cm/yr, subtropics 25–35° ≈ +83,
   midlat 40–60° ≈ +101, polar ≈ +28). Only the **extratropics (poleward ~40°) are right**. So the
   benchmark test asserts **only** equatorial export + extratropical `P>E` convergence, and **one test
   pins the subtropical mislocation** as the honest limitation (guarding against a silent "fix"). Same
   "trade, not a win" the staircase keeps banking; rung-0 `precip.py` stays the default map.
   **The sub-grid WALL = the prescribed `R_atm` slope** (`R_ATM_SLOPE = 2 W m⁻² K⁻¹`, cited Held & Soden;
   the rate is a *cited-closure* result, **not derived** — and explicitly **not** `B_OLR`, a different
   2 W m⁻² K⁻¹ quantity). **Triad (re-classed):** *tight* — `q_sat` exact C–C (textbook values + ~7 %/K
   log-slope) + the operator reproduces the **P₂ eigenvalue `−6`** (it *is* the EBM transport operator,
   the `transport.py` anchor); *real-but-loose (unlock)* — the ~2–3 %/K rate, slower than C–C, doubling
   with the closure slope; *plumbing* — `∫(P−E)=0` machine-exact + `q→0` reduction + unity-at-reference;
   *benchmark (loose)* — equatorial export + extratropical convergence (the named extratropical-only
   trade). **Deferred:** ITCZ/Hadley (mean circulation); the resolved storm-track precip *pattern* (rung
   3, the vertical, spike-confirmed); the full MSE-diffusing moist EBM where `T` responds (**rung 2.5** —
   re-opens the `(A,B,D)` calibration — **now BUILT, see below; refined to D-only**). Tests:
   `planet/tests/test_moist.py` (15, all **fast**); full planet gate **179 passed, 1 skip**. No engine
   edit; `uses` unchanged.

**Rung 2.5 — BUILT (the MSE-diffusing moist EBM where `T` responds; `planet/moist_ebm.py`, 2026-06-12).**
The named step up from rung-2 Phase A: rung 2 added moisture as a **pure diagnostic** that never touched
the `T`-equation; rung 2.5's whole point is that **`T` itself responds**, so the headline is an *emergent
climate response* the diagnostic could not produce — **polar-amplified warming** (Hwang & Frierson 2010;
Flannery 1984; Siler–Roe–Armour 2018). It is a **separate model alongside rung-0, NOT a replacement** (the
dry EBM stays the default; this is the opt-in moist sibling, as circ-informed precip and the
energy-constrained rate were). Built **spike-first** (`outputs/rung25_moist_ebm_spike.py`, gitignored) and
advisor-pressure-tested — **the advisor caught a load-bearing math error and added the attribution
backbone** (recorded below).
- **The mechanism.** A moist atmosphere diffuses **moist static energy** `m = c_pT + L·q`, `q=RH·q_sat(T)`;
  in temperature-equivalent units this is a **temperature** diffusion with a *moisture-amplified*
  effective coefficient `D_eff(T) = D_s·(1 + β(T))`, `β = (L/c_p)·RH·dq_sat/dT`
  (`moisture_amplification`). Because C–C `q_sat` is steep, **β is large in the warm tropics and ≈0 at
  the cold pole** (Earth: `D_eff` ~1.3 equator → ~0.35 pole). As the climate warms the *tropical* β grows
  fastest → tropics export more heat poleward → **poles warm more than the tropics**, emergent from
  moisture alone (no ice-albedo feedback, no change in `D_s`). The dry EBM (constant `D`) warms **exactly
  uniformly** under a uniform forcing — the clean null.
- **THE ADVISOR'S MATH CATCH (load-bearing).** `D_eff` must sit **INSIDE** the flux divergence
  `∂/∂x[(1−x²)·D_eff·∂T/∂x]`, not outside (`D_eff·∂/∂x[…]`) — the cross-term `(1−x²)(∂D_eff/∂x)(∂T/∂x)`
  is the *same order* as the main term and is *literally where part of the PA mechanism lives*; the
  outside form also breaks energy conservation. **`ebm.py` already places the callable `D(x)` inside** (the
  rung-1 array-`D` path), so the design passes `D_eff` as that callable and lets the engine place it
  (spike-confirmed at machine precision via conservation + a direct operator match).
- **Design (5th reuse of the diffusion spine).** The moist climate is **one nonlinear relaxation**, not
  nested solves: each Strang substep **freezes `D_eff` at the current `T`** and rebuilds the conservative
  transport operator — the **identical idiom the ice-albedo `α(T)` already uses** in
  `ebm.equilibrate` (a state-dependent coefficient re-frozen each substep). **Self-contained — does NOT
  modify `ebm.py`** (the validated rung-0 hot path untouched; the ~one duplicated radiation helper is the
  correct price). `face="harmonic"` pinned (the bit-for-bit reduction relies on the per-step cells
  matching the dry model's once-built harmonic cells; `face="exact"` would silently break it).
- **THE HEADLINE FRAMING (advisor) — redistribution around a PINNED mean, not added pole warmth.** With
  constant albedo the global-mean response to a uniform OLR forcing `ΔA` is **pinned**: `⟨δT⟩ = ΔA/B` to
  machine precision for *any* `D` (the diffusion conserves `∫T dx`, so transport cannot change the mean).
  **Moisture REDISTRIBUTES that fixed `⟨δT⟩` poleward.** Asserted **tight** (conservation). The PA factor
  is **reported two ways (name the metric, advisor):** the **single-endpoint** ratio `δT(pole)/δT(equator)`
  ≈ **2.05** (Earth, RH 0.8 — the *most generous*, polar cell on the harmonic-face bias), and the
  **area-band** ratio `mean(δT|φ≥60°)/mean(δT|φ≤30°)` ≈ **1.80** (less generous) — both honest. **Direction
  banked** (PA>1 robustly), **magnitude loose** (the observed ~2–3× also needs the ice-albedo + lapse-rate
  feedbacks held out of scope; and ~2.05 is the model's *converged* number, not a *validated* one — the
  formulation sets it).
- **THE dt-FREE RE-BANK (2026-06-14, a splitting-error finding).** The headline was originally banked as
  endpoint **≈1.5** / band **≈1.4** — those were a **first-order operator-splitting artifact**, not the
  physical answer. The PA had been read off the Strang relaxation at the default `n_tau=0.5`, whose **shape**
  carries an **O(Δt) splitting bias** (backward-Euler transport split against the exact radiation half-step;
  the global *mean* stays exact — `pa_dry≡1.0` exactly, the clean control), and that bias *suppresses* the
  amplification. Recomputing **dt-free** lifts it to **endpoint ≈2.05 / band ≈1.80** (Richardson dt→0 and a
  direct solve agree to 3 digits). The fix: `moist_steady_direct` — a **Picard iteration on the
  frozen-`D_eff` linear solve `(L_T−B·I)T=A−S(1−α)`**, the *nonlinear generalisation of
  `ebm.steady_linear`* (exactly as rung-4 went Newton), **no time-stepping ⟹ no splitting bias**, ~20
  iterations / ~10 ms. The relaxation `moist_equilibrium` is kept (bit-for-bit rung-0 reduction +
  animations); the **headline magnitude is the direct solve's**. *Direction was never in doubt* — only the
  magnitude moved, and it moved the same way the under-converged number always under-reported. Spike +
  findings: `outputs/rung25_splitting_dt_spike.py`, `outputs/rung25_picard_prototype.py` (gitignored).
- **THE ATTRIBUTION NULL (advisor — the backbone).** The moist model differs from dry in *two* ways (a
  recalibrated `D`-shape AND the T-dependent `D_eff`); which causes PA? **Freeze `D_eff` at its present
  profile and warm → PA = 1.0 *exactly*** (uniform `δT=ΔA/B` solves the perturbation for *any* frozen
  `D(x)`), proving the PA is **100 % the `dD_eff/dT` feedback, 0 % the `D`-shape** (spike + test:
  PA=1.0000, spread ~3e-10, via the genuine array-`D` `EnergyBalanceModel`).
- **The recalibration = the named wall (the double-count).** Rung-0's `D=0.555` is an *effective*
  diffusivity already absorbing latent transport; explicit MSE diffusion would double-count it, so
  `recalibrate_sensible_D` re-derives a **smaller sensible `D_s≈0.28`** matching the dry present-day
  **equator-pole contrast** (`⟨T⟩` is automatically equal — energy balance fixes it from `A,B,ᾱ`
  independent of `D` — so contrast is the natural single scalar). **The TRADE (not a win):** a *single
  scalar* can't reproduce all of dry `T(x)` — same mean+contrast but a **higher-moment shape residual**
  (matched-contrast moist profile flatter in the interior, curvature toward the edges). **Target is a
  modeling choice (named):** matching the P₂ amplitude `T₂` instead moves PA **< 5 %** (PA set by the
  *shape* of β, not the scaling).
- **Scope edges (named).** Refines the plan's "rung 2.5 re-opens `(A,B,D)`" → **only `D` re-opened**: `A`
  is the forcing knob (uniform `ΔA` = CO₂ proxy), **`B` held FIXED** — re-deriving B's water-vapour
  content is *local radiation = the rung-4 wall*, not opened. Forcing is **uniform `ΔA`, not `ΔS₀`** (`S₀`
  is equator-weighted by `(1+s₂P₂)`, imposing tropical structure that *fights* PA). Fixed RH; **constant
  albedo** for the clean Hwang–Frierson experiment (ice/lapse-rate feedbacks out of scope). **Separate
  from the rung-2 `P−E` budget** — a `T`-response model, doesn't touch `moist.py`; rung-2's
  `test_subtropical_evaporative_belt_is_not_reproduced` did **not** flip (a different deliverable).
- **Triad.** *Tight* — exact analytic `dq_sat/dT` (vs finite-diff + ~7 %/K log-slope); `D_eff`-inside
  conservation; pinned `⟨δT⟩=ΔA/B`; the frozen-`D_eff` PA=1 null. *Real-but-loose (unlock)* — the PA
  itself (dt-free endpoint ratio ~2.05 / area-band ~1.80, direction banked / magnitude loose). *Plumbing* —
  RH=0 ∧ `D_s=0.555` ⟹ the genuine `EnergyBalanceModel` rung-0 solve **bit-for-bit** for the relaxation,
  and to **machine precision** for the dt-free `moist_steady_direct` → `steady_linear`. A test also pins the
  splitting artifact (the `n_tau=0.5` relaxation reads ~1.5, climbing toward the dt-free value as `n_tau→0`,
  the mean exact throughout). *Named choices* — recalibration to present contrast (`D_s<0.555`) + its
  target-invariance + the **named PA metric** (endpoint vs area-band, the advisor's "don't let the
  most-generous number read as *the* number"). Tests: `planet/tests/test_moist_ebm.py` (14, all **fast**,
  incl. a `.converged` guard); full gate **443 passed, 1 skip**. No engine edit; `uses` unchanged.

**Rung 2.x — BUILT (the full-sphere EBM + the energetic ITCZ; `planet/sphere_ebm.py`, 2026-06-14).** The
chosen slice of the user's "ITCZ/Hadley + moist precip pattern" deferral. Rung-0's EBM is **hemisphere-only**
(`x=sinφ∈[0,1]`, equatorial-symmetry BC) — so an ITCZ that **migrates off the equator breaks that symmetry
by construction** and *cannot* be represented; this rung lifts rung 0 to the **full sphere** `x∈[−1,1]`
(two real poles, equator interior) and adds the **energetic ITCZ**: the latitude where the emergent
atmospheric energy transport `H(x)=−2πa²·D(1−x²)∂ₓT` crosses zero — the **energy-flux equator (EFE)**, which
migrates toward the warmer hemisphere under an interhemispheric imbalance (Kang+2008; Bischoff–Schneider
2014; Schneider–Bischoff–Haug 2014; Donohoe+2013). **A SIBLING model — `ebm.py`/`moist.py` UNTOUCHED** (the
`moist_ebm`/`baroclinic_qg` discipline): "re-validate Phase-1" = a **cross-model reduction check** (sibling
reproduces hemisphere `ebm.py` to **1e-9** under symmetric forcing), the SW↔QG-bridge pattern. Built
**spike-first** (`outputs/rung_itcz_fullsphere_spike.py`, gitignored) + **advisor-pressure-tested twice**,
which set the honest altitude (below). Sources pinned at build (the `[[…-source]]` discipline).
- **Banked (tight):** the full-sphere sibling; reduction to `ebm.py` (1e-9); North two-mode via the direct
  `steady_linear` (constant albedo, **harmonic-face polar floor ~0.16 °C — NOT clean 2nd order**, same as
  `ebm.py`); the **closed form `δ/AHT = 1/(2πa²·D·T̄ₓₓ(0))` reproduced by the engine**; EFE=0 exactly for
  symmetric forcing; global energy balance machine-exact. The EFE diagnostic + the opt-in **`itcz_center_deg`
  precip seam** (`precip.precip_pattern`/`precipitation` gain it, default 0 → rung-0 **bit-for-bit**; the
  ITCZ band uses *signed* latitude so it can migrate; `sphere_ebm.itcz_informed_precip` feeds it `φ_EFE`,
  the rung-1 `circ_precip` pattern for the ITCZ).
- **Banked at the LOWER altitude (the unlock — advisor's load-bearing catch):** the ITCZ sensitivity is a
  **closed-form consequence of the already-calibrated `D` and the mean-state curvature**, **NOT an emergent
  prediction**. In a dry EBM the EFE is just the **temperature maximum** (`∂ₓT=0`), so the sensitivity is the
  algebraic `1/(2πa²·D·T̄ₓₓ(0))`; the forcing-independence (Q-flux and albedo-asymmetry give the **same**
  number) is a **linear-operator identity, not robustness**, and the shift **direction is by-construction**
  (the "guaranteed result" trap the QG rung carried). Values: **≈ −6.3 deg/PW (no-ice, splitting-free) and
  ~ −5 deg/PW (present-day ice)** — the **same ORDER** as the observed **~3 deg/PW** but a **factor ~1.5–2
  high**, and `∝ (6+B/D)` (curvature `∝1/(6D+B)`, so a *pure function of D*, **not** the naive "∝1/D"). The
  honest one-liner: ***"corroborates that `D` is realistic," NOT "predicts the ITCZ migration."*** (The
  spike's first −3.78, suspiciously near observed-3, was an **operator-splitting artifact** at the default
  `n_tau=0.5`; the EFE/sensitivity must be read off a **converged** profile — `steady_linear` or small
  `n_tau` — a real gotcha, code- and test-enforced.)
- **Named edges (every one carried):** **dry EBM** → the EFE is *identified* with the ITCZ via external moist
  theory and the wire **relocates a prescribed band, NOT emergent rainfall**; the asymmetry is **imposed**
  (Q-flux/albedo, the coupler's synthetic-gradient precedent), not from an ocean model; **`moist.py`'s
  `moisture_convergence` stays backwards in the deep tropics** — this rung adds ITCZ **position**, the Hadley
  moisture-convergence fix (the literal backwards-`P−E`) is **still deferred**; annual-mean (no seasonal
  migration). **Demo banked + CI-guarded** (`planet/demo_sphere_itcz.py` → `docs/figures/planet-sphere-itcz.png`;
  the `slow` `test_demo_reproduces_the_banked_headline` pins engine==closed-form, observed-order, right sign,
  band relocates). Tests: `planet/tests/test_sphere_ebm.py` (15: tight/unlock/plumbing fast + 2 slow); full
  gate **399 passed, 1 skip**. `ebm.py`/`moist.py` untouched; `precip.py` additive (default-preserving).

**Rung 2 — the Hadley moisture-convergence fix BUILT (the literal deep-tropical backwards-`P−E`;
`planet/moist.py`, 2026-06-14).** The last named rung-2 deferral (the one rung 2.x explicitly *left*
deferred — it fixed ITCZ *position*, this fixes the deep-tropical *sign*). The eddy-only
`moisture_convergence` (`P−E=(D/c_p)·∂/∂x[(1−x²)∂q/∂x]`) is **structurally** backwards at the moist equator:
down-gradient diffusion *exports* moisture from a maximum — there is **no diffusive way to converge moisture
at the ITCZ** (advisor-confirmed: not gettable by tuning `D`; the mean-circulation term is genuinely
needed). Built **spike-first** (`outputs/rung2_hadley_moisture_spike.py`, gitignored) + advisor-pressure-
tested **before** the build (the honesty classification was set up front — the load-bearing move). The fix
adds the **mean Hadley circulation** as an **opt-in** term (`hadley_moisture_convergence`; `moisture_budget(
..., hadley=True)`), eddy-only stays the **default** (every existing benchmark test stays green — the
opt-in/independent-diff discipline of `energy_constrained_precip_field` / `circ_precip` / `itcz_center_deg`).
- **The model.** A **prescribed** tropical overturning: the northward MMC moisture flux `F(x)=−strength·
  ψ(x)·q(x)` is **equatorward** in the tropics (the low-level moist branch flows toward the ascent; dry-
  upper-branch `Δq≈q_surface`), with ψ the normalized cell (`hadley_streamfunction`) and `q=RH·q_sat(T)`.
  Convergence `P−E=−∂F/∂x` in **conservative face form** (`_mean_flux_convergence`, the MMC analogue of the
  diffusive `_spherical_flux_divergence`) ⟹ `∫(P−E)=0` **machine-exact** (a conserving *budget*, ITCZ
  convergence paid for by the descending dry belt — not a painted band). Ascent **pinned at the equator**
  (hemisphere model; rung 2.x owns *migration*). The cell vanishes poleward of its edge ⟹ the extratropical
  eddy budget is **untouched** (asserted identical with/without).
- **THE HONESTY CLASSIFICATION (advisor, load-bearing — set up front).** Convergence-at-ITCZ / divergence-
  under-descent is **GUARANTEED BY CONSTRUCTION** for any prescribed equatorward flux ⟹ that is **plumbing,
  NOT a benchmark win** (the same "guaranteed result" trap as QG down-gradient+irr~1 and the 2.x warm-ward
  shift-direction). `HADLEY_STRENGTH≈4.2e-4 kg/m²/s` is the named, **prescribed WALL** (calibrated to
  observed *order* ~1–2 m/yr, **NOT** derived; calibrated **transparently** — not tuned-then-cited). **The
  genuinely emergent, non-vacuous nugget = the AMPLITUDE:** `q(T)` is carried from the EBM, so the ITCZ
  convergence **intensifies at the ~C–C moisture rate (~6.6 %/K)** under warming — *faster* than the energy-
  constrained global mean (~2.5 %/K) — the observed **"rich-get-richer" P−E scaling** (Held & Soden 2006).
  That is a prediction, not a prescription, and it is the bankable physics.
- **THE TRADE (advisor warned, the cubic profile REVEALED — record it):** the fix flips the ITCZ **sign**
  robustly but does **NOT relocate the desert**. The emergent dry belt comes out **equatorward (~12°) of the
  canonical 25–35° subtropics**: the hyper-peaked fixed-RH C–C `q` pulls the moisture flux `ψ·q` equatorward
  — the **same** mislocation the eddy budget has (`test_subtropical_evaporative_belt_is_not_reproduced`), so
  25–35° stays `P>E` on **both** paths. **A half-sine ψ first *appeared* to fix 25–35° — but only via an
  edge-discontinuity artifact** (ψ′(edge)≠0 ⟹ a ~210 cm/yr jump in `P−E` at 30°); switching to the **cubic
  `ψ=(27/4)u(1−u)²`** (ψ′(0)>0 strong ITCZ, ψ′(edge)=0 smooth merge) removed the jump **and** exposed that
  the canonical-subtropics flip was artifactual. Relocating the desert needs a realistic (less peaked) `q` =
  moist dynamics / the resolved vertical, **rung 3+** (where the **gross-moist-stability / overturning**
  route — the fully *emergent* cell, not an imposed ψ — is the honest framework; GMS in a column model just
  moves the prescription onto vertical-structure quantities it doesn't have, and deep-tropics is where GMS→0
  — named, NOT built). The pinned subtropical test's NOTE updated to say the mislocation **persists** past
  this fix.
- **Triad.** *Plumbing (by-construction)* — the ψ cell shape (0 at ascent & edge, smooth merge); `∫=0`
  machine-exact; `strength=0` ⟹ eddy-only **bit-for-bit**; tropics-confined ⟹ extratropics identical; the
  ITCZ **sign flip** itself. *Real-but-loose (the emergent unlock)* — the ~C–C-rate intensification (faster
  than the energy-constrained mean); equatorial convergence of observed *order* (~1–2 m/yr, loose band, not
  a tuned match). *Named trade* — the desert is NOT relocated (equatorward dry belt; canonical 25–35° still
  `P>E`). **Demo banked + CI-guarded** (`planet/demo_hadley_moisture.py` → `docs/figures/planet-hadley-
  moisture.png`: before/after `P−E` + the rate-bar; the `slow` `test_demo_reproduces_the_hadley_fix_headline`
  pins sign-flip + dry-belt-equatorward + conservation + emergent-rate). Tests: `planet/tests/test_moist.py`
  (+9 fast Hadley + 1 slow demo guard; full file 29 incl slow). No engine edit; `sphere_ebm.py` breadcrumb
  updated; `precip.py` untouched. Sources extend [[moist-ebm-source]] (Held & Soden 2006 rich-get-richer;
  the diffusive-moist-EBM + mean-circulation moisture budget — Hartmann GPC; Hwang & Frierson 2010).

**Rung 2.x — Emergent ITCZ rain BUILT (the full-sphere moisture budget co-located with the EFE;
`planet/sphere_moist.py`, 2026-06-14).** The "meatiest new finding" slice: rung 2.x's precip wire only
**relocated a prescribed Gaussian band** to the EFE (a *dry* model painting a belt); this rung carries the
rung-2 **column moisture budget** onto the full sphere so the ITCZ rain **emerges from a conserving `P−E`
budget** whose convergence maximum **sits on the EFE** — *rained*, not *painted*. **A SIBLING** —
`moist.py`/`sphere_ebm.py` UNTOUCHED (reuses `moist.specific_humidity`/`q_sat`/`HADLEY_STRENGTH`, re-derives
the conservative operators on the doubled grid). Built **spike-first** (`outputs/rung2x_sphere_moist_spike.py`,
gitignored) + **advisor-pressure-tested twice** — the second round **overturned the advisor's own predicted
headline** (record it, it's load-bearing). The model = full-sphere eddy convergence `(D/c_p)·∂ₓ[(1−x²)∂ₓq]`
+ a **two-cell Hadley** circulation whose ascent is **anchored at the EFE** (`hadley_streamfunction`,
asymmetric widths — descent edges pinned at ±30° ⟹ the physical cross-equatorial cell widens), both in
**conservative face form** with two real polar Neumann-0 ends ⟹ `∫(P−E)=0` **machine-exact** for *any*
asymmetric cell.
- **Banked (tight/structural):** the **cross-model reduction** to the hemisphere `moist.moisture_budget` on
  the NH at `φ_EFE=0` — **machine-exact** (eddy 0.0, full ~4e-11: a symmetric `q` zeroes the equatorial
  face-flux, so the full-sphere stencil collapses to the hemisphere's equatorial-symmetry boundary);
  `∫(P−E)=0` machine-exact (symmetric **and** displaced); symmetric climate ⟹ even `P−E` peaking at the
  equator. The two-cell `Ψ` (0 at the EFE + at the ±edge descents, sign-flips across the EFE).
- **The two real nuggets (advisor's final altitude):** (1) **co-location of the NET `P−E` on the EFE = a
  FALSIFIABLE CHECK, not a given** — the down-gradient eddy term *exports* moisture from the warm EFE (the
  rung-2 ITCZ trade), so the prescribed cell must **beat that export at the displaced latitude**; it does, by
  a **~2.6× margin** (eddy ≈ −113 cm/yr vs Hadley ≈ +287 at the EFE) for the calibrated strength ⟹ the net
  rain max lands on the EFE to **<1°** (checked, not assumed). (2) **The displaced-ITCZ peak intensification
  is GEOMETRIC, NOT emergent `q` — a CLEAN NEGATIVE RESULT** (pinned so it is not silently re-read as a win):
  the peak grows because the pinned-edge near cell **narrows**, not because the warm hemisphere is moister —
  replacing `q(T)` with a hemispherically-symmetric `q` leaves the peak unchanged (180.8 ≈ 180.5 cm/yr), and
  symmetric cell widths remove the intensification entirely. The only clean emergent-`q` signature is the
  **~C–C warming response** (the Hadley-fix nugget, re-confirmed on the full sphere).
- **THE ADVISOR OVERTURN (record it).** The advisor first predicted the headline = "the displaced ITCZ is
  more intense because it sits in the warmer/moister hemisphere; magnitude tracks the interhemispheric
  q-contrast = emergent." A decomposition spike **refuted it**: under the (constant-albedo) Q-flux the
  q-contrast is tiny (~2–11 %), and **neither** the peak **nor** the wet-NH/dry-SH dipole tracks it — the
  peak is geometric (above), the **dipole is displacement-driven** (present at full strength with a
  *symmetric* `q`; its *direction*, toward the warm hemisphere, is by-construction). Surfaced as a
  reconcile-call (the "don't silently switch on primary evidence" discipline); the advisor **conceded** and
  re-pitched the altitude to: architectural win (budget-not-band) = the headline "what"; the two nuggets
  above; **drop** any claim that `q` shapes the meridional asymmetry in either direction (the realq<symq
  "damping" was confounded — symmetrizing also relocates the q-peak). The genuinely-new emergent content is
  **modest** — largely 2.x's EFE displacement × the rung-2 Hadley convergence recombined into a conserving
  budget.
- **Named edges / walls (carried):** the cell is **prescribed** (`HADLEY_STRENGTH` the wall; the fully
  emergent `Ψ∝H/GMS` cell = **rung 3+** — it double-counts the mean transport already lumped into the EBM's
  `D` and needs a GMS closure); anchoring the ascent at the EFE is a **placement**, not a derivation (it is
  what makes "rain co-locates with the EFE" true — by-construction); the asymmetry is **imposed** (Q-flux/
  albedo); the **subtropical desert stays mislocated** (hyper-peaked C–C `q`, the rung-2 wall, unfixed);
  `P−E` not a full `P` (no honest zonal `E` keeps `P≥0` — `moist.py`'s discipline). **Demo banked +
  CI-guarded** (`planet/demo_sphere_moist.py` → `docs/figures/planet-sphere-moist.png`: rained-vs-painted +
  the co-location margin + the conserving dipole/geometric-negative; the `slow`
  `test_demo_reproduces_the_banked_headline` pins co-location + the margin + conservation + dipole signs +
  geometric-not-q). Tests: `planet/tests/test_sphere_moist.py` (10 tight/loose/plumbing fast + 1 slow guard);
  full gate **507 passed, 1 skip**. `moist.py`/`sphere_ebm.py` untouched (sphere_ebm breadcrumb added);
  `precip.py` untouched. Extends [[planet-rung2x-itcz]], [[planet-rung2-hadley-fix]], [[moist-ebm-source]].

**Rung 3 — SCOPED + spike-validated (vertical structure → baroclinic instability; 2026-06-12).** The
biggest jump on the §5 staircase: the **first *structural* edit to a shared engine** since Phase 3 (so it
triggers the full-repo gate + the import-drift guard, ADR 0003) and the **first compute wall**. A
single-layer model *categorically cannot* be baroclinically unstable (no available potential energy, no
vertical shear), so the rung adds the **minimal honest vertical structure: a two-layer free-surface
shallow-water model (Phillips 1954)** — two stacked SW layers of slightly different density, coupled only
through pressure. The table's "multi-level / 3-D" is realized minimally as **2 layers**; **N-layer is the
within-rung upgrade**, not a separate rung (no cleaner anchor above 2, just more cost). The interface
displacement *is* the dynamical temperature/buoyancy — so heat transport becomes **intrinsic to the
dynamics**, the qualitative leap past the passive rung-1 tracer (the CONTRACT already names "an active
buoyancy tracer feeding `h` = a two-layer model"). Scoped advisor-pressure-tested **and de-risked
spike-first** (`outputs/rung3_baroclinic_spike.py`, gitignored) — the spike **answered all three of the
advisor's load-bearing risks** before any build:
- **The tight anchor = the two-layer SW equations' OWN linear stability, NOT QG Phillips (advisor).** Making
  QG Phillips primary would import an approximation gap (ageostrophy, free-surface coupling) into what
  should be the tight leg — so the rung-appropriate move (matching Phase 3's Poincaré/Rossby anchor and rung
  2.5's exact `dq_sat/dT`) is to **linearize the two-layer SW equations and root the 6×6 dispersion matrix
  for the growing mode** (built from the equations, **not a recalled quartic**). The spike's solver is
  **first-principles-validated**: at zero shear it is **neutral to machine precision** (`max|Im ω|≈3e-20`)
  and **recovers the two-layer Poincaré dispersions** (external `ω²=f₀²+gH_tot k²`, internal
  `ω²=f₀²+g'H_e k²`). It gives σ(k) with a **short-wave cutoff**, a most-unstable wavelength **~6.8× the
  layer deformation radius** (`λ*≈680 km` at idealized params). **The baroclinic `G_k` terms are
  externally anchored** (zero-shear Poincaré leaves them untested → otherwise they'd rest only on the
  self-consistent Phillips convergence = two hand-derivations that could share a compensating error): the
  max-growth coefficient `σ_max≈0.304·U_s/L_d` matches the **Eady** model (continuous stratification, a
  *wholly independent* derivation) `0.310` to **2 %**, and the critical-shear *formula* `U_s,crit=β·g'H/f₀²`
  matches the literature `β/k²_int` (Pedlosky/Vallis, web-confirmed). The solver is **f-plane** (β is not in
  the perturbation operator → it is correctly **Eady-like: no critical shear, unstable for all shear**); a
  *finite* critical shear needs β = a Phase-A (β-capable, QG/PV-gradient) item, not a spike claim.
- **QG Phillips = the LOOSE cross-check that confirms the tight leg.** Derived as a 2×2 and compared to the
  SW solver: in the **rigid-lid limit (`g→∞`, the external mode infinitely fast) they converge to <0.5%**
  (`σ_SW/σ_Phillips→1.004`), mutually validating both; at the actual free-surface params the gap is a clean
  **~4%** (SW slightly *more* unstable — the extra free-surface degree of freedom). **The spike caught +
  fixed exactly the recalled-coefficient trap the advisor flagged** (the β=0 cutoff is `K²=2F`,
  `F=f₀²/(g'H)` — *not* `2F`, which was √2 too short and is what first surfaced the gap).
- **The external-mode CFL is the real cost (advisor's #1 risk) — quantified, and AFFORDABLE.** A
  free-surface two-layer model carries a fast **barotropic** gravity wave at `√(gH_tot)` that sets the
  explicit RK3 step while the **slow baroclinic mode** is all we care about. The spike measures the penalty
  `c_ext/c_int≈14×` at idealized params, but **4 e-folds cost only ~3–16 s wall** (nx=32→128) on the C-grid
  route → **GO on the free-surface engine extension**; the **rigid-lid fork** (a barotropic *elliptic solve*
  = a structural change to the explicit engine — "a different animal") is the **named within-rung upgrade**
  if Phase-B's saturated, higher-resolution, many-life-cycle runs make the penalty bite. **Idealized
  `(g, H, Δρ/ρ)`** chosen so `√(gH)` is modest and the internal `L_d` is resolvable — **honest at rung 3**
  (validates the *mechanism + growth rate*, not Earth jet speeds; the same "config-tuned, not Earth-
  calibrated" honesty banked at rungs 1–2).
- **The route is sound: the C-grid engine reproduces the analytic σ.** A linearized two-layer C-grid SSP-RK3
  solver (mirroring `engines/fluid`'s discretization) integrating a single growing mode lands the measured
  growth rate **within ~4 % of the analytic σ, monotonically converging with resolution** (5.2→4.1→3.8 % at
  nx=32→64→128). The stacked-layer advection + coupled-pressure + Coriolis + SSP-RK3 all behave.

**Calibrate the "GO" (advisor done-check — the hedge is load-bearing).** The spike de-risks **buildability**
at **high confidence** — linear baroclinic instability exists, the engine route reproduces it, the CFL is
affordable — but that was the **textbook-guaranteed** part (70-year-old physics + an implementation/CFL
check). **The Phase-B *quantitative* payoff remains an OPEN BET:** that the *saturated* eddy flux comes out
**irreversible** (vs rung-1's ~90 %-reversible barotropic finding) at a **realistic magnitude** (near rung-0's
`κ~2.2e6`, vs rung-1's ~1000× off) — i.e. that the reduction-to-EBM is genuinely **non-vacuous** — is
**untouched by a linear spike** and is *exactly* the claim class this project **downgraded at rung 1**
("named, not banked") and **overturned at rung 2** (the extratropical-only trade). Phase A succeeding moves
that risk **zero**. So "spike-validated / GO" means *the route is sound to build*, **not** *the rung is
de-risked* — the saturation reality-check is owed before the payoff is banked.

**The A/B split (advisor — the same machinery-vs-emergent split as rungs 1 & 2):**
- **Phase A = the linear baroclinic growth rate (tight, the current energy-conserving engine).** Extend
  `engines/fluid` to **N layers** (the CONTRACT's promised leading-axis seam): layers couple *only* through
  the Montgomery/Bernoulli pressure term (`P₁=g(h₁+h₂)`, `P₂=g(h₁+h₂)+g'h₂`); the vector-invariant momentum
  generalizes per-layer; **single-layer `tracer=None` stays bit-for-bit** (the by-construction rung-0
  reduction held every rung). Anchor: a small perturbation on a thermal-wind-balanced two-layer basic state
  grows at the **analytic two-layer SW linear rate** (the spike's eigenvalue solver), within a few %. Clean
  on the current engine (small amplitude, early time, before saturation).
- **Phase B = the nonlinear saturated eddy field → THE HEADLINE PAYOFF (loose magnitude; needs
  dissipation).** This is the staircase's long-promised closure: rung 1 found the eddy-flux→EBM "tight
  reduction" **near-vacuous** (barotropic eddies ~1000× too weak, ~90 % reversible) and named that it
  becomes "**non-vacuous only at rung 3 (a strong baroclinic flux)**." So Phase B integrates the unstable
  two-layer flow to saturation and measures the **emergent, irreversible, down-gradient baroclinic eddy heat
  flux** `⟨v'·(interface)'⟩` — now genuine — and re-runs the rung-1 reduction-to-EBM, expecting a
  **non-vacuous `D_eff` near a realistic magnitude**. The energy-conserving engine has **no hyperviscosity /
  no limiter** (the CONTRACT warns turbulent runs cascade enstrophy to the grid), so Phase B needs **(a)** a
  long post-saturation run, **(b)** the named-but-unbuilt **hyperviscosity** dissipation operator, and
  **(c)** the external-mode cost over many life cycles. Phase A's engine extension is designed to **leave
  room for the Phase-B dissipation operator**.

**Anchor classification (the triad, projected):** *tight* — the two-layer SW linear growth rate (eigenvalue
solver, Poincaré-validated) reproduced by the engine within a few %; the single-layer bit-for-bit
reduction. *real-but-loose (the unlock)* — the **emergent baroclinic eddy heat flux** + the now-non-vacuous
reduction-to-EBM (Phase B; magnitude loose / config-tuned). *loose cross-check* — QG Phillips (most-unstable
λ ~ L_d, cutoff, critical shear), rigid-lid-convergence-validated. **Held–Suarez is DEFERRED, not the
rung-3 anchor** — it is a **3-D sphere primitive-equation** benchmark (Newtonian relaxation + Rayleigh
friction → a statistically-steady storm track) that belongs at **rung 5** (the idealized GCM); naming it at
rung 3 would over-reach a 2-layer β-plane laptop model. **Named walls / deferrals:** the **hyperviscosity**
dissipation operator (Phase B's prerequisite; would break the energy-conserving symmetry — a deliberate
CONTRACT non-goal now forced); the **rigid-lid / external-mode fork**; **rigid channel walls in y** (the
classic baroclinic-lifecycle geometry — the spike used the doubly-periodic `l=0` mode to sidestep them, but
the saturated nonlinear field may want a channel = the named BC extension); the **sphere** (rung 5).
**Sources to pin at build** (the `[[…-source]]` discipline — extending `[[shallow-water-source]]`):
baroclinic instability → **Phillips 1954 / Eady 1949 / Charney 1947**; the two-layer SW formulation →
**Vallis 2017 (*AOFD*) / Cushman-Roisin & Beckers**; Held–Suarez (named-deferred) → **Held & Suarez 1994**.
**Rung 3 Phase A — BUILT (2026-06-12; the linear growth rate).** `engines/fluid/layered.py`
(`LayeredShallowWater` + `LayeredState` + `ThermalWindBackground`) + `engines/fluid/stability.py`
(`TwoLayerStability`), 22 tests: the unstable two-layer mode grows at the analytic rate within ~4 %,
converging with resolution; the single-layer `nl=1` reduction is **byte-identical** to `ShallowWater`;
the basic state enters as constant background coefficients (the periodicity resolution). First
*structural* `engines/fluid` edit → the full-repo gate. See `[[planet-rung3-scoped]]`.

**Rung 3 Phase B — SPIKE FINDING: the free-surface SW engine OUTCROPS at saturation; the payoff needs
two-layer QG (2026-06-13).** Phase B was de-risked **spike-first**
(`outputs/rung3_baroclinic_phaseB_spike.py`, gitignored) *before* any production code — and the spike
returned a clear **negative result** that re-routes the build (spike-first doing exactly its job). The
dissipation the plan named was built in the spike (**hyperviscosity on momentum *and* thickness +
linear bottom Ekman drag**, both default-off → the Phase-A engine bit-for-bit; per-layer mass stays
machine-exact with dissipation on), and the advisor's under-specification was confirmed: **a fixed-`G`
background is an infinite APE reservoir, so hyperviscosity (a small-scale enstrophy sink) cannot
saturate it — linear bottom drag (Held–Larichev 1996) is required to arrest the inverse cascade.** But
the deeper blocker is the **free surface**: at the idealized parameters the saturated/overshoot
interface displacement reaches the layer depth and the layer **outcrops** (`h→0` ⟹ `PV=(f+ζ)/h`
detonates). The control parameter is the **Froude ratio `Fr = U_s/√(g'H) ≈ η_sat/H`**; empirically the
first-saturation **overshoot** drives the peak displacement to **`η/H ≈ 12·Fr`** (far above the ~4–5×
RMS vortex-tail factor — the first baroclinic life cycle dumps mean APE into eddies before
equilibrating). The outcropping is **robust** across `g'∈{0.2,0.8,1,2}`, `H∈{400,500}`, `U_s∈{2,4}`,
drag `r∈{0.5…3σ}`: stronger drag only **delays** the overshoot (`r=3σ` still `η/H=0.91`); stronger
stratification gives a **bigger** overshoot (more stored APE). Avoiding it needs `Fr≲0.04` ⟹ `U_s≲1` ⟹
e-fold `≳370 h` — full **QG-regime cost in the tool worst-suited to it**. The linear growth rate is
reproduced at every config tried (growth-err 2.2–2.6 %), so this is **not** leaving the validated regime
— it is the **finite-amplitude free-surface wall**. **Interpretation:** this is the empirical reason
two-layer turbulence is done in **QG** (Held–Larichev 1996 *is* two-layer QG); the free-surface explicit
SW model has a hard thickness floor QG does not, so **the saturated payoff does not live in this
engine**. The plan **pre-named** the rigid-lid/QG fork as the within-rung upgrade "if the saturated runs
make the penalty bite" — they bit, **via outcropping** (not just the CFL). The pre-outcrop
finite-amplitude flux was considered and **rejected as the headline** (near-vacuous: an unstable mode
fluxes down-gradient *by definition* and `irr~1` is trivial without sloshing to overcome — no meaningful
contrast with rung-1's reversible barotropic finding). **Status:** Phase A banked; **Phase B remains the
OPEN BET**, now with a sharpened obstacle — it requires the **two-layer QG (rigid-lid) solver**: an
FFT-based 2×2 spectral PV-inversion model (~200 lines, canonical Held–Larichev; the linear anchor is the
spike's QG-Phillips cross-check), a **new model outside `engines/fluid`** (no C-grid reuse, no
bit-for-bit reduction — Phase A and B would validate *different* models). QG makes the experiment
**possible**; it does **not** pre-guarantee the κ comes out irreversible / well-scaled — that stays the
open bet, finally testable. See `[[planet-rung3-phaseB-outcropping]]`.

**Rung 3 Phase B — next build: two-layer QG, SCOPED (2026-06-13; user picked "bank + scope," not "build
now").** The Phase-B payoff is re-routed to the **canonical Held–Larichev 1996 doubly-periodic two-layer
QG turbulence model** — the standard tool for exactly this experiment. This is the **plan of record for the
next increment**; no code shipped this session (the spike finding above is what was banked).
- **The model.** Two layers of QG potential vorticity `q_k = ∇²ψ_k + (−1)^k F_k(ψ_1−ψ_2) + βy`,
  `F_k = f₀²/(g'H_k)`, on a doubly-periodic plane. Prognostic = the PV anomaly advected by its own flow;
  **ψ is recovered from q by a 2×2 *spectral* inversion** — at each wavenumber the coupled
  `(−K²−F_k)ψ̂_k + F_kψ̂_{3−k} = q̂_k` is a 2×2 solve (FFT-diagonalized, vectorized). Pseudospectral
  Jacobian advection with the 2/3 dealias rule; RK3/4 in time. ~200 lines. The fixed thermal-wind shear
  `(U_1−U_2)` is the background APE reservoir (as in the SW engine's `background`); **β returns** (β-plane
  QG, vs the SW f-plane) — it sets the **Rhines arrest** of the inverse cascade and restores a **finite
  critical shear** `U_crit=β/F` (the named β-item from Phase A, natural in QG).
- **Why it works where SW failed.** QG is **linearized in the thickness** → the interface displacement has
  **no layer-depth floor**, so saturation is **well-posed (no outcropping)**; and the **rigid lid filters
  the fast barotropic gravity wave** → **no external-mode CFL** (timestep set by the slow advective speed),
  so long post-saturation runs are cheap. Exactly the two stiffnesses the free-surface model paid.
- **Reuse + the dissipation already de-risked.** Port the spike's **hyperviscosity `−ν₄∇⁴q` + linear
  bottom Ekman drag `−r∇²ψ_2`** (built, mass/PV-clean, CFL-aware) — both default-off. Reuse the **κ→D
  bridge** (`planet/transport.py`) and the **pre-registered discriminators** unchanged: emergent thickness
  diffusivity `κ = −⟨v'τ'⟩/(dτ̄/dy)` (`τ ∝ ψ_1−ψ_2` the interface/buoyancy), the **homogeneous-turbulence
  κ** (the clean GM/thickness-diffusivity definition — **no operator-shape test**, that needs a meridional
  channel = the named BC extension), **irr fraction O(1)** (vs rung-1 ~0.1) and **κ_eff/(v'_rms·L_d)
  O(0.1–10)** (vs rung-1 ~1e-3), **validated dimensionless** (idealized `κ_ML` is intrinsically 15–60×
  below Earth's `κ₀=2.2e6`).
- **The triad (projected).** *tight* — the QG **linear stability = the spike's QG-Phillips cross-check**
  (most-unstable `λ~L_d`, the `K²=2F` short-wave cutoff, `U_crit=β/F`); the spectral inversion is exact;
  zero-shear → neutral. *real-but-loose (the unlock)* — the emergent **saturated, irreversible,
  down-gradient eddy thickness diffusivity** → the **now-non-vacuous reduction-to-EBM** (direction +
  irreversibility banked; magnitude dimensionless / config-tuned). *plumbing* — zero shear ⟹ no eddies
  (decay); `q↔ψ` round-trips.
- **Day-one build order + the K=0 trap (advisor de-risk).** **Sequence the linear anchor FIRST:** before
  any turbulence run, reproduce the linear growth rate against the spike's QG-Phillips cross-check — it
  catches the two things most likely wrong on day one cheaply (the **mean-PV-gradient sign**
  `∂q̄_k/∂y = β ∓ F_k U_s`, a convention-dependent flip, and the inversion). **Pin the K=0 inversion
  singularity:** the 2×2 determinant is `K²(K²+F_1+F_2)` — well-conditioned for all `K>0` but **zero at
  `K=0`**, so the **domain-mean ψ is undetermined from q** and must be handled explicitly (set domain-mean
  `ψ=0`, or carry the barotropic domain-mean as a separately drag-damped quantity) — the #1 day-one
  spectral-QG trap.
- **Honest edges (named).** A **new model OUTSIDE `engines/fluid`** (pseudospectral, not the C-grid) →
  **no bit-for-bit single-layer reduction; Phase A and Phase B validate *different* models** — the bridge
  between them is the **shared two-layer linear instability** (SW solver ↔ QG agree to <0.5 % in the
  rigid-lid limit, per the spike). Home is a **scope decision for next session**: a planet module
  (`planet/baroclinic_qg.py`, single-consumer) **vs** a third shared engine (`engines/spectral/`,
  reusable) — lean planet-module first (rule-of-three not yet met). **QG makes the experiment *possible*,
  not pre-guaranteed:** the saturated κ coming out irreversible / well-scaled is **still the open bet**
  (the class downgraded@rung1 / overturned@rung2), now finally *testable*. Held–Suarez (sphere
  primitive-eq) stays **rung 5**. Sources to pin at build: **Held & Larichev 1996** (two-layer QG
  turbulence + the diffusivity scaling), Phillips 1954 / Eady 1949, Vallis 2017 *AOFD*; extends
  `[[shallow-water-source]]`. See `[[planet-rung3-phaseB-outcropping]]`, `[[planet-rung3-scoped]]`.

**Rung 3 Phase B — BUILT 2026-06-13: two-layer QG turbulence, and THE OPEN BET IS WON.** The
re-routed Phase-B engine is built (`planet/baroclinic_qg.py` — the **home scope-decision settled to the
lean single-consumer planet module**, not `engines/spectral/`; rule-of-three unmet) and the saturated,
irreversible, down-gradient baroclinic eddy thickness flux comes out at an **order-unity dimensionless
mixing efficiency** → the rung-1 reduction-to-diffusive-EBM is **finally non-vacuous**. This is the claim
class downgraded@rung1 (barotropic flux ~1000× too weak, ~90 % reversible) and overturned@rung2 — *won*
at rung 3, as the staircase predicted. Tests: `planet/tests/test_baroclinic_qg.py` (tight + plumbing
fast; the linear-reduction + saturated-flux checks slow). The advisor done-check caught that the flux
test pinned only **necessary-not-sufficient** conditions (κ>0 + irr~1 hold for *any* sustained baroclinic
state) → added the **sufficient** assertion `k_peak < k*` (the inverse-cascade spectral peak, via
`TwoLayerQG.ke_spectrum`) so CI guards *turbulence*, and made the evidence **durable**: a committed demo
(`planet/demo_baroclinic_qg.py` → `docs/figures/planet-baroclinic-qg-turbulence.png` — the PV-vortex
fields + inverse-cascade KE spectrum) reproducible from a fresh clone, not a gitignored spike. **Full
pytest (all slow, `-n auto`): 386 passed, 1 skipped.**
- **The engine.** A doubly-periodic two-layer QG model: the PV anomaly `q_k` (a `QGState`, leading layer
  axis) advected by its own flow, **ψ recovered by a 2×2 *spectral* PV inversion** (`det = K²(K²+F₁+F₂)`,
  the `K=0` domain-mean gauge set to `ψ=0`), a **pseudospectral 2/3-dealiased Jacobian** for the
  eddy–eddy nonlinearity, **SSP-RK3** (mirroring `engines/fluid`), the mean shear `(U_k, ∂q̄_k/∂y)` as
  **background coefficients** (eddies = the prognostic fields, same design as the SW engine's
  `background`), **β-plane**, and the spike's **hyperviscosity `−ν₄∇⁴q` + lower-layer bottom Ekman drag
  `−r∇²ψ_2`**, all **default-off**. ~430 lines incl. docstrings.
- **Why QG won where the free-surface SW engine outcropped.** QG is linearized in the thickness → the
  interface has **no layer-depth floor** (saturation well-posed, **never outcrops** — pinned by a test),
  and the rigid lid **filters the fast external gravity wave** (no external-mode CFL — the advective step
  is far larger), so the long saturated runs are affordable. Both stiffnesses the SW engine paid are gone.
- **The tight leg (linear anchor first — advisor's day-one sequencing).** The rooted 2×2 QG dispersion
  equals the **analytic Phillips closed form to 2e-15** (equal layers — *the same equations*; advisor:
  the bar is ~1e-9 **not** "a few %" — the ~4 % I'd carried forward was the free-surface-SW-vs-Phillips
  gap, which does not exist here, and a loose tolerance would hide a partially-compensated PV-gradient
  sign bug). The `K²=2F` short-wave cutoff, **neutral at zero shear** (machine ε), and — **β's re-entry,
  the thing the f-plane SW solver could not test** — a finite **Charney–Stern critical shear `U_crit=β/F`**
  (sub-critical neutral, super-critical growing; the lower-layer mean PV gradient `β−F·U_s` reverses sign).
  The **cross-model bridge** to Phase A: the SW 6×6 solver (`TwoLayerStability`) in the rigid-lid limit
  (`g→∞`) converges to the QG rate, **σ_SW/σ_QG → 1 to <0.5 %** — *asserted*, since it is the only tie
  between the Phase-A (SW) and Phase-B (QG) models (no bit-for-bit reduction across the model boundary).
  The full **nonlinear** engine reduces to this linear operator: a single growing eigenmode grows at the
  analytic σ to **0.1 %**.
- **The WIN (the open bet) — and what actually decides it (advisor, load-bearing).** A drag sweep
  `r∈{0.5,1,2}σ` at the saturated state gives `κ>0` (down-gradient), `irr=|⟨F⟩_t|/⟨|F|⟩_t = 0.96–1.00`
  (vs rung-1 ~0.1), and **`κ/(v'_rms·L_d) = 0.71–1.27`** (vs rung-1 ~1e-3) — robust across drag. **BUT the
  advisor's catch is the headline:** down-gradient + `irr≈1` are *guaranteed* for **any** sustained
  baroclinic state (the flux *is* the APE→EKE conversion that powers the eddies; sign-pinned + spatially
  averaged ⟹ `irr≈1` automatically) — they are **necessary, not sufficient**, and cannot distinguish
  developed turbulence from a quasi-steady wave (exactly the **P2 "irr~1 trivial *without sloshing*"**
  rejection carried from the spike). So the bet is won **only by showing genuine turbulent mixing**, which
  the **weak-drag (`r=0.5σ`) condensate candidate** (`v'_rms≈16 ≫ U_s=4`) does on three independent
  diagnostics: **EKE(t) irregular + drifting** (`std/mean≈0.25`, not flat), an **isotropic KE spectrum
  that is an inverse-cascade condensate** (84 % of the energy below the injection band, the peak migrated
  to the box scale `0.33k*`, broadband-continuous with a clean dissipation tail — *not* spikes at
  `k*,2k*,3k*`), and a **PV snapshot of coherent vortices + rolled-up filaments across scales** (not a
  wave train). The dimensionless `κ/(v'L_d)~O(1)` is the discriminating quantity; the qualitative
  inverse-cascade is what makes it turbulence, not a wave.
- **Honest edges (advisor, banked — do NOT overclaim).** (1) The banked claim is **dimensionless +
  qualitative**; the **dimensional `κ ≈ 0.7–4×10⁶ m²/s`** happens to land in Earth's observed band
  (1–5×10⁶) but that is **coincidental and box/drag-dependent** (varies 5× across the sweep) — *not* a
  reproduction of Earth's κ. (2) The large-scale **condensate dominates `v'_rms`**, so read
  `κ/(v'L_d)~1` as "mixing length ~ L_d-ish," not a precise efficiency. (3) **`κ₁=κ₂` is an estimator
  identity** (`F₁−F₂=⟨½∂_x(τ²)⟩=0` on a periodic domain), *not* an independent cross-layer check (it was
  for the SW spike's separate per-layer thicknesses). (4) A **new model outside `engines/fluid`**
  (pseudospectral, not the C-grid) → no bit-for-bit reduction; Phase A and B validate *different* models,
  bridged by the <0.5 % rigid-lid linear cross-check. (5) **Homogeneous box → a domain-bulk κ**, not
  `κ(y)`; a meridional channel (the operator-shape test) is the named BC extension. (6) Resolution
  **nx=96 ≈ 3.5 pts/L_d** is marginal at the deformation scale (the condensate lives at well-resolved
  large scales; the spectrum's clean dissipation tail shows no grid-scale pileup, and the dimensionless
  ratio is already robust across the *drag* sweep); the **nx=128 firm-up (4.7 pts/L_d) holds the
  dimensionless ratio** (`κ/(v'L_d)=1.10`, `irr=1.00`) while the *dimensional* κ shifts ~45 % with
  resolution (`1.44→2.10×10⁶`) — exactly why only the dimensionless ratio is banked. Held–Suarez
  (sphere primitive-eq) stays **rung 5**. Sources pinned: **Held &
  Larichev 1996**, Phillips 1954 / Eady 1949, Vallis 2017 *AOFD*; extends `[[shallow-water-source]]`. See
  `[[planet-rung3-qg-built]]`, `[[planet-rung3-phaseB-outcropping]]`, `[[planet-rung3-scoped]]`.

**Visualization rungs A/B/C — DECIDED to build all three (animated eddy flow; 2026-06-11;
rungs A+B+C BUILT — A 2026-06-11, B 2026-06-12, C 2026-06-13; build detail in §9.5).** A forward decision (user): animate the emergent eddy life cycle across three
rising-cost **visualization** rungs (distinct from the §5 GCM staircase) — **A** a matplotlib
two-panel mechanism animation (the repo's first time-animation primitive, finally the §9.4
rule-of-three third consumer), **B** a Plotly-globe animation (existing `[webviz]` stack), and
**C** a WebGL particle globe (`mapbox/webgl-wind` ISC + `cambecc/earth` MIT — the
Ventusky / *Perpetual-Ocean* showcase). Shared prerequisite: bank `(h,u,v,θ)` frames from the
`eddy_flux` release loop as a **diagnostic-pure** opt-in side-channel (`n_frames=0`
bit-for-bit unchanged — the inert-seam discipline). Honesty edges carried through all three
(hardest at C): the flow is a **doubly-periodic midlatitude band**, not a global field, and
the instantaneous flux is **~90 % reversible** — so rung A's second panel (cumulative `∫F̄dt`
vs `|F̄|` throughput) makes the small-net-residual finding *visible* rather than letting a
stirring movie overclaim transport. Recommended build order is **A first** (honest by
construction, in-repo, CI-testable), then judge B-vs-C after seeing real frames move. Full
roadmap, the rungs table, and the named scope edges: **§9.5**.

**Rung C — build approach decided, then amended to a true 3-D sphere; BUILT 2026-06-13 (full record in §9.5).** A
forward decision (user): build the showcase as a **general-purpose flow-on-a-globe renderer** aimed at
*one day* visualizing a full **GCM / ESM** field — rendering lesser models (today, the one eddy band) in
the meantime. The renderer was first locked as an *original Canvas2D orthographic* globe, then **amended the
same day (user) to a true 3-D sphere — a `three.js` / WebGL perspective globe with particles streaming on a
real, rotatable sphere.** The pivot **touches only the renderer**: the **renderer-agnostic data contract**
(grid + `(u,v)` + scalar + frames + **coverage-extent** + **provenance/honesty label**, nothing about
projection/particles) was built to absorb exactly this, so the contract, the honesty carve-out, §9.3, and
the ADR 0002 note are unchanged — the locality of the change is the proof the boundary was right. **B and C
are now both 3-D globes; their *roles* separate them** — B is the faithful **scalar field** on the one true
band (honest-by-construction), C the **immersive particle-streaming** showcase (honest-by-disclosure). The
choice **reverses the §6 stance**: three.js is a *library* (vendored, not reimplemented), so its MIT licence
now **owes a `NOTICE`/attribution file** (a named build deliverable), and it is **vendored inline** (not
CDN) to keep the artifact offline-self-contained off `file://`; the original Canvas2D globe demotes to a
lighter §6-free fallback. The **GPU-advection seam shifts down a level** (CPU `BufferGeometry` advection v1
→ GPU ping-pong upgrade *within* WebGL; same <~30 fps @ GCM-resolution trigger). The one **policy change**
is unchanged by the amendment: rungs A/B stay *honest-by-construction*; **Rung C is *honest-by-disclosure***
(user, 2026-06-12: illustrate freely *"if documented … currents carry heat, when they do not"*).
Asymmetric verification: **physics-fidelity relaxes** (approximate is fine, no byte-golden — the figure was
never in the correctness path, ADR 0002 #2) but **documentation tightens** — a test machine-checks the
on-screen disclaimer (now a DOM overlay over the WebGL canvas) is present, *because the disclaimer is the
entire license*. Carve-out is narrow (showcase renderer only; science + A/B untouched). Doctrine recorded in
**§9.5** + **§9.3** + ADR 0002 status note; build deliverables listed in §9.5. **BUILT 2026-06-13 exactly as
locked** (`planet/flow_globe.py` generic renderer + `FlowField` contract, `demo_eddy_particles.py`,
`docs/figures/planet-eddy-particles.html`, `NOTICE` with three.js' full MIT body, `planet/vendor/three.min.js`
r137 UMD, structural + disclaimer tests; the browser play-through handed to the user to eyeball) — see §9.5.

**Rung 4 — BUILT (gray radiative transfer: where the prescribed OLR ``A + B·T`` comes from;
`planet/radiation.py`, 2026-06-14).** The named wall every rung from 2.5 on cites — **``B`` held fixed** —
retired: the linear OLR is made **emergent** by a **gray radiative–convective column** (Schwarzschild
two-stream over an optical depth set by the greenhouse-gas amount), a **sibling** module (``ebm.py``
untouched, the `moist_ebm`/`sphere_ebm`/`baroclinic_qg` discipline). Built **spike-first**
(`outputs/rung4_radiation_spike.py`, gitignored) + advisor-pressure-tested **before** the build, which
**sharpened the headline into a decomposition** (below). Sources pinned at build: the gray-RE closed form
+ the two-stream/Eddington closure → **Pierrehumbert *PoPC* §4 / Goody & Yung**; the present operating
point (OLR 239, Ts 288, the 33 K greenhouse) → **Trenberth–Fasullo–Kiehl 2009**; the feedback orders the
emergent slopes are validated against (Planck/λ_wv/λ_LR) → **Soden & Held 2006** (read off the paper's own
text, local PDF); the logarithmic-CO₂ contrast → **Myhre+ 1998**.
- **The tight anchor (DERIVED, not recalled — the rung-3 ``K²=2F`` lesson).** The gray-RE closed form is
  derived in-module from the two-stream equations: ``σT⁴(τ)=½σTe⁴(1+τ)``, skin ``Te/2^¼``, **ground
  ``σTg⁴=½σTe⁴(2+τ_s)``** (the surface–air discontinuity, *exactly* where a recalled coefficient goes
  wrong). A numerical two-stream RE solver (`solve_gray_equilibrium`, an independent radiosity relaxation
  using *no* analytic input) reproduces the profile + skin + **ground** to **~2nd order** in layer
  thickness with **``OLR=σTe⁴`` machine-exact** (energy conservation) — the FV-engine/C-grid "reproduces
  analytic" pattern; the ``Tg`` → derived-coefficient convergence is the recalled-coefficient guard.
- **THE HEADLINE — a DECOMPOSITION, not a trade (advisor's load-bearing reframe).** Present-day energy
  balance *forces* ``OLR≈ASR≈239``; calibrating ``τ_s`` to the 33 K greenhouse makes the emergent OLR pass
  through that operating point **by construction**, so gray and ``A+B·T`` agree in *value* at present and
  **the slope ``B`` is the finding** (``A,B`` linked through the point, ``A=239−B·T̄`` — the *point* is
  forced, the slope is open, not both recovered independently). The slope **decomposes**: **no-WV
  ``B≈3.41``** sits at the ``4σTe³≈3.75`` **emission-level Planck touchstone** (the advisor's tight-ish bug
  guard: ≈3.8 validates the radiative core, ~1.5 = a bug) — **above climlab's 2**; turning **water vapour**
  on (``τ(Ts)`` via Clausius–Clapeyron, reusing `moist.saturation_specific_humidity`) lifts the emission
  level to colder air and **subtracts ≈2** (to ~1.33 at the nominal 50% WV loading). So **climlab's ``B=2``
  ≈ Planck − water-vapour + the lapse-rate feedback the gray column omits** — and every term is
  **ORDER-VALIDATED against Soden & Held 2006, NOT tuned** (advisor's WV-hardening, sources pinned from the
  paper itself via the local PDF): no-WV ``3.41`` ≈ Planck ``|λ₀| 3.1–3.2``; the WV feedback ``2.08`` ≈
  ``λ_wv 1.8`` (clear-sky order); the gray net ``1.33`` ≈ the clear-sky Planck+WV ``3.2−1.8``; and the gap
  from ``1.33`` to climlab's 2 ≈ the lapse-rate feedback ``|λ_LR| 0.84`` a **fixed** lapse rate (uniform
  warming) **cannot produce** (climlab's obs-tuned ``B`` folds it + clouds in). The decomposition is thus
  **non-circular** (the rung-2.5 frozen-``D_eff`` / ITCZ closed-form attribution flavour); the residual
  tuning lives only in the *exact* magnitudes via `WATER_VAPOUR_FRACTION` (the **wall**).
- **The named WALL + edges.** The wall = the **gray (band-independent) absorption** + the prescribed
  **τ↔GHG-column** mapping (calibrated to *order*, not line-by-line — the ``R_ATM_SLOPE``/``HADLEY_STRENGTH``
  cited-closure status). Edges, each pinned: **CO₂ forcing is SATURATING, not logarithmic** — a gray band
  gives a concave ``OLR(τ)`` (per-doubling ``ΔF`` 48→53→41→25→20 W/m², *decreasing*, not the constant-per-
  doubling Myhre log) at an **unrealistic whole-band magnitude** → the log law was **band physics = a
  named within-rung upgrade, now BUILT** (`SpectralCO2Band`: exponential band wings → constant
  per-doubling `ΔF` ~4.5 W/m² in the Myhre band; form banked, magnitude the wall; reduction-to-gray the
  machine-precision anchor — see §12.2 + [[planet-rung4-radiation]]); **clear-sky only** (clouds out of scope); **no
  lapse-rate feedback** — the fixed convective Γ means warming is a uniform profile shift (zero LR
  feedback), which is *exactly* why the gray net sits below climlab's 2 by ``≈λ_LR`` (a moist-adiabatic Γ
  is the named upgrade that would supply it); **single column** — wiring ``OLR(Ts,τ)`` *per-latitude* into
  ``ebm.py`` (real radiation **driving** the climate as an opt-in sibling EBM) is the natural rung-4
  completion, **left to a user call, NOT foreclosed** — deferred not because the feedback is wrong (at the
  climlab-matched loading the global-mean ``B`` is 2 exactly) but because the emergent ``OLR(Ts)`` is
  **nonlinear**, so the per-latitude slope differs (cold pole vs warm equator) and the wire re-opens the
  meridional profile (a **feature** — emergent latitudinal radiative structure — as much as a risk);
  **linearization breaks far from present** (a steep WV loading → a Komabayashi–Ingersoll runaway, the hot
  analogue of the snowball — ``B`` not linearized across it). **Reduction:** near present ``OLR(Ts)`` is **locally affine** (rung-0's line is its tangent, residual
  <0.5 W/m² over ±3 K; the wide-range curvature *is* the feedback), and feeding climlab's ``B=2`` through
  the forced point recovers climlab's ``A≈210`` (the ``A/B`` linkage, asserted vs ``ebm.A_OLR``). **Triad:**
  *tight* — the numerical solver↔analytic gray RE (2nd-order, OLR machine-exact, ground coefficient) + the
  ``4σTe³`` Planck touchstone; *real-but-loose (the unlock)* — the ``B=Planck−WV`` decomposition (direction
  banked, magnitude on the WV-loading wall); *plumbing* — present operating point matched by construction +
  local affinity + the climlab-``(A,B)`` consistency. **Demo banked + CI-guarded** (`planet/demo_radiation.py`
  → `docs/figures/planet-radiation.png`: emergent OLR + the decomposition **waterfall** (Planck − WV +
  lapse-rate = 2.17 ≈ climlab's 2) + the saturating forcing; the `slow` `test_demo_reproduces_the_radiation_headline`
  pins Te, the operating point, the decomposition bracket, the WV-loading recovery, and saturation, and the
  fast `test_the_decomposition_is_order_validated_against_soden_held` pins the Soden–Held orders). Tests:
  `planet/tests/test_radiation.py` (14 fast + 1 slow). No engine edit; `uses` unchanged. **Rung 4 core
  COMPLETE** (the spectral-band log law + a moist-adiabatic lapse-rate feedback + clouds = named within-rung
  upgrades; the per-latitude EBM wire is now BUILT, below).

**Rung 4 completion — the per-latitude EBM wire BUILT** (2026-06-14, `planet/radiative_ebm.py`, 11 fast +
  2 demo tests; the natural rung-4 completion, was *[left to user call]*). Wires the emergent gray
  ``OLR(Ts,τ)`` into the EBM **per latitude** so *real radiation drives the climate* — a **separate sibling
  alongside rung-0** (`ebm.py` untouched, the `moist_ebm`/`sphere_ebm`/`baroclinic_qg` discipline). Built
  **spike-first** (`outputs/rung4_radiative_ebm_spike.py`, gitignored) + advisor-pressure-tested. **§12's
  scoping guess was OVERTURNED by the spike: the amplification is TROPICAL, not "radiative polar."** The
  headline = the OLR slope is **not one number**: its *local* value ``B_loc(Ts)=dOLR/dTs``
  (:func:`local_radiative_slope`) collapses to ~1.0 at the warm equator (water-vapour feedback) and rises to
  ~2.4 by the cold pole — so under a uniform forcing warming **concentrates in the tropics** (endpoint
  ``δT(pole)/δT(equator)≈0.68``, band 0.73), the **mirror image** of rung-2.5's moisture-*transport* polar
  amplification (dt-free ~1.8–2.05). **The SIGN was measured, not assumed** (the advisor's discriminator: smallest
  ``B_loc`` warms most; WV pulls ``B_loc`` *down* at the equator, the Planck ``4σT³`` pulls it *up* — WV wins
  decisively). A clean "two mechanisms pull opposite ways" pair with rung 2.5: the WV *radiative* feedback
  alone favours the tropics; *transport* (+ lapse-rate + ice, out of scope) make Earth's poles amplify.
  **Numerical core = a coupled Newton steady solve** of ``L_T·T + S(1−α) − OLR(T)=0`` (Jacobian
  ``L_T − diag(B_loc)``), the **nonlinear generalisation of ``ebm.steady_linear``'s direct mode** — reusing
  the engine-pinned transport tridiagonal so transport cannot drift. **NOT the Strang relaxation** (two
  findings forced this): (a) the relaxation half-step that is *analytically exact* for rung-0's linear OLR
  carries an O(Δt²) **splitting error that does not vanish at equilibrium** once OLR is nonlinear (the relaxed
  steady ⟨T⟩ drifts with the step, converging onto Newton; **rung-0's *own* relax default carries a sizeable
  contrast error** — 47.9 vs ``steady_linear``'s dt-free 38.4 — so all rung-0 comparisons here use the direct
  reference); (b) near the warm-equator runaway the local half-step goes unstable. **Runaway finding:** the
  per-latitude wire **exposes the local Komabayashi–Ingersoll edge the global column hid** — at rung-4's
  *default* WV loading 0.5 the equatorial column is *past* its local runaway (``B_loc<0`` for ``Ts≳32 °C``,
  no stable local equilibrium), so the wire runs at the **climlab-matched loading**
  (:func:`climlab_matched_column`: the WV fraction ≈ **0.348** giving global-mean ``B=2``, where ``B_loc>0``
  everywhere and Newton converges in ~6 iters as transport stabilises the equator). **Triad:** *tight* — a
  **linear ``olr_fn`` reproduces ``ebm.steady_linear`` bit-for-bit** (4.5e-13; whole gray departure
  attributable to OLR **curvature** alone) + **net-TOA=0 machine-exact** at convergence (conservation);
  *the discriminator (tight, pure column)* — ``B_loc`` minimised at the warmest latitude + the WV feedback
  **flips the Planck ordering**; *unlock (real but loose)* — tropical amplification **at the Earth-calibrated
  loading**, with **both the sign and the magnitude riding the WV loading** (the wall): a drier,
  Planck-dominated planet is **polar** (``amp>1`` below a crossover loading ≈0.15; the bare ``B_noWV`` rises
  with ``Ts`` so the dry equator damps most) — *unlike* rung-2.5's polar direction, which IS robust to its RH
  wall, so the "mirror" holds at Earth's loading but is not a symmetry of equally-robust mechanisms (advisor's
  catch — the single-loading blind spot, now pinned by a low-loading sign-flip test). Part of the ~0.68
  magnitude is **runaway-proximity** (the warmed equator hits ~39 °C / ``B_loc≈0.5``, stable but near the
  hot edge; a test asserts ``B_loc>0`` at the warmed state). **Advisor's Jensen catch (load-bearing): the global mean is NOT pinned at
  ``ΔA/B``** — OLR is concave so ``⟨OLR(T)⟩≠OLR(⟨T⟩)`` and the WV feedback **amplifies the mean response too**
  (``⟨δT⟩=6.99 > ΔA/B_tan=5.94``); the moist-EBM "redistribution around a pinned mean" framing **deliberately
  does not transfer**. Present mean-state = a **Jensen warm shift** (gray ⟨T⟩ ~2 °C above rung-0) with
  **contrast ≈ unchanged** (the loading-matched *average* slope ≈ 2), so the latitudinal signal lives in the
  *warming response*, not the present climate. **``D`` NOT recalibrated** (§12 expected "recalibrate ``D`` as
  rung 2.5 did" as the cost — but the present contrast is already ≈ rung-0's at this loading, so there is
  nothing to recalibrate *for*; ``D`` sets the amplification magnitude not its sign). The clean **null** = the
  present-day *tangent* (uniform ``B``), which warms exactly uniformly (``amp_null=1``), isolating the
  nonlinearity. Demo banked + CI-guarded (`planet/demo_radiative_ebm.py` → `docs/figures/planet-radiative-ebm.png`:
  the ``B_loc(φ)`` discriminator, the Jensen warm shift, the tropical-amplification warming). Scope edges
  (each named): constant albedo / fixed lapse rate / clear-sky (the within-rung upgrades), uniform ``ΔA``
  forcing not ``ΔS₀``. [[planet-rung4-radiation]]; extends [[moist-ebm-source]] [[ebm-radiation-source]].

**Rung 4 — the emergent lapse-rate feedback BUILT** (2026-06-14, `planet/radiation.py`:
  `GrayRadiationColumn(moist_adiabat=True)` + `moist_adiabat_temperature` + `feedback_kernel`; 9 fast + 1
  slow tests; the §12 "moist-adiabatic lapse-rate feedback" slice — the named within-rung upgrade that
  *retires the omitted feedback* gray's fixed ``Γ`` left out). Swaps the constant convective ``Γ`` for a
  **moist adiabat** (:func:`moist_adiabatic_lapse_rate`, derived by its limits — dry ``→ g/c_p ≈ 9.8 K/km``,
  flattening to ``~4`` when warm) that **flattens as it warms**, so surface warming amplifies in the upper
  troposphere (``ΔT_aloft/ΔTs ≈ 2.8`` peak) and ``OLR(Ts)`` steepens — making the lapse-rate feedback
  **emergent** where rung 4 had *imported* it from Soden & Held. A **default-off flag** (`moist_adiabat=False`
  is bit-for-bit the rung-4-core column). Spike-first (`outputs/rung4_lapse_rate_spike.py`, gitignored) +
  advisor-pressure-tested before + after.

  **THE §12 SCOPED ANCHOR WAS OVERTURNED** (the rung-4-wire "radiative→tropical OVERTURNED" pattern). §12
  scoped it as "supplies ``λ_LR≈0.84``, closing the gap from the gray net ``1.33`` up to climlab's ``2``." It
  does **not**: the emergent value is **``≈+1.5``** and the moist-adiabat column **OVERSHOOTS** — its
  with-water-vapour ``B≈3.1`` sits *above* climlab's 2, not at it. **Banked tight:** the *sign* (``λ_LR>0`` —
  it adds to ``B``), the *kind* (upper-troposphere amplification, measured not assumed via ``ΔT_aloft/ΔTs>1``),
  the **kernel closure** (advisor's load-bearing design: a one-column Soden–Held split — Planck = uniform
  warming, LR = the profile's *departure* from uniform, WV = the ``τ(Ts)`` change — sums to ``B_total`` at
  ~9e-4, *cleaner* than a two-column difference which conflates LR with the Planck-base shift), and
  **resolution convergence** (``λ_LR`` 1.5105→1.5132, n=100→800). **Magnitude LOOSE for two named reasons**
  (advisor): (a) a single *global* moist-adiabat column applies the **tropical** mechanism everywhere, missing
  the extratropical branch (bottom-heavy warming, opposing sign) that pulls the global mean down to ``0.84`` —
  so the column recovers the *tropical* feedback, not the global mean; and (b) it rides the prescribed vertical
  ``τ`` shape + :data:`WATER_VAPOUR_FRACTION` (the rung-4 wall), which set where the emission level sits — the
  same loading the column's *null* is not perfectly clean about (fixed ``Γ`` itself shows a small ``≈−0.25``
  **tropopause-migration residual**, not a true lapse-rate feedback). **Reconciliation with the existing
  decomposition** (mandatory, advisor — else the record self-contradicts): the docstring's ``λ_LR≈0.84`` is a
  **global-mean touchstone** for what the *fixed*-``Γ`` default omits; this emergent single-column value is
  the **tropical branch** — both true about different things. **Two-column cross-check demoted** to a
  consistency check: the ``B_WV(moist)−B_WV(fixed)=+1.79`` reconciles as ``ΔLR(1.77) + ΔPlanck(−0.01) +
  ΔWV(+0.03)`` where ``ΔLR = λ_LR(moist 1.51) − λ_LR(fixed −0.25)`` — proving the kernel cleanly isolated LR
  (note ``B_noWV`` is **Planck+LR**, not the kernel Planck, the subtlety that makes the two numbers differ).
  **Honesty edge (advisor):** the clean WV/LR *separation* is partly a model artifact — here ``τ_wv`` tracks
  the **surface** ``Ts``, not the profile, so the upper-troposphere moisture–temperature coupling that links
  the two feedbacks in reality is absent. The remaining within-rung upgrades (a moist adiabat *with*
  latitudinal structure on the per-latitude wire, recovering the extratropical branch and the global mean;
  the spectral-band log law; clouds) are named. Demo `planet/demo_lapse_rate.py` →
  `docs/figures/planet-lapse-rate.png` (the warming-amplification profile, the emergent kernel waterfall,
  the overturn bars). [[planet-rung4-radiation]]; extends [[moist-ebm-source]].

**Reference sources — pin at build (the `[[…-source]]` discipline, not carried from
memory).** Phase 1 pinned `[[ebm-radiation-source]]` (`A, B, D, α, T_freeze` — Budyko
1969 / North 1975 / climlab defaults). Phase 2 pins **`[[whittaker-biome-source]]`** +
**`[[precip-parameterization-source]]`**; Phase 3 pins **`[[shallow-water-source]]`**.
The §9.1 knobs pin their own when built: a **`[[stellar-spectrum-ice-albedo-source]]`** for the
spectrum-as-albedo-modifier knob, a **`[[obliquity-insolation-source]]`** for the obliquity knob (the
daily-mean-insolation formula — Hartmann *GPC* §2.7 / Berger 1978 / Rose's climlab notes — and the
mean-annual Legendre context, Nadeau & McGehee 2017 / North 1975). **Rung 4 pins
[[planet-rung4-radiation]]** (the gray radiative-transfer build — Pierrehumbert *PoPC* §4 / Goody & Yung;
Trenberth–Fasullo–Kiehl 2009; Soden & Held 2006; Myhre+ 1998); a line-by-line *spectral*-radiation source
is still future (the named band upgrade).

## 11. Spin-out roadmap — the editable-ocean GPU project (born here, across a contract seam)

**The decision (recorded, not yet acted on).** A *separate* future project — an **editable land/ocean
world with 3-D GPU visualization**, driven by a real ocean engine (**ClimaOcean.jl**) — is **born from
this repo across a documented contract seam**, *not* built inside it and *not* cold-started. Rejected:
*stay-within* (planet-sim stays Python/pedagogical/atmosphere-only — a different language, stack, and
audience) and *cold-start* (throws away the emergent atmosphere this repo already produces). The spin-out
is a **consumer** of planet-sim's output, the same "consume-don't-vendor" relationship planet-sim has with
its shared `engines/` — the seam is the product boundary, the language boundary (Python ↔ Julia), **and**
the physical boundary (atmosphere ↔ ocean), all the same line.

**The collapsed-seam insight (why this is cheap).** ClimaOcean runs the ocean + sea-ice with the
**atmosphere *prescribed* as forcing** (its default is the JRA55 reanalysis). planet-sim **is** an
atmosphere — EBM temperature + the emergent shallow-water jet + rung-2 moisture/`P−E`. So the spin-out's
one-way coupling is a **substitution**: *planet-sim's emergent atmosphere replaces JRA55.* One-way coupling
this way needs **neither ClimaAtmos nor ClimaCoupler** — only ClimaOcean. ClimaCoupler (the full coupled-ESM
orchestrator) is owed **only** for the deferred two-way loop (S5). The GPU-native 3-D path on the Julia side
is **Makie.jl** (GL/WGL) — the synergistic renderer the new repo gets for free in its own language; on the
planet-sim side the same fields render through the §9.5 globe stack already built.

> **Verified against primary sources (2026-06-12).** The four load-bearing architecture claims above were
> checked against the CliMA documentation, not assumed: (1) ClimaOcean is *"a framework for realistic
> ocean-only and coupled ocean + sea-ice simulations **driven by prescribed atmospheres**"* — atmosphere
> prescribed, not prognostic, confirmed; (2) `JRA55PrescribedAtmosphere` is the canonical built-in, and the
> headline example *"initialized from the ECCO state estimate, coupled to a prescribed atmosphere derived from
> the JRA55-do reanalysis"* — note the README frames JRA55 as **the example built-in, not a blessed sole
> default**, which *strengthens* the substitution argument: the prescribed-atmosphere slot is pluggable by
> design, so dropping planet-sim's atmosphere in where JRA55 sat is the interface working as intended, not a
> hack; (3) ClimaOcean is **standalone**, *"built on Oceananigans and ClimaSeaIce"* — no ClimaAtmos / no
> ClimaCoupler dependency, confirming one-way needs ClimaOcean only; (4) ClimaCoupler is the orchestration
> layer *"bringing atmosphere, land, ocean and sea ice together"* with feedback (CMIP/AMIP) — exactly the
> two-way loop, correctly deferred to S5. ECCO confirmed as a **dataset/state estimate** (ECCO2Daily,
> `ECCORestoring`), the init + restoring anchor. **One live wrinkle the living-staircase rule already covers:**
> the CliMA docs note that *generic coupling machinery is migrating to a separate package* (`NumericalEarth.jl`),
> so the exact package/type names S2–S3 bind against may shift before we get there — which is precisely why S2
> is *"design against the **seen** API, not the guessed one"* and sits after S3. Re-confirm at S1.

**The honesty ceiling (carried from the staircase, hardened here).** A custom world's ocean output **cannot
be validated** — there is no observation of a planet that doesn't exist. So Earth is the **only** anchor:
ECCO (below) is the ocean's `[[ebm-radiation-source]]`-class ground truth. The discipline is **anchor on
Earth, then trust the delta** for edited worlds — *honest-by-construction* on Earth, *honest-by-disclosure*
off it (the §9.5 Rung C carve-out, now load-bearing for an entire project). And it is **not real-time** for
a custom world (the expensive engine is the truth; a cheap learned emulator for instant editing is a named,
deferred layer, not in this roadmap).

### 11.1 The two seams — do not fuse them

The single most important correction baked into this roadmap (advisor-caught): "the seam" is **two** seams,
in opposite directions. Conflating them silently drops the atmosphere→ocean forcing direction, which is the
actual engine-coupling value of the whole spin-out.

| | **Viz / output seam** | **Forcing / input seam** |
|---|---|---|
| **Direction** | model → screen | atmosphere → ocean engine |
| **Payload** | grid + `(u,v)` **velocity** + scalar + frames + **coverage-extent** + **provenance/honesty label** | wind **stress** + surface heat-flux components + `P−E` (+ optional SST-restoring) |
| **Purpose** | *display* a field | *drive* an ocean |
| **Consumer** | a renderer (§9.5 globe; Makie; the Rung C showcase) | ClimaOcean's forcing API |
| **Where designed** | **in planet-sim, now** (§9.3/§9.5 — already mostly built) | **at the boundary, later** — once ClimaOcean's input API is *seen*, not guessed |
| **Status** | the near-term planet-sim work (R1) | the spin-out's first physics deliverable (S2/S4) |

They share a backbone — *pin a schema, not a format* (§9.3): one schema, two encodings (**JSON + `.npz`**
browser-friendly, **NetCDF** Julia/ClimaOcean-friendly), with `import(export(s)) == s` round-trip-identity
as the real correctness test — but they carry different physical quantities for different ends, and the
forcing seam is **deliberately not designed until ClimaOcean's API is known** (designing an input contract
against an unseen API is the classic over-fit). Note: **neither ECCO-ingest nor ClimaOcean-vs-ECCO validation
needs the forcing seam at all** — those run on JRA55. The forcing seam only goes live at S4.

### 11.2 planet-sim rungs — before & at the seam (R1–R3, with Rung C in parallel)

These finish *in this repo*, and bank the viz/output seam so the spin-out has something stable to consume.
planet-sim stays **atmosphere-only** throughout — it never ships an ocean visual.

> **Status note (2026-07-06): the *visual* half of that sentence is amended by §9.6.** planet-sim now ships
> a real-data ocean **visual** (the O-rungs: OSCAR/ECCO ingest → viz/output seam → Rung-C globe) — the
> §11.4 named alternative, exercised under the living-staircase rule. The **engine** half stands: no ocean
> *physics* in this repo, ever; ClimaOcean stays the spin-out, and the forcing seam stays S2's.

- **R1 — materialize & serialize the viz contract (this is what the spin-out binds on).** Today the
  circulation field is *computed-then-viewed* (§9.1 `vector_overlay`, the Phase-4 jet) — R1 **serializes**
  it: pin the §9.3 planet-spec schema to carry a **vector-field layer** (grid + `(u,v)` + scalar + frames +
  coverage-extent + provenance), in both encodings, with the **round-trip-identity test** extended to cover
  it. The decisive move: **add a second, *synthetic global-coverage* producer** and read **both** the real
  eddy band **and** the synthetic field back through the **already-built Rung B renderer** (`eddy_globe.py`)
  — proving **producer-agnosticism** (the renderer/serialization does not care *what* generated the field).
  *Producer-agnosticism is the exact property ClimaOcean later relies on* — an ocean engine's `(u,v)` + SST
  flows through the same contract as the eddy band. **Not** Rung C — this rung — is what the spin-out depends
  on. *Retarget-when-done:* the serialized schema's shape is the first input to S1's ECCO ingest; revisit it
  after seeing a real ocean field's dimensions.

  > **BUILT 2026-06-14** (`planet/flow_serialize.py`, `planet/demo_flow_serialize.py`,
  > `planet/tests/test_flow_serialize.py`; the committed globe `docs/figures/planet-flow-serialize.html`).
  > **The seam reuses the existing schema — no new format** (ADR 0004 #3 "one structure serialized, not a
  > second one invented"): a `FlowField` is expressed as a `VECTOR_OVERLAY` `Layer` in a
  > `PlanetView`/`PlanetSpec`, so `planet_spec.save`/`load` round-trips it unchanged. Coverage-extent +
  > provenance + `radius_m` + honesty ride in the JSON-safe `style` dict, **cast native** (the win32
  > landmine: numpy ints don't `json.dumps`, and a surviving np scalar breaks the round-trip `==`). **The
  > load-bearing proof is the round-trip identity on *both* producers** (`load(save(spec)) == spec` — the
  > real eddy band *and* the synthetic global field), pinned in `test_flow_serialize.py` + the crown-jewel
  > `_synthetic_spec`. Two design calls the advisor settled:
  > - **Renderer = `planetmap.render` (the generic cone overlay), not `eddy_globe.py`.** The plan's
  >   "already-built Rung B renderer (`eddy_globe.py`)" was a **stale reference** — `eddy_globe` consumes an
  >   `eddy` *object* (it is not field-generic). The producer-agnostic vehicle is `planetmap.render`, which
  >   dispatches on `LayerKind` and paints any `VECTOR_OVERLAY` as cones with no edit. (`render` gained one
  >   backward-compatible `caption=` param so a longitudinally-structured field can carry an honest caption
  >   instead of the biome family's zonal-mean one — the geometry renderer is otherwise untouched.) A common
  >   **speed = |(u,v)|** scalar layer is the surface for both producers (`render` needs a scalar surface;
  >   speed is native to any vector field → the two views are structurally parallel).
  > - **The band is embedded on a *full-globe* grid, zeros outside its coverage box** — not a patch grid
  >   (whose poleward `_polecapped` padding would smear the band across the sphere) and **not**
  >   `circulation_layer`'s mirror-and-wrap (valid only for the zonal-mean jet). The emergent eddy band is
  >   the project's only longitudinally-structured field, so it is laid **NH-only, true ~55° sector only,
  >   zeros elsewhere** (honest-by-construction; pinned by a test). `is_global=False` on a full-globe grid is
  >   not a contradiction — coverage records *where the data is*. This is also exactly the **ECCO target
  >   shape** (a full-globe `(u,v)`+SST with data everywhere) — *the grid you serialize is the grid you
  >   render*. **`frames`** (a time axis) is a **named, deferred** schema increment (orthogonal to
  >   producer-agnosticism; trivially a stacked `.npz` later) — R1 serializes a single saturated snapshot.
  >   **The deferral is now scheduled: §9.6 O4** (built after the real ocean field's dimensions are seen
  >   at O2 — exactly the "revisit against a real field" this note reserved).
- **R2 — toolkit promotion (§9.4 rule-of-three, the natural co-rung).** With the frame side-channel, the
  flow-globe renderer, and the serialization now serving a **third** consumer (the synthetic producer / the
  spin-out), the rule-of-three is met: promote the viz+serialization machinery to a documented, shared
  contract. *Retarget-when-done:* whatever the promotion reveals as "still planet-specific" is a candidate
  cut before the seam freezes.

  > **BUILT (documentation) 2026-06-14 — the promotion is *documentation, not extraction*.** R1 created the
  > third consumer for the **serialization** machinery (`planet_spec` now serves the biome export + the
  > two-world diff + the vector-field interchange), so that contract's rule-of-three is met — but it is
  > *already* a clean shared module, so "promote to a documented, shared contract" = **document it**
  > (§9.4 post-split status note + this record), not move code. **The globe-geometry viz *helpers* stay
  > two-consumer** (`_sphere_xyz`: planetmap→eddy_globe; `_band_geometry`/`_earth_radius`:
  > eddy_globe→flow_globe — R1 added a *serialization* consumer, not a *geometry* one), so extracting them
  > to a `planet/viz/` now would be the **pre-emptive promotion §9.4 forbids** — named + held, to re-trip on
  > a genuine third geometry consumer. The other half of R2: the **stale cross-repo `viz/` /
  > ARCHITECTURE.md §6 references** (dead since the 2026-06-10 split) were refreshed to the post-split
  > reality (§9.4, §12.3).
- **R-parallel — Rung C (the three.js/WebGL particle showcase), OFF the critical path.** Rung C proves
  **renderer-agnosticism** (a *different* axis from R1's producer-agnosticism) and is the immersive
  *honest-by-disclosure* showcase. It is a planet-sim viz deliverable (§9.5) that **does not gate the
  spin-out** — it can land during or after R1–R3. If it slips, the split still proceeds.
- **R3 — bank planet-sim (the clean hand-off).** The atmosphere-only capstone: the viz/output seam is
  documented + tested, the schema is round-trip-pinned, Rung C ideally landed. This is the point at which the
  new repo is worth starting — *not before R1 is banked.*

### 11.3 The spin-out repo — after the split (S1–S5)

Numbered by *logical grouping*; the **execution order** has S3 teach the API that S2 needs, so it runs
**S1, then S3 → S2 → S4**, then decide S5.

- **S1 — ECCO ingest + viz, pure Python (the new repo's first rung; no Julia engine yet).** Ingest the
  **ECCO** ocean state estimate and render it as a globe (the *Perpetual-Ocean* visual) through the viz
  seam. **Dual role, which is why it earns its own rung:** (i) a striking **real-data** deliverable that
  stands up the data pipeline cheaply, and (ii) the **validation anchor** for S3. Independent of all the
  Julia work. *Retarget-when-done:* the ECCO field's real dimensions/coverage retarget the viz-seam schema
  inherited from R1.

  > **Status note (2026-07-06): S1 narrows — role (i) moves to planet-sim as §9.6 O2.** The real-data
  > ingest + the Perpetual-Ocean visual now land in planet-sim (the §11.4 named alternative, taken for the
  > viz half). What remains of S1 in the new repo: role (ii) — **ECCO as the validation anchor for S3**
  > (standing up the *Julia-side* data path, re-confirming the CliMA package landscape) — inheriting O2's
  > pipeline and schema knowledge instead of rebuilding them. The R1-schema retarget duty moves to O2's
  > own retarget note.
- **S2 — design the forcing/input seam (NetCDF encoding, at the boundary).** The atmosphere→ocean contract:
  wind stress + heat-flux components + `P−E` (+ optional SST-restoring), encoded as NetCDF (the Julia-friendly
  encoding of the §9.3 schema). **Designed against ClimaOcean's *seen* input API** — which is why it follows,
  not precedes, S3. *Retarget-when-done:* the seam shape is provisional until a real ClimaOcean run accepts it.
- **S3 — ocean-1: ClimaOcean + JRA55, validated against ECCO (no planet-sim seam).** Stand up ClimaOcean on
  its default JRA55 forcing and **reproduce ECCO** — the Earth anchor, and the rung where ClimaOcean's input
  API is *learned* (the knowledge S2 consumes). Runs entirely on JRA55; the planet-sim forcing seam is **not
  involved** here. *Retarget-when-done:* the validation gap (where ClimaOcean+JRA55 misses ECCO) sets the
  honesty floor for every custom-world claim downstream.
- **S4 — swap JRA55 → planet-sim's emergent atmosphere (the forcing seam goes live; the payoff).** Replace
  JRA55 with planet-sim's EBM+jet+moisture forcing through the S2 seam; the ocean response flows **back**
  through the viz seam to the renderer. This is the actual engine-coupling win: *an edited world's emergent
  atmosphere drives a GPU ocean, rendered in 3-D.* Honest-by-disclosure off Earth. *Retarget-when-done:* the
  realism of the driven ocean retargets whether S5 is worth its weight.
- **S5 — DEFERRED decision point: ocean-2, two-way loop closure (ClimaCoupler territory).** Feed ocean SST
  **back** into planet-sim — which closes rung-2's *faked* evaporation `E` (the prescribed `P−E` whose honest
  `E` the staircase never had). This is where **ClimaCoupler** finally earns its place. **Reassess here**,
  with real ClimaOcean experience in hand, before committing — the roadmap deliberately *stops* and re-plans
  at this gate.

### 11.4 The settled fork & the living-staircase rule (the user's caveat, made the framing principle)

**One fork, settled: ECCO-ingest lives in the *new* repo as S1** (planet-sim stays atmosphere-only). The
alternative — ECCO as planet-sim's last rung, to de-risk against real data *before* the split — stays named,
because (per the rule below) even settled forks can be retargeted.

> **Retargeted (2026-07-06) — the named alternative taken, proving the rule.** Real-data ocean
> ingest-for-*display* now lives in planet-sim as §9.6 O2 (user-ordered: the renderer and the seam live
> there, and the deliverable is wanted now). Only the **viz half** of the fork moved — the engine/forcing
> half (ClimaOcean, S2–S5) is unmoved, and S1 keeps its validation-anchor role (§11.3 status note).

**The rule (the user's caveat, 2026-06-12): every rung is provisional until the previous one lands.** This
is not a frozen waterfall — it is the **same living-contract discipline the rest of this plan runs on**
(engines are living contracts / ADR 0005; spike-first de-risking; "a trade, not a win"; anchor-then-trust-
delta). Each rung above carries a **`Retarget-when-done`** note precisely because its successor's *target*
is expected to move once the predecessor's real output is in hand: R1's serialized schema is re-judged
against S1's real ECCO dimensions; S2's forcing seam is provisional until an S3 run accepts it; S5 is
re-planned from scratch after S4. **Plan the next rung concretely, hold the rung after it loosely**, and
**revalidate the whole chain at each landing** — the staircase is climbed one validated step at a time, with
the banister redrawn after every step.

## 12. Upgrade backlog — every named growth axis, consolidated

**What this is.** A single index of every upgrade *named but not built* across the project, pulled
together from the per-rung records in §10, the macro-staircase in §5, and the memory files — so "what
could be built next" is one scan, not a hunt. This is an **index, not a re-derivation**: each line is
*what it is · what it unlocks · the wall-or-cost · → the detailed record*. The physics lives in §10 / the
`[[…]]` memory; do not duplicate it here.

**Maintenance rule (or this rots).** When an upgrade ships, **strike its line here and record it in its
§10 rung block** — the same "plan §10 + memory is the record" discipline the rest of the file runs on. No
separate backlog file (a third record to hand-sync) and no MEMORY.md line for this section.

**Status tags** preserve the project's honesty gradient — do not flatten them to "TODO":
**[named upgrade]** buildable on the current engine, no new anchor · **[the wall]** a prescribed
closure / cited constant a future rung would *derive* · **[left to user call]** scoped, explicitly not
foreclosed, awaiting a go · **[scoped]** designed + pressure-tested, not built · **[feasibility sketch]**
explored + costed, not committed · **[decided — separate repo]** belongs to the spin-out · **[deferred —
~free]** cheap docs/notebook work, offered not asked.

### 12.1 The macro-staircase — rungs 5–6 (see §5, not re-tabulated)

- **Rung 5 — idealized GCM** · sphere core + convection/PBL/cloud params + slab ocean + sea ice →
  aquaplanet GCM. The deferred **Held–Suarez 3-D sphere primitive-eq** test lives here (over-reaches the
  rung-3 2-layer β-plane). **Leaves the laptop.** → §5 table; [[planet-rung3-scoped]].
- **Rung 6 — full GCM / ESM** · topography + full moist physics + clouds + ocean GCM + carbon cycle. The
  infeasible end (compute + context-coherence + validation walls). → §5 table.

### 12.2 Within-rung upgrades on the landed rungs (1–4)

Split per the honesty gradient (the §12 intro's "do not flatten to TODO"): **buildable slices** are
pick-up-able now on the current engine — written as `- [ ]` work items carrying their **done-condition
(the anchor)**, the field that makes a slice more than a wish. **Walls / awaits a higher rung** are
prescribed closures or caveats a *future* rung must clear — left as descriptive bullets, deliberately
*not* checkboxed. Slice fields: *deliverable · anchor (the triad test that pins it) · cost · tag · → record.*

**Rung 1 — two-way coupler (complete).**
*Buildable slices:*
- [x] ~~**TVD/WENO flux limiter**~~ **BUILT 2026-06-14** (`engines/fluid/shallowwater.py`,
  `tracer_limiter=` ∈ {minmod, vanleer, mc, superbee}; `test_tracer_limiter.py`, 10 tests; full-repo gate
  **432 fast**). A **TVD** flux limiter (van Leer default) on the passive-tracer face value, written as
  `θ_face = θ_up + ½·ψ(r)·(θ_down − θ_up)` — and **ψ≡1 IS the existing centered average**, so the unlimited
  path is literally the ψ≡1 special case. **Opt-in, default-off → centered scheme byte-for-bit** (the
  `None` branch is the original expression verbatim, not routed through the ψ form — float-reorder would
  break `array_equal`; advisor #1). **Anchor met:** a uniform x-flow advecting a step develops **no new
  extrema** (θ∈[IC] to 1e-9, positive) under all four limiters where centered overshoots — *rigorous 1-D
  TVD* (grid-aligned flow, h frozen by an exact steady state); `∫hθ` stays machine-exact (conservative
  flux form untouched); variance/passivity seals stay green. **Honest edge (advisor #3, sharper than my
  draft):** strict TVD is 1-D — in 2-D the dimension-split limiting (Goodman–LeVeque) gives no maximum
  principle, so the 2-D test asserts only that the limiter *reduces* the overshoot vs centered (+ it is
  dissipative ⇒ variance one-sided-decreases; gentle van Leer clips smooth extrema only mildly). Sweby
  region `0≤ψ≤min(2,2r)` + **raise on unknown limiter name** (advisor #4: ψ≡2 is bounded yet non-monotone).
  No reuse ripple (`LayeredShallowWater`/`baroclinic_qg` have their own `_rhs`; kwarg last). WENO not
  needed — TVD meets the anchor at a fraction of the code. **The advisor done-check then caught a real
  sign error** in the `U<0`/`V<0` upwind smoothness ratio (must be `θ−θ_plus`, not `θ_plus−θ` — pinned by
  the linear-ramp `r=+1` criterion in both flow directions); the `+x`-only anchor missed it. Its empirical
  signature is **direction-asymmetric over-diffusion** (negative-flow peak 0.83 vs 0.96), **not** the
  predicted overshoot — a pure step stays monotone under the buggy limiter too (no interior extremum to
  over-amplify) — so the regression is a **reflection-symmetry** test (peak retention is *bit-identical*
  for ±x/±y on the correct scheme) plus the ±x/±y-parametrized monotone step. **[named upgrade — BUILT]** →
  [[planet-rung1-two-way-coupler]].
- [x] ~~**Wet-get-wetter precip shape/amplitude**~~ **BUILT 2026-06-14** — shipped in **`planet/moist.py`**
  (`wet_get_wetter_precip_field` + `_amplify_contrast` + `WetGetWetter`/`wet_get_wetter`), **NOT**
  `circ_precip.py`: the §12 pointer was **stale** (it predated the rung-2 build that added
  `energy_constrained_factor`); the advisor sited it next to its ingredients, where it is literally the
  **generalization** of the one-line `energy_constrained_precip_field` (that scales mean+anomaly together
  at one rate; this **splits** them — mean at the energy rate, the wet−dry anomaly at C–C). Held & Soden
  2006 "rich-get-richer", on the precip pattern: warming **wettens the ITCZ/storm-track while DRYING the
  deserts** — fixing the rung-0 uniform-`CC` flaw (it wettens deserts too). **Honesty (advisor):** a
  better *prescribed* parameterization, **not** derived — the split *direction* and the two cited rates
  are calibrated; what is structurally exact is the **mean-zero anomaly split** (so ⟨P⟩ scales at the
  energy rate — plumbing, not a finding) + the **reduction to BOTH** existing fields when the rates
  coincide. Opt-in/default-off (rung-0 `precip.py` unchanged); **deliberately not fused** with the
  storm-track position seam (the non-composition rule). Edges: global-`T̄` (local-`q_sat(T(φ))` = richer
  named upgrade), `P≥0` floor (deep warming → total aridity), thermodynamic-only (dynamic shift = rung 3+).
  Demo `planet/demo_wet_get_wetter.py` → `docs/figures/planet-wet-get-wetter.png`; tests in `test_moist.py`
  (7 fast + 1 slow guard). **[named upgrade — BUILT]** → [[planet-rung1-two-way-coupler]].

*Walls / awaits a higher rung:*
- **Realistic-magnitude eddy κ** · the emergent barotropic `κ~10³` is ~1000× below rung-0's `2.2e6` and
  config/window-tuned (config can't be separated from physics). The *dimensionless* win moved to rung-3
  QG; a realistic *dimensional* magnitude on this single-layer engine stays unbankable. **[the wall]** →
  [[planet-rung1-two-way-coupler]].

**Rung 2 — moist dynamics (complete, incl. 2.5 / 2.x / Hadley fix).**
*Buildable slices:*
- [x] ~~**Emergent ITCZ rain**~~ **BUILT 2026-06-14** (`planet/sphere_moist.py`) — the full-sphere moisture
  budget (eddy + a two-cell Hadley cell anchored on the EFE), conservative `∫(P−E)=0` machine-exact, rain
  *rained* not *painted*. **The advisor's predicted "q-contrast" headline was OVERTURNED** by a decomposition
  spike: the displaced-ITCZ intensity is **geometric, not emergent `q`** (clean negative result), and the
  wet/dry dipole is **displacement-driven** (direction by-construction). The two real nuggets: co-location of
  the **net** `P−E` on the EFE is a **falsifiable check** (the prescribed cell beats the eddy export ~2.6×),
  and the geometric-not-`q` negative. Genuinely-new emergent content is **modest** (2.x × Hadley recombined
  into a conserving budget); the architectural win (budget-not-band) is the "what". **[named upgrade —
  BUILT]** → [[planet-rung2x-itcz]].
- [ ] **Tighten the ITCZ-migration sensitivity** · *deliverable:* re-derive `D` in `sphere_ebm.py` ·
  *anchor:* the closed-form `δ/AHT` lands within ~factor-1 of observed (currently ~2× high, rides the
  calibrated `D`) · *cost:* one re-derive. **[named refinement]** → [[planet-rung2x-itcz]].
- [x] ~~**Recalibrate `D_s≈0.28` (rung 2.5)**~~ **ALREADY BUILT** — this backlog line was **stale**: the
  recalibration shipped *inside* the parent rung-2.5 build (`recalibrate_sensible_D` in `moist_ebm.py`,
  commit `e65f33c`) and was upgraded to the dt-free `moist_steady_direct` path in `be61dad`. Live values
  (RH 0.8): `D_s = 0.2784` (≈ 0.28, < the dry 0.555), present-contrast matched, endpoint PA 2.05;
  **target-invariance is 0.81%** (contrast vs P₂ amplitude — well under the <5% anchor, confirming the
  docstring's "<2%"). All three anchors pinned by passing tests
  (`test_recalibration_matches_the_present_contrast`, `test_recalibrated_Ds_is_below_the_dry_default`,
  `test_pa_is_invariant_to_the_recalibration_target`; 14/14 green). **No work was outstanding** — only this
  checkbox. The genuine remaining rung-2.5 residual is the *structured `D(x)` shape* wall below.
  **[named upgrade — was already BUILT]** → [[planet-rung25-mse-diffusion]].

*Walls / awaits a higher rung:*
- **Derive `R_ATM_SLOPE` (=2 W/m²/K)** · the prescribed atmospheric radiative response (cited-closure;
  explicitly **not** `B_OLR`, so rung 4 did *not* retire it). **[the wall]** → [[planet-rung2-scoped]].
- **Derive `HADLEY_STRENGTH`** · the prescribed overturning amplitude that flips the ITCZ sign
  (convergence-by-construction is plumbing). Needs **gross moist stability (GMS) = rung 3+**. **[the
  wall]** → [[planet-rung2-hadley-fix]].
- **Relocate the subtropical desert** · both eddy and Hadley paths leave the dry belt ~12° equatorward of
  the canonical 25–35° (hyper-peaked C–C `q`); a real fix needs GMS-resolved moisture transport.
  **[the wall / rung 3+]** → [[planet-rung2-hadley-fix]].
- **Structured `D(x)` shape residual (rung 2.5)** · the `D_s` recalibration slice fixes the *level*; a
  single scalar still leaves a higher-moment **shape residual** only a structured `D(x)` closes.
  **[the wall]** → [[planet-rung25-mse-diffusion]].

**Rung 3 — baroclinic instability (complete: Phase A linear + Phase B QG flux won).**
*Buildable slices:*
- [ ] **N-layer QG** · *deliverable:* extend `planet/baroclinic_qg.py` from 2 layers to N (multi-level
  vertical structure) · *anchor:* the single-layer/low-N reduction stays exact + the multi-level
  dispersion matches the analytic Phillips matrix (the rung-3 tight leg) · *cost:* more layers, no cleaner
  anchor above it. **[named upgrade]** → [[planet-rung3-qg-built]].

*Walls / awaits a higher rung:*
- **Realistic dimensional κ** · the QG win is **dimensionless + qualitative**; the dimensional κ landing
  in Earth's band is coincidental + box/drag-dependent. A non-coincidental magnitude needs realistic
  forcing/geometry = rung 5. **[caveat → rung 5]** → [[planet-rung3-qg-built]].

**Rung 4 — gray radiative transfer (complete).**
*Buildable slices:*
- [x] ~~**Per-latitude EBM wire**~~ **BUILT 2026-06-14** (`planet/radiative_ebm.py`) — gray `OLR(Ts,τ)`
  drives each band via a coupled Newton solve. §12's guess of "radiative polar amplification" was
  **OVERTURNED → TROPICAL** amplification (the mirror of rung-2.5): `B_loc` smallest at the warm equator
  (WV feedback beats Planck). `D` was **not** recalibrated (present contrast already ≈ rung-0's at the
  climlab-matched loading); the global mean is **not** pinned (Jensen). → §10 rung-4 completion;
  [[planet-rung4-radiation]].
- [x] ~~**Spectral-band log law**~~ **BUILT 2026-06-14** (`planet/radiation.py`,
  `SpectralCO2Band`) — the CO₂ 15-µm band resolved into spectral bins whose absorption falls off
  **exponentially in the wings**, each bin a gray sub-problem solved with the *same* two-stream
  emission kernel. The §12 anchor **held** (not overturned, the rung-4-wire/lapse-rate pattern this
  time confirms): per-doubling `ΔF` becomes **constant** (~4.5 W/m²/doubling, the Myhre band) where
  gray's *decreases* 48→…→20 — the exponential wing's `τ=1` level spreads a constant spectral width
  per doubling. **Triad:** the independent anchor is **reduction-to-gray** (the band kernel
  `_transmission_emission`, written independently, reproduces the gray `_olr_from` to machine
  precision; a *uniform*-`k` band saturates like gray — the wing is the whole ingredient); the unlock
  is loose (the **magnitude is the wall** — rides the wing scale `l`, band-centre τ, half-width;
  calibrated to *order*, "wings ≈ exponential" is itself empirical — so the *functional form* is the
  win, not the ~3.7 coefficient); the derivation `dF/dlnC=2l·π[B(Ts)−B_strat]` matches the τ=1 sharp
  limit ~1% (the column smear realizes ~20–30% more). **Range-limited (named edges):** linear/√ below
  band-centre saturation, saturates again above where the wings exhaust — `0.5×–8×` sits in the flat
  middle. Demo `planet/demo_spectral_band.py` → `docs/figures/planet-spectral-band.png`.
  **[named upgrade]** → [[planet-rung4-radiation]].
- [x] ~~**Moist-adiabatic lapse-rate feedback**~~ **BUILT 2026-06-14** (`planet/radiation.py`,
  `moist_adiabat=True` + `feedback_kernel`) — variable-Γ moist adiabat makes the lapse-rate feedback
  *emergent*. **The §12 anchor was OVERTURNED:** it does **not** supply `λ_LR≈0.84` to land at climlab's 2
  — it supplies `≈+1.5` and *overshoots* (with-WV `B≈3.1`, *above* 2). Banked tight: sign + kind
  (upper-trop amplification), kernel closure (Planck+LR+WV=B, ~1e-3), resolution convergence,
  bit-for-bit reduction at `moist_adiabat=False`. Loose magnitude (two reasons): single *global* column =
  the *tropical* branch only (misses the extratropical branch that pulls the global mean to 0.84) + rides
  the τ shape & WV loading (the wall). → §10 rung-4 lapse-rate block; [[planet-rung4-radiation]].
- [ ] **Clouds** · *deliverable:* a cloud layer in the column (clear-sky today) · *anchor:* clear-sky
  reduction stays exact + cloud feedback order-validated · *cost:* **large** — cloud radiative properties
  + fraction are their own modelling problem, not a one-sitting slice. **[named upgrade — but big]** →
  [[planet-rung4-radiation]].

*Walls / awaits a higher rung:*
- **Komabayashi–Ingersoll runaway** · the hot analogue of the snowball — a steep WV loading where the
  linearization breaks (`B` not linearized across it). An exploration, not a build. **[named edge]** →
  [[planet-rung4-radiation]].

### 12.3 Visualization & interactivity

- [x] ~~**Rung-C GPU advection**~~ **BUILT 2026-06-13** (`planet/flow_globe.py`) — the named seam shipped
  **ahead of** its <~30 fps trigger (user-requested): advection now runs entirely on the GPU by default
  (RGBA32F state texture, off-screen `UPDATE_FS` ping-pong, `DRAW_VS` reads positions from the state tex),
  the CPU `BufferGeometry` loop demoted to a runtime fallback (WebGL can't run in CI → feature-detect +
  GLSL compile-validate + console diagnostics, so a GPU failure degrades to a working globe). Renderer-only:
  `FlowField`/`_build_data`/disclaimer/carve-out all unchanged. Full build record in §9.5; 7th structural
  test pins both pipelines. **[shipped — seam closed]** → [[planet-viz-animation-rungs]].
- [x] ~~**§9.4 toolkit promotion (rule-of-three)**~~ **LANDED (as documentation) 2026-06-14 — R1+R2.** R1
  (`flow_serialize`) gave the **serialization** machinery its third consumer (export + diff + vector-field
  interchange), meeting rule-of-three — but it is already a clean shared module, so the promotion is
  **documentation** (§9.4 post-split note + §11.2 R2 record), not extraction. The globe-geometry *helpers*
  stay **two-consumer** (deliberately not promoted — pre-emption is what §9.4 forbids). Also refreshed the
  stale cross-repo `viz/` / ARCHITECTURE.md §6 references (dead since the 2026-06-10 monorepo split).
  **[shipped — documented, not extracted]** → §9.4, [[planet-spinout-roadmap]].
- **Live notebook widgets** · §2 snowball (live hysteresis / two stable states), §4 winds, §5 jet.
  Advisor flag: time the continuation/coupler sims first → `continuous_update=False` or precompute before
  promising "live". **[deferred]** → [[interactive-what-if]].
- [x] ~~**Browser what-if 3rd axis (obliquity)**~~ **BUILT 2026-06-14** (`planet/interactive.py`) — the
  precomputed grid gained a third axis: a tilt slider over 0…45° (9 values, capped at
  `OBLIQUITY_FAITHFUL_MAX`, the exact `OBLIQUITY_EARTH` float included so the detent stays bit-for-bit),
  wired through `obliquity_params` → `EBMParams.s2` and narrated by the *existing* `explain.py` `obliquity_deg`
  rules (no prose work — it already handled the knob). As predicted it touched only the data axis: S0/CO2
  untouched (S0 keeps the snowball cliff + detent), `_LAT_STRIDE` 3→6 (free, 30 lats indistinguishable on a
  300px canvas) keeps the page ~4 MB (≈ the eddy-globe precedent). `cells[(i·nCo2 + j)·nObl + k]`; cells
  ≈3.7 k, the slow byte-golden stays CI-skipped. **[shipped — ~free, as scoped]** → [[interactive-what-if]].
- **Ocean-currents showcase, rungs O1–O5** · real global surface currents through the R1 seam onto the
  Rung-C globe (§9.6, scoped 2026-07-06 — the §11.4 fork retargeted, viz half only): **O1** per-cell
  validity mask (`FlowField.mask` — the first contract growth past R1) · **O2** OSCAR/ECCO producer +
  provenance disclaimer clause (the deliverable; narrows spin-out S1) · **O3** renderer beauty pass (land
  base layer, accumulate-and-fade trails, speed-weighted seeding; unlocks the §9.5 control-surface seam) ·
  **O4** frames time axis (the R1 deferral; OSCAR monthly climatology → the Somali-Current reversal) ·
  **O5** QG producer (independent; re-trips §9.4 for the geometry helpers → **re-affirmed HOLD**). Wall-or-cost:
  O2's data acquisition/auth (Earthdata login) was the one external unknown — spiked + settled.
  **[BUILT + BANKED 2026-07-06 — all five rungs shipped]** → §9.6.

### 12.4 Pedagogy (the notebook)

- **Bucket A — predict-then-check** · one hypothesis prompt before each section's slider (concrete prompts
  already drafted: pole-vs-equator warming, the hysteresis gotcha, desert migration, jet direction).
  **[deferred — ~free, markdown-only]** → [[pedagogy-novice-intermediate]].
- **Bucket B — mission cards** · goal-directed challenges read straight off the figure (freeze the planet,
  make a desert world, find the habitable-zone edges, build an M-dwarf world); missions drafted. **[deferred
  — cheap]** → [[pedagogy-novice-intermediate]].
- **Two-world diff** · load two saved specs → Earth-vs-your-world side by side (the layer registry already
  supports multiple views). A follow-on to the built bucket C. **[deferred — stretch]** →
  [[pedagogy-novice-intermediate]].

### 12.5 Editable geography & seasonality (§5 / §9.3)

- **Cheap tier (rides rungs 0–1)** · elevation → a lapse-rate map diagnostic; land/ocean → an albedo
  difference; fraction-per-band → continentality-lite. **No engine change — buildable now.** → §9.3.
- **Seasonal cycle / ocean heat capacity** · annual-mean v1 drops `C` at equilibrium, so thermal lag +
  continentality need a seasonal cycle first. **[named, not built — the §3 scope edge]** → §9.3.
- **True 2-D longitudinal geography** · regional climate, orographic precip, rain shadows (the north
  star) — new transport that **leaves the 1-D engine** = the rung-5 exit from the zonal-mean planet.
  Literal rung 5 (a full idealized GCM) is infeasible-tier ("leaves the laptop", §5), so it is climbed
  as a **spike-first sub-ladder of reduced laptop models**, each banking an analytic anchor — the same
  decomposition the gas-giant sketch uses beneath its infeasible deep-shell ceiling (§5).
  **[→ rung 5]** → §5.
  - **Rung 5A — linear orographic precipitation · BUILT 2026-07-09** (`planet/orographic.py`,
    `test_orographic.py`). Smith & Barstad (2004) linear theory: a **diagnostic** on a *prescribed*
    uniform wind over a 2-D terrain → a wavenumber-space transfer function (one FFT) → windward rain +
    a lee **rain shadow**. **Wakes the dormant elevation seam** (§9.3 — carried inert since v1). Tight
    anchor = convergence to the **closed-form triangle-ridge solution** — but honestly scoped: in that
    limit `(1 − i m H_w) → 1`, so it pins only the *reduced* transfer function (`C_w`, upslope `iσ`,
    fallout `τ_f`); the `sgn(σ)` vertical-wavenumber **branch** is guarded *solely* by the rain-shadow
    direction test (a branch flip reddens exactly that one). Plus the upslope limit `C_w·max(0,U·∇h)`,
    a wind-reversal *reflection*-symmetry check (not a branch test), flat-ground null. Constants cited
    ([[smith-barstad-orographic-source]]). **Honest scope: a *trade*, not the engine leaving the zonal
    mean** — the *precipitation* becomes 2-D, the *temperature* climate stays zonal-mean (the Phase-2
    diagnostic-precip precedent, one rung out). **[→ 5A.2]** = sphere placement of the patch, where the
    cross-mountain wind comes from on a zonal-jet globe (prescribed, not emergent — the named caveat),
    the cm/yr↔mm/hr + lat×lon integration, serialization, and the demo/figure. See [[planet-rung5a-orographic]].

### 12.6 Spin-outs — separate repos, not upgrades of this one

Two destinations now decided to live *outside* this repo. They differ in **maturity and kind** — keep them
from reading as equally committed:

- **Editable-ocean GPU project** · a Julia / ClimaOcean + Makie repo born here across a contract seam — a
  *downstream consumer* of this repo's serialized output (R1's producer-agnostic schema is what it binds
  on). **Roadmapped:** spin-out steps S1–S5, ECCO ingest at S1; provisional per the living-staircase rule
  (§11.4). **[decided — separate repo · roadmapped S1–S5]** → §11; [[planet-spinout-roadmap]].
- **Gas-giant atmosphere** · the *same* shallow-water / two-layer-QG engines pointed at a new planet —
  destination now decided (a sibling repo, **not** a within-this-repo rung), but still a feasibility
  sketch, *not* scoped to the ocean's S1–S5 depth. Three tiers, costed in
  [`docs/explorations/gas-giant-atmosphere.md`](../explorations/gas-giant-atmosphere.md):
  - **Tier 1 — β-plane banded jets** · ≈ one rung on `baroclinic_qg.py` (the Williams-1978 Jovian-jet
    model; pyqg is its published twin). Not "turn up β" — needs scale separation (a Rhines window),
    anisotropic diagnostics (zonal-mean `ū(y)` bands; today's azimuthal spectrum is jet-blind), and a
    gas-giant forcing (small-scale stochastic + large-scale drag).
  - **Tier 2 — sphere-correct globe** · global jet count / polar polygons / equatorial superrotation — a
    new spherical-geometry engine (Dedalus / EPIC); both current engines are doubly-periodic Cartesian
    β-planes.
  - **Tier 3 — deep convective interior** · Busse-annulus QG + rotating Rayleigh–Bénard are *reduced*
    laptop-scale entries (a steeper reach, **not** out of scope); only the realistic anelastic-deep-shell +
    MHD dynamo regime is the wall.
  - **[decided — separate repo · still feasibility sketch]** → [[gas-giant-feasibility]];
    [[planet-spinout-roadmap]].
