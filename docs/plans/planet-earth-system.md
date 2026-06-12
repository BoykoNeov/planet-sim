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
    plots.py                  # planet-local static figures (→ promote to viz/ by rule-of-three)
    planetmap.py              # the deep-end INTERACTIVE map: a LAYER REGISTRY painted by Plotly+ipywidgets (ADR 0004 #1, §9.1)
    planet_spec.py            # the planet-spec interchange schema: export/import the layer stack; round-trip-identity tested (ADR 0004 #3-4, §9.3)
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
  phases progress" requirement made structural. The registry stays planet-local; it is
  the **third consumer** that will eventually promote the 2-D-field / animation
  primitives to shared `viz/` by rule-of-three (it does not pre-empt that — ADR 0004).
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

### 9.5 Animated flow — the visualization rungs (decided 2026-06-11; **rungs A+B BUILT** (A 2026-06-11, B 2026-06-12), C pending)

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
| **C** | three.js / WebGL perspective globe | self-contained inline HTML (three.js vendored, no CDN) | the **showcase** — a general flow-on-a-globe **particle-streaming** renderer; *honest-by-disclosure* (decided 2026-06-12, amended same day to a true 3-D sphere, below) |

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
perspective; *honest-by-disclosure*).** The build approach is locked. Rung C is what §9.5 already calls
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
  stays responsive, self-contained (inline three.js, no CDN), and golden-friendly.

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
  `pole/equator δT` ≈ **1.5** (Earth, RH 0.8); **direction banked** (PA>1 robustly), **magnitude loose**
  (1.43→1.50 across RH 0.6→0.8 — RH-dependent; the observed ~2–3× also needs the ice-albedo + lapse-rate
  feedbacks held out of scope).
- **THE ATTRIBUTION NULL (advisor — the backbone).** The moist model differs from dry in *two* ways (a
  recalibrated `D`-shape AND the T-dependent `D_eff`); which causes PA? **Freeze `D_eff` at its present
  profile and warm → PA = 1.0 *exactly*** (uniform `δT=ΔA/B` solves the perturbation for *any* frozen
  `D(x)`), proving the PA is **100 % the `dD_eff/dT` feedback, 0 % the `D`-shape** (spike + test:
  PA=1.0000, spread ~3e-10, via the genuine array-`D` `EnergyBalanceModel`).
- **The recalibration = the named wall (the double-count).** Rung-0's `D=0.555` is an *effective*
  diffusivity already absorbing latent transport; explicit MSE diffusion would double-count it, so
  `recalibrate_sensible_D` re-derives a **smaller sensible `D_s≈0.30`** matching the dry present-day
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
  itself (~1.5, direction banked / magnitude loose). *Plumbing* — RH=0 ∧ `D_s=0.555` ⟹ the genuine
  `EnergyBalanceModel` rung-0 solve **bit-for-bit**. *Named choices* — recalibration to present contrast
  (`D_s<0.555`) + its target-invariance. Tests: `planet/tests/test_moist_ebm.py` (12, all **fast**); full
  planet gate **261 passed, 1 skip**. No engine edit; `uses` unchanged.

**Visualization rungs A/B/C — DECIDED to build all three (animated eddy flow; 2026-06-11;
rungs A+B BUILT — A 2026-06-11, B 2026-06-12; C pending — build detail in §9.5).** A forward decision (user): animate the emergent eddy life cycle across three
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

**Rung C — build approach decided, then amended to a true 3-D sphere; not yet built (2026-06-12).** A
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
**§9.5** + **§9.3** + ADR 0002 status note; build deliverables listed in §9.5. **Plan only this session —
not executed.**

**Reference sources — pin at build (the `[[…-source]]` discipline, not carried from
memory).** Phase 1 pinned `[[ebm-radiation-source]]` (`A, B, D, α, T_freeze` — Budyko
1969 / North 1975 / climlab defaults). Phase 2 pins **`[[whittaker-biome-source]]`** +
**`[[precip-parameterization-source]]`**; Phase 3 pins **`[[shallow-water-source]]`**.
The §9.1 knobs pin their own when built: a **`[[stellar-spectrum-ice-albedo-source]]`** for the
spectrum-as-albedo-modifier knob, a **`[[obliquity-insolation-source]]`** for the obliquity knob (the
daily-mean-insolation formula — Hartmann *GPC* §2.7 / Berger 1978 / Rose's climlab notes — and the
mean-annual Legendre context, Nadeau & McGehee 2017 / North 1975), and (rung 4) a radiative-transfer
source if spectral radiation is ever computed rather than parameterized.

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
- **R2 — toolkit promotion (§9.4 rule-of-three, the natural co-rung).** With the frame side-channel, the
  flow-globe renderer, and the serialization now serving a **third** consumer (the synthetic producer / the
  spin-out), the rule-of-three is met: promote the viz+serialization machinery to a documented, shared
  contract. *Retarget-when-done:* whatever the promotion reveals as "still planet-specific" is a candidate
  cut before the seam freezes.
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

**The rule (the user's caveat, 2026-06-12): every rung is provisional until the previous one lands.** This
is not a frozen waterfall — it is the **same living-contract discipline the rest of this plan runs on**
(engines are living contracts / ADR 0005; spike-first de-risking; "a trade, not a win"; anchor-then-trust-
delta). Each rung above carries a **`Retarget-when-done`** note precisely because its successor's *target*
is expected to move once the predecessor's real output is in hand: R1's serialized schema is re-judged
against S1's real ECCO dimensions; S2's forcing seam is provisional until an S3 run accepts it; S5 is
re-planned from scratch after S4. **Plan the next rung concretely, hold the rung after it loosely**, and
**revalidate the whole chain at each landing** — the staircase is climbed one validated step at a time, with
the banister redrawn after every step.
