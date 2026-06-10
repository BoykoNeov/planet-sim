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

### 9.4 Toolkit promotion

Plot primitives start planet-local (`plots.py` static floor; the `planetmap.py` layer
registry for the deep end). The **2-D field / heatmap** and **time-animation** primitives
are the ADR 0002 §3 candidates whose third reuse (after steel/chip) would **promote to
the shared `viz/`** by rule-of-three (ARCHITECTURE.md §6); the layer registry is that
third consumer-in-waiting. Promotion is **not** done pre-emptively (the existing three
`plots.py` share conventions, not copy-pasted code — the thin-extraction finding,
2026-06-09).

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
diagnostic). **THE anchor — reduction-to-EBM:** the closure `⟨v'θ'⟩=−D_eff·∂θ̄/∂y` has the *same form*
as the EBM transport term, so the two-way model with a constant flow-diagnosed `D_eff` *is* a rung-0
diffusive EBM with that `D`; tested by (a) re-equilibrating at an *independently*-chosen κ recovering a
rung-0 EBM at the bridge-implied `D`, and (b) **rung-0 being a fixed point** of the map (the EBM's own
diffusive flux → `D_eff=D`, climate unchanged). Plus the bridge round-trip/magnitude and the
**right-signed response** (stronger flux ⇒ flatter contrast). **De-risked the A/B split** (advisor):
Phase A drives the machinery with a **synthetic** down-gradient flux (the Phase-4 synthetic-gradient
playbook), so it lands *independent* of the (tuned) eddy sim; the `flux_fn` argument is the seam Phase B
plugs into. Planet gate **151 passed, 1 skip** (no shared-engine edit → planet gate, not full-repo).
**Next (step 2, Phase B):** the *emergent* eddy flux — advect θ (relaxed to the EBM target) on the
barotropically-unstable jet, diagnose `⟨v'θ'⟩` **post-saturation** via a life-cycle integral (release
mode), magnitude named tuned; the loop-closes claim scoped to **one feedback pass** (a converged
fixed-point is a `slow` demo if it converges cleanly). Then **step 3** circulation-informed precip.

**Reference sources — pin at build (the `[[…-source]]` discipline, not carried from
memory).** Phase 1 pinned `[[ebm-radiation-source]]` (`A, B, D, α, T_freeze` — Budyko
1969 / North 1975 / climlab defaults). Phase 2 pins **`[[whittaker-biome-source]]`** +
**`[[precip-parameterization-source]]`**; Phase 3 pins **`[[shallow-water-source]]`**.
The §9.1 knobs pin their own when built: a **`[[stellar-spectrum-ice-albedo-source]]`** for the
spectrum-as-albedo-modifier knob, a **`[[obliquity-insolation-source]]`** for the obliquity knob (the
daily-mean-insolation formula — Hartmann *GPC* §2.7 / Berger 1978 / Rose's climlab notes — and the
mean-annual Legendre context, Nadeau & McGehee 2017 / North 1975), and (rung 4) a radiative-transfer
source if spectral radiation is ever computed rather than parameterized.
