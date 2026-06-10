# `planet` — the Earth-system / Planet simulator (the capstone)

*Planetary knobs in, climate & habitability out.* Project #3 of the program and its **capstone** —
the first project to reuse the diffusion/heat spine (`engines/diffusion`) a **third** time
(as a sphere's latitudinal heat transport) and, in later phases, to build the program's one
remaining shared engine (`engines/fluid`, the shallow-water solver). Full plan:
[`docs/plans/planet-earth-system.md`](../../docs/plans/planet-earth-system.md).

> **Units — SI / climlab-conventional** (the deliberate contrast with Chip's per-module native units):
> **W m⁻²** (S₀, insolation, the OLR offset A), **W m⁻² K⁻¹** (the OLR slope B, the transport D),
> **°C** (T, the freeze isotherm Tf — the climlab convention `A+B·T` and `Tf` assume °C), and the
> dimensionless **area coordinate `x = sin φ`** on `[0, 1]` (equator → pole; equal Δx = equal area on
> the sphere, so the global mean is a plain `∫₀¹ T dx`). Latitudes are reported in **degrees**. The
> engine is fed the latitudinal transport in these units directly.

## Load pointer (per-session working set, ARCHITECTURE.md §11)

- **To work on the EBM machinery (Phase 1):** `ebm.py` + `tests/test_ebm.py`. It loads
  `engines/diffusion/CONTRACT.md` (**heat mode**: array diffusivity `D_eng(x) = (D/C)(1−x²)`,
  insulated Neumann(0) at both ends) and **Strang-splits the radiation around it** — the
  **Jominy-2a idiom reused** (`projects/steel/jominy.py`): the linear `−B·T` relaxation is an exact
  exponential half-step, the albedo threshold makes the local step nonlinear (what creates the
  bistability). Public API: `EnergyBalanceModel` (the transport + split-radiation solver),
  `equilibrium_temperature_0d` / `two_mode_solution` (the analytic anchors), `insolation`,
  `legendre_P2`, `ice_line_latitude`, `ClimateState`. The module docstring is its contract.
- **Three interchangeable steady-state modes (two orthogonal knobs)** — the design call from the
  Phase-1 review, useful for accuracy/speed *and* as a mutual cross-check web:
  - `face=` on the model — `"harmonic"` (plain `(1−x²)`; the engine's harmonic-mean faces, a named
    ~0.1 °C polar bias) vs `"exact"` (cell values **pre-distorted** so the harmonic mean reproduces
    the true face coefficient, no bias). *The engine is never modified either way.*
  - `method=` on `EnergyBalanceModel.equilibrium` — `"relax"` (the Strang-split relaxation, the
    general / **only** nonlinear-capable path, used by the Snowball sweep) vs `"direct"` (a dt-free
    linear solve, the splitting-error-free **reference** for the constant-albedo North check; it
    *raises* on the ice feedback). The direct path's operator is **pinned to the engine** by a test.

    The default `face="harmonic"`, `method="relax"` is the simple, general, snowball-capable combo.
- **To work on the ice-albedo feedback & the Snowball hysteresis (Phase 1):** `albedo.py` +
  `tests/test_albedo.py`. The step-function albedo (`planetary_albedo`, `absorbed_shortwave`), the
  parameter bundle (`EBMParams`), the present-day finite-cap branch (`present_day_climate`), and the
  **continuation-sweep hysteresis** (`snowball_hysteresis` → `HysteresisLoop`). The module docstring
  is its contract.
- **To work on the banked artifact (Phase 1):** `demo_snowball.py` + `tests/test_demo_snowball.py`
  (the end-to-end integration test, `slow`-marked) and `plots.py` (the figure — `[viz]` extra). The
  demo wires `present_day_climate` + `snowball_hysteresis` → `plots.snowball_figure` and saves
  `docs/figures/planet-snowball.png` (the hysteresis loop + ice-line loop + present-day `T(φ)` profile).
- **To work on the biome map (Phase 2 — the payoff, banked early):** `biomes.py` + `precip.py` +
  `tests/test_biomes.py` + `tests/test_precip.py`. `precip.py` is the **diagnostic precipitation**
  parameterization — a circulation-set Gaussian latitude pattern (ITCZ-wet / subtropics-dry /
  midlat-wet / poles-dry, in **cm/yr**) times a **Clausius–Clapeyron** global-moisture amplitude
  (`precipitation`, `precip_field`). `biomes.py` is the **Whittaker `(T,P)→biome` classifier** — an
  *original*, total, sloped-boundary partition of the (T,P) plane (the **Irvin precedent**: the
  copyrighted diagram is reproduced by an independent computation calibrated to it, not digitized),
  with cold biomes temperature-limited and warm biomes moisture-limited on diagonal thresholds
  (`Biome`, `classify`, `classify_field`, `biome_area_fractions`). Both reuse only the EBM — no new
  engine. The module docstrings are their contracts.
- **To work on the Phase-2 banked artifact:** `demo_biomes.py` + `tests/test_demo_biomes.py`
  (`slow`) and `plots.biomes_figure` (`[viz]`). The demo composes `present_day_climate` → `precip` →
  `biomes` and saves `docs/figures/planet-biomes.png` — the **biome-band map** + the **Whittaker
  (T,P) plane** shaded by biome with the planet's climate trajectory drawn through it (the mechanism:
  *why* deserts sit at ~30°) + the `T(φ)`/`P(φ)` profiles. It also warms the planet (a CO₂ proxy) to
  show the bands migrate poleward.
- **To work on the benchmark (Phase 1):** `climate_reference.py` + `tests/test_climate_reference.py`.
  A **frozen reference table** of the climlab/North benchmark facts (present ice line ~70°, the
  Snowball threshold, the hysteresis) keeps the triad green without the `[climate]` extra; the live
  climlab cross-check (`climlab_present_day`) is a `slow` / `importorskip` test (the pycalphad pattern).
- **To work on the deep-end interactive map (§9 / ADR 0004):** `planetmap.py` + `tests/test_planetmap.py`.
  The **layer registry** — `LayerKind` (`scalar_field` / `vector_overlay` / `annotation`) + `Layer` +
  `Grid` + `PlanetView` — and the builders that turn a climate result into the v1 **biome-map** layer
  stack (`build_view`, `climate_view` = the rung-0 live recompute over S₀ / CO₂→A / D + the two §9.1
  exoplanet knobs star `T_star` / planet `size`). `render` paints
  it as a Plotly globe (3-D `go.Surface` + the ice-line annotation; lazy import, `[webviz]` extra);
  `save_html` banks the standalone globe (`docs/figures/planet-map.html`); `interactive_map` is the
  in-notebook live-slider loop (the `main()` analogue — ipywidgets, not unit-tested). The renderer is
  **generic over `LayerKind`**: later phases *register* layers (circulation = `vector_overlay`, Phase 4;
  it raises `NotImplementedError` until then), never edit the renderer. Built matplotlib-free; the
  module docstring is its contract.
- **To work on the state-interchange schema (§9.3 / ADR 0004 #3–4):** `planet_spec.py` +
  `tests/test_planet_spec.py`. The **planet-spec** — grid + explicit units + the layer list (the
  registry *is* the manifest) + the knob values + a `schema_version` — with `from_view` / `save` /
  `load` (the v1 lean encoding = a JSON manifest + a `.npz`). **Round-trip identity
  `load(save(spec)) == spec` is the one *real* correctness property of the deep end** (`__eq__` is
  array-aware, `equal_nan` for the annotation NaN gaps), so it is an always-green test, not a smoke
  test. Carries an **inert elevation** layer (the geography seam — round-tripped, not yet consumed by
  the climate; §9.3). NumPy/JSON only — no render deps.
- **To work on the exoplanet knobs (§9.1 — stellar spectrum & planet size):** `exoplanet.py` +
  `tests/test_exoplanet.py`. Two **parameter-deriving** knobs (no engine, no new EBM physics) the
  interactive map wires: **stellar spectrum → ice albedo** (`stellar_ice_albedo` — a two-band
  blackbody-weighted snow/ice albedo, applied as the *ratio* to the solar value so the Sun recovers
  the climlab `ai = 0.62` exactly; a redder star lowers it → harder to snowball; pinned empirical bands
  [[stellar-spectrum-ice-albedo-source]]) and **planet size → transport** (`transport_for_size` —
  `D ∝ 1/size²`, **derived** from the spherical Laplacian in `x = sin φ`; a bigger planet sharpens the
  equator-pole gradient, the 0-D mean size-invariant). `exoplanet_params(T_star, size, base)` composes
  both onto an `EBMParams`. The module docstring is its contract.
- **To work on the §9.1 banked artifact:** `demo_exoplanet.py` + `tests/test_demo_exoplanet.py`
  (`slow`) and `plots.exoplanet_figure` (`[viz]`). The demo traces a Sun-vs-M-dwarf Snowball loop pair
  + size-scaled `T(φ)` profiles → `docs/figures/planet-exoplanet.png` (a redder star → a narrower loop;
  a bigger planet → a steeper gradient).
- **To work on the obliquity knob (§9.1 — axial tilt → the insolation gradient):** `obliquity.py` +
  `tests/test_obliquity.py`. A third **parameter-deriving** knob (no engine, no new EBM physics): the
  axial tilt `ε` sets the insolation P₂ coefficient `s₂` the EBM already accepts. `s₂(ε)` is **computed**
  by integrating the pinned daily-mean-insolation formula over a circular-orbit year and projecting onto
  P₂ (`insolation_p2_coefficient` — *not* a memorized coefficient), then applied as the *ratio* to the
  Earth value (`insolation_s2`) so Earth's tilt recovers the climlab `s₂ = −0.48` exactly. Validated by
  the **exact `s₂(0) = −5/8`** analytic limit + the independent `≈−0.48`-at-23.44° climlab cross-check
  ([[obliquity-insolation-source]]); scope edge = the single-P₂-mode truncation degrades at high tilt
  (the ≈55° sign reversal is real but surfaced as a loose bracket). `obliquity_params(ε, base)` composes
  it onto an `EBMParams`. Banked: `demo_obliquity.py` + `tests/test_demo_obliquity.py` (`slow`) +
  `plots.obliquity_figure` (`[viz]`) → `docs/figures/planet-obliquity.png` (the `s₂(ε)` curve + the
  relaxed `T(φ)` flattening with tilt). The module docstring is its contract.
- **To work on the circulation / shallow-water engine (Phase 3):** `circulation.py` +
  `tests/test_circulation.py`. It loads `engines/fluid/CONTRACT.md` (the program's
  **second shared engine** — a rotating shallow-water solver, hyperbolic/explicit, sharing no
  machinery with the parabolic-implicit diffusion spine) and pins the **planetary** numbers the
  engine leaves to its consumer: Earth's `f₀ = 2Ω sinφ`, `β = 2Ω cosφ/a`, and an equivalent depth
  giving the cited extratropical deformation radius `L_R ≈ 1000 km` ([[shallow-water-source]]).
  Public API: `midlatitude_beta_plane` (the configured β-plane), `geostrophic_adjustment` /
  `rossby_wave` (the two banked demos returning plain-array records), `coriolis_f0` / `coriolis_beta`.
  Phase 3 is **dry one-layer dynamics in isolation** — no EBM coupling yet (that is Phase 4's
  `coupler.py`, where a forced jet emerges and the map registers its `vector_overlay` layer).
- **To work on the Phase-3 banked artifact:** `demo_shallowwater.py` + `tests/test_demo_shallowwater.py`
  (`slow`) and `plots.shallowwater_figure` (`[viz]`). The demo composes `geostrophic_adjustment` +
  `rossby_wave` → `docs/figures/planet-shallowwater.png`: the adjustment (bump → balanced vortex over
  `L_R`, with the conservation diagnostics holding flat) beside a **westward** Rossby wave.
- **To work on the one-way coupler (Phase 4 — the capstone payoff):** `coupler.py` +
  `tests/test_coupler.py`. It **couples the two engines**: the EBM's meridional temperature
  gradient (a `ClimateState`) is mapped to a doubly-periodic, zero-mean **height target**
  (`height_target` — windowed for the wall-less engine), and the flow is forced toward it by
  **operator splitting around the bare engine** (`couple_jet`: exact-exponential thermal relaxation +
  weak Rayleigh drag — the third reuse of the EBM/Jominy idiom) until a **geostrophically-balanced
  westerly jet emerges** (`geostrophic_balance` is the anchor). The jet sits at the EBM gradient
  maximum (emergent, *not* the channel centre — `gradient_peak_latitude`); the periodic channel exacts
  a **flanking easterly return** (named). Conservation is **reframed** (forced–dissipative): mass
  forced-exact + a **release test** (forcing off → the bare engine conserves + the jet persists).
  One-way, dry single layer (two-way = rung 1, the `tracer` seam). The module docstring is its contract.
- **To work on the Phase-4 banked artifact:** `demo_coupler.py` + `tests/test_demo_coupler.py` (`slow`)
  and `plots.coupler_figure` (`[viz]`). The demo forces the present-day climate → emergent jet →
  `docs/figures/planet-coupler.png` (the jet on the geostrophic estimate + the forcing chain + the 2-D
  jet + the forced/release conservation diagnostics); with `[webviz]` it also banks the **interactive
  globe** `docs/figures/planet-coupler-map.html` — the jet drawn as a `circulation` `vector_overlay`
  (Plotly cones) over the temperature field (`planetmap.circulation_layer` / `_vector_overlay_trace`,
  `build_view(jet=…)`; computed-then-viewed, not in the live-slider loop).
- **To use the shallow-water engine directly:** load `engines/fluid/CONTRACT.md` only — the
  one-page API contract (`ShallowWater`, `SWState`, `uniform_grid`; mass + tracer-mass machine-exact,
  energy/PV/variance bounded; the `tracer` slot is now advected as a **passive** scalar — rung 1
  built, ADR 0005). Never the engine internals.
- **To use the diffusion/heat spine:** load `engines/diffusion/CONTRACT.md` only — never Steel's or
  Chip's internals. Planet instantiates the same contract in **heat mode**, with the radiation
  composed around it by operator splitting (the Jominy precedent).

## Status

- **Phase 1 — the latitudinal EBM & the Snowball bifurcation: BUILT** (2026-06-09). `ebm.py` (the
  engine transport + Strang-split radiation, the three interchangeable A/B/C modes) +
  `albedo.py` (ice-albedo feedback + the continuation-sweep hysteresis) + the banked Snowball demo
  (`docs/figures/planet-snowball.png`). **Validation web:** the **North two-mode** analytic profile
  reproduced to ~1e-4 °C by the exact-face direct solve (the harmonic floor named at ~0.1 °C), the
  Strang relaxation shown to **converge** to it as dt→0, the **0-D mean** matching the discrete
  energy balance (net-TOA machine-exact) and the continuous `T̄` to the grid's O(1/n²) limit, and the
  direct operator **pinned to the engine**. **Banked numbers:** present-day global mean
  ≈ 14.7 °C / ice line ≈ 73° (the finite-cap branch — Earth's, in the bistable zone); the planet
  **freezes over** (Snowball) when the sun dims ~8 % and **re-melts only ~580 W/m² brighter** (a wide
  hysteresis loop). 26-test triad green (+1 skipped live-climlab cross-check).
- **Phase 2 — biomes & habitability (the consequence payoff, banked early): BUILT** (2026-06-09).
  `precip.py` (the diagnostic precipitation parameterization — circulation pattern × C–C amplitude,
  cm/yr) + `biomes.py` (the original total Whittaker `(T,P)→biome` partition) + the banked biome-map
  demo (`docs/figures/planet-biomes.png`). **Validation:** the classifier is an **exact, total
  partition** (every (T,P) → one biome, area fractions sum to 1), nine canonical **probe points** land
  in their textbook biome, the present-day **band ordering** (equator→pole: rain forest → savanna →
  desert/grassland → temperate forest → boreal → tundra) is reproduced, and the C–C global-water
  budget is **monotone** in T̄ (a *consistency* check, not a conservation law — named). **No new
  engine** (project-local reuse of the EBM only); 20-test triad added (3 `slow`). Scope edges named:
  Whittaker not Köppen (annual `T,P`, no seasonal precip), prescribed precip not a water cycle, fixed
  band centres (migration is the rung-1/2 circulation enhancement), the C–C 7 %/K is moisture-capacity
  not the energy-constrained ~2–3 %/K global precip rate.
- **The deep-end interactive map — `planetmap.py` v1 + `planet_spec.py`: BUILT** (2026-06-09, plan §9 /
  ADR 0004). The interactive map's **first version is the biome map**: a Plotly 3-D globe painted from a
  **layer registry** (temperature / precipitation / biome scalar fields + the ice-line annotation +
  an inert elevation seam), with **S₀ / CO₂→A / D** knob-sliders driving an instant recompute-and-remap
  (the rung-0 live loop). Pure consumer of `demo_biomes.compute` — no new physics. *(The §9.1 exoplanet
  knobs — star `T_star` / planet `size` — and the **obliquity** knob were since wired in; see the §9.1
  status entries below.)* The **planet-spec** schema
  exports/imports the registry (JSON + `.npz`) with a **real round-trip-identity test** (the deep end's
  one genuine correctness property, vs the render's smoke-tests). Banked globe:
  `docs/figures/planet-map.html`. 31-test pair added (the Plotly render smoke-tests `importorskip` on
  the `[webviz]` extra — fast, not `slow`). **No new engine / no gate-manifest change.**
- **Phase 3 — the shallow-water engine (`engines/fluid`): BUILT** (2026-06-09). The program's
  **second shared engine** — a rotating shallow-water solver on a doubly-periodic β-plane (Arakawa
  C-grid, vector-invariant, explicit SSP-RK3), built standalone and sealed behind its suite
  (`engines/fluid/CONTRACT.md`) before any coupling. **Validation triad:** gravity-wave `√(gH)` &
  **Poincaré dispersion `ω²=f₀²+gHk²`** to ~1e-3 (the rotation check), **Rossby waves** westward &
  dispersive (loose, converging to analytic with resolution), **geostrophic balance steady** +
  **geostrophic adjustment** to the analytic Helmholtz state over `L_R` (~1%), **mass conserved to
  machine precision**, **energy bounded & dt³-convergent**, and — the discriminating Coriolis leg —
  **potential vorticity / enstrophy bounded at FINITE amplitude** (a vortex at Rossby ~0.5). The
  symmetric scheme conserves energy semi-discretely (not enstrophy); claims are stated as measured
  (honest, not aspirational). `circulation.py` pins the planetary numbers (`L_R ≈ 960 km` at 45°) and
  banks the artifact (`docs/figures/planet-shallowwater.png`). **First multi-engine gate row**
  (`planet` uses `{engines/diffusion, engines/fluid}`); the import-drift guard is now live.
- **Phase 4 — the one-way EBM→circulation coupler: BUILT** (2026-06-09 — the capstone complete). The
  EBM's meridional temperature gradient forces the dry shallow-water flow (thermal relaxation + weak
  drag split around the bare engine) and a **geostrophically-balanced westerly jet emerges** (~16.5 m/s
  @ ~42°, core residual ~0.6%) at the **EBM gradient maximum** — emergent, not channel-placed — flanked
  by the easterly return the doubly-periodic channel requires. **Conservation reframed**
  (forced–dissipative): mass forced-exact + a **release test** re-confirming the engine's invariants
  while the jet persists. The interactive map now **paints** the `vector_overlay` (the deferred
  Phase-4 seam): a `circulation` layer over the temperature field. One-way / dry single layer; two-way
  is rung 1 (the `tracer` seam). No engine modified; `planet` uses `{engines/diffusion, engines/fluid}`
  unchanged.
- **§9.1 — the exoplanet knobs (stellar spectrum & planet size): BUILT** (2026-06-09, a growth-axis
  batch on the complete capstone). `exoplanet.py` adds two **parameter-deriving** knobs (no engine, no
  new EBM physics — both compute an `EBMParams` the existing machinery consumes): **(1) stellar
  spectrum → ice albedo** — a two-band blackbody-weighted snow/ice albedo (bright visible / dark
  near-IR, [[stellar-spectrum-ice-albedo-source]]), applied as the *ratio* to the solar value so the
  Sun recovers `ai = 0.62` exactly; a redder host star lowers the ice albedo → the ice-albedo feedback
  weakens → **harder to snowball** (banked: an M-dwarf's Snowball loop ~83% narrower, freeze threshold
  lower; matches Joshi & Haberle 2012). **(2) planet size → transport** — `D ∝ 1/size²`, **derived**
  from the spherical Laplacian in `x = sin φ` (the `D` value [[ebm-radiation-source]]); a bigger planet
  transports heat less per unit area → a **sharper equator-pole gradient** (ice line 90°→71°→42° over
  0.5→1→2 R⊕), the 0-D mean size-invariant (the relaxed mean drifts only via the ice feedback — named).
  Both wired into the interactive map (`climate_view` / `interactive_map` sliders) — defaults (Sun,
  Earth-size) recover the present-day map **bit-for-bit**. Banked `docs/figures/planet-exoplanet.png`;
  **no new engine, no gate-manifest change** (numpy-only, planet-local). Scope edge: only the bright-ice
  albedo responds to the star (ocean/land left unchanged); size is transport-only (rotation effects route
  through the fluid engine — a different rung).
- **§9.1 — the obliquity knob (axial tilt → the insolation gradient): BUILT** (2026-06-10, the
  growth-axis batch's final knob — *obliquity was the lone deferred slider, now wired*). `obliquity.py`
  adds a third **parameter-deriving** knob (no engine, no new EBM physics): the axial tilt `ε` sets the
  insolation P₂ coefficient `s₂`. `s₂(ε)` is **computed** by integrating the pinned daily-mean-insolation
  formula over a circular-orbit year and projecting onto P₂ ([[obliquity-insolation-source]] — *not* a
  memorized coefficient), then applied as the *ratio* to the Earth value so Earth's tilt recovers the
  climlab `s₂ = −0.48` **bit-for-bit**. Validated by the **exact `s₂(0) = −5/8`** analytic limit + the
  independent `≈−0.48`-at-23.44° climlab cross-check (the numerical projection also reproduces the known
  closed form `−(5/8)(1−1.5·sin²ε)` across the range). More tilt spreads sunlight poleward → a flatter
  planet, the ice cap retreats (ice line 60°→71°→ice-free over 0°→23.44°→40°); past ≈55° the gradient
  reverses (poles warmer — surfaced as a loose bracket). Wired into the interactive map (the obliquity
  slider, formerly disabled, is now live). Banked `docs/figures/planet-obliquity.png`; **no new engine,
  no gate-manifest change**. Scope edge: the single-P₂-mode insolation truncates at high tilt (the real
  annual-mean grows `s₄`), and the model is annual-mean (no seasonal extremes); eccentricity/precession
  are a separate deferred Milankovitch axis.
- **The teaching notebook now covers the full arc (`planet.ipynb`, §4–6 added): BUILT** (2026-06-10).
  The chip-style thin skin (built 2026-06-09 for Phases 1–2 + the deep-end globe) gained the two
  dynamics sections it was missing — **§4 the winds** (the second shared engine: geostrophic adjustment
  + a westward Rossby wave, `demo_shallowwater`) and **§5 the emergent jet** (the one-way coupler's
  capstone payoff, `demo_coupler`) — both static banked renders at `nx=48` (the demos still bank the
  full-resolution figures; the notebook is *reach*, so the Helmholtz match prints the honest ~2 % and
  names it as the converges-with-resolution beat). The deep-end globe section (now §6) was de-staled:
  obliquity + the `vector_overlay` are now-built, not deferred, and "what's next" points up the GCM
  staircase (rung-1 two-way coupler / editable elevation). The notebook stays matplotlib-only; the
  globe's obliquity/exoplanet sliders run in `planetmap.interactive_map` behind `[webviz]`, described
  not embedded. Executes clean top-to-bottom (the `slow` smoke test green, 29 s); **no module code, no
  new engine, no gate-manifest change.** **Pedagogy tiered (2026-06-10):** each of §1–§5 gained an
  expert **`<details>` "Going deeper"** collapsible — the convention chip & steel already use, the one
  notebook that lacked it — so the three reading depths are explicit (narrative = novice, sliders =
  intermediate, the collapsible = expert: derivations, the Strang-split / C-grid / energy-not-enstrophy
  machinery, and the named scope edges). Markdown-only, grounded in the cited module docstrings; banked
  code-cell outputs byte-identical, smoke test still green.
- **Engines un-frozen — ADR 0005: DONE** (2026-06-10). The freeze-before-reuse ceremony is dropped;
  `engines/*` are now *living, versioned* contracts (extend directly + test + Changelog; the suite +
  the full-repo gate on an engine edit are the guardrails). Both `CONTRACT.md` status headers flipped.
- **Rung 1 (the two-way coupler) — STARTED, step 1 BUILT** (2026-06-10). `engines/fluid` now advects
  its long-declared `tracer` slot: a scalar `θ` carried in **flux form** through the same SSP-RK3,
  **strictly passive** (no back-reaction on `h,u,v` → dry dynamics bit-for-bit unchanged). New
  `tracer_mass` (`∫hθ`, machine-exact — the anchor) + `tracer_variance` (`∫½hθ²`, bounded) diagnostics
  and `engines/fluid/tests/test_tracer.py` (translation analytic limit, conservation, consistency, the
  not-monotone scope edge). **Step-0 finding:** the Phase-4 jet is *barotropically unstable* (so a
  passive tracer gets an emergent `⟨v'θ'⟩` flux — no imposed wave needed). **Anchor (for step 2):**
  reduction-to-EBM (`D_eff`), not a tuned PW number. **Next:** step 2 the two-way coupler, step 3
  circulation-informed precip.

## Test runner (tiered gate, ADR 0003 — the per-project successor)

```powershell
python -m tools.gate planet -m "not slow"   # routine commit gate: planet's tests + the frozen engine's
python -m tools.gate planet                  # full gate for planet (incl. the slow demo sweep) — EXCEPTIONAL
./run_tests.ps1 -m "not slow"               # whole-repo fast lane (release / CI / shared-engine edit)
```

`pyproject.toml`'s `testpaths` already carries `projects` and `engines`, so both
`planet/tests/` and `engines/fluid/tests/` are collected with no config change;
`pythonpath = ["."]` lets planet import the frozen engines as `engines.diffusion…` / `engines.fluid…`.
The full-resolution Snowball sweep (`test_demo_snowball`), the shallow-water demos
(`test_demo_shallowwater`, `test_circulation` integration), the figure renders, and the live climlab
cross-check are `slow`-marked / extra-gated, so the fast lane deselects them. Planet's gate `uses` is
now `{engines/diffusion, engines/fluid}` — the manifest's **first genuinely multi-engine row** (a
`planet` commit runs both engines' seals), kept honest by the import-drift guard in
`tools/tests/test_gate.py`.
