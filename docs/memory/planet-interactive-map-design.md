---
name: planet-interactive-map-design
description: "Planet's deep-end interactive map — the converged design, ADR 0004, and the viz/ thin-extraction finding (design done 2026-06-09; planetmap.py v1 + planet_spec.py BUILT 2026-06-09)"
metadata: 
  node_type: memory
  type: project
  originSessionId: 5c6cf01d-ca3b-414c-9680-b83ef02898dd
---

The "visual engine" session (2026-06-09): user asked to work on the visual
engine/representation. Two tracks emerged.

**Track A — shared `viz/` toolkit (the literal "visual engine", ADR 0002 §3): NOT
built, deliberately.** Read all three `plots.py` bodies (steel 589 / chip 384 /
planet 109): rule-of-three is satisfied by *count* but **not by substance** — they
share **conventions + styling, not copy-pasted code** (near-identical ADR-doctrine
docstrings, color-constant blocks, the `axvline/axhline + marker + annotate-arrow`
idiom, grid/legend defaults). Heatmap = steel-only (one use); stacked-bars =
steel-only; **time-animation exists nowhere**. So a `viz/` extraction today = thin
styling shim = the premature abstraction ADR 0002 §3 itself warns against. User
**leans "thin extraction + doc changes only," deferred**. The `planetmap.py` layer
registry (Track B) is the eventual **third consumer** that will earn the
field/animation promotion to `viz/` by rule-of-three — don't pre-promote. **Finding
recorded in docs 2026-06-09** (ADR 0002 status-note + ARCHITECTURE §12 #3 status
clause) — the doc-pass is DONE, don't redo it. The *thin extraction itself* (the
styling/doctrine helpers) is leaned-to but **still unbuilt** — an available future task.

**Track B — the deep-end interactive planet map (the priority).** Design **converged
with the user** and written up; **build pending** (after Phase-2 biomes). New ADR
**`docs/decisions/0004-interactive-maps-and-state-interchange.md`** + planet plan §9
rewritten (9.1–9.4) + §5 two staircase consequences + §4 module map + §10. Locked
decisions:

- **D4 layer registry** — renderer is a generic painter over a stack of
  `(name, kind, array, style, z-order)` layers (`kind ∈ {scalar field, vector/line
  overlay, annotation}`); phases **register** layers, never edit the renderer. = the
  user's "show more features as phases progress" made structural.
- **D5 tech = Plotly + ipywidgets** behind a **`[webviz]` extra** (the
  `[viz]`/`[climate]` pattern); matplotlib too weak, full web app = the editable future.
- **Slider lifecycle (ADR 0004 #2)** — renderer is array-consumer ⇒ **invariant up
  the GCM staircase**; only the *trigger* changes: live-slider→instant-remap (rung 0)
  → set-params→launch-run→view (heavy rungs). User raised this; named so nobody keeps
  the live loop alive past where compute allows.
- **Editable planet = preplanned, NOT built (D1/D2).** Two distinct seams: (1)
  *renderer-input* built 2-D-ready now (v1 broadcasts `T(φ)` across longitude → globe
  paints honest **latitude bands**, since v1 is zonal-mean); (2) *geography-physics* =
  a documented **"geography spec"** (elevation+bathymetry+land/ocean mask as plain
  INPUT arrays) — contract written, **no machinery**. **Honesty flag the user
  explicitly accepts:** imported/edited geography is **INERT** at v1 (carried,
  displayed, round-tripped — does NOT change climate). North star = regional climate /
  orographic precip / rain shadows = the rung-5 2-D exit from zonal-mean. Cheap 1-D
  tier (elevation→lapse-rate diagnostic, land/ocean→albedo diff) rides rungs 0–1;
  ocean-heat-capacity→seasonality needs the seasonal cycle annual-mean v1 lacks.
- **State interchange (ADR 0004 #3–4): pin a SCHEMA not a format.** A versioned
  **planet-spec** (grid geometry + explicit units + the layer list [the registry IS
  the manifest] + knob values + `schema_version`). Encoding per consumer: v1 lean =
  JSON manifest + `.npz`; editable-geography heightmaps = **16-bit grayscale PNG**
  (browser/paint-tool native, 8-bit too coarse); NetCDF = future climate-interop
  encoding, **deliberately not v1** (browser-hostile, and the future *web editor* is
  the consumer). **Round-trip identity `import(export(s))==s` is a REAL test** (unlike
  the map's smoke-tests).
- **Exoplanet knobs (§9.1):** *amount* of radiation = `S₀` (already the Snowball
  lever); *stellar spectrum/type* → an **albedo modifier** (redder star → near-IR ice
  albedo ↓ → feedback weakens → harder to snowball; citeable, modest); *planet size* →
  **leaves 0-D global-mean T untouched** (that's S₀,α,A,B only), enters **only
  transport** (bigger → steeper equator-pole gradient); richer size effects route
  through **rotation** = Phase-3 fluid engine. All deepen up the staircase, not faked.

**Build order unchanged:** next = **Phase-2 biomes** (small, spine-reuse, no new
engine) so the map's FIRST version *is* the biome map (the §9 centerpiece), **then**
`planetmap.py` v1 + `planet_spec.py`. See [[planet-plan]], [[bigsim-program]],
[[end-of-batch-ritual]].

**BUILT 2026-06-09** — `planetmap.py` v1 + `planet_spec.py` (`projects/planet/`, 31-test
pair, planet gate 65→96, full suite **470**; opt-in behind a new `[webviz]` extra =
plotly+ipywidgets). First version *is* the biome map: a Plotly 3-D `go.Surface` globe
painted from the **layer registry** — `LayerKind{scalar_field,vector_overlay,annotation}`
+ `Layer(name,kind,data,units,style,z_order,inert)` + `Grid` + `PlanetView`; builders
`build_view(BiomeResult)`/`climate_view(S0,A,D)` = the rung-0 **live recompute = pure
consumer of `demo_biomes.compute`** (added an `n_tau` kwarg so the live loop keeps the
fine present-day step — coarsening shifts the bands via the O(Δt) splitting bias). v1
layers = temperature/precipitation/biome (scalar) + ice-line (annotation, **±φ_ice circles
w/ a NaN-gap separator**) + an **inert elevation** geography seam. **Three-layer thin-skin
discipline** (Steel app.py): registry/builders NumPy-only (always-green) → `render()`
plotly-lazy (smoke-test `importorskip`, **NOT slow** — no kernel) → `interactive_map()`
ipywidgets-lazy (the untested `main()` analogue); guard test = subprocess import-isolation
(in-process `sys.modules` check is session-fragile once render tests import plotly).
`planet_spec` = **pin-the-schema** (grid+units+layer-list[=the registry IS the manifest]+
knobs[=`asdict(EBMParams)`]+`schema_version`), v1 lean encoding = **JSON manifest +
`.npz`**; **round-trip-identity `load(save(s))==s` = the deep end's ONE real (non-smoke)
test** (`__eq__` array-aware, `np.array_equal(equal_nan=True)` for floats only — int biome
codes compare exact; negative-control test guards against a trivially-true `__eq__`).
**TWO advisor "do-less" calls (do not re-litigate):** (1) **obliquity DEFERRED** = a
named/disabled slider — wiring it needs `s₂(obliquity)` pinned to a source (the
`[[…-source]]` discipline), so ship only the validated **S₀/CO₂→A/D** knobs; (2)
**`vector_overlay` declared-but-UNPAINTED** — renderer raises `NotImplementedError` naming
Phase 4 (build the seam not the machinery; a synthetic vector overlay = "machinery without
a consumer"). The extensibility proof = the **inert elevation scalar layer** §9.3 requires
(existing scalar renderer paints it free, planet_spec round-trips it) — better than a
synthetic vector. **NO gate.py edit** (no new engine; `GATES["planet"]` stays
`{engines/diffusion}`). Banked artifact = `docs/figures/planet-map.html` (write_html,
plotly-only, no kaleido). **`planet.ipynb` BUILT 2026-06-09** (chip-style thin skin, planet gate
96→97): the matplotlib teaching notebook — arc = EBM `T(φ)` (live, overlaying the **same climate
with the ice feedback ON vs OFF** so the gap IS the feedback and → 0 when ice-free) → Snowball
hysteresis (static banked, a 120-equilibria sweep not a live point) → biome map (live S₀/CO₂
sliders) → §4 markdown pointer to *this* `[webviz]` globe. **DURABLE GOTCHA (advisor-caught by
eyeballing the banked figure, commit 2b57cee):** the first cut overlaid North's **constant-α
two-mode** as the "no-ice anchor" and claimed gap=ice-feedback — WRONG: the constant-0.30 analytic
and the real poleward-brightening ice-free albedo (a₀+a₂·P₂ = 0.26 at equator) differ *even with no
ice* → curves split at the equator, never merge. Right anchor = the **ice-feedback-OFF climate**
(`present_day_climate(replace(params, T_freeze=-1000.0))` → freeze isotherm below any T → α never
ices, SAME solver) → equator-coincident / poleward-split, merges exactly when ice-free; present-day
max gap 4.92°C (not the misleading 7.73 vs two-mode). North's two-mode stays credited as the tight
validation leg in the "where the numbers come from" table, just not plotted. Thin-skin/flicker discipline
reused (every validated call in a *direct* cell before the `interact` sugar; figure-before-print);
`slow` subprocess smoke-test mirrors `test_chip_notebook.py`. **NEXT = Phase 3** (`engines/fluid` →
the map registers a circulation `vector_overlay`, no renderer edit).

**§9.1 exoplanet knobs BUILT** (5b8d598 + c797ce0 advisor fix, merged 2026-06-10):
`exoplanet.py` (`stellar_ice_albedo` two-band ratio-to-solar + `transport_for_size`
D∝1/size² derivation + `exoplanet_params` composing both onto a base `EBMParams` —
parameter derivations only, no engine/EBM change; Sun + Earth-size defaults recover v1
**bit-for-bit**), `demo_exoplanet.py` + banked `docs/figures/planet-exoplanet.png`
(M-dwarf ~83 % narrower Snowball loop; ice line 90°→71°→42° over 0.5→1→2 R⊕), both knobs
wired into `climate_view` + `interactive_map` sliders. Source = [[stellar-spectrum-ice-albedo-source]].
**Advisor catch (c797ce0): two-level mean honesty** — the *analytic* 0-D mean (≈14.5 °C)
is the size-invariant quantity; the *relaxed* mean is NOT (drops sharply at 2 R⊕ via the
enlarged ice cap's albedo feedback); `print_summary` shows BOTH columns and tests assert
the two-level structure, not a false "mean ~stable" band. Planet fast gate 121→136.

**§9.1 obliquity knob BUILT** (d760731 + advisor docstring fix, merged 2026-06-10): the
**last deferred knob, now wired** — `obliquity.py` derives the insolation `s₂(ε)` from first
principles (integrate the daily-mean-insolation formula over a circular-orbit year, project
onto P₂), applied as the ratio to Earth so 23.44° recovers climlab `s₂=−0.48` bit-for-bit;
pure parameter-feed into `EBMParams.s2` (no engine/EBM/gate change). `OBLIQUITY_FAITHFUL_MAX
= 45°` = a named UI/scope cap (the faithful pre-reversal regime), NOT a pinned physics number.
Source = [[obliquity-insolation-source]]. **All three §9.1 exoplanet knobs (spectrum, size,
obliquity) are now wired into `climate_view` + the interactive map's sliders.** Banked
`docs/figures/planet-obliquity.png`. Planet fast gate 136→146 (+10), +2 slow demo tests.

**`planet_spec.build_spec` + `planetmap.climate_params` BUILT 2026-06-10** (for the notebook §7
design-a-world sandbox, [[pedagogy-novice-intermediate]]): extracted the knobs→`EBMParams` composition
(obliquity `s2` + exoplanet `ai`/`D`) out of `climate_view` into `climate_params` — the **single source
of truth** `climate_view` / `build_spec` / the notebook bench all ride on (behavior-preserving refactor;
the present-day round-trip canary stayed green untouched). `build_spec(**knobs) -> PlanetSpec` is the
"knobs → saveable spec" entry point and **fixes a real trap**: `climate_view` builds the composed params
*and discards them*, so `from_view(climate_view(...), EBMParams(S0,A,D))` would store the **un-perturbed**
knobs for any non-default star/size/tilt → the exported world re-runs to the WRONG climate; `build_spec`
serializes the *same* composed params it solved, so `to_params()` reconstructs the exact world. Stays
headless (the subprocess import guard holds). +4 tests (3 in `test_planet_spec`: defaults-recover-
`climate_params` + round-trip, the composed-knobs trap guard, equals-the-explicit-`from_view` path; 1 in
`test_planetmap`: `climate_params` composes onto exactly `ai`/`D`/`s2`). Notebook §7 = matplotlib design
bench (save fires only on the button → Run-All writes nothing; worlds in gitignored `outputs/worlds/`).

**Two-world diff/compare BUILT 2026-06-10** (the "follow" on the §7 bench; user split it into two
committed batches: 18411a8 = diff+matplotlib cell, 10559a9 = webviz globes). **Batch 1 — `planet_spec.diff(a,b)
-> SpecDiff`** (the *compare* sibling of `PlanetSpec.__eq__`): changed `knobs` (`name→(a,b)`, with units —
an exoplanet reports the *composed* `ai`/`D`/`s2`, never the raw star/size/tilt), a per-cell `FieldDelta`
per scalar field (signed `b−a` for numeric; an `a≠b` **changed-mask** for the categorical biome — code
subtraction is meaningless; `mean_delta` is honestly area-weighted on the uniform-in-sinφ grid), layer-set
edges (`only_in_a`/`only_in_b`/`other_changed`), and a `grids_compatible` guard (skips per-cell deltas
across different `n_cells`, keeps the knob diff). **Data-scoped** (knob values + layer *data*, not layer
metadata) ⇒ the invariant `bool(diff(a,b)) == (a!=b)` holds for built specs (asserted in tests); reuses the
schema's `_arrays_equal`. Notebook §7 `compare_worlds` (knob dicts OR saved specs) = matplotlib: knob table
+ per-field Δ summary + ΔT/Δprecip band profiles + side-by-side biome bands; surfaces `grids_compatible`.
**Batch 2 — the diff on the globe (webviz triptych A·B·Δ).** `planet_spec.delta_view(a,b,active)` = the
**headless** Δ-globe data (the active field's `FieldDelta` as a 1-layer `PlanetView` — diverging `b−a` +
`cmid=0` for numeric, the 2-tone changed-mask for biome; raises on grid-mismatch / non-scalar `active`;
JSON-safe style keeps `planet_spec` Plotly-free — **delta_view lives in `planet_spec` not `planetmap`
because `planetmap→planet_spec` is a real import cycle**). `planetmap.render_comparison(view_a,view_b,
view_delta,active)` = a 1×3 **scene**-subplot triptych that **reuses `render()` unchanged** — refactored out
`_scalar_surface` (now honoring a `cmid` style hint) + `_overlay_traces`, render() externally identical;
takes *views* not specs (no cycle); A shares B's colorbar (A hidden), Δ its own, positioned so the 3 bars
don't stack (verified via distinct scene x-domains [0,.32]/[.34,.66]/[.68,1], the thing a build-only smoke
test can't catch). `save_comparison_html` = the `save_html` analogue (NOT wired into `main()`). Notebook §7
`compare_globes` is **defined-not-run** on Run-All (needs `[webviz]`, writes to gitignored `outputs/` behind
an explicit call — mirrors §6 "describe, don't embed"). +10 tests (delta_view always-green + render_comparison
/save_comparison_html `importorskip` smoke). Full suite green.
