# `projects/planet` — the Earth-system / Planet simulator (the capstone)

*Planetary knobs in, climate & habitability out.* Project #3 of the program and its **capstone** —
the first project to reuse the frozen diffusion/heat spine (`engines/diffusion`) a **third** time
(as a sphere's latitudinal heat transport) and, in later phases, to build the program's one
remaining shared engine (`engines/fluid`, the shallow-water solver). Full plan:
[`docs/plans/planet-earth-system.md`](../../docs/plans/planet-earth-system.md).

> **Units — SI / climlab-conventional** (the deliberate contrast with Chip's per-module native units):
> **W m⁻²** (S₀, insolation, the OLR offset A), **W m⁻² K⁻¹** (the OLR slope B, the transport D),
> **°C** (T, the freeze isotherm Tf — the climlab convention `A+B·T` and `Tf` assume °C), and the
> dimensionless **area coordinate `x = sin φ`** on `[0, 1]` (equator → pole; equal Δx = equal area on
> the sphere, so the global mean is a plain `∫₀¹ T dx`). Latitudes are reported in **degrees**. The
> frozen engine is fed the latitudinal transport in these units directly.

## Load pointer (per-session working set, ARCHITECTURE.md §11)

- **To work on the EBM machinery (Phase 1):** `ebm.py` + `tests/test_ebm.py`. It loads the frozen
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
    the true face coefficient, no bias). *The frozen engine is never modified either way.*
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
  stack (`build_view`, `climate_view` = the rung-0 live recompute over S₀ / CO₂→A / D). `render` paints
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
- **To use the diffusion/heat spine:** load `engines/diffusion/CONTRACT.md` only — never Steel's or
  Chip's internals. Planet instantiates the same contract in **heat mode**, with the radiation
  composed around it by operator splitting (the Jominy precedent).

## Status

- **Phase 1 — the latitudinal EBM & the Snowball bifurcation: BUILT** (2026-06-09). `ebm.py` (the
  frozen-engine transport + Strang-split radiation, the three interchangeable A/B/C modes) +
  `albedo.py` (ice-albedo feedback + the continuation-sweep hysteresis) + the banked Snowball demo
  (`docs/figures/planet-snowball.png`). **Validation web:** the **North two-mode** analytic profile
  reproduced to ~1e-4 °C by the exact-face direct solve (the harmonic floor named at ~0.1 °C), the
  Strang relaxation shown to **converge** to it as dt→0, the **0-D mean** matching the discrete
  energy balance (net-TOA machine-exact) and the continuous `T̄` to the grid's O(1/n²) limit, and the
  direct operator **pinned to the frozen engine**. **Banked numbers:** present-day global mean
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
  (the rung-0 live loop; **obliquity** is a *named, deferred* slider awaiting its pinned `s₂(obliquity)`
  source). Pure consumer of `demo_biomes.compute` — no new physics. The **planet-spec** schema
  exports/imports the registry (JSON + `.npz`) with a **real round-trip-identity test** (the deep end's
  one genuine correctness property, vs the render's smoke-tests). Banked globe:
  `docs/figures/planet-map.html`. 31-test pair added (the Plotly render smoke-tests `importorskip` on
  the `[webviz]` extra — fast, not `slow`). **No new engine / no gate-manifest change.**
- **Phases 3–4 — shallow-water engine / coupler: PENDING** (plan §3). Phase 3 builds & freezes
  `engines/fluid` (the map then *registers* a `vector_overlay` circulation layer — no renderer edit);
  Phase 4 is the one-way coupler. The chip-style teaching notebook `planet.ipynb` is the remaining §9
  thin-skin surface (pending).

## Test runner (tiered gate, ADR 0003 — the per-project successor)

```powershell
python -m tools.gate planet -m "not slow"   # routine commit gate: planet's tests + the frozen engine's
python -m tools.gate planet                  # full gate for planet (incl. the slow demo sweep) — EXCEPTIONAL
./run_tests.ps1 -m "not slow"               # whole-repo fast lane (release / CI / shared-engine edit)
```

`pyproject.toml`'s `testpaths` already carries `projects`, so `projects/planet/tests/` is collected
with no config change; `pythonpath = ["."]` lets planet import the frozen engine as
`engines.diffusion…`. The full-resolution Snowball sweep (`test_demo_snowball`), the figure render,
and the live climlab cross-check are `slow`-marked / extra-gated, so the fast lane deselects them.
Planet's gate `uses` is `{engines/diffusion}` today; `engines/fluid` joins it in Phase 3 (making it
the manifest's first genuinely multi-engine row).
