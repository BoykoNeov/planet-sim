# planet-sim — a planetary climate simulator

[![full-gate](https://github.com/BoykoNeov/planet-sim/actions/workflows/full-gate.yml/badge.svg)](https://github.com/BoykoNeov/planet-sim/actions/workflows/full-gate.yml)

*Stellar flux + planet parameters in, climate, circulation, and biomes out.* An educational
simulator for planetary climate: a latitudinal energy-balance model (EBM) with the Snowball-Earth
bifurcation, a Whittaker biome classifier, rotating shallow-water atmospheric circulation, a
one-way EBM→circulation coupler that grows an emergent jet, exoplanet knobs (stellar spectrum,
planet size, obliquity), and an interactive Plotly globe — each leg validated against cited
climate references. Beyond that rung-0 core a **research staircase** climbs toward a GCM one validated
rung at a time: a moist water cycle, an energetic ITCZ, baroclinic eddy turbulence, a spectral greenhouse,
orographic rain shadows, a seasonal land–sea map with winter snow, and the **complete equilibrium
diagram** of the ice-albedo climate — every branch, both folds, the small-ice-cap instability.

It is the program's capstone: the **only simulator built on two** separately-validated solver
engines — the 1-D diffusion/heat solver (`engines/diffusion`, the meridional heat transport) and
a rotating shallow-water solver (`engines/fluid`, the circulation). Each engine is a *living,
versioned* contract ([ADR 0005](docs/decisions/0005-engines-are-living-contracts.md)) with its own
`CONTRACT.md` and test suite — extended directly as the simulator grows, the suite as the guardrail.

## Layout

```
engines/diffusion/   # the 1-D diffusion/heat solver (+ its own tests)
engines/fluid/       # the rotating shallow-water solver (C-grid / SSP-RK3) (+ its own tests)
planet/              # the simulator: ebm, albedo, climate_reference, biomes, precip,
                     #   circulation, coupler, exoplanet, obliquity, planetmap + planet_spec,
                     #   plots, demos, planet.ipynb
ARCHITECTURE.md      # the repo map: layers, the one-way dependency rule, how a rung is added
docs/decisions/      # ADRs 0001–0005 (incl. 0004 interactive maps, 0005 engines are living contracts)
docs/plans/          # planet-earth-system.md — the full build plan
docs/figures/        # banked figures (planet-*.png / .gif) + the interactive globes
                     #   (planet-map.html, planet-coupler-map.html)
```

## Quickstart

```powershell
pip install -e ".[viz,webviz,notebook]"   # compute + figures + interactive globes + notebook
python -m planet                          # the front door: a menu of every demo, globe & the notebook
```

`python -m planet` is the **one command to remember** — it lists every demonstration and runs
or opens them for you (including launching the notebook straight into your browser). Want to go
straight to one?

```powershell
python -m planet snowball     # run one demo (prints its validation table + banks a figure)
python -m planet bifurcation  # every climate the sun allows — the S-curve, both folds, the second cliff
python -m planet seasonal_ice_map   # the seasons, the continents and the snow — a monthly GIF + a month-slider globe
python -m planet list         # print the full catalogue of demos
python -m planet interactive  # drag four knobs in your browser — a live what-if + plain-language "why"
python -m planet notebook     # open the teaching notebook in JupyterLab (opens your browser for you)
python -m planet globes       # just open a saved interactive globe — no compute
python -m planet site         # build & open the landing page — a clickable gallery of everything
```

**Turn a knob, build a climate — no install.** `python -m planet interactive` opens
[`docs/interactive/index.html`](docs/interactive/index.html): drag the Sun, the greenhouse, the axial
tilt, and how much of the world is ocean, and the planet's temperature, polar ice, and bands of life
respond instantly, with a plain-language explanation of *what changed and why* (`planet/explain.py`).
You can also flip a world's *starting climate*, warm or frozen, to meet its bistable twin — at today's
Sun a frozen start stays a Snowball where a warm start is temperate. It's a lookup over a precomputed
grid of the real `planet.demo_biomes.compute` runs — instant, deterministic, and **self-contained** (the
data is inlined, so it opens straight off disk and serves from GitHub Pages alike). Want continuous
knobs, live re-runs, and the *full* hysteresis loop — the catastrophic freeze and the late re-melt?
That lives in the notebook.

**Browse it as a webpage.** `python -m planet site` (re)generates [`docs/index.html`](docs/index.html)
— a self-contained, clickable gallery linking every demo figure, the three interactive globes, and
the notebook. It's generated from the demo catalogue (`planet/catalog.py`), so it never drifts: a
test fails the build if the committed page is out of date. Open it straight off disk, or publish the
`docs/` folder via **GitHub Pages** (Settings → Pages → *Deploy from branch* → `main` / `/docs`) for a
shareable site — the figure and globe links are relative, so they resolve either way.

Each demo banks its figure under `docs/figures/` (and `outputs/`). The render stacks are opt-in —
`.[viz]` for the PNG/GIF figures, `.[webviz]` for the interactive Plotly globes, `.[notebook]` for
JupyterLab — and a demo whose stack is missing still prints its physics summary and tells you which
extra to install rather than erroring. The notebook also hosts the live-slider globe
(`planet.planetmap.interactive_map()`), the one view that needs a running kernel.

**Run the tests** (the tiered gate — [ADR 0003](docs/decisions/0003-test-execution-policy.md)):

```powershell
./run_tests.ps1 -m "not slow"     # routine fast lane — 661 tests
./run_tests.ps1                   # full suite — 729 tests (adds slow live-solver + notebook)
```

The suite is **729 tests** (661 in the fast lane), all green. The one **live-climlab** cross-check needs the
`[climate]` extra and otherwise skips — it is an opt-in bonus on top of the EBM's analytic +
frozen-table validation, so it skips in CI by design. The Plotly map render smoke-tests need
`[webviz]`; the planet-spec round-trip-identity test (the deep end's one real correctness
property) is NumPy-only and always runs. Optional stacks are importorskip-gated, so a headless
checkout skips rather than errors.

## Provenance

planet-sim was developed inside the **BigSim** monorepo — an educational program of three
simulators (steel, microchip, planet) sharing two separately-validated solver engines — then
extracted into a standalone repo with its history. The diffusion/heat engine was first built &
validated by the steel simulator; planet added the shallow-water engine and is the capstone that
couples the two. The sibling simulators live in their own repos. The archive:
[github.com/BoykoNeov/BigSim](https://github.com/BoykoNeov/BigSim).
