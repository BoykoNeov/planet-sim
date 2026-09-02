# planet-sim — architecture

*How the repository is put together, why it is shaped that way, and how to grow it.* This is the
standalone repo's map. It replaces the BigSim monorepo's `ARCHITECTURE.md` that the plan, the engine
contracts and several docstrings still cite by section number (those references were written before
the 2026-06-10 split; where one says "ARCHITECTURE.md §N", the matching rule is restated here).

If you are looking for *what the simulator does*, start at the [README](README.md); for *the build
record and the growth axis*, [`docs/plans/planet-earth-system.md`](docs/plans/planet-earth-system.md);
for *the decisions*, [`docs/decisions/`](docs/decisions/). This file is the layer in between: the
structure and the rules of the road.

## 1. The shape in one paragraph

Two headless, separately-validated **solver engines** (`engines/diffusion`, `engines/fluid`) sit under
one **simulator package** (`planet/`) that composes them into a planetary climate. The simulator is
organised as a **staircase of rungs**: rung 0 is a complete, banked reduced climate (EBM + biomes +
shallow-water circulation + a one-way coupler + an interactive globe); every later rung is a **sibling
module** that adds one piece of physics beside the earlier ones, never by editing them, and banks one
demo figure with a validation anchor. The **demo catalogue** (`planet/catalog.py`) is the single source of
truth that the CLI launcher, the generated landing page and the drift-guard tests all read. Rendering is
strictly downstream of the numbers.

## 2. Directory layout

```
engines/                 the shared solver toolkit — headless, array-in/array-out, no planet imports
  diffusion/             1-D diffusion / heat solver (implicit FV, tridiagonal)   + CONTRACT.md + tests/
  fluid/                 rotating shallow-water (C-grid, SSP-RK3) + two-layer     + CONTRACT.md + tests/
planet/                  the simulator (flat package — one module per physical piece, one demo per rung)
  tests/                 the simulator's triads + demo integration tests + surface drift guards
  vendor/three.min.js    the one vendored asset (the WebGL particle globe), MIT
  planet.ipynb           the teaching notebook (a thin skin over the modules; executed by a slow test)
docs/
  plans/                 the build plan — the running record (§10 rung blocks, §12 backlog)
  decisions/             ADRs 0001–0005
  memory/                per-topic memory files + MEMORY.md index (the agent's project memory, vendored)
  explorations/          feasibility sketches (not built)
  figures/               banked artifacts: planet-*.png / .gif stills, planet-*.html interactive globes
  interactive/index.html the no-install browser what-if (precomputed grid, self-contained)
  index.html             the GENERATED landing page (python -m planet site) — byte-pinned by a test
```

`planet/` is deliberately **flat**: ~60 modules, each named for its physical content, so a session can
load exactly the working set it needs (the README's *load pointer* names, for each rung, the module +
test pair to open). Sub-packaging was considered and rejected: the module names already carry the
grouping, the notebook and catalogue import them by name, and moving files would cost a rename churn
across the plan, the memory and the notebook for no reduction in what a session must read.

## 3. The layers (and the one-way dependency rule)

Dependencies point **downward only**. A module may import from its own layer or below; never above.

| Layer | What lives there | Imports allowed from |
|---|---|---|
| **0 · engines** | `engines/diffusion`, `engines/fluid` — solvers behind a one-page `CONTRACT.md` | numpy/scipy only |
| **1 · physics** | the climate modules: `ebm`, `albedo`, `bifurcation`, `sphere_ebm`, `moist*`, `seasonal*`, `radiation*`, `orographic*`, `transport`, `eddy_flux`, `baroclinic_qg`, `circulation`, `coupler`, `precip`, `biomes`, `exoplanet`, `obliquity`, `ocean` | engines + layer 1 |
| **2 · state & interchange** | `planetmap` (layer registry), `planet_spec` (schema, round-trip identity), `flow_serialize`, `explain`, `interactive` (the precomputed grid) | layers 0–1 |
| **3 · render** | `plots` (matplotlib), the Plotly globe in `planetmap.render`, `flow_globe` / `eddy_globe` (three.js) | layers 0–2, plotting libs (opt-in extras) |
| **4 · surfaces** | `demo_*` (compute → print → bank a figure), `catalog`, `site`, `__main__`, `planet.ipynb`, `docs/interactive` | everything |

Consequences that are enforced by tests, not convention:

- **Engines never import a plotting library** and never import `planet` (ADR 0002 §1). The compute
  core and the whole test suite run on numpy + scipy alone; matplotlib / plotly / jupyter are
  `importorskip`-gated extras.
- **A figure is never evidence** (ADR 0002 §2). Every claim in a demo's printed summary is asserted by a
  test on the arrays first; the figure only draws them.
- **Every demo is a pure function chain** `compute() → print_summary() → save_figure()` so the slow
  integration test can run `compute()` and assert the story without rendering, and the render test can
  write to `tmp_path`.

## 4. The rung discipline — how a piece of physics is added

Every rung on the staircase (plan §5) lands the same way. Follow it and the repository stays coherent
for the next session; skip a step and the drift guards will tell you.

1. **Spike first, in `outputs/`** (git-ignored). Establish that the mechanism exists and that an anchor
   is reachable before writing the module.
2. **A sibling module, never an edit.** `sphere_ebm.py` did not modify `ebm.py`; `seasonal_map.py`
   *wraps* `seasonal.py` as its source of truth; `bifurcation.py` reads `ebm.py`'s operator. Reuse the
   engine's operator **assembly** (the pinned tridiagonal) rather than re-deriving it, so the new module
   cannot drift from the marcher. Default-off / opt-in flags keep the parent bit-identical.
3. **The validation triad** (the invariant every module carries — plan §3):
   - *tight*: an analytic limit or a reduction to the parent module, to machine precision or the
     scheme's order (and a **convergence** check when there is a step or a grid);
   - *conservation*: a budget that closes;
   - *loose (calibrated)*: the payoff number, stated with its cited source and an honest band.
   Say which is which in the test-file docstring. Do not flatten "loose" into "validated".
4. **Name the walls.** Each module docstring ends with its scope edges — what is prescribed, what a
   higher rung would derive. The plan's §12 backlog indexes them with status tags; strike a line there
   when you ship it.
5. **Bank one demo** (`demo_<rung>.py`) with a `slow` integration test that asserts the story, a figure
   under `docs/figures/planet-<rung>.png` (and a GIF/HTML when the mechanism is a motion), and a
   `Demo(...)` row in `planet/catalog.py`. Then `python -m planet site` and commit `docs/index.html`
   (the golden test fails otherwise).
6. **Record it** in three places, once each: the plan's §10 rung block (the full record) and §12 line
   (struck), a `docs/memory/<topic>.md` file plus a one-line hook in `MEMORY.md`, and the README's
   status list. Extend the notebook's §8 when the rung is cheap enough to run live or has a figure to
   embed (the notebook is executed by a slow test; embed-only cells cost nothing).

## 5. Numerical conventions worth knowing before touching a solver

- **Coordinates.** Latitude is carried as the area coordinate `x = sin φ` (equal Δx = equal area, so
  the global mean is a plain cell mean). Rung 0 uses the hemisphere `[0, 1]` with symmetry at the
  equator; the seasonal and full-sphere siblings use `[−1, 1]` because the seasonal cycle is
  antisymmetric. Longitude is cell-centred on `[0, 2π)` with no cell at the seam or at a pole.
- **Radiation is split around the engine.** The engine carries only linear transport; the radiative
  source (state-dependent through the albedo step) is composed by Strang splitting. This is exact for
  the linear sink but leaves the marcher's *fixed point* with an O(Δt) bias in the profile shape (the
  transport substep is backward-Euler). Where an equilibrium is the target, prefer the dt-free direct
  solvers (`ebm.steady_linear`, `sphere_ebm.steady_linear`, `seasonal.spectral`,
  `bifurcation.equilibrium_curve`) as the reference and treat the marcher as the *method*.
- **The ice edge is cell-quantised in the marcher.** A whole cell flips albedo; the exact diagram places
  the edge at a face by interpolation. They agree to within one cell width, and the gap halves with
  resolution — a Δx effect, not a Δt one. Near a fold the relaxation also suffers critical slowing, so a
  per-step tolerance can declare convergence early; use a tight `tol` and many iterations, or the exact
  curve.
- **Conservation is structural.** Neumann(0) / periodic finite-volume operators telescope to zero, so
  the transport conserves `∫C T dA` exactly; global-annual net TOA ≈ 0 is the test that the split did
  not break it.
- **Stated as measured.** The shallow-water scheme conserves energy semi-discretely, not enstrophy;
  claims say so. Do not upgrade a measured property to an aspirational one in prose.

## 6. Surfaces and their drift guards

| Surface | Source of truth | Guard |
|---|---|---|
| `python -m planet` menu / CLI | `planet/catalog.py` | `test_launcher` — every key resolves to a module with `main()` |
| `docs/index.html` landing page | generated from the catalogue by `planet/site.py` | `test_site` — committed page == fresh build |
| `docs/interactive/index.html` | `planet/interactive.py` precomputed grid | `test_interactive` (a slow byte-golden, CI-skipped) |
| `planet.ipynb` | thin skin over the modules | `test_planet_notebook` — executes clean top to bottom (slow) |
| Interactive globes `docs/figures/*.html` | `planetmap.save_html`, `flow_globe`, `eddy_globe` | render smoke-tests, `[webviz]`-gated |

## 7. Tests — the tiered gate (ADR 0003)

```
./run_tests.sh -m "not slow"    # the inner loop: every structural anchor, seconds to a couple of minutes
./run_tests.sh                  # the full gate: + demo integrations, notebook, render smoke-tests (CI)
```

`slow` marks anything that drives a live external solver, kernel or subprocess, or marches a 2-D field
to a limit cycle. The suite runs in parallel under `pytest-xdist` (`-n auto`, capped at half the
cores by the root `conftest.py`); every file-writing test uses `tmp_path`. CI (`.github/workflows/full-gate.yml`)
installs the optional stacks and fails loudly if one is missing, so a silent-skip cannot masquerade as a
pass.

## 8. Where the rules came from

- **ADR 0001** — Python + NumPy/SciPy, arrays as the boundary, compiled cores only when a scope ceiling
  demands it.
- **ADR 0002** — compute/render separation; the matplotlib floor is universal, interactive surfaces are
  per-need, the deep end (Plotly / WebGL) is opt-in.
- **ADR 0003** — the tiered test gate.
- **ADR 0004** — the interactive map's layer registry + the pin-a-schema interchange
  (`load(save(spec)) == spec` is the deep end's one real correctness property).
- **ADR 0005** — engines are living contracts: extend directly, test, changelog; no re-sealing.
