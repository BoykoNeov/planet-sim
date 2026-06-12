---
name: interactive-what-if
description: The no-install browser what-if + explain.py shared prose engine + notebook enrichment (built 2026-06-12); design rationale and the deferred §2/§4/§5 live widgets
metadata: 
  node_type: memory
  type: project
  originSessionId: 14b69ce7-9ac6-49bc-b843-e078d2bad808
---

User asked (2026-06-12) for "more interactiveness — change parameters, run the computation, see
what changes + explanation"; audience = **users**, richer prose at **both** depths, in **browser AND
notebook**. They had not run the notebook (so never saw its existing widgets) → the browser page is
the real user front door.

**BUILT (commits db899d0 / 0b27556 / 2694750):**
- `planet/explain.py` — the SINGLE rule-based source of "what changed + why" prose, keyed to the
  *computed* deltas (never overclaims). Two depths: `oneline` causal chain + `paragraph` mechanism.
  `Knobs`/`Diagnostics`/`diagnose()`/`explain()`. Used live in the notebook AND baked into the
  browser page → can't drift.
- `planet/interactive.py` → `docs/interactive/index.html` — no-install browser what-if (2 knobs:
  Sun S0 × Greenhouse CO2). Self-contained: data **inlined** (file:// blocks fetch), vanilla canvas
  (latitude-banded planet disk + T(lat) curve), no CDN. `python -m planet interactive`; hero card on
  the landing page. Wired in `__main__`/`site.py`/README.
- Notebook: §3 + §7 what-ifs render `explain_panel()` (live on every drag) in place of terse prints;
  `KNOB_STYLE` (gold non-white handle + `description_width:"initial"` so labels stop truncating) on
  §1/§3/§7 sliders; **all code cells `source_hidden`** (collapsed by default). Re-banked via a clean
  headless Run-All (`store_widget_state=True` keeps widgets rendering).

**KEY DESIGN DECISIONS (non-obvious):**
- Browser = **precomputed grid** of real `demo_biomes.compute` runs (NOT Pyodide/Flask) — instant,
  deterministic, golden-testable, shareable; ~493 cells / ~70s / ~0.5MB inlined. Earth detent exact
  on S0=1365, CO2=0. Snowball cliff lives between S0≈1245→1265 (drag freezes the planet).
- Hysteresis can't be a lookup → browser **names** the path-dependence and defers the live two-branch
  demo to the notebook §2 (use `present_day_climate(ic_equator, ic_pole)` for warm/cold starts).
- Explanations baked in Python (no JS climate rules) → no drift. 2-knob v1 is structured so a 3rd
  axis (tilt) is a trivial add; D/star/size stay in the notebook bench.
- The slow byte-exact golden (`test_committed_page_is_up_to_date`) is **CI-gated** (`_SKIP_IN_CI`)
  like the notebook test — 408 live EBM solves compared byte-for-byte is fragile cross-platform
  (LAPACK last-bit near a Whittaker threshold flips a biome-string digit). Fast structural tests
  cover CI. Gate is `M:\claud_projects\planet-sim\planet\tests\test_interactive.py`.

**DEFERRED (offered, user hadn't explicitly asked):** NEW live widgets on the still-static §2
snowball (the genuine live hysteresis / two stable states), §4 winds, §5 jet. Advisor flagged: time
the §2 continuation (~3s) / coupler sims before promising "live" — lean on `continuous_update=False`
or precompute. Minor open: browser snowball stat shows "Tundra 100%" while prose says "bands
collapsed" (honest model output, mild tension). See [[planet-plan]], [[pedagogy-novice-intermediate]],
[[viz-prose-novice-intermediate]].
