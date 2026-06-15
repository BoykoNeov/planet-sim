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
- `planet/interactive.py` → `docs/interactive/index.html` — no-install browser what-if (originally 2 knobs:
  Sun S0 × Greenhouse CO2; **3rd axis obliquity/tilt BUILT 2026-06-14**, see below). Self-contained: data
  **inlined** (file:// blocks fetch), vanilla canvas (latitude-banded planet disk + T(lat) curve), no CDN.
  `python -m planet interactive`; hero card on the landing page. Wired in `__main__`/`site.py`/README.
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
- Explanations baked in Python (no JS climate rules) → no drift. 2-knob v1 was structured so a 3rd
  axis (tilt) is a trivial add; D/star/size stay in the notebook bench.

**3rd axis — obliquity/tilt — BUILT 2026-06-14** (was the "[deferred — ~free]" §12.3 item; built it +
reconciled the *stale* §12.3 "Rung-C GPU advection" line in the same pass — GPU ping-pong had already
shipped 2026-06-13, the backlog read pending). The trivial-add prediction **held**: only the data axis
moved. A tilt slider over **0…45°, 9 values, capped at `OBLIQUITY_FAITHFUL_MAX`** (reuse the constant so
it can't drift from the plots.py faithful-band shade), **including the exact `OBLIQUITY_EARTH` float** so
the obliquity factor is exactly 1 there → `s2` bit-identical → the (1365, 0, 23.44°) cell stays the
baseline. Wired `obliquity_params(obl, EBMParams(...))` → `EBMParams.s2`; narrated by the **existing**
`explain.py` `obliquity_deg` rules (it already carried the knob → **zero prose work**, confirmed before
building). Advisor's load-bearing call = **axis SHAPE, not size, is the constraint**: S0 has the snowball
cliff + Earth detent → **left S0/CO2 untouched**; obliquity is smooth/cliff-free → coarse 9 is plenty.
Size held ~4 MB (≈ the eddy-globe precedent) via **`_LAT_STRIDE` 3→6** (free — 30 lats indistinguishable on
a 300px canvas), *not* axis-coarsening (regression risk) or prose-dedup (complexity, exact-number prose has
few true dupes). JS decode `cells[(i·nCo2 + j)·nObl + k]` must track the `for s0: for co2: for obl` loop
order. Grid ≈3.7 k solves (~10 min regen); slow byte-golden stays CI-skipped; +1 fast test
(`test_obliquity_axis_moves_the_climate`: flat 0° world → more ice + colder than Earth's tilt; prose names
the knob). D/star/size still notebook-only.

**4th axis — ocean fraction — BUILT 2026-06-15** (`3612377`; the user's "water knob" question, resolved to
**option B = water that affects TEMPERATURE**, over option A = a biome-only precip-amplitude multiplier).
**NOT a free reuse** like obliquity (which already had its `explain.py` clause): ocean needed a NEW cited
`planet/ocean.py` param-derivation, a NEW `Knobs.ocean_fraction` field + cause/gradient clause + an honesty
caveat. `ocean_params(w, base)` maps sea fraction → **two** EBM params it already accepts — `a0` (FIRM leg:
wetter=darker=warmer, Donohoe&Battisti 2011, full-range planetary-albedo swing 0.07) and `D` (LOOSE/flagged
leg: more ocean spreads heat poleward, Trenberth&Caron 2001, all-land −25%/all-ocean +10%); **commutes with
`obliquity_params`** (disjoint: s2 vs a0/D). Earth's **0.71 = exact identity** → baseline bit-for-bit. **No
`Knobs.D` double-count** (ocean drives D internally; the ocean clause narrates it). Honesty ceilings surfaced
in prose + slider hint: heat-capacity/seasons INVISIBLE (equilibrium annual-mean), rain PATTERN fixed (precip
= prescribed pattern × C–C, no ocean source) → a wetter world rains more *everywhere*, not somewhere new.
[[ocean-albedo-transport-source]] pins the citations. **The obliquity "only the data axis moves" pattern did
NOT fully hold**: ocean moves TEMPERATURE, so every ΔT/ice/biome/prose field refreshes per ocean value →
**multiplicative on the page** (a biome-only knob would have shared the temperature fields, far cheaper). And
the obliquity build's size lever was already spent: prose is **77% of the page and onelines/paragraphs are
100% UNIQUE** (they embed exact ΔT/ice/%), so **string-interning is dead** and `_LAT_STRIDE` can't drop
again. So the budget was bought elsewhere: **CO2 trimmed 17→9** (2 W/m², advisor: over-resolved at 1). Axis
= **5 ocean steps {0, 35, 71, 85, 100}%**, including the exact `OCEAN_FRACTION_EARTH`. JS decode is now 4-D
`cells[((i·nCo2+j)·nObl+k)·nOcean+l]` (track the `for s0: for co2: for obl: for ocean` order). **Size: 18.4
MB** — my ~14 MB estimate **undershot by 30%** because ocean prose runs ~1890 vs ~1450 bytes/cell; user
**confirmed keep-5 after I flagged it** (revealed preference = resolution over size, since they'd already
picked the richest option). ~9.7 k solves (~regen minutes); slow byte-golden still CI-skipped; +`test_ocean.py`
+ ocean clauses in `test_explain.py` + `test_ocean_axis_moves_the_climate`. Browser eyeball owed →
[[pending-interactive-hover-check]]. D/star/size still notebook-only.
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
