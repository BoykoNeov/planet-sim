---
name: planet-rung5b4-seasonal-sici
description: Rung 5B.4 BUILT 2026-09-04 — the seasonal small-ice-cap sweep (`seasonal_sici.py`): the SEASONS DISSOLVE the fold rung 0+ found (perennial cap grows one cell at a time through θ_c, zero hysteresis at 720 AND 1440 cells, where the annual-mean parent's own loop is ~7 W/m²), and DEEPENING the mixed layer brings it back; obliquity REJECTED as the sweep axis (it moves the parent too — no SICI at all at ε=0)
metadata:
  type: project
---

**Rung 5B.4 — does the small-ice-cap instability survive the seasons? · BUILT 2026-09-04**
(`planet/seasonal_sici.py`, `test_seasonal_sici.py` (13 fast + 3 slow), `demo_seasonal_sici.py` +
`test_demo_seasonal_sici.py` (2 slow),
`plots.seasonal_sici_figure` → `docs/figures/planet-seasonal-sici.png`; catalogue key `seasonal_sici`).
Closes the **"the *seasonal* SICI sweep stays named"** line that [[planet-bifurcation-sici]] left open when
it resolved the annual-mean SICI on rung 0, and the bifurcation study 5B.1+ named-deferred
([[planet-rung5b-seasonal]]). A **sibling**: `seasonal.py` is untouched except one `T_init` branch.

**THE FINDING — the fold is an artifact of the annual-mean idealisation.** The annual-mean parent of this
very model says a polar cap smaller than `θ_c ≈ 9.9°` is held by **no** sun (hysteresis loop
**6.8 W/m²**). March the *same* model with the seasons on at Earth tilt / 50 m mixed layer and none of it
survives: the perennial cap grows **one grid cell at a time** (0→1→2→…, `max_cell_jump == 1` at ΔS₀ = 0.5)
straight **through** `θ_c`, and the dimming and brightening legs retrace each other — largest down↔up cap
gap **0.18° at 720 cells, 0.08° at 1440**, both far under the ½-polar-cell detection threshold, so
**loop width = 0** where the parent's is ~7. Stable caps exist at 3.1°, 4.7°, 5.9° … — sizes the
annual-mean model forbids outright. **And it comes back:** plant a cap of `θ_c` at the parent's fold sun
and deepen the mixed layer (damping the seasonal swing without touching the annual-mean sunlight) —
gap 0.05° at 50 m (ONE climate), 0.13° at 200 m, **8.46° at 800 m** (TWO climates: the warm start finds
**no** cap, the planted one **survives** as a ~4-cell 8.5° cap). At **720 cells on purpose** — the first
version ran this control at 360, where a `θ_c` cap is a single cell, so "a planted cap survived" would have
meant one frozen cell versus none, on the same grid the resolution table calls *too coarse to tell*
(advisor-caught, pre-commit). Guard: `cells_in_cap(θ_c) >= 4` **and** the surviving cap must exceed one
polar cell. Direction banked; magnitudes ride the
calibrated `C`, `T_f`, `a_ice`. Agrees in direction with Wagner & Eisenman 2015 ("stability ... vastly
increases with ... a seasonal cycle") and with Huang & Bowman 1992's mechanism (the **amplitude of the
seasonal cycle**) — but is **not a replication** of either: step-function albedo on temperature alone (no
sea-ice thermodynamics), one tilt, uniform two-tile mix (no real geography). Banked as the **mechanism**.

**THE AXIS DECISION — obliquity REJECTED, mixed-layer depth chosen (the load-bearing call).** ε looks like
the natural "no seasons" axis and is **confounded**: tilt sets the *annual-mean* insolation gradient too,
so the **parent moves with it**. Measured: at **ε = 0 the annual-mean model has NO SICI at all**; at
ε = 15° `θ_c` is still shrinking at 2880 cells (2.45°, loop 0.93 and falling = resolution-limited, not
converged); only at ε = 23.44° does it converge (θ_c → 9.87°, loop → 6.76 at 2880). Depth moves **only**
the seasonal amplitude ⇒ one fixed converged reference for the whole sweep. *(Advisor first recommended
the obliquity axis; the ε-sweep spike overturned it — record the evidence, not the recommendation.)*

**READING IT WITHOUT BEING FOOLED BY THE GRID (the trap this rung is really about).** The polar cell's
*latitude* width is `≈ √(2Δx)` — **6.0° at 360 cells, 4.3° at 720, 3.0° at 1440** — and falls only as
`√Δx`, so 4× the cells halves it. Therefore *"the cap shrank smoothly to nothing"* and *"θ_c fell below the
polar cell"* are **the same data** on a coarse grid. Three guards, weakest → strongest:
1. **loop width in S₀** — resolution-robust (a real fold's interval converges; quantization's shrinks),
   defined against an explicit ½-polar-cell threshold so *"no loop"* ≠ *"a loop below detection"*;
2. **seed dependence at one sun** (`plant_cap` + march) — independent of sweep direction, because a
   warm-started continuation *can* walk continuously through a bistable band while pinned to one branch;
3. **the perennial ice-CELL COUNT** — the load-bearing one. The cap *radius* is **interpolated**
   (`ice_line_latitude` reads the T crossing → sub-cell values, e.g. a 3.10° cap inside a 4.27° cell) but
   the albedo feedback only ever flips **whole cells**, so a radius curve can look smooth over a stepping
   state. `+1 per S₀ step` is what makes "continuous" a fact rather than an interpolation — **but only
   against a scale.** A fold switches a whole `θ_c` cap on in ONE step: `cells_in_cap(θ_c)` = **5 at 720,
   11 at 1440, and exactly 1 at 360**. So at 360 cells "the count grew by one cell" and "a fold fired" are
   the *same observation* and the claim asserts nothing — caught while writing the payoff test, which had
   been making it at 360. The test now carries `assert fold_cells >= 4` as a standing guard, and the demo
   prints the fold scale in every row of the resolution table.

**Two traps found.** (1) **Seeding the seasonal marcher with the annual-mean parent's profile does NOT
work** — winter at high latitude runs tens of K below the annual mean, so it is effectively a *cold* start
and runs away to a **snowball** (first `seed_dependence` attempt: capped seed → 90° cap). Fix: perturb an
**already-converged seasonal limit cycle** (`plant_cap` cools only the cells poleward of the radius to
`T_f − 5 K`), so everything but the cap stays on the attractor. (2) **Critical slowing near a fold** —
a march that quits early looks like a *smooth shrink*, i.e. the one systematic error that biases this
experiment toward its expected answer; hence `tol=1e-7`, `max_years` up to 8000, and a per-point
`converged` flag that a fold must never be read across.

**The reduction that ties child to parent (tight).** `annual_mean_curve(cfg)` transplants rung 0+'s inverse
solve onto the *seasonal* model's **own** full-sphere `L_T` and its **own** annual-mean insolation (the
time-mean of `insolation_series()`, **not** the P₂-truncated `insolation()` — the 5B.1 lesson; it equals the
pinned `obliquity.annual_mean_insolation` to machine precision at `n_steps=720`), and returns a
`bifurcation.EquilibriumCurve` so the fold algebra (slope-stability, parabola-refined turning points,
`small_ice_cap_fold`, `finite_cap_window`) is **reused, not re-derived**. Anchors: the curve satisfies its
own steady EBM to ~1e-8 and hits `T_f` at the prescribed face to 1e-9; hemispherically symmetric to 1e-9;
**cross-model** — it finds the same fold as rung 0+'s independent hemisphere solve, differing by
**θ_c 9.85° vs 10.86°** and S₀ 1372.0 vs 1367.3 (0.34 %), which is the **P₄+ moments the P₂ truncation
drops — stated in advance, not discovered**; and the **ε=0 marched ice line converges onto the exact curve
at first order in dt** (gap 0.149 → 0.076 → 0.039 → 0.019 over n_steps 45→360, ratios 1.95/1.98/1.99) — the
Strang-splitting rate, the same shape as 5B.1+'s ε=0 self-consistency check. Conservation: a swept point
still closes global+annual net TOA with the realized co-albedo.

**THE INSTRUMENT'S POSITIVE CONTROL** (advisor-caught pre-commit; now
`test_the_loop_detector_reads_a_real_loop_when_one_exists`). Every *other* marched `hysteresis_loop` call in
this rung returns 0, so the headline would otherwise rest on a detector shown to read non-zero only on
**synthetic** arrays (the unit test). Pointed at an 800 m mixed layer the **same** detector and the same
sweep find the fold: **width 8.00 W/m² against the parent's 8.86** (7.50 at the cheap test settings, ~90 s),
in the right band of suns, with a sun inside the loop where the dimming leg carries **no** year-round ice
while the brightening leg still does. *A detector that always reads zero is indistinguishable from a broken
one* — reuse that framing on any future null result.

**Upstream edit (one, guarded).** `SeasonalEBM.march`'s `T_init` grew array and per-tile `(T_land, T_ocean)`
branches (the continuation seed) **beside** the untouched `None`/scalar paths; 5B.1+'s two bit-identical
reductions are the regression guard. Sources: North 1984 JAS 41 (the annual-mean SICI), Huang & Bowman 1992
*Clim. Dyn.* 7:205 (the SICI in seasonal EBMs — SH yes / NH no, via the seasonal amplitude), Wagner &
Eisenman 2015 *J. Climate* 28:3998, Cahalan & North 1979. → [[planet-rung5b-seasonal]],
[[planet-bifurcation-sici]], [[seasonal-ebm-source]], [[ebm-radiation-source]].
