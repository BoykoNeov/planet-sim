---
name: planet-phase2
description: "Planet Phase 2 — climate→biome map; diagnostic precip parameterization + Whittaker (T,P) biome classifier as an original sloped-boundary TOTAL partition (the Irvin graphical→computation precedent, NOT plotbiomes polygons), cm/yr units, Clausius-Clapeyron by global-mean T̄, band ordering equator→pole, warming migrates bands poleward"
metadata: 
  node_type: memory
  type: project
  originSessionId: 4a878190-f97f-4c5b-84c9-1e0d35fd121f
---

**Planet Phase 2 BUILT 2026-06-09 — the climate→biome map** (the payoff, banked EARLY per
the locked scope decision): `projects/planet/` (precip + biomes + demo_biomes +
plots.biomes_figure, 20-test triad incl 3 slow; planet gate **65 +1 skip**, full suite
**412 green** at build, 11 figures).

`precip.py` = a **diagnostic** precip parameterization (**NOT** a water cycle):
`P(φ,T̄) = pattern(φ)·CC(T̄)` — a Gaussian latitude pattern (ITCZ-wet / subtropics-dry /
midlat-wet / poles-dry, **cm/yr**) × a Clausius–Clapeyron **global-T̄** amplitude.
`biomes.py` = a **Whittaker (T,P) → biome classifier** = an **original, TOTAL,
sloped-boundary rule partition** (9 biomes), **NOT** embedded polygons. **NO new engine**
(project-local EBM reuse only — gate `uses` stays `{diffusion}`).

**Durable advisor calls:** (1) the **classifier = an original partition, NOT plotbiomes
polygons** = the **[[irvin-sheet-resistance-source]] precedent** (a copyrighted GRAPHICAL
diagram → independent calibrated computation, not digitized) → also total-by-construction
(no nearest-polygon fallback → area-fractions-sum-to-1 holds). (2) **UNITS = cm/yr**
(Whittaker's axis; **mm/yr was the trap** → 10× too wet → all rainforest) — both modules
pin cm/yr. (3) **C-C scales by GLOBAL-mean T̄, NOT local T** (a deliberate refinement of the
advisor's "local T": local-T-vs-fixed-ref over-amplifies the warm equator ×2.6 & breaks
calibration; pattern = circulation = FIXED, only the global moisture amplitude responds → an
honest pattern/amplitude split). (4) **boundaries DIAGONAL** (precip thresholds linear in T)
not axis-aligned — a forest needs more rain when warmer; cold biomes (tundra < −5, boreal
−5..3) = vertical T-limited cuts.

**Non-circularity:** the thresholds are *calibrated to* the diagram (loose); the **probe
points = independent canonical facts** (rainforest = warm+wet, tundra = cold, desert =
warm+dry). **Validation:** exact total partition + area-fractions = 1 (a consistency leg,
**NOT** a conservation law — named), 9 probes land, present-day **band ordering equator→pole
= rainforest → savanna → desert/grassland → temperate-forest → boreal → tundra** reproduced,
⟨P⟩(T̄) **monotone**; banked figure = the biome-band map + the **Whittaker (T,P) plane with
the climate trajectory through it** (mechanism: why deserts @ ~30°) + warming **migrates
bands poleward** (CO₂ knob ≈ −8 W/m² in A → +5 °C; rainforest 7→14 %, tundra 10→1 %;
CO₂-doubling = 4 W/m² so −8 ≈ a strong scenario).

**Scope edges NAMED:** Whittaker-not-Köppen (annual T,P, no seasonal precip),
prescribed-not-derived precip, **fixed band centres** (migration = circulation =
rung-1/2 enhancement), **7 %/K = moisture-capacity NOT the energy-constrained ~2–3 %/K
global precip rate** (rung-2). Sources → [[whittaker-biome-source]] +
[[precip-parameterization-source]]. See [[planet-plan]]; siblings [[planet-phase1]] /
[[planet-phase3-engine]] / [[planet-phase4-coupler]].
