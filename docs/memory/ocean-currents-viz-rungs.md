---
name: ocean-currents-viz-rungs
description: Ocean-currents showcase rungs O1–O5 (plan §9.6) — O1 mask + O2 OSCAR producer + O3 beauty pass BUILT 2026-07-06 (real currents banked on the Rung-C globe, land/ocean base + trails + speed styling, FlowField untouched a 4th time); O4 frames next; §11.4 fork retargeted viz-half-only
metadata:
  type: project
---

**Scoped 2026-07-06; O1 BUILT 2026-07-06 (see below) — plan §9.6 is the record.** The user's "visualize beautiful ocean
currents" ask exercised §11.4's *named alternative* under the living-staircase rule: real-data ocean
ingest-for-**display** moves INTO planet-sim (the renderer + seam live here), while the engine/forcing
half of the spin-out boundary is unmoved (no ocean physics here; ClimaOcean stays S2–S5; spin-out **S1
narrows** to ECCO-as-validation-anchor-for-S3). The §11.2 "never ships an ocean visual" line is amended
(visual yes, engine never); status notes threaded at §11.2/§11.3-S1/§11.4 + backlog §12.3.

**The ladder (order O1→O2→O3→O4; O5 independent):**
- **O1 — mask increment: BUILT 2026-07-06** (all three spike retargets in). `FlowField.mask` (True=valid,
  `None`=all-valid=pre-O1-unchanged); Coverage stays the bounding box. Renderer: mask on the velocity
  texture's free 4th channel `(u,v,θ,mask)` — GPU recycles land-drifting/land-respawned particles via
  the *invisible-retry* idiom (shader can't loop a rejection sample; ~2 frames at 44% land), CPU
  rejection-samples spawns + recycles in step(); both share the bilinear+0.5-threshold coastline rule.
  Serialize: `_nearest_mask` = nearest-neighbor (boolean, coastline-crisp) + lat-beyond-source→False
  (pole honesty; lon edge-clamps, periodic) + mask **applied** (zero u,v / NaN scalar on invalid cells);
  rides as an additive categorical 0/1 `mask` layer in the same `.npz` (`render(active="mask")` = a free
  coverage globe). Round-trip `==` extended to an OSCAR-shaped masked probe. NaN→0 fill = O2's job.
- **O2 — real-ocean producer: BUILT 2026-07-06** (the deliverable + S1 de-risk; spike settled auth+format
  earlier same day). `planet/ocean_currents.py` NumPy-only at import: `load_oscar` (lazy `h5netcdf`, new
  `[ocean]` extra; owns `(time,lon,lat)`→`(lat,lon)` transpose, `_FillValue`→NaN, `stride`,
  `geostrophic_only`) → `OceanSnapshot` (convention-raw: 0–360 lon, NaN land, provenance from granule
  attrs) → `flow_field_from_ocean` (±180 rewrap + argsort re-sort → monotone axis; mask-from-finiteness
  THEN NaN→0 fill — mask carries "no data", never a filled zero; scalar=speed; `is_global=True`, coverage
  honestly ±89.75°). **Both conscious calls made:** render = **0.5°** (stride 2; 2° interchange stays
  proof-only) → banked HTML **5.2 MB**; pace = additive default-off `crossing_seconds` knob on
  `flow_globe_html` (eddy artifact byte-unchanged; ocean passes 45 s/360°). Provenance clause ("REAL
  data … OSCAR L4 v2.0 PO.DAAC DOI 10.5067/OSCAR-25F20 … NOT computed by planet-sim's models")
  machine-checked in the DOM **and on the committed artifact**. Demo `demo_ocean_currents.py`
  (catalogued): `EARTHDATA_TOKEN` bearer download to gitignored `outputs/`, asserts the **R1 round-trip
  on the real field** before banking `docs/figures/planet-ocean-currents.html`. Tests (17) off a
  committed **14 KB 5° fixture** `planet/tests/fixtures/oscar_subsample.npz` (raw conventions kept so
  rewrap/mask/poles are exercised; dir named `fixtures/` because `.gitignore` ignores any `data/`) +
  synthetic-granule loader test (`importorskip`). Gate 531 green. Spike facts (auth = one
  `Authorization: Bearer` header, no `.netrc`; `cftime`/`julian` calendar relevant only at O4) recorded
  in plan §9.6. **O3 forward flags observed:** particles recycle (not wrap) at ±180°; uniform-in-lat
  spawning over-densifies high lat on a globe → cos-weighted seeding belongs in O3 with speed styling.
- **O3 — beauty pass: BUILT 2026-07-06** (renderer-only, `FlowField` untouched a 4th time; `flow_globe.py`
  only, 3 commits + eddy re-bank). **(a)** land/ocean base is **honest-by-construction** NOT a CanvasTexture
  — a base fragment shader inverts each point to `(lat,lon)` with the SAME `sph()` mapping (`atan(n.z,n.x)`)
  + samples the SAME mask on the SAME 0.5 rule as `validAt()`, so coast-under-particles can't drift from
  coast-under-base; no-mask/compile-miss → solid sphere. **(b)** trails **default-OFF** behind a kwarg (no
  WebGL CI + blind hand-off ⇒ ocean globe exercises them first, eddy can't silently regress); the
  **depth-only occluder prepass is load-bearing** (kills back-hemisphere particles BEFORE the accum buffer
  ⇒ nothing bleeds through the planet); **additive** One+One accumulation (not alpha-over) sidesteps
  premultiplied fringing + IS the glow; **rotation-smear fix = `decay=0 while dragging`** (screen-space
  buffer smears when projection moves → history pauses mid-drag, resumes still — a NAMED property, not a
  bug). New shaders `compileOK`-gated, RTs try/catch + resize-realloc, any miss → plain single-pass (never
  the CPU fade-only fallback). **(c)** colour path was ALREADY producer-driven (eddy=θ/ocean=speed) so
  nothing to recolour — the new part is **speed-weighted seeding in the RESPAWN path** (seed-only relaxes to
  uniform as particles age), composing w/ the mask reject via the invisible-retry idiom, floored for calm
  water; + a **sequential** speed ramp opt-in default (RdBu_r diverging bleaches 0→max; stays default for θ,
  ocean opts in). Both §9.5 knobs shipped: **density** (GPU rank-cut/CPU tail-hide) + **trail length**
  (decay). Verified on the committed **5° OSCAR fixture** (masked pipeline, node-`--check` clean), gate
  **538 green**. The 0.5° banked ocean artifact re-bank needs `EARTHDATA_TOKEN` (user hand-off; code +
  fixture prove it). Advisor's 3 load-bearing calls (occluder-prepass / additive / respawn-not-seed) all in.
- **O4 — frames time axis**: the R1 deferral, built after O2's real dims are seen; `(nt,ny,nx)` stacked
  npz + crossfaded `uVelA`/`uVelB` textures; OSCAR monthly climatology → Somali Current monsoon reversal.
- **O5 — QG producer** `flow_field_from_qg` (independent): second EMERGENT producer, box coverage no mask;
  a third producer re-trips §9.4 rule-of-three for the two-consumer geometry helpers.

Out of scope, restated: no ocean engine, no forcing seam (S2 designs against the SEEN ClimaOcean API), no
real-time what-if on real data. [[planet-spinout-roadmap]] [[planet-viz-animation-rungs]]
[[planet-rung3-qg-built]]
