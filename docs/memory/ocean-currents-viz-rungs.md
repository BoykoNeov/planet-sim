---
name: ocean-currents-viz-rungs
description: "Ocean-currents showcase rungs O1–O5 (plan §9.6) — ALL FIVE BUILT + BANKED 2026-07-06. O1 mask + O2 OSCAR producer + O3 beauty pass + O4 seasonal frames (browser-verified) + O5 QG emergent producer (flow_field_from_qg; rule-of-three re-affirmed HOLD, no transpose, explicit 45° display latitude not derived, fresh irreversible-not-reversible honesty); §11.4 fork retargeted viz-half-only"
metadata: 
  node_type: memory
  type: project
  originSessionId: 2c3400c6-1bf1-4351-9bb8-4b7944907744
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
- **O4 — frames time axis: BUILT + VERIFIED + BANKED 2026-07-06.**
  All four pieces built, fast gate **552 pass/1 skip**, **default-off so pre-O4 is bit-for-bit** (single-snapshot
  shaders byte-untouched). **Contract** `FlowFrames` on `FlowField` (`u`/`v` `(nt,ny,nx)`+`labels`; `frames=None`
  = exact pre-O4 path — O1-mask/O3-trails discipline). **Producer** `flow_field_from_ocean_series`: rewrap→mask→fill
  over N snaps, stacked; **static mask = finite-in-EVERY-frame** (advisor: conservative, never blinks; sea-ice
  folded into "no data in any frame → bare"); **acquisition-agnostic** (advisor: producer never says "climatology"
  — caller passes the honest `period` phrase). **Renderer** (GPU-only like trails): **separate** crossfade shaders
  `UPDATE_FS_F`/`DRAW_VS_F` (`velAt()=mix(velA,velB,uMix)` substitution — single-path `UPDATE_FS`/`DRAW_VS`
  untouched), N frame textures once, `stepSeason` cyclic **Dec→Jan `(k+1)%NT`** wrap, `seconds_per_year` pace,
  **live month time badge** (the showpiece). **Advisor payload calls taken:** colour = in-shader mixed speed ⇒
  **per-frame scalar dropped** (−⅓); animation grid **coarser than the O2 still** (demo `STRIDE=6`=1.5°, motion
  hides res) ⇒ 12-frame HTML ~5 MB (12×`(u,v)`@0.5° would be ~30 MB). **Serializer**: the R1 frames deferral acted
  on — `(nt,2,ny,nx)` stack rides as ONE additive `FRAMES_LAYER` VECTOR_OVERLAY (labels in style), round-trip `==`
  free; interactive map **skips the 4-D stack** (one-line `planetmap._overlay_traces` `ndim>=4` guard — paints the
  primary snapshot). **Demo** `demo_ocean_seasonal.py` (bearer token, 12 mid-month-day granules of 2020). **Frame
  data (user 2026-07-06): 12 monthly SNAPSHOTS = one day/month of 2020, NOT means/climatology** (Somali reversal
  reads in a day-per-month series; label says exactly that). **The honesty gap (advisor):** `node --check` ≠ GLSL
  compile + no WebGL CI ⇒ the whole frames GPU path (crossfade shaders, stepSeason, badge, trails+frames) ran
  NOWHERE; a frames-shader compile error degrades **silently** to CPU fallback (static frame-0 = reads as "not
  animating"). ⇒ **that gap is CLOSED by the owed browser play-through (user 2026-07-06, PASS): month badge
  cycles Jan→Dec, particles stream, Somali reversal reads — GPU crossfade compiled + ran, no CPU fallback.**
  Artifact **banked** `docs/figures/planet-ocean-currents-seasonal.html` (**4.23 MB**, real 2020 OSCAR series:
  52% valid-in-every-frame ocean, |current| max 2.83 m/s, round-trip identity OK); demo **catalogued**
  (`ocean_seasonal`, "Interactive globes") + on the landing page (`python -m planet site` regen). Rung-B/C
  eddy-band animation back-port = deferred-not-owed.
- **O5 — QG producer: BUILT + BANKED 2026-07-06** (independent bonus, did not gate O1–O4).
  `flow_field_from_qg` (`flow_globe.py`, beside `flow_field_from_eddy`) = the **second EMERGENT** producer:
  the rung-3 two-layer QG condensate (coherent vortices + PV filaments, [[planet-rung3-qg-built]]) as
  particles. `(u,v)` from the model's own `invert`→`velocities`; **advisor: axes already match the contract
  → NO transpose** (unlike OSCAR); colour = **upper-layer PV anomaly `q₁`** (signed → diverging RdBu_r, the
  eddy-θ twin) **÷f₀** (advisor caught PRE-COMMIT: raw PV is `O(1e-4/s)` → `_build_data`'s 3-dp `flat(scalar,3)`
  rounds it to a constant 0 = every particle one flat colour, erasing the vortex field; ÷f₀ = an `O(1)`
  Rossby-like field that survives the rounding, monotone so RdBu_r still centres on 0; fix in the PRODUCER
  not the shared renderer; a payload-fidelity regression test now pins a non-flat rendered scalar — the class
  the Python-side `np.allclose`+disclaimer-grep tests missed); **box coverage, no mask, no frames = the plain
  pre-O1 contract shape.** Two advisor content
  calls, both taken: (1) display latitude **explicit `center_lat_deg=45°`, NOT derived from f₀/β** — the
  idealized `(f₀,β)` are independent knobs, not a consistent `(sinφ,cosφ)` pair (demo's f₀→~43°, β→~44°), so
  `atan(f₀/βa)` would *manufacture* a latitude never put in; box maps by `Δlon=Δx/(a cosφ_c)`, `Δlat=Δy/a`,
  centred, its honest ~box-width sector, never 360°-wrapped. (2) honesty string **FRESH, NOT the eddy's
  reversibility clause** — the rung-3 win is this saturated flux is *irreversible* (persistent turbulence),
  so "~90%-reversible/mostly-sloshes" is *false* here; carries the QG edges (idealized model / not real data
  / box-not-planet-wide / inverse-cascade condensate) + names the colour (PV). QG-specific honesty test pins
  the fresh clauses + asserts the eddy's absent (the generic disclaimer machine-check runs only on the eddy
  field, doesn't misfire). **§9.4 rule-of-three (the architectural side-effect): third *geometry* consumer
  arrived → re-affirmed HOLD** — `flow_field_from_qg` **cannot** call `_band_geometry` (takes a frames obj
  with `.phi/.x/.y` the metre-space box lacks) + never touches `_sphere_xyz` (renderer-side), so the sector
  formula is *inlined*; extracting a shared helper would force the banked eddy path to recompute lat from `y`
  (ULP-risk on banked art) = the pre-emptive promotion R2 forbids. Resolution = **documentation, not
  extraction** (mirrors R2). **R1 round-trip = named deliverable, PASSES** (plain box shape → `save`/`load ==`
  free; a `test_flow_serialize` case pins it as the third producer). Demo `demo_qg_particles.py` (catalogued
  `qg_particles`, "runs a short sim") reuses `demo_baroclinic_qg` params but keeps the **model** (renderer
  needs the grid metric), runs to a saturated `v'≫U_s` condensate (nx=96), banks
  `docs/figures/planet-qg-particles.html`. Producer tests fast-lane (cheap `random_state`, no spin-up); bank
  = slow smoke-test. **Verification LIGHTER than O4** (advisor): reuses the already-browser-verified
  single-snapshot GPU path — no new shaders, so O4's "GLSL-never-runs-in-CI" hand-off gap does not reapply.

Out of scope, restated: no ocean engine, no forcing seam (S2 designs against the SEEN ClimaOcean API), no
real-time what-if on real data. [[planet-spinout-roadmap]] [[planet-viz-animation-rungs]]
[[planet-rung3-qg-built]]
