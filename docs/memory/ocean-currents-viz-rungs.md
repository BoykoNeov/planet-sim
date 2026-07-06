---
name: ocean-currents-viz-rungs
description: Ocean-currents showcase rungs O1–O5 (plan §9.6) — O1 mask BUILT 2026-07-06 (contract+renderer+serialize, all O2-spike retargets in); O2 producer next (data spike DONE, auth+format settled); real OSCAR/ECCO currents through the R1 seam onto the Rung-C globe; §11.4 fork retargeted viz-half-only
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
- **O2 — real-ocean producer** (the deliverable + S1 de-risk): one OSCAR 0.25° snapshot (PO.DAAC) →
  `flow_field_from_ocean` → full-globe `is_global=True` (first true-global consumer) → serialize (round-trip
  on REAL data) → unchanged Rung-C globe → `docs/figures/planet-ocean-currents.html`. Honesty flips class:
  a new **provenance clause** ("real reanalysis-class currents, NOT this model's output"), machine-checked
  DOM like Rung C's. Raw netCDF never committed; reader = optional demo dep (NumPy-only-at-import holds).
  **Spike-first: Earthdata login/data acquisition = the one external unknown.**
  **SPIKE DONE 2026-07-06** (`OSCAR_L4_OC_FINAL_V2.0`, one granule, 33 MB, fetched to local temp workspace
  only — never committed): **auth** = an EDL bearer token as `Authorization: Bearer <token>` against
  `archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-protected/...` works directly, no `.netrc`/URS
  redirect dance; **format** = netCDF4/HDF5 via `h5netcdf`+`h5py` (`cftime` also needed — time uses a
  **`julian`** calendar `pandas` can't decode, relevant at O4); `lat`(-89.75…89.75,719)/`lon`(0…359.75,1440,
  **0–360 not ±180**) both already ascending (no N→S flip), but the lon rewrap to ±180 un-sorts the array
  (needs re-sort/roll); **dim order is `(time, lon, lat)`** — opposite of `FlowField`'s `(n_lat,n_lon)`,
  transpose required; `u`/`v`=total(geo+Ekman, the fuller signal) vs `ug`/`vg`=geostrophic-only(named future
  knob); units already m/s, no conversion. **This makes O1's mask retarget concrete, not hypothetical**
  (see O1 above): 44% NaN=land/missing poisons `_bilinear`'s plain `np.interp` at every coastline →
  fill-NaN→0-before-resample + nearest-neighbor mask-resample + explicit pole-masking (OSCAR doesn't reach
  the true poles; a global field must mask them, not edge-clamp-extrapolate). **Forward flag for O3, not
  yet acted on:** interchange round-trip grid is 2° but OSCAR is 0.25° — confirm the *render* path uses
  native/finer res (2° spec should stay proof-only) before calling it "beautiful"; native-res global texture
  costs much more than the eddy showcase's 758 KB, a size tradeoff needing a conscious call.
- **O3 — beauty pass** (renderer-only, contract untouched): (a) land/ocean two-tone base texture from the
  O1 mask; (b) **accumulate-and-fade trails** (third ping-pong render-target pair; CPU fallback keeps
  fade-only); (c) speed colormap + speed-weighted seeding (boundary currents dominate). Unlocks the §9.5
  control-surface seam (trail length + density = first knobs).
- **O4 — frames time axis**: the R1 deferral, built after O2's real dims are seen; `(nt,ny,nx)` stacked
  npz + crossfaded `uVelA`/`uVelB` textures; OSCAR monthly climatology → Somali Current monsoon reversal.
- **O5 — QG producer** `flow_field_from_qg` (independent): second EMERGENT producer, box coverage no mask;
  a third producer re-trips §9.4 rule-of-three for the two-consumer geometry helpers.

Out of scope, restated: no ocean engine, no forcing seam (S2 designs against the SEEN ClimaOcean API), no
real-time what-if on real data. [[planet-spinout-roadmap]] [[planet-viz-animation-rungs]]
[[planet-rung3-qg-built]]
