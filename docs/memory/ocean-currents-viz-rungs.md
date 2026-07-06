---
name: ocean-currents-viz-rungs
description: Ocean-currents showcase rungs O1–O5 scoped 2026-07-06 (plan §9.6, NOT built) — real OSCAR/ECCO currents through the R1 seam onto the Rung-C globe; §11.4 fork retargeted viz-half-only
metadata:
  type: project
---

**Scoped 2026-07-06 (NOT built) — plan §9.6 is the record.** The user's "visualize beautiful ocean
currents" ask exercised §11.4's *named alternative* under the living-staircase rule: real-data ocean
ingest-for-**display** moves INTO planet-sim (the renderer + seam live here), while the engine/forcing
half of the spin-out boundary is unmoved (no ocean physics here; ClimaOcean stays S2–S5; spin-out **S1
narrows** to ECCO-as-validation-anchor-for-S3). The §11.2 "never ships an ocean visual" line is amended
(visual yes, engine never); status notes threaded at §11.2/§11.3-S1/§11.4 + backlog §12.3.

**The ladder (order O1→O2→O3→O4; O5 independent):**
- **O1 — mask increment**: `FlowField.mask` per-cell validity (land), `None`=all-valid=bit-for-bit;
  Coverage stays the bounding box; seed only masked cells (band-zeros honesty style); mask rides the
  velocity texture's free 4th channel `(u,v,θ,mask)`; round-trip `==` extended.
- **O2 — real-ocean producer** (the deliverable + S1 de-risk): one OSCAR 0.25° snapshot (PO.DAAC) →
  `flow_field_from_ocean` → full-globe `is_global=True` (first true-global consumer) → serialize (round-trip
  on REAL data) → unchanged Rung-C globe → `docs/figures/planet-ocean-currents.html`. Honesty flips class:
  a new **provenance clause** ("real reanalysis-class currents, NOT this model's output"), machine-checked
  DOM like Rung C's. Raw netCDF never committed; reader = optional demo dep (NumPy-only-at-import holds).
  **Spike-first: Earthdata login/data acquisition = the one external unknown.**
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
