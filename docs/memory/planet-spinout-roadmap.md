---
name: planet-spinout-roadmap
description: "The editable-ocean GPU spin-out project — DECIDED + roadmapped 2026-06-12 (NOT built): a separate Julia/ClimaOcean+Makie repo born FROM planet-sim across a contract seam; the TWO-seams split (viz/output vs forcing/input), planet-sim rungs R1–R3+Rung-C, spinoff rungs S1–S5, ECCO=the ocean anchor at S1, + the living-staircase retarget rule; full record plan §11"
metadata:
  node_type: memory
  type: project
  originSessionId: 499fdfc2-b79a-4a22-a0c4-6b2bb365d651
---

A **separate future project** — an **editable land/ocean world with 3-D GPU viz, driven by ClimaOcean.jl**
— **DECIDED + fully roadmapped 2026-06-12, NOT built** (this is forward planning, no code). Full record:
`docs/plans/planet-earth-system.md` **§11**. Continues [[planet-plan]] / [[planet-viz-animation-rungs]].

**The shape:** a separate **Julia/GPU repo born FROM planet-sim across a documented contract seam** —
*consume-don't-vendor* (the same relationship planet-sim has with `engines/`). Rejected: stay-within (wrong
language/audience) + cold-start (throws away the emergent atmosphere). **The collapsed-seam insight:**
ClimaOcean runs ocean+sea-ice with the **atmosphere PRESCRIBED as forcing** (default JRA55) — and
planet-sim **IS** an atmosphere (EBM + jet + rung-2 moisture/`P−E`) → one-way coupling is a **substitution**
(planet-sim's emergent atmosphere replaces JRA55), needing **ClimaOcean ONLY** — NOT ClimaAtmos, NOT
ClimaCoupler. ClimaCoupler is owed **only** for the deferred two-way loop (S5).
**VERIFIED against CliMA docs 2026-06-12** (all 4 architecture claims hold): ClimaOcean = *"ocean-only and
coupled ocean+sea-ice driven by prescribed atmospheres"*, standalone *"built on Oceananigans and ClimaSeaIce"*
(no ClimaAtmos/ClimaCoupler dep); `JRA55PrescribedAtmosphere` is the canonical **built-in (pluggable, not a
blessed sole default → the substitution slot is the interface working as intended)**; ClimaCoupler = the
atmosphere+land+ocean+ice orchestration layer = the S5 two-way loop; ECCO = `ECCO2Daily`/`ECCORestoring`
dataset. **Live wrinkle (living-staircase already covers it):** generic coupling machinery is migrating to
`NumericalEarth.jl`, so exact package/type names S2–S3 bind against may shift — re-confirm at S1. Julia-side GPU renderer =
**Makie.jl** (GL/WGL), free in-language; planet-sim side reuses the §9.5 globe stack. Python↔Julia =
atmosphere↔ocean = project boundary, all the same line.

**THE KEY CORRECTION (advisor-caught — do not re-fuse): "the seam" is TWO seams, opposite directions.**
- **Viz/output seam** (planet-sim, NOW): grid + `(u,v)` **velocity** + scalar + frames + coverage-extent +
  provenance/honesty-label → *display* → a renderer (§9.5 globe / Makie / Rung C). This is R1.
- **Forcing/input seam** (spinoff, LATER, at the boundary): wind **STRESS** + heat-flux components + `P−E`
  (+ optional SST-restoring) → *drive* → ClimaOcean's input API. **Designed only once that API is SEEN, not
  guessed.** This is S2/S4. **Neither ECCO-ingest nor ClimaOcean-vs-ECCO validation needs it** (they run on
  JRA55) — it goes live only at S4.
- Shared backbone = **pin-a-schema-not-a-format** (§9.3): one schema, two encodings (JSON+npz browser /
  NetCDF Julia), `import(export(s))==s` round-trip-identity = the real test.

**planet-sim rungs (finish in THIS repo; stays ATMOSPHERE-ONLY, never ships an ocean visual):**
- **R1 = materialize+serialize the viz contract — THIS is what the spinoff binds on, NOT Rung C.** Serialize
  the today-computed-then-viewed `vector_overlay` jet into the §9.3 schema (vector-field layer, both
  encodings, round-trip test extended). Decisive move: **add a synthetic global-coverage 2nd producer** and
  read BOTH the real eddy band AND the synthetic field through the **already-built Rung B renderer**
  (`eddy_globe.py`) → proves **producer-agnosticism** = the exact property ClimaOcean later relies on.
- **R2 = §9.4 toolkit promotion** (rule-of-three now met: frame side-channel + flow-globe + serialization
  have a 3rd consumer). Natural co-rung.
- **R-parallel = Rung C** (three.js/WebGL particle showcase) **OFF the critical path** — proves
  *renderer*-agnosticism (a DIFFERENT axis from R1's *producer*-agnosticism); spinoff does NOT wait on it.
- **R3 = bank planet-sim** (atmosphere-only capstone; seam documented+tested). Start the new repo only AFTER
  R1 is banked.

**spinoff rungs (new repo; numbered by logic, EXECUTION order = S1 → S3 → S2 → S4 → decide S5):**
- **S1 = ECCO ingest + viz, pure Python** (the new repo's first rung, NO Julia engine). Dual role:
  cheap *Perpetual-Ocean* real-data globe (stands up the pipeline) **+** the validation anchor for S3.
- **S2 = design the forcing/input seam** (NetCDF; at the boundary, after S3 reveals the API).
- **S3 = ocean-1: ClimaOcean+JRA55 validated vs ECCO** (no planet-sim seam — runs on JRA55; this is where
  ClimaOcean's input API is *learned*, which S2 then consumes).
- **S4 = swap JRA55 → planet-sim's emergent atmosphere** (forcing seam goes LIVE = the payoff: an edited
  world's atmosphere drives a GPU ocean, rendered back through the viz seam). Honest-by-disclosure off Earth.
- **S5 = DEFERRED gate: ocean-2 two-way loop** (SST back into planet-sim closes rung-2's *faked* `E` — the
  prescribed `P−E` that never had an honest evaporation; **ClimaCoupler** finally earns its place).
  Re-plan from scratch here with real ClimaOcean experience.

**The settled fork:** ECCO lives in the **new repo as S1** (planet-sim stays atmosphere-only); the
alternative (ECCO as planet-sim's last rung, de-risk earlier) stays named. **ECCO = a DATASET (ocean state
estimate / reanalysis), the ocean's `[[ebm-radiation-source]]`-class ground truth** — the Earth anchor.
Honesty ceiling carried from the staircase: a custom world's ocean **cannot be validated** →
*honest-by-construction on Earth (ECCO checks you), honest-by-disclosure off it* (the §9.5 Rung C carve-out,
now load-bearing for a whole project); and **not real-time** for custom worlds (a cheap emulator layer is
named+deferred, not in this roadmap).

**THE FRAMING RULE (user caveat 2026-06-12): every rung is PROVISIONAL until the previous lands** — not a
frozen waterfall but the **same living-contract discipline** as the rest of the plan ([[engines-living-contracts]]
ADR 0005, spike-first, "a trade not a win", anchor-then-trust-delta). Each rung carries a
`Retarget-when-done` note; plan the next rung concretely, hold the one after loosely, revalidate the chain
at each landing. See also [[planet-interactive-map-design]] (the §9.3 schema), [[always-push-commits]].
