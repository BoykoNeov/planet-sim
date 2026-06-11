---
name: ebm-radiation-source
description: Planet P1 — cited climlab/North EBM radiation+albedo+transport constants and the present-ice-line / snowball-threshold benchmark numbers that ebm.py + climate_reference.py pin
metadata: 
  node_type: memory
  type: reference
  originSessionId: 488c0c78-5c50-4c1e-823a-1c7eea857983
---

**Planet Phase 1 ([[planet-plan]], [[bigsim-program]]) pinned reference** — the energy-balance-model
radiation/albedo/transport **constants** and the **benchmark outputs** that `projects/planet/ebm.py`
and `climate_reference.py` use. Pinned at build (not carried from memory), the `[[…-source]]`
discipline. Two halves — the constants (the *inputs*) and the benchmark facts (the *outputs*):

**Constants = the climlab `EBM` class defaults** (North 1975 / Budyko 1969 lineage), fetched from the
climlab source (`climlab.model.ebm`, readthedocs `_modules/climlab/model/ebm.html`):
- `S0 = 1365.2 W/m²`, `s2 = −0.48` (insolation `(S0/4)(1 + s2·P2(sinφ))`, more sun at the equator)
- OLR `A + B·T` with `A = 210 W/m²`, `B = 2 W/m²/°C` (T in **°C** — the climlab convention)
- transport `D = 0.555 W/m²/K` in `D·d/dx[(1−x²)dT/dx]`, x = sinφ
- ice-line isotherm `Tf = −10 °C`; **albedo** = `a0 + a2·P2(x)` ice-free (`a0=0.30, a2=0.078`) /
  `ai = 0.62` where `T < Tf` (climlab `StepFunctionAlbedo`); `water_depth = 10 m` (heat capacity =
  timescale only). These ARE `climate_reference.REFERENCE.climlab_*` (single source of truth, tested).

**Benchmark outputs (independent of the constants — the non-circular targets):**
- **Present-day ice line ≈ 70–72°** = the *observed* perennial snow/ice edge (North 1975 / Sellers),
  NOT this model's output (model gives ~73° on the finite-cap branch → genuine independent target).
- **Snowball threshold = the large-ice-cap instability:** Budyko/North put it at a *few-percent*
  solar dimming; **Voigt & Marotzke 2010 (a modern coupled GCM) find 6–9 %** (via Isaac Held's GFDL
  diffusive-EBM blog post #40). This model's **~8.3 % dimming** sits inside that band.
- **Small- AND large-ice-cap instabilities** (North 1975, *J. Atmos. Sci.* 32:2033): a cap below a
  critical size recedes unstably to ice-free; beyond a critical size it runs away to the **Snowball**
  (ice line at the equator). The ice-age ice line corresponds to a ~1.3 % solar reduction (North 1975).
- The **bistability/hysteresis is wide** (the white planet re-melts only far brighter than it froze —
  the Snowball is a deep trap); model gives a ~580 W/m² loop.

**Live cross-check** = climlab `EBM` itself, consumed as a *reference tool* (the pycalphad pattern):
opt-in `[climate]` extra, `slow`/`importorskip` test, never copied. `climlab_present_day()` seeds a
capped IC so it lands on the same finite-cap branch (the system is bistable at present S₀).

Sources (pinned 2026-06-09): climlab EBM source
`climlab.readthedocs.io/en/latest/_modules/climlab/model/ebm.html`; North 1975 *J. Atmos. Sci.*
32(11):2033 (`journals.ametsoc.org`); Voigt & Marotzke 2010 via I. Held GFDL blog #40
`gfdl.noaa.gov/blog_held/40-playing-with-a-diffusive-energy-balance-model/`.
