---
name: planet-rung5b-seasonal
description: Rung 5B.1 seasonal zonal EBM BUILT 2026-07-10 — heat capacity woken, continentality from the land/ocean C contrast; 5B.2 (2-D map) named
metadata:
  type: project
---

**Rung 5B — the seasonal 2-D EBM (seasons × continents).** The "needs seasonality" half of the geography
seam (plan §12.5 / §9.3 row 611). *"2-D EBM"* is a **term of art = two spatial dims (lat × lon)** — North,
Mengel & Short 1983, *"…resolving the seasons and the continents"* (near-verbatim the task). Advisor caught
my lean toward the small build; WebSearch confirmed the term. **Climbed staged (user's call, the 5A
pattern):** 5B.1 zonal core now, 5B.2 the full map later. → [[seasonal-ebm-source]], [[planet-rung5a-orographic]].

**The core insight (why this rung has to exist):** every EBM before it solved for an **equilibrium**, where
`C` **cancels** — so a land and an ocean column at the same latitude reached the **identical** temperature
and **continentality was exactly zero**. Seasonality is the prerequisite to *any* continentality.

**Rung 5B.1 — seasonal zonal EBM · BUILT 2026-07-10** (`planet/seasonal.py`, `test_seasonal.py`,
`demo_seasonal.py`, `plots.seasonal_figure`; also refactored `obliquity.daily_mean_insolation` out of
`annual_mean_insolation`). North & Coakley 1979. A **sibling** full-sphere EBM (`x ∈ [−1, 1]` — the
seasonal cycle is hemispherically *anti*-symmetric, NH summer = SH winter, so the rung-0 hemisphere grid
can't carry it; spine EBMs untouched), same diffusive transport + linear radiation **marched under
axial-tilt insolation `S(x,t)`** to a converged **annual limit cycle** → `C` finally **load-bearing**.

- **Continentality mechanism:** two heat-capacity tiles per latitude (land ≈ atmospheric column `c_p p_s/g`
  + ~2 m soil; ocean + ~50 m mixed layer, ~12× bigger), sharing meridional transport on the zonal mean via
  an **effective transport heat capacity `C_a = (f_L/C_L + f_O/C_O)⁻¹`** — chosen so the per-tile
  energy-flux redistribution `ΔT_i = (C_a/C_i)ΔT̄` **exactly reproduces** the engine's `T̄` step and
  conserves column energy (uniform `f_L` → constant `C_a`). Advisor verified the algebra.
- **Two solvers:** `march()` (engine-reuse method, Strang-split to limit cycle — carries ice-albedo later)
  and `spectral()` (the tight reference — one complex banded solve per temporal harmonic; **`n=0` harmonic
  *is* the annual-mean EBM**, so reduction-to-parent falls out as the DC term).
- **PAYOFF:** at 45° land tile swings ~30 K (range ~59 K) nearly in step with the sun; ocean tile ~3 K,
  lags ~2.7 months → **~10× continentality contrast from `C` alone**. Beautiful antisymmetric pole-to-pole
  Hovmöller in the figure.

**Anchors banked (tight):** 0-D slab `amp=F₁/√(B²+ω²C²)`, `lag=arctan(ωC/B)` both solvers, transport off
(the mechanism); reduction spectral `n=0` == `SphereEBM.steady_linear` to 1e-11 **against the true `⟨S⟩`
not P₂-truncated `insolation()`** (advisor-caught, else ~1e-2) **and `⟨T_L⟩=⟨T_O⟩`** (continentality is
*entirely* in the amplitude, zero in the mean); **marcher → spectral at 1st order in `dt`** = the
anti-damping cross-check (advisor-caught: every other tight anchor is blind to time-accuracy in the
backward-Euler transport substep — slab turns transport off, reduction/conservation time-average — so a
quiet amplitude-damping bug would ship un-caught); hemispheric antisymmetry `T(x,t)=T(−x,t+½yr)`;
annual+global energy balance. **Loose (calibrated):** amplitudes/lag ride the calibrated `C`; land runs
high-end because damping is `B`-only (amp & lag tied through the single `C`).

**Named scope:** same albedo both tiles (continentality = *pure* heat capacity); fixed albedo (exact
reduction; seasonal ice-albedo + small-ice-cap instability are the marcher's future); uniform land fraction
(exact conservation); `D` = atmospheric transport shared over both tiles (well-mixed-atmosphere closure).

**Rung 5B.2 — the true `T(φ,λ,t)` land–sea map · NAMED, not built.** Longitude axis + 2-D transport
(ADI/operator-split **reuses** the 1-D tridiagonal solver per sweep — bounded engine step, not a rewrite;
wrinkle = polar meridian convergence) + a real land mask → continental interiors as a *map* (NMS83). The
5B.1 zonal two-tile reduction becomes resolved 2-D diffusion of one field with spatially-varying `C(φ,λ)`.
