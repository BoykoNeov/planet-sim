---
name: planet-rung5b-seasonal
description: Rung 5B BUILT 2026-07-10 — 5B.1 seasonal zonal EBM (continentality from land/ocean C) + 5B.1+ seasonal ice-albedo (migrating ice edge, ice asymmetry, bistability; SICI → resolved on rung 0, see planet-bifurcation-sici) + 5B.2 the full 2-D lat×lon map (NMS83) + 5B.3 seasonal ice ON the map + albedo maps (2026-09-02: the snow map, the annual mean SEES the mask by rectification); rung 5B COMPLETE
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

**Named scope:** same *ice-free* albedo both tiles (continentality = *pure* heat capacity); fixed albedo
*for the spectral solve + tight reductions* (seasonal ice-albedo now BUILT as a marcher-only extension —
5B.1+ below; small-ice-cap instability still deferred); uniform land fraction (exact conservation); `D` =
atmospheric transport shared over both tiles (well-mixed-atmosphere closure).

**Rung 5B.1+ — the seasonal ice-albedo feedback · BUILT 2026-07-10** (`seasonal.py` `ice_coalbedo` +
`march(coalbedo_fn=…)`, `test_seasonal_ice.py`, `demo_seasonal_ice.py`, `plots.seasonal_ice_figure`). The
fixed-albedo edge lifted **on the marcher only**: the rung-0 step-function ice-albedo (`albedo.planetary_albedo`,
REUSED — Phase 1's Snowball nonlinearity) as `march`'s `coalbedo_fn`, re-frozen on **each tile's own T** every
half-step (`absorbed_i = S(x,t)·(1−α(x,T_i))`) → land/ocean freeze *independently*. Marcher-only (state-dependent
α ⟹ spectral solve can't carry it); opt-in (`coalbedo_fn=None` ⟹ fixed path **bit-identical**). **PAYOFF:** a
**migrating seasonal ice edge** (forms in winter, melts in summer — impossible in any equilibrium EBM below);
**continentality → an ice asymmetry** (Earth insol: land seasonal ice reaches ~33° vs ocean ~66°; land at 45°
frozen ~60% of year, ocean 0%); and `C` decides whether ice MELTS — **land ice purely seasonal** (tiny-`C`
climbs above freezing every summer, no year-round land ice) while **ocean ice ~perennial once formed** (huge-`C`
barely warms). Restores **Phase-1 bistability inside the seasonal cycle**: warm seed ⟨T⟩≈+11 °C (finite ice) vs
cold seed ⟨T⟩≈−40 °C (snowball, ice to equator), one sun. **Tight anchors:** warm-limit reduction (unreachable
`Tf` never trips → bit-identical to fixed ice-free march, <1e-12) + frozen-limit (high `Tf` → const co-albedo →
bit-identical to fixed-`a_ice` march); **ε=0 in-model self-consistency** — seasons off ⟹ steady fixed point of
`L_T·T̄ + S(1−α(T)) − A − B·T = 0`, residual **O(dt)** (reads model vs ITSELF — the advisor's fix for two ghosts:
no external insolation-truncation-match, no ice-branch-select). **Conservation:** global-annual net TOA w/
*realized* co-albedo ≈0 (transport untouched). **Loose:** ice edges/fractions ride calibrated `C`+`Tf`/`a_ice`;
direction banked. **Named-DEFERRED:** small-ice-cap instability (North&Coakley critical-cap-size / ice-edge jump)
— the sweep this feedback *enables* but this build didn't attempt (plain reading = mechanism, not bifurcation).
Advisor: design + scope-defer-SICI + all three anchor traps. → [[ebm-radiation-source]] (the reused α), [[planet-phase1]].

**Rung 5B.2 — the true `T(φ,λ,t)` land–sea map · BUILT 2026-07-10** (`planet/seasonal_map.py`,
`test_seasonal_map.py`, `demo_seasonal_map.py`, `plots.seasonal_map_figure`). North, Mengel & Short 1983.
The longitude axis added: a **single** field `T(φ,λ,t)` on `x=sinφ∈[−1,1] × λ∈[0,2π)`, a **mask-set** heat
capacity `C(φ,λ)` (5B.1's `land_/ocean_heat_capacity` reused; the two zonal *tiles* become one resolved
field), marched to a limit cycle by an **ADI operator-split**: `½rad → meridional sweep → zonal sweep →
½rad` (ordered to match 5B.1 so the reductions are exact).

- **Method = reuses the engine's tridiagonal ASSEMBLY, hand-rolls the batched solves** (the honest reading
  of "ADI reuses the 1-D solver"). The meridional sweep **reconstructs** the engine's exact harmonic-mean-
  face `Lₓ` (as `sphere_ebm`/`seasonal` do) + a per-cell `C/Δt` diagonal → batched Thomas; the zonal sweep
  is the **periodic** counterpart (wrap couples last↔first lon → **cyclic** tridiagonal, `solve_banded`'s
  band broken) → batched **Sherman–Morrison** of the same Thomas kernel. **Why not call the engine `step()`:**
  varying `C` (mask) can't ride the engine's *uniform-*`C` step, and periodicity breaks its banded solve —
  so reconstruct + hand-roll, the sibling pattern (ADR 0005), **not** spine surgery. The **meridian-
  convergence** wrinkle (zonal coeff `D/(1−x²)`→∞ at poles) is tamed by the cell-centered grid (no cell
  *at* a pole) + unconditionally-stable backward-Euler (relaxes each polar ring to its zonal mean).
- **`SeasonalMapEBM` wraps a `SeasonalEBM`** as its source-of-truth for grid/operator/insolation/co-albedo/
  `dt` → the reduction to 5B.1 is *bit-identical*, drift-proof. Insolation & albedo are lon-independent
  (named scope) — **every drop of zonal structure comes from the mask.**
- **PAYOFF — continentality is a MAP.** Coarse idealized-Earth mask (land fraction 0.29): at 45°N the broad-
  continent **interior swings ~40 K**, its **coast ~21 K** (adjacent ocean moderates by diffusion), the
  **open ocean ~7 K** — continentality now varies *within* a latitude. Beautiful map figure: range map
  (interiors blaze), annual-mean map (visibly zonally FLAT), lon cross-section, three-point cycle.
- **The NMS headline (banked tight):** average over the year and `C` cancels, so `⟨T⟩` solves the annual-mean
  EBM with **lon-independent** forcing → **zonally flat for *any* mask** and == the 1-D parent
  `SphereEBM.steady_linear` (the `(D∇²−B)` operator has trivial kernel for `B>0` → unique = flat parent;
  advisor-proved). This is 5B.1's `⟨T_L⟩=⟨T_O⟩` generalized: **the land/sea contrast is entirely in the
  seasonal amplitude, the annual mean is blind to the mask.** Machine-tight in the all-land case, ~0.2 K
  split residual with a real mask.

**Anchors banked:** per-cell 0-D slab (`D=0`, machine); **zonal invariance** (zonal-uniform mask → λ-flat to
1e-11) + **bit-for-bit reduction to 5B.1** (all-land≡`f=1`, all-ocean≡`f=0`, 1e-9 → inherits 5B.1's
spectral-validated anti-damping); **cyclic solver == circulant eigenmodes** `cos mλ` (machine — the one new
numerical object) + `_thomas_columns` vs `solve_banded` at **varying** diagonal (machine, advisor-flagged
gap closed); annual-mean reduction (above); hemispheric antisymmetry (symmetric mask); global+annual energy
`∫C T dA` conserved (net TOA ~6e-7). **Loose:** interior/coast/ocean magnitudes ride the 5B.1 `C`; direction
banked. **Named scope (from 5B.1):** fixed ice-free albedo same on land/sea (continentality = *pure* `C`);
**diffusive** continentality only (interior extremes + coastal moderation, no wind-driven downwind tilt);
prescribed geography. **Deferred:** a 2-D frequency-domain solver (would make the annual-mean reduction
machine-tight like 5B.1's spectral; the marcher's structural anchors already pin every piece). **Rung 5B
COMPLETE.**

**Rung 5B.3 — seasonal ice ON the map + albedo maps · BUILT 2026-09-02** (`seasonal_map.py`:
`march(coalbedo_fn=…)`, `albedo=` now accepts a `[n_x, n_lon]` **map**, `T_init` a field; `ice_free_albedo_map`
(the cheap-tier land/sea contrast knob, §12.5 "land/ocean → an albedo difference", default offsets 0),
`masked_ice_coalbedo` (the step feedback on top of a map), `SeasonalMapClimate.ice_fraction/frozen/zonal_anomaly`;
`test_seasonal_ice_map.py` (8 fast + 4 slow); `demo_seasonal_ice_map.py` → `planet-seasonal-ice-map.png` + a
**month-by-month GIF** `planet-seasonal-ice-map.gif` (`plots.seasonal_ice_map_figure/_animation`) + a **month-slider
Plotly globe** `seasonal_globe.py` → `planet-seasonal-ice-globe.html` (one `go.Surface`, 12 `surfacecolor` frames,
frozen cells pinned to an ice-white top colour stop; reuses `planetmap._sphere_xyz/_polecapped`; `showspikes=False`);
catalogue `seasonal_ice_map` (+`interactive=`); notebook §8.8). The 5B.1+ tile feedback broadcast over longitude — `sea.ice_coalbedo` works
as-is on `x[:,None]`; the marcher's `forcing_at(s, T)` picks one of three forcings (fixed 1-D = the 5B.2
arithmetic verbatim / fixed map / live `coalbedo_fn`).

- **PAYOFF 1 — the seasonal-ice MAP** (45×90×360, 51 yr, 37 s): winter snow over the continental interiors
  (at 53°N the interior frozen 34 % of the year, the ocean at the same latitude 0 %; NH seasonal-ice reach
  land 45° vs ocean 58°), **land ice purely seasonal (0 % of land cells perennial), polar sea ice lingers
  (5 % of ocean cells perennial)**. The 2-D interior freezes LESS than the 5B.1 zonal tile (34 % vs 60 % at
  midlatitude) — the zonal sweep moderates it; 45° was too marginal (19 %), the story reads at 55°.
- **PAYOFF 2 — the annual mean now SEES the mask (5B.2's theorem broken by design).** Zonal anomaly of
  the annual mean: interior **−0.6 K**, open ocean **+0.2 K** (0.9 K east–west spread; **exactly 0** for
  the fixed-albedo march — re-verified in the same test). Mechanism = **rectification**: winter snow
  reflects sun the ocean keeps absorbing → a nonlinear (albedo-step) effect no linear model can show
  whatever its `C`. Also the linear route to a visible mask: a fixed albedo **map** (brighter land).
- **Anchors (tight):** warm-limit reduction **bit-identical** (`array_equal`) to the fixed 5B.2 march; a
  zonally-uniform albedo map ≡ per-latitude albedo bit-for-bit; `masked_ice_coalbedo(offset-free map)` ≡
  `sea.ice_coalbedo` bit-for-bit; all-land/all-ocean under ice → the 5B.1+ tile marcher (<1e-6, WITH ice
  present); zonal invariance under ice (1e-11); **fixed albedo map is linear** → `D=0` per-cell annual mean
  = its own radiative equilibrium (2e-6), and with transport the ZONAL MEAN of the annual-mean map = the
  1-D parent driven by the ZONAL-MEAN co-albedo (<0.3 K); hemispheric antisymmetry under ice; global-annual
  net TOA with the realized co-albedo ≈ 0. **Loose:** the ice fractions / the −0.6 K ride the calibrated
  `C`, `T_f`, `a_ice`; direction banked. **Named:** the land/sea ice-free albedo *offset* is a knob
  (default 0 — planetary contrast is muted by clouds vs the surface 0.06–0.10 ocean / 0.15–0.35 land,
  Hartmann GPC Table 4.2; pick ~0.05, not the surface value); no snow-depth / melt physics (the step
  albedo); idealized blocky mask. **SICI** for this rung: resolved on rung 0 ([[planet-bifurcation-sici]]);
  the *seasonal* SICI sweep (Huang & Bowman 1992) is **BUILT 2026-09-04 as rung 5B.4**
  ([[planet-rung5b4-seasonal-sici]]) — on THIS 5B.1+ marcher, and it finds the fold **gone** at Earth's
  50 m mixed layer (returning only as the layer deepens). Note for reuse: `march`'s `T_init` now also
  accepts an array / a per-tile `(T_L, T_O)` continuation seed.
- **Globe render traps (5B.3):** a cell-centred longitude grid leaves a **gash at the 0/360° seam** on a
  `go.Surface` — append a wrap column (`_wrap_lon`, the lon-analogue of `planetmap._polecapped`); the banked HTML
  loads Plotly from the **CDN**, so a sandboxed/offline browser screenshot renders blank — inline
  `plotly.offline.get_plotlyjs()` into a scratch copy to verify (Playwright + `/opt/pw-browsers/chromium-1194/…/chrome`,
  drive frames via `Plotly.animate(gd, ['Jul'], …)`). Verified Jan (snow over Eurasia, blue Arctic) / Jul (bare) 2026-09-02.
