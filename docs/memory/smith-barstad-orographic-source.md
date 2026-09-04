---
name: smith-barstad-orographic-source
description: Rung 5A cited source pin — the Smith & Barstad (2004) linear-theory orographic-precip constants + transfer function + the closed-form triangle-ridge exact solution that planet/orographic.py pins and test_orographic.py anchors on; + the 5A.3 drying-ratio band (Roe 2005 / Smith & Evans 2007) and the 5A.4 freezing-level + polar-inversion pins (Harris Bowman & Shin 2000; Serreze et al. 1992)
metadata: 
  node_type: memory
  type: reference
  originSessionId: fdb095c9-cdad-43a3-ac05-2e4a17bc3da3
---

The `[[...-source]]` pin for **`planet/orographic.py`** (Rung 5A, built 2026-07-09). Numbers cited at
build, cross-checked against the **PISM/QGIS reference implementation**
(`pism/LinearTheoryOrographicPrecipitation`, `linear_orog_precip.py`), NOT from memory.

**Paper:** R. B. Smith & I. Barstad (2004), "A Linear Theory of Orographic Precipitation",
*J. Atmos. Sci.* **61**, 1377–1391.

**The transfer function** (applied to the FFT of the terrain `ĥ`):
`P̂(k,l) = C_w·iσ·ĥ / [(1 − i·m·H_w)(1 + iσ·τ_c)(1 + iσ·τ_f)]`, with
`σ = Uk + Vl` (intrinsic frequency) and the vertical wavenumber
`m² = (N_m² − σ²)(k² + l²)/(σ² − f²)`.

**Three numerics that are each an anchor's failure mode:**
- **Branch of `m` (the rain-shadow sign):** propagating modes (`m² ≥ 0`) take the root signed by
  `sgn(σ)` so the wave tilts *upwind* → drying in the *lee*. Wrong branch flips the shadow to the
  windward side. Evanescent modes (`m² < 0`) take the decaying/bounded branch (numpy complex sqrt).
- **`σ = 0` locus** (modes ⊥ wind, incl. DC): numerator `iσĥ = 0` there → set `P̂ = 0` (else 0/∞).
  Also regularize `σ² = f²` (the wide-mountain resonance) to ±eps.
- **FFT periodicity:** zero-pad by a full domain width so lee drying can't wrap into the windward edge.

**Pinned constants (the loose-magnitude knobs):** `τ_c = τ_f = 1000 s`, `N_m = 0.005 s⁻¹`,
`H_w = 2500 m`, `ρ_Sref = 7.4e-3 kg/m³`, moist lapse `Θ_m = −6.5 K/km`, env lapse `γ = −5.8 K/km`,
so **`C_w = ρ_Sref·Θ_m/γ ≈ 8.29e-3 kg/m³`** (both lapses negative → C_w > 0). Reference wind
`U = 15 m/s`; `direction` is meteorological (270° = westerly = wind FROM west, u = +speed).
Output is a condensation mass flux kg/(m²·s); with ρ_water = 1000 (1 kg/m² ≡ 1 mm), **×3600 → mm/hr**.
Truncate `P < 0 → 0` (the linear model's unphysical "anti-rain" at a downwind kink).

**THE TIGHT ANCHOR — the closed-form triangle-ridge solution** (reduced limit `H_w = τ_c = 0`, `f = 0`,
wind → +x): for a triangle ridge `h = A·max(0, 1 − |x|/d)`, with `C = C_w·U·A/d` and `Uτ = U·τ_f`,
- windward `−d ≤ x < 0`: `P = C·(1 − e^{−(x+d)/Uτ})`,
- **lee `0 ≤ x ≤ x_c`: `P = C·(e^{−x/Uτ}(2 − e^{−d/Uτ}) − 1)`** — the rain shadow, decaying to 0 at
  `x_c = Uτ·ln(2 − e^{−d/Uτ})`; zero beyond.
The FFT model converges to this ~O(dx²) **once truncated** (untruncated it has a −C-scale "anti-rain"
spike at the downwind kink `x = d` that dominates max-error until clipped). At dx = 1 km on A = 500 m,
d = 50 km, U = 15 m/s: peak ≈ 4.32 mm/hr = `C·(1 − e^{−d/Uτ})` (NOT the asymptotic `C ≈ 4.48` — finite
fallout advection reduces the crest), rel error < 0.1 %. **Scope of this anchor (do not overclaim):**
in the limit `H_w = 0` the factor `(1 − i·m·H_w) → 1`, so `m` and its `sgn(σ)` branch DROP OUT — the
triangle anchor pins only the *reduced* response `C_w·iσ/(1 + iσ·τ_f)` (the `C_w` scaling, upslope `iσ`,
fallout `τ_f`), NOT the branch or `H_w`. The branch's sole validator is the rain-shadow *direction*
check (a branch flip reddens only that). See [[planet-rung5a-orographic]].

**Rung 5A.3 addendum — the drying-ratio calibration** (`planet/orographic_depletion.py`, built
2026-07-10). The lee-depletion moisture budget is calibrated through the **orographic drying ratio**
`DR ≡ ∫P_oro dx/(U·W₀)` (the fraction of the passing moisture flux that rains out crossing a range) — a
standard *observed* quantity. Cited band **DR ≈ 0.3–0.5** for mid-latitude ranges — the Sierra Nevada /
Cascades constrained to **0.48 ± 0.02**, the Southern Alps of NZ and the southern Andes similar with more
scatter. **Verified sources (WebSearch-checked 2026-07-10 — my first-pass citations were WRONG from
memory, advisor-caught):** **G. H. Roe (2005)**, "Orographic Precipitation", *Annu. Rev. Earth Planet.
Sci.* **33**, 645–671 (the review; drying-ratio concept + observed bands, incl. the upslope form
`DR = 1 − e^{−h/H}`); **R. B. Smith & J. P. Evans (2007)**, "Orographic Precipitation and Water Vapor
Fractionation over the Southern Andes", *J. Hydrometeorol.* **8**(1). (**Do NOT reuse the retracted
first-pass cite** "Smith 2003 / Kirshbaum & Smith 2008 *Tellus A* 60, 543–561" — fabricated; the real
Kirshbaum & Smith 2008 is *Q. J. R. Meteorol. Soc.* **134**, 1183–1199 on moist-stability effects, not a
drying-ratio pin.) The one new knob is the incoming column precipitable water **`PWV_IN_MM ≈ 30 mm`**
(≈ 3 cm, a typical mid-latitude column), tuned so the Cascades demo `DR ≈ 0.47` lands in-band. Also
derived (not pinned): the evaporative refill length `L = U·W₀/P_base` (~16 000 km at Earthlike numbers) —
the scale that makes the no-refill (`L→∞`) budget honest on a ~450 km patch. NOT from memory — cite at
any reuse.


**Rung 5A.4 addendum — the freezing-level benchmark + the polar-inversion caveat**
(`planet/elevation_temperature.py`, built 2026-09-04). Two numbers, both **WebSearch-verified at build,
not from memory** (the 5A.3 lesson above).

1. **`OBSERVED_TROPICAL_FREEZING_LEVEL_M = (4500, 5000)`** — the observed altitude of the **0 °C
   isotherm** (freezing level) in the deep tropics, the single observational check that decides which
   lapse rate this module defaults to. **Source:** **A. R. Harris, K. P. Bowman & D.-B. Shin (2000)**,
   "Comparison of Freezing-Level Altitudes from the NCEP Reanalysis with TRMM Precipitation Radar
   Brightband Data", *J. Climate* **13**(23), 4137–4148 — a 20-year global climatology of the 0 °C
   isotherm from 6-hourly NCEP reanalysis, cross-checked against TRMM radar brightband heights:
   **tropical freezing levels are the highest on the planet at ≈ 5000 m**, with the lowest intramonth
   and interannual variability; they fall and grow more variable through the subtropics and
   midlatitudes. Corroborated by **R. S. Bradley et al. (2009)**, "Recent changes in freezing level
   heights in the Tropics…", *Geophys. Res. Lett.* **36**, L17701 (doi:10.1029/2009GL037712) — same
   ~5000 m tropical level, rising ~45 m over the last ~30 years.
   **Read it honestly:** the cited central value is ≈ **5.0 km**, so the module's band is a generous
   read of "the tropics" (its warm edge is the pinned number). The 6.5 K/km constant puts the level at
   **4.38 km — just BELOW the band**, not inside it; the moist adiabat puts it at **7.09 km, ~45 %
   above**. The verdict that keeps the constant as the default therefore rests on the **ordering**
   (the constant is close, the adiabat is far), which is robust, *not* on the constant landing inside
   the band. Say "just below" — do not upgrade it to "in band".
2. **The polar surface inversion** — why the module declines to bank its own predicted latitude
   contrast (`Γ_m` steepening toward the dry adiabat in cold air ⟹ mountains cooling *more* at high
   latitude). **Source:** **M. C. Serreze, J. D. Kahl & R. C. Schnell (1992)**, "Low-Level Temperature
   Inversions of the Eurasian Arctic and Comparisons with Soviet Drifting Station Data", *J. Climate*
   **5**, 615–629: from rawinsonde climatology, low-level inversions occur in **over 95 % of winter
   soundings** east of Novaya Zemlya and are typically **surface-based**, strongest in winter under the
   large surface radiative deficit. So the real high-latitude lower troposphere is stably stratified —
   often *warming* with height — and a saturated-adiabat lapse rate is the wrong idealisation there.
   NOT from memory — cite at any reuse.