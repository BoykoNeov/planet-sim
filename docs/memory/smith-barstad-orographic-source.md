---
name: smith-barstad-orographic-source
description: Rung 5A cited source pin — the Smith & Barstad (2004) linear-theory orographic-precip constants + transfer function + the closed-form triangle-ridge exact solution that planet/orographic.py pins and test_orographic.py anchors on
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
