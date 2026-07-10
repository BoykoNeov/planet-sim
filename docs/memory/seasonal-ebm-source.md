---
name: seasonal-ebm-source
description: Rung 5B.1 cited sources seasonal.py pins — North&Coakley1979 seasonal EBM, NMS83 seasons+continents, Hartmann heat capacities
metadata:
  type: reference
---

The `[[…-source]]` pins for **rung 5B.1** (`planet/seasonal.py`, the seasonal zonal EBM +
continentality, BUILT 2026-07-10) — see [[planet-rung5b-seasonal]]. Extends [[ebm-radiation-source]]
and [[obliquity-insolation-source]].

**Model structure (cited):**
- **North & Coakley 1979** (*J. Atmos. Sci.* 36, 1189–1204, "Differences between seasonal and mean
  annual energy-balance model calculations…") — the seasonal EBM, land/ocean distinguished **strictly
  by thermal inertia** (heat capacity). The 1-D-in-latitude, two-surface reduction 5B.1 builds.
- **North, Mengel & Short 1983** (*JGR* 88, C11, 6576, "Simple energy balance model resolving the
  seasons and the continents") — the full **2-D (lat × lon)** seasonal EBM; the title is near-verbatim
  the task ("seasonality for continentality"). This is what **"2-D EBM"** means (two *spatial* dims) and
  the target of the deferred **5B.2**. Confirmed by WebSearch at build (I could not close-read the 1983
  tabulated `C` values from search — so `C` is BUILT from textbook pieces, not memorized).

**Heat capacities — BUILT from textbook constants, not memorized (the non-fabrication choice):**
- Land `C_L = C_atmos + ρ_w c_w · soil_depth`, ocean `C_O = C_atmos + ρ_w c_w · mixed_depth`.
- `C_atmos = c_p·p_s/g` (the atmospheric column heat capacity, Hartmann *Global Physical Climatology*):
  `c_p=1004`, `p_s=1.013e5`, `g=9.81` → **≈ 1.04e7 J/m²/K** (a pure textbook constant).
- Calibrated (loose) depths: **soil ≈ 2 m** (the seasonal soil thermal-penetration depth), **mixed layer
  ≈ 50 m** (Hartmann; the seasonal ocean mixed layer). → `C_O/C_L ≈ 12×`.
- **Why loose:** the model damps only through `B` (no evaporative damping over land), so amplitude and lag
  are tied through the single `C` (slab law). Land amplitude therefore runs **high-end** (~30 K at 45°,
  Siberian) — direction banked, magnitude calibrated. Deeper soil lowers amplitude *and* grows lag toward
  the ocean's (shrinks the contrast both ways) — 2 m is the balance.

**The tight anchor is the 0-D slab, not a benchmark number:** forced `C dT/dt = F₀ + F₁cos ωt − A − BT`
→ amplitude `F₁/√(B²+ω²C²)`, lag `arctan(ωC/B)` (exact; textbook forced-oscillator). The spectral
`n=0` harmonic == the annual-mean EBM ([[ebm-radiation-source]] constants), reproduced to 1e-11.

Seasonal forcing = the **pinned daily-insolation kernel** ([[obliquity-insolation-source]]:
`H₀ sinφ sinδ + cosφ cosδ sinH₀`, `S₀/π` factor for absolute W/m², global-annual mean `S₀/4`),
factored out of `annual_mean_insolation` into `obliquity.daily_mean_insolation` so the seasonal model
keeps the time axis the obliquity knob averages away — one kernel, two consumers.
