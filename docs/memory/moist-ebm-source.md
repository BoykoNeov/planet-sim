---
name: moist-ebm-source
description: Planet rung-2 Phase A — cited sources moist.py pins (Clausius-Clapeyron q_sat constants; the energy-constrained ~2-3%/K precip rate; the fixed-RH diffusive moist EBM formulation)
metadata: 
  node_type: memory
  type: reference
  originSessionId: 0369b58c-6c5c-4ddf-a80e-44cc29351f78
---

**Pinned sources for Planet rung-2 Phase A** (`planet/moist.py`, the column moist-EBM diagnostic), the
`[[…-source]]` discipline. Built 2026-06-11; see [[planet-rung2-scoped]], extends
[[precip-parameterization-source]] / [[ebm-radiation-source]].

**The cited facts (three distinct source sets).**
- **Clausius–Clapeyron `q_sat`** — the integrated, constant-`L` form `e_sat(T)=e₀·exp[(L_v/R_v)(1/T₀−
  1/T)]` and `q_sat=ε·e_sat/(p−(1−ε)e_sat)`, with pinned `e₀=611.2 Pa @ T₀=273.15 K`, `L_v=2.5e6 J/kg`,
  `R_v=461.5 J/kg/K`, `ε=0.622`, `p_s=1.013e5 Pa`. Source: **Hartmann *Global Physical Climatology* §1.5
  / Bohren & Albrecht *Atmospheric Thermodynamics***. Validated tight as an EXACT function (textbook
  values ~3.8 g/kg @ 0 °C, ~14.7 @ 20 °C; ~7 %/K log-slope) — the Whittaker-partition precedent (exact
  testable function, not a fit).
- **The energy-constrained global precip rate** — global-mean `⟨P⟩` is limited by the atmospheric ENERGY
  budget `L⟨P⟩≈R_atm−SH` (atmospheric radiative cooling), giving **~2–3 %/K**, much slower than the C–C
  moisture-capacity **7 %/K**. Source: **Held & Soden 2006 / Allen & Ingram 2002** (the gap
  [[precip-parameterization-source]] already names). Pinned closure: `R_ATM_SLOPE=2 W/m²/K` (the
  atmospheric-column cooling sensitivity `dR_atm/dT̄`) + `⟨P⟩₀≈100 cm/yr` → rate `=slope/(L⟨P⟩₀)≈2.5%/K`.
  **The honesty:** a *cited-closure* result, NOT derived; the slope is the **named sub-grid WALL**. TRAP
  (advisor): `R_ATM_SLOPE` is coincidentally 2 W/m²/K like `B_OLR` but a DIFFERENT quantity (column
  cooling vs the TOA longwave feedback) — never set `r=B` and call it derived.
- **The fixed-RH diffusive moist EBM** — moisture slaved as `q=RH·q_sat(T)` (fixed RH) with down-gradient
  latent diffusion sharing the heat diffusivity (so the latent heat `L` cancels → no new `D_q`). Source:
  **Flannery 1984 / Hwang & Frierson 2010 / Siler, Roe & Armour 2018**. The emergent `P−E` pattern from
  this is the **extratropical-only trade** (see [[planet-rung2-scoped]] for the overturned-subtropics
  finding); the full MSE-diffusing version where `T` responds is the deferred rung 2.5.
