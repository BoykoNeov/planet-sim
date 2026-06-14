---
name: planet-rung4-radiation
description: Planet rung 4 (gray radiative transfer) BUILT — the emergent OLR decomposes climlab's B=2 into the Planck slope minus the water-vapour feedback
metadata:
  type: project
---

**Rung 4 (gray radiative transfer) BUILT 2026-06-14** (`planet/radiation.py`, 13 fast + 1 slow tests,
full planet gate green). Retires the wall every rung from 2.5 on cites — **`B` (the OLR slope) held
fixed**: the linear `OLR = A + B·T` is made **emergent** by a **gray radiative–convective column**
(Schwarzschild two-stream over an optical depth set by the greenhouse-gas amount). A **separate sibling
model alongside rung-0** — `ebm.py` UNTOUCHED (the `moist_ebm`/`sphere_ebm`/`baroclinic_qg` discipline).
Built spike-first (`outputs/rung4_radiation_spike.py`, gitignored) + advisor-pressure-tested before the
build. Retires the wall named in [[planet-rung25-mse-diffusion]]; extends [[ebm-radiation-source]].

**THE HEADLINE = a DECOMPOSITION, not a trade (advisor's load-bearing reframe).** Present-day energy
balance *forces* `OLR≈ASR≈239`; calibrating `τ_s` to the 33 K greenhouse makes the emergent OLR pass
through that operating point **by construction**, so gray and `A+B·T` agree in *value* at present and
**the slope `B` is the finding** (`A,B` linked via `A=239−B·T̄` — the *point* forced, the slope open, NOT
both recovered independently). The slope decomposes: **no-WV `B≈3.41`** sits at the `4σTe³≈3.75`
**emission-level Planck touchstone** (advisor bug-guard: ≈3.8 validates the radiative core, ~1.5 = a bug),
**above climlab's 2**; turning **water vapour** on (`τ(Ts)` via C–C, reusing
`moist.saturation_specific_humidity`) lifts the emission level to colder air → **`B` drops THROUGH 2** (to
~1.33 at the nominal 50% WV loading). So **climlab's `B=2` ≈ Planck − water-vapour + the lapse-rate
feedback the gray column OMITS**, and every term is **ORDER-VALIDATED against Soden & Held 2006, NOT
tuned** (advisor's WV-hardening; pinned off the paper's own text via the local PDF): no-WV `3.41`≈Planck
`|λ₀|3.1–3.2`; WV feedback `2.08`≈`λ_wv 1.8`; gray net `1.33`≈clear-sky Planck+WV `3.2−1.8`; the gap from
`1.33`→2 ≈ `|λ_LR|0.84`, the lapse-rate feedback a **fixed** Γ (uniform warming) **cannot produce**
(climlab's obs-tuned B folds it+clouds in). **Non-circular** (rung-2.5 frozen-`D_eff` / ITCZ closed-form
flavour); only the *exact* magnitudes ride on `WATER_VAPOUR_FRACTION` (the WALL; Planck+WV alone=2 at
loading 0.35). The demo's panel-2 **waterfall** shows Planck−WV+LR=2.17≈climlab's 2 at a glance.

**Tight (DERIVED, not recalled — the rung-3 `K²=2F` lesson):** the gray-RE closed form is derived
in-module from the two-stream eqs — `σT⁴(τ)=½σTe⁴(1+τ)`, skin `Te/2^¼`, **ground `σTg⁴=½σTe⁴(2+τ_s)`**
(the surface–air discontinuity, exactly where a recalled coefficient goes wrong). `solve_gray_equilibrium`
(independent radiosity relaxation, NO analytic input) reproduces profile+skin+**ground** at **~2nd order**
with **`OLR=σTe⁴` machine-exact** (conservation); the `Tg→`derived-coefficient convergence IS the
recalled-coefficient guard. **Wall** = gray band-independent absorption + prescribed `τ↔GHG-column` mapping
(calibrated to *order*, the `R_ATM_SLOPE`/`HADLEY_STRENGTH` cited-closure status).

**Named edges (each pinned):** (1) **CO₂ forcing is SATURATING, not logarithmic** — gray gives a concave
`OLR(τ)` (per-doubling `ΔF` 48→53→41→25→20 W/m², *decreasing*, not the constant-per-doubling Myhre log) at
an **unrealistic whole-band magnitude** → log law + realistic magnitude = **band physics = named
within-rung upgrade**; (2) **clear-sky only** (clouds out of scope); (3) **no lapse-rate feedback** —
fixed Γ ⟹ uniform-shift warming ⟹ zero LR feedback (= why gray net sits below climlab's 2 by ~λ_LR; a
moist-adiabatic Γ is the named upgrade that supplies it); (4) **single column** — wiring `OLR(Ts,τ)`
*per-latitude* into `ebm.py` (real radiation DRIVING the climate, opt-in sibling EBM) = the natural rung-4
completion, **LEFT TO A USER CALL, not foreclosed** (advisor corrected my deferral reasoning: NOT "B_eff
off" — at the climlab-matched loading global-mean B=2 exactly; the real reason = emergent `OLR(Ts)` is
**nonlinear** ⟹ per-latitude slope differs (cold pole vs warm equator) ⟹ re-opens the meridional profile,
a **feature** as much as a risk); (5) **linearization breaks far from present** (steep WV →
Komabayashi–Ingersoll runaway — `B` not linearized across it). **Reduction:** `OLR(Ts)` locally affine near present (rung-0's line = its
tangent, residual <0.5 W/m² over ±3 K; wide-range curvature IS the feedback) + climlab-`(A,B)` consistency
(`B=2` through the forced point recovers `A≈210`, asserted vs `ebm.A_OLR`). **Rung 4 core COMPLETE**
(per-latitude EBM wire + spectral-band log law = named within-rung upgrades). Sources pinned at build:
Pierrehumbert *PoPC* §4 / Goody & Yung (gray RE + two-stream closure); Trenberth–Fasullo–Kiehl 2009
(operating point); Myhre+ 1998 (the log law). Demo `planet/demo_radiation.py` → `docs/figures/planet-
radiation.png` (3 panels: emergent OLR + slope-decomposition/WV-loading sweep + saturating forcing).
[[moist-ebm-source]]; plan §10.
