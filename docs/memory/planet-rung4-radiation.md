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

---

**Rung-4 completion — the per-latitude EBM wire BUILT 2026-06-14** (`planet/radiative_ebm.py`, 11 fast +
2 demo tests; the natural completion, was *left-to-user-call*). Wires the emergent gray `OLR(Ts,τ)` into the
EBM **per latitude** (real radiation drives the climate) — a **separate sibling alongside rung-0** (`ebm.py`
untouched). Spike-first (`outputs/rung4_radiative_ebm_spike.py`, gitignored) + advisor-pressure-tested.

**§12's scoping guess "radiative polar amplification" was OVERTURNED by the spike → it is TROPICAL.** The
headline: the OLR slope is **not one number** — its *local* value `B_loc(Ts)=dOLR/dTs` collapses to ~1.0 at
the warm equator (water-vapour feedback) and rises to ~2.4 by the cold pole, so under a uniform forcing
warming **concentrates in the tropics** (endpoint `δT(pole)/δT(equator)≈0.68`, band 0.73) — the **mirror** of
rung-2.5's moisture-*transport* polar amplification (dt-free ~1.8–2.05). **The SIGN was measured, not assumed** (advisor's
discriminator: smallest `B_loc` warms most; WV pulls it *down* at the equator, Planck `4σT³` pulls it *up* —
WV wins). Clean "two mechanisms pull opposite ways" pair with rung 2.5.

**Core = a coupled Newton steady solve** of `L_T·T + S(1−α) − OLR(T)=0` (Jacobian `L_T − diag(B_loc)`) — the
**nonlinear generalisation of `ebm.steady_linear`'s direct mode**, reusing the engine-pinned transport
tridiagonal. **NOT the Strang relaxation**, for two reasons found in the spike: (a) the half-step that is
*exact* for rung-0's linear OLR carries an **O(Δt²) splitting error that does not vanish at equilibrium** once
OLR is nonlinear (relaxed steady ⟨T⟩ drifts with the step → Newton; and **rung-0's own relax default has a
sizeable contrast error**, 47.9 vs `steady_linear` 38.4 — so all rung-0 comparisons here use the *direct*
reference); (b) it goes unstable at the warm-equator runaway. **Runaway finding:** the per-latitude wire
**exposes the local Komabayashi–Ingersoll edge the global column hid** — at rung-4's *default* WV loading 0.5
the equatorial column is *past* local runaway (`B_loc<0` for `Ts≳32°C`), so the wire runs at the
**climlab-matched loading** (`climlab_matched_column`: WV fraction ≈ **0.348** giving global-mean `B=2`, where
`B_loc>0` everywhere; Newton converges ~6 iters as transport stabilises the equator).

**Advisor's Jensen catch (load-bearing): the global mean is NOT pinned at `ΔA/B`** (unlike rung 2.5) — OLR is
concave so `⟨OLR(T)⟩≠OLR(⟨T⟩)` and the WV feedback **amplifies the mean too** (`⟨δT⟩=6.99 > ΔA/B_tan=5.94`);
the moist-EBM "redistribution around a pinned mean" framing **does not transfer**. Present mean-state = a
**Jensen warm shift** (~2°C above rung-0) with **contrast ≈ unchanged** (loading-matched mean slope ≈2) → the
signal is in the *warming response*, not the present climate. **`D` NOT recalibrated** (§12 expected
"recalibrate `D` as rung 2.5" as the cost; but the present contrast is already ≈ rung-0's, nothing to
recalibrate *for*; `D` sets magnitude not sign). **Triad:** *tight* — linear `olr_fn` → `ebm.steady_linear`
**bit-for-bit** (4.5e-13, all departure = OLR curvature) + **net-TOA=0 machine** at convergence; *discriminator
(pure column)* — `B_loc` min at the warmest latitude + WV flips the Planck ordering; *unlock (loose)* —
tropical amp **at the Earth loading**, **both sign AND magnitude ride the WV loading** (advisor's catch — the
single-loading blind spot): a drier Planck-dominated planet is **POLAR** (`amp>1` below crossover ~0.15;
`B_noWV` rises with Ts) — *unlike* rung-2.5's polar direction which IS RH-robust, so the "mirror" is not a
symmetry of equally-robust mechanisms. Part of the 0.68 is **runaway-proximity** (warmed eq ~39°C, `B_loc≈0.5`,
stable but near the hot edge). Null = the present *tangent* (uniform `B`, warms uniformly, `amp=1`). Demo `planet/demo_radiative_ebm.py` → `docs/figures/planet-radiative-ebm.png`.
Remaining within-rung upgrades: spectral-band log law + moist-adiabatic lapse-rate feedback + clouds.
[[planet-rung25-mse-diffusion]] [[moist-ebm-source]]; plan §10.
