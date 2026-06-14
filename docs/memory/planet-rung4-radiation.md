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

---

**Rung-4 lapse-rate feedback BUILT 2026-06-14** (`planet/radiation.py`: `GrayRadiationColumn(moist_adiabat=True)`
+ `moist_adiabat_temperature` + `feedback_kernel` + `LapseRateFeedback`; 9 fast + 1 slow tests
`test_radiation_lapse_rate.py`, full planet gate green). The §12 within-rung slice that makes the lapse-rate
feedback **EMERGENT** (rung-4 core had *imported* `λ_LR=0.84` from Soden & Held). Swaps the fixed convective
`Γ` for a **moist adiabat** (derived by its limits — dry `→g/c_p≈9.8 K/km`, flattens to ~4 when warm) that
**flattens as it warms** → surface warming amplifies in the upper troposphere (`ΔT_aloft/ΔTs≈2.8` peak) →
`OLR(Ts)` steepens. **Default-off flag** (`moist_adiabat=False` = bit-for-bit the rung-4-core column).
Spike-first (`outputs/rung4_lapse_rate_spike.py`, gitignored) + advisor-pressure-tested before AND after.

**THE §12 SCOPED ANCHOR WAS OVERTURNED** (the rung-4-wire "radiative→tropical OVERTURNED" pattern). §12 said
"supplies `λ_LR≈0.84`, closing the gap from gray net 1.33 to climlab's 2." It does **not**: emergent value
**≈+1.5**, the moist-adiabat column **OVERSHOOTS** — with-WV `B≈3.1` sits *above* climlab's 2, not at it
(fixed-Γ default = 1.33, below 2). **Banked TIGHT:** sign (`λ_LR>0`, adds to B); kind (upper-trop
amplification, **measured** via `ΔT_aloft/ΔTs>1`, not assumed); **kernel closure** (advisor's load-bearing
design — a ONE-column Soden–Held split: Planck=uniform warming τ-fixed, LR=profile's *departure* from uniform
τ-fixed, WV=`τ(Ts)` change profile-fixed; sums to `B_total` ~9e-4, **cleaner than a two-column diff** which
conflates LR with the Planck-base shift); resolution convergence (1.5105→1.5132, n=100→800). **Magnitude
LOOSE, TWO named reasons (advisor):** (a) single *global* moist-adiabat column = the **TROPICAL** branch only
(extratropics not moist-adiabatic; the extratropical bottom-heavy-warming branch pulls the global mean to
0.84) — recovers the tropical feedback NOT the global mean; (b) rides the prescribed vertical `τ` shape +
`WATER_VAPOUR_FRACTION` (the rung-4 wall, headline at default 0.5), which set the emission level. **Null isn't
perfectly clean:** fixed-Γ itself shows `≈−0.25` **tropopause-migration residual** (the strat floor doesn't
warm) — not a true LR feedback. **Reconciliation (mandatory, advisor — docstring updated):** the existing
`λ_LR≈0.84` is a **global-mean touchstone** for what *fixed*-Γ omits; this emergent value is the **tropical
branch** — both true about different things. **Two-column cross-check DEMOTED:** `B_WV(moist)−B_WV(fixed)=+1.79`
= `ΔLR(1.77)+ΔPlanck(−0.01)+ΔWV(+0.03)`, `ΔLR=λ_LR(moist 1.51)−λ_LR(fixed −0.25)` — proves the kernel
isolated LR (`B_noWV`=**Planck+LR**, NOT kernel Planck — the subtlety that made the two estimates differ;
advisor's arithmetic correction: their `ΔPlanck=0.24` was really `Δ(Planck+LR)`). **Honesty edge (advisor):**
clean WV/LR *separation* is partly a model artifact — `τ_wv` tracks SURFACE `Ts` not the profile, so the
upper-trop moisture–temp coupling that links the two feedbacks in reality is absent. Within-rung upgrades
left: moist adiabat WITH latitudinal structure on the per-latitude wire (recovers the extratropical branch +
global mean); spectral-band log law; clouds. Demo `planet/demo_lapse_rate.py` →
`docs/figures/planet-lapse-rate.png`. [[moist-ebm-source]]; plan §10.

---

**Rung-4 spectral-band log law BUILT 2026-06-14** (`planet/radiation.py`: `SpectralCO2Band` +
`planck_flux_per_wavenumber` + `_transmission_emission` + module band/Myhre constants; 9 fast + 1 slow
tests in `test_radiation.py`, full planet gate green). The §12 within-rung slice that fixes
the rung-4-core CO₂-forcing edge: gray's band-independent absorption SATURATES (per-doubling `ΔF`
48→53→41→25→20, decreasing — adding CO₂ pushes the WHOLE Planck spectrum to the cold upper atmosphere), but
the observed law is LOGARITHMIC (Myhre 5.35·ln(C/C₀) ≈ 3.7 W/m² *per doubling, constant*). **Separate
opt-in construct; gray's `co2_forcing` + its saturation test UNTOUCHED** (advisor #5). Spike-first
(`outputs/rung4_spectral_spike.py`, gitignored) + advisor-pressure-tested before the build.

**THE §12 ANCHOR HELD (not overturned — the rung-4-wire/lapse-rate slices both overturned; this one
confirms):** per-doubling `ΔF` becomes CONSTANT (~4.5 W/m²/doubling, the Myhre band) vs gray's decreasing
48→20. **Model = a band-RESOLVED COLUMN** (advisor #2 let the spike pick it over the sharp τ=1 fallback):
the CO₂ 15-µm band split into `n_bins` spectral bins with absorption `k(ν)=k_c·e^(−|ν−ν₀|/l)` (EXPONENTIAL
WINGS = the whole ingredient), each bin a gray sub-problem solved with the SAME two-stream
transmission-weighted emission kernel over the column's fixed-Γ profile + `(p/p_s)` CO₂ shape, spectral
Planck `πB_ν` source. Band centre deeply saturated (`k_c=1000`) so forcing comes only from the wings; an
exp wing's `τ=1` level spreads a CONSTANT spectral width `2l·d(lnC)` per doubling ⟹ constant `ΔF`.

**Triad (advisor #3 — the tight leg here is WEAKER than the gray core's, framed honestly):** (1)
**independent anchor = REDUCTION-TO-GRAY** — `_transmission_emission` (written independently of
`GrayRadiationColumn._olr_from`) with the gray whole-spectrum `σT⁴` source reproduces `_olr_from` to
MACHINE PRECISION (`rel<1e-12`, residual = float mul-order ULP — exactly what two independent
implementations agreeing looks like); + an END-TO-END check (the full `band_olr` path, full-spectrum
uniform-`k` band on a `wv_fraction=0` column, reproduces `outgoing_longwave` to ~0.1% = Planck-grid
truncation — advisor's belt-and-suspenders, exercises the per-bin weighting directly); + a UNIFORM-`k` band
(no wings) SATURATES like gray (the wing is the ingredient); + `π∫B_ν dν=σT⁴` pins the spectral Planck. (2) **the unlock (loose) = constant per-doubling in
the Myhre band** (2–6 W/m²) vs gray's 20–53 — but the **MAGNITUDE IS THE WALL** (advisor #4): rides the
wing scale `l`, band-centre τ, half-width — calibrated to ORDER, and "CO₂ wings ≈ exponential" is itself an
empirical input; the **FUNCTIONAL FORM (logarithmic) is the win**, NOT the ~3.7 coefficient (the column
realizes ~4.5, ~20% above Myhre — layer smear). (3) **derivation/consistency (NOT independent)** —
`log_law_coefficient` = `dF/dlnC = 2l·π[B_ν(Ts)−B_ν(T_strat)]` (cold-to-space τ=1 limit) matches the SHARP
emission-level model to ~1% (3.54 vs 3.50); the column's finite-layer emission sits ~20–30% above it (4.5).

**RANGE-LIMITED — both named edges spike-confirmed (advisor #1, the #1 spike risk):** constant only in the
flat middle (band centre saturated AND wings not exhausted). **Low-CO₂ edge** (`C≲1/k_c`): per-doubling
GROWS (linear/√, band centre un-saturated) — pinned with a weak `k_c=4` band. **High-CO₂ edge**: per-doubling
FALLS again once the active wing reaches the finite band edge — spike-confirmed with a narrow band. The
realistic `0.5×–8×` test window sits solidly in the flat middle (both edges far outside at the build params).
Demo `planet/demo_spectral_band.py` → `docs/figures/planet-spectral-band.png` (2 panels: per-doubling gray
vs spectral vs Myhre; cumulative `F(C)` on log-x = straight line for spectral, gray bends over). Sources
pinned: Myhre+ 1998 (the log law); Pierrehumbert *PoPC* §4 / Wilson–Gea-Banacloche 2012 (the exp-wing
emission-level mechanism). Rung-4 within-rung upgrades left: clouds; moist adiabat with latitudinal
structure. [[moist-ebm-source]]; plan §12.2 + §10.
