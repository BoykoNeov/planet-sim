---
name: planet-rung2x-itcz
description: "Planet rung 2.x (full-sphere EBM + the energetic ITCZ) BUILT 2026-06-14 (planet/sphere_ebm.py): pole-to-pole sibling, ebm.py/moist.py UNTOUCHED; ITCZ = energy-flux-equator migrates toward warm hemisphere; sensitivity is a CLOSED-FORM consequence of the calibrated D (~5–6 deg/PW, factor ~2 above observed ~3), NOT an emergent prediction; the chosen slice of the ITCZ/Hadley deferral"
metadata: 
  node_type: memory
  type: project
  originSessionId: d1abcdcc-5ec6-4b7a-8c87-7f05451c89e9
---

Planet **rung 2.x** — the **full-sphere EBM + the energetic ITCZ** — **BUILT 2026-06-14**
(`planet/sphere_ebm.py`, 15 tests; full record plan §10 "Rung 2.x — BUILT"). The user's chosen slice of the
**"ITCZ/Hadley + moist precip pattern"** rung-2 deferral ([[planet-rung2-scoped]]): of three slices I
surfaced (wet-get-wetter / symmetric-Hadley-fix / **full-sphere + ITCZ migration**) the user picked the
**structural** one. Built **spike-first** (`outputs/rung_itcz_fullsphere_spike.py`, gitignored) +
**advisor-pressure-tested twice** — the second round reset the headline (record it, it's the load-bearing
correction).

**Why a new model:** rung-0's `ebm.py` is **hemisphere-only** (`x=sinφ∈[0,1]`, equatorial-symmetry BC), so
an ITCZ that **migrates off the equator breaks that symmetry by construction** and can't be represented.
This rung lifts to the **full sphere `x∈[−1,1]`** (two real poles, equator interior) and locates the ITCZ at
the **energy-flux equator (EFE)** — the zero of the emergent atmospheric energy transport
`H(x)=−2πa²·D(1−x²)∂ₓT` — which migrates toward the warmer hemisphere under an interhemispheric imbalance.
**A SIBLING — `ebm.py`/`moist.py` UNTOUCHED** (the [[planet-rung25-mse-diffusion]] / [[planet-rung3-qg-built]]
discipline): "re-validate Phase-1" = a **cross-model reduction check** (sibling reproduces hemisphere
`ebm.py` to **1e-9** under symmetric forcing), the SW↔QG-bridge pattern.

**THE LOAD-BEARING ADVISOR CATCH (the honest altitude — do NOT regress to "EBM predicts the observed ITCZ"):**
In a *dry* EBM the EFE is just the **temperature maximum** (`∂ₓT=0`), so the migration sensitivity is the
**closed form `δ/AHT_eq = 1/(2πa²·D·T̄ₓₓ(0))`** — a near-algebraic consequence of the **already-calibrated
transport `D`** and the mean-state curvature, **NOT an emergent prediction**. The forcing-independence
(cross-equatorial Q-flux and antisymmetric albedo give the **same** number) is a **linear-operator identity,
not robustness**; the shift **direction (toward the warm hemisphere) is by-construction** (the same
"guaranteed result" trap as down-gradient+irr~1 on the QG rung). Values **≈ −6.3 deg/PW (no-ice,
splitting-free) and ~ −5 deg/PW (present-day ice)** — the **same ORDER** as observed **~3 deg/PW** (Donohoe
2013) but a **factor ~1.5–2 high**, and `∝ (6+B/D)` (curvature `∝1/(6D+B)` ⟹ a *pure function of D*, **not**
the naive "∝1/D" — a test premise the build caught). Honest one-liner: **"corroborates `D` is realistic,"
NOT "predicts the ITCZ."** The spike's first **−3.78** (suspiciously near observed-3) was an
**operator-splitting artifact** at the default `n_tau=0.5`; the EFE/sensitivity MUST be read off a
**converged** profile (`steady_linear` or small `n_tau`) — a real gotcha, code+test-enforced (`ebm.py`
itself anchors its North check with `method="direct"`, never the relaxation).

**Triad.** *Tight* — reduction to `ebm.py` (1e-9); North two-mode via direct `steady_linear` (constant
albedo, **harmonic-face polar floor ~0.16 °C, NOT clean 2nd order** — same as `ebm.py`); the **closed form
reproduced by the engine**; EFE=0 exactly for symmetric; global energy balance machine-exact. *Real-but-loose
(unlock, lower altitude)* — the sensitivity above. *Plumbing* — symmetric⟹EFE=0; the precip wire reduces to
rung-0 at `φ_EFE=0` **bit-for-bit**.

**The precip wire (the "moist precip pattern" half):** `precip.precip_pattern`/`precipitation` gain an
**`itcz_center_deg`** seam (default 0 → rung-0 **bit-for-bit**; the ITCZ band now uses *signed* latitude so
it can migrate; midlat bands still symmetric); `sphere_ebm.itcz_informed_precip` feeds it `φ_EFE`. **Honest
scope: a DRY model RELOCATES A PRESCRIBED BAND — NOT emergent rainfall**, opt-in like [[planet-rung1-two-way-coupler]]'s
`circ_precip`. **`moist.py`'s `moisture_convergence` stays backwards in the deep tropics by DEFAULT** — this
rung adds ITCZ **position**; the Hadley moisture-convergence fix (the backwards-`P−E` *sign*) was then BUILT
2026-06-14 as an opt-in mean-circulation term → [[planet-rung2-hadley-fix]].
Other edges: asymmetry **imposed** (Q-flux/albedo, not an ocean); annual-mean (no seasonal migration).

**THE ITCZ-SENSITIVITY "TIGHTENING" — RESOLVED 2026-07-10 as a RADIATION LIMIT (a negative + an identity, NOT a tightening).**
The §12.2 backlog line ("re-derive `D` to land the sensitivity within factor-1 of observed; it *rides* the calibrated `D`")
rested on the naive `∝1/D` premise **this very build already overturned** to `∝(6+B/D)` — so it was self-contradictory with
the code's own tests. Chasing it (advisor-caught) surfaced the real structure: at the symmetric steady state the **equatorial
energy balance pins `D·T̄ₓₓ(0) = −NEI(0)`** (the transport divergence at `x=0`, where `∂ₓT̄=0`, collapses to `D·T̄ₓₓ(0)`, and
steady state sets it to minus the net radiative input). So the sensitivity is **identically `δ/AHT = −1/(2π a² NEI(0))`** —
the **Bischoff & Schneider 2014** `δ ≈ −AHT/NEI` (already in the sources; now the EBM's sensitivity *is* that formula). It is a
**radiation** quantity that `D` **cancels out of** → **no transport tightens it**. Turning `D` up only slides the equatorial
`T` (hence `NEI(0)`) along one curve: `D→0` radiative-equilibrium limit (`NEI→0`, infinite sens) up to the **isothermal
ceiling** `D→∞` (`NEI≈57 W/m²`, floor **`−3.9 deg/PW`**; the `6` in `∝(6+B/D)` = the untunable `P₂` eigenvalue `n(n+1)=2·3`).
Observed `−3` needs `NEI(0)≈75 W/m²` — **above the ceiling** ⟹ unreachable by *any* transport; the lever is a stronger
equatorial radiative surplus (**rung 4**) or GMS dynamics (rung 3+). **Upgrades loose→tight**: the functional form is now the
cited machine-checkable identity (the `NEI` form matches the *measured* migration tighter than the curvature fit), not "factor
~2 corroborates `D`". `∝1/D` doc-fix in `sphere_ebm.py`/tests (added `itcz_sensitivity_from_nei`, `net_radiative_input_equator`,
`efe_from_transport`; +2 tight tests).

**Then the user chose to ALSO try the moist/MSE build** (hypothesis: observed `AHT` is *moist* static energy, ~`(1+β)`× the
sensible energy/degree at the warm equator → smaller deg/PW). BUILT **`planet/sphere_moist_ebm.py`** (`SphereMoistEBM` — a
SIBLING composing rung-2.5's `D_eff(T)` MSE diffusion, [[planet-rung25-mse-diffusion]], onto the rung-2.x full sphere; dt-free
Picard; `RH=0`+`D_s=0.555` ⟹ `steady_linear` **bit-for-bit** 2.8e-13; `test_sphere_moist_ebm.py` 10 tests). **ALSO a negative,
and the identity says why:** recalibrated to the dry contrast (`D_s≈0.28`), moist moves `−6.3`→`−5.7` only ~10 %, **saturating
across RH** — because `D_eff(0)` rises (×1.67) but `T̄ₓₓ(0)` flattens in lockstep (×0.66), so the product `D_eff·T̄ₓₓ = −NEI(0)`
is pinned (two faces of the *same* equatorial moisture amplification). The **entire** moist effect is the ~1.7 K cooler moist
equator shifting `NEI(0)` by *exactly* `−B·ΔT_eq` (machine-tight `1e-6`) — a radiation nudge, not transport. Demo
**`planet/demo_itcz_radiation_limit.py`** → `docs/figures/planet-itcz-radiation-limit.png` (master curve `deg/PW` vs `NEI(0)`
with the transport-reachable band + the dry/moist/floor/observed markers; + the cancellation bars), slow-guarded.

**Demo banked + CI-guarded:** `planet/demo_sphere_itcz.py` → `docs/figures/planet-sphere-itcz.png` (H(x)+EFE;
φ_EFE vs imposed AHT on the closed-form line beside observed ~3; the relocated rain band); the `slow`
`test_demo_reproduces_the_banked_headline` pins it. Tests `planet/tests/test_sphere_ebm.py` (15: tight/unlock/
plumbing fast + 2 slow). Full gate **399 passed, 1 skip**. Sources: Kang+2008 / Bischoff–Schneider 2014 /
Schneider–Bischoff–Haug 2014 / Donohoe+2013; North 1975 (two-mode, already in `ebm.py`). Extends
[[ebm-radiation-source]], [[precip-parameterization-source]].
