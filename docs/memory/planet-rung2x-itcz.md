---
name: planet-rung2x-itcz
description: "Planet rung 2.x (full-sphere EBM + the energetic ITCZ) BUILT 2026-06-14 (planet/sphere_ebm.py): pole-to-pole sibling, ebm.py/moist.py UNTOUCHED; ITCZ = energy-flux-equator migrates toward warm hemisphere; sensitivity is a CLOSED-FORM consequence of the calibrated D (~5–6 deg/PW, factor ~2 above observed ~3), NOT an emergent prediction; the chosen slice of the ITCZ/Hadley deferral"
metadata:
  type: project
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
`circ_precip`. **`moist.py`'s `moisture_convergence` stays backwards in the deep tropics** — this rung adds
ITCZ **position**, the Hadley moisture-convergence fix (the literal backwards-`P−E`) is **STILL DEFERRED**.
Other edges: asymmetry **imposed** (Q-flux/albedo, not an ocean); annual-mean (no seasonal migration).

**Demo banked + CI-guarded:** `planet/demo_sphere_itcz.py` → `docs/figures/planet-sphere-itcz.png` (H(x)+EFE;
φ_EFE vs imposed AHT on the closed-form line beside observed ~3; the relocated rain band); the `slow`
`test_demo_reproduces_the_banked_headline` pins it. Tests `planet/tests/test_sphere_ebm.py` (15: tight/unlock/
plumbing fast + 2 slow). Full gate **399 passed, 1 skip**. Sources: Kang+2008 / Bischoff–Schneider 2014 /
Schneider–Bischoff–Haug 2014 / Donohoe+2013; North 1975 (two-mode, already in `ebm.py`). Extends
[[ebm-radiation-source]], [[precip-parameterization-source]].
