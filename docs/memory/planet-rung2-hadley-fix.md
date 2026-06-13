---
name: planet-rung2-hadley-fix
description: "Planet rung-2 Hadley moisture-convergence fix BUILT 2026-06-14 (planet/moist.py): the deferred deep-tropical backwards-P−E — an opt-in prescribed Hadley cell flips the ITCZ SIGN; eddy-only stays default; the emergent nugget is the ~C–C amplitude (q(T)); the TRADE is it does NOT relocate the desert (hyper-peaked q); the cubic ψ both removed a half-sine edge artifact AND exposed the canonical-subtropics fix was artifactual"
metadata: 
  node_type: memory
  type: project
  originSessionId: d1abcdcc-5ec6-4b7a-8c87-7f05451c89e9
---

Planet **rung-2 Hadley moisture-convergence fix** — **BUILT 2026-06-14** (`planet/moist.py`; full record
plan §10 "Rung 2 — the Hadley moisture-convergence fix BUILT"). The **last named rung-2 deferral**: the
literal **backwards-`P−E` in the deep tropics** that [[planet-rung2-scoped]] / [[planet-rung2x-itcz]] both
left deferred (2.x fixed ITCZ *position*; this fixes the deep-tropical *sign*). Built **spike-first**
(`outputs/rung2_hadley_moisture_spike.py`, gitignored) + **advisor-pressure-tested BEFORE the build** — the
load-bearing move was setting the **honesty classification up front**, not the code.

**Why it's needed (advisor-confirmed):** the eddy-only `moisture_convergence` (`(D/c_p)·∂/∂x[(1−x²)∂q/∂x]`)
is **structurally** backwards — down-gradient diffusion *exports* moisture from the moist equator; there is
**no diffusive way to converge moisture at a maximum**, so it is **not** gettable by tuning `D`. The mean
circulation is a genuinely-needed *separate* term.

**The model (opt-in; eddy-only stays DEFAULT — every existing benchmark test stays green):**
`hadley_moisture_convergence` / `moisture_budget(..., hadley=True)`. A **prescribed** tropical overturning:
northward MMC flux `F(x)=−strength·ψ(x)·q(x)` is **equatorward** in the tropics (low-level moist branch →
ascent; dry-upper-branch `Δq≈q_surface`), `q=RH·q_sat(T)`. Convergence `P−E=−∂F/∂x` in **conservative face
form** (`_mean_flux_convergence`) ⟹ `∫(P−E)=0` **machine-exact** (a conserving *budget*, not a painted band).
Ascent **pinned at the equator** (hemisphere model; 2.x owns migration); cell vanishes poleward of its edge
⟹ extratropical eddy budget **untouched**.

**THE HONESTY (advisor, load-bearing):** convergence-at-ITCZ / divergence-under-descent is **GUARANTEED BY
CONSTRUCTION** ⟹ **plumbing, NOT a win** (same "guaranteed result" trap as QG down-gradient+irr~1 and the
2.x warm-ward shift). `HADLEY_STRENGTH≈4.2e-4 kg/m²/s` = the named **WALL** (calibrated to observed *order*
~1–2 m/yr, **transparently** — not tuned-then-cited; **NOT** derived). **The emergent, non-vacuous nugget =
the AMPLITUDE:** `q(T)` carried from the EBM ⟹ ITCZ convergence **intensifies at ~C–C (~6.6 %/K)** under
warming, *faster* than the energy-constrained global mean (~2.5 %/K) = the observed **"rich-get-richer" P−E
scaling** (Held & Soden 2006). That's the bankable physics.

**THE TRADE (advisor WARNED it was contingent; the cubic profile REVEALED it):** the fix flips the ITCZ
**sign** robustly but does **NOT relocate the desert**. The emergent dry belt sits **equatorward (~12°) of
the canonical 25–35° subtropics** — the hyper-peaked fixed-RH C–C `q` pulls the flux `ψ·q` equatorward, the
**same** mislocation the eddy budget has, so 25–35° stays `P>E` on **both** paths. **A half-sine ψ first
*appeared* to fix 25–35° — but only via an edge-discontinuity artifact** (ψ′(edge)≠0 → a ~210 cm/yr jump in
`P−E` at 30°). Switching to the **cubic `ψ=(27/4)u(1−u)²`** (ψ′(0)>0 strong ITCZ ascent, ψ′(edge)=0 smooth
merge) **both** removed the jump **and** exposed that the canonical-subtropics flip was artifactual.
Relocating the desert needs a realistic (less peaked) `q` = moist dynamics / resolved vertical, **rung 3+**
(where the **gross-moist-stability / overturning** route — the fully *emergent* cell — is honest; column GMS
just moves the prescription onto vertical quantities it lacks, and deep-tropics is where GMS→0 — named, NOT
built). The pinned `test_subtropical_evaporative_belt_is_not_reproduced` NOTE updated: mislocation
**persists** past this fix (stays green; guards the default).

**Banked:** demo+CI-guard (`planet/demo_hadley_moisture.py` → `docs/figures/planet-hadley-moisture.png`;
`slow` `test_demo_reproduces_the_hadley_fix_headline`). Tests `planet/tests/test_moist.py` (+9 fast +1 slow).
No engine edit; `sphere_ebm.py` breadcrumb updated; `precip.py` untouched. Extends [[moist-ebm-source]],
[[planet-rung2-scoped]], [[planet-rung2x-itcz]].
