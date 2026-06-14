---
name: planet-rung2x-emergent-itcz
description: "Planet rung-2.x Emergent ITCZ rain BUILT 2026-06-14 (planet/sphere_moist.py): full-sphere moisture budget (eddy + two-cell Hadley anchored on the EFE), conservative ∫(P−E)=0 machine-exact, rain RAINED not PAINTED; the advisor's predicted q-contrast headline was OVERTURNED — displaced-ITCZ intensity is GEOMETRIC not emergent q, the wet/dry dipole is DISPLACEMENT-driven; the two real nuggets are co-location-as-a-falsifiable-check + the clean negative result"
metadata:
  node_type: memory
  type: project
  originSessionId: d1abcdcc-5ec6-4b7a-8c87-7f05451c89e9
---

Planet **rung 2.x — Emergent ITCZ rain** — **BUILT 2026-06-14** (`planet/sphere_moist.py`, 11 tests;
full record plan §10 "Rung 2.x — Emergent ITCZ rain BUILT"). The **"meatiest new finding"** slice of the
ITCZ/Hadley deferral: rung 2.x's precip wire ([[planet-rung2x-itcz]]) only **relocated a prescribed Gaussian
band** to the EFE (a *dry* model painting a belt); this rung carries the rung-2 **column moisture budget**
([[planet-rung2-hadley-fix]]) onto the **full sphere** so the ITCZ rain **emerges from a conserving `P−E`
budget** whose convergence maximum **sits on the EFE** — **RAINED, not PAINTED**. **A SIBLING** —
`moist.py`/`sphere_ebm.py` UNTOUCHED (reuses `moist.specific_humidity`/`q_sat`/`HADLEY_STRENGTH`, re-derives
the conservative operators on the doubled grid). Built **spike-first** (`outputs/rung2x_sphere_moist_spike.py`,
gitignored) + **advisor-pressure-tested twice**.

**The model.** Full-sphere eddy convergence `(D/c_p)·∂ₓ[(1−x²)∂ₓq]` + a **two-cell Hadley** circulation
whose ascent is **anchored at the EFE** (`hadley_streamfunction`: signed Ψ, +south/−north of the EFE,
asymmetric widths — descent edges **pinned at ±30°** ⟹ the physical cross-equatorial cell widens). Both in
**conservative face form** with **two real polar Neumann-0 ends** (NOT moist.py's hardcoded equator-symmetry
BC) ⟹ `∫(P−E)=0` **machine-exact** for *any* asymmetric cell.

**THE ADVISOR OVERTURN (the load-bearing correction — the second round reset the headline).** The advisor
first predicted the headline = *"the displaced ITCZ is more intense because it sits in the warmer/moister
hemisphere; magnitude tracks the interhemispheric q-contrast = emergent."* A **decomposition spike refuted
it** and I surfaced the conflict as a reconcile-call (the "don't silently switch on primary evidence"
discipline); the advisor **conceded**. Under the (constant-albedo) Q-flux the q-contrast is **tiny (~2–11 %)**
and **neither the peak nor the dipole tracks it**:
- **Displaced-ITCZ peak intensification is GEOMETRIC, NOT emergent `q` — a CLEAN NEGATIVE RESULT** (pinned by
  a test so it is not silently re-read as a win): the peak grows because the **pinned-edge near cell narrows**,
  not because the warm side is moister — replacing `q(T)` with a hemispherically-symmetric `q` leaves the peak
  unchanged (**180.8 ≈ 180.5 cm/yr**); symmetric cell widths remove the intensification entirely.
- **The wet-NH/dry-SH dipole is DISPLACEMENT-driven** (present at full strength with a *symmetric* `q`); its
  *direction* (toward the warm hemisphere) is **by-construction**. **DROPPED** (advisor): any claim that `q`
  shapes the meridional asymmetry in either direction — the realq<symq "damping" I first measured was
  **confounded** (symmetrizing `q` also relocates the q-peak from EFE back to equator).

**The two real nuggets (advisor's final altitude).** (1) **Co-location of the NET `P−E` on the EFE = a
FALSIFIABLE CHECK, not a given:** the down-gradient eddy term *exports* moisture from the warm EFE (the rung-2
ITCZ trade), so the prescribed cell must **beat that export at the displaced latitude** — it does by a
**~2.6× margin** (eddy ≈ −113 vs Hadley ≈ +287 cm/yr at the EFE) ⟹ net rain max lands on the EFE to **<1°**
(checked). (2) The **geometric-not-`q` negative** above. The only clean emergent-`q` signature is the **~C–C
warming response** (the Hadley-fix nugget, re-confirmed on the full sphere). Genuinely-new emergent content is
**MODEST** — largely **2.x's EFE displacement × the rung-2 Hadley convergence recombined into a conserving
budget**; the architectural win (**budget-not-band**) is the headline "what".

**Triad.** *Tight* — cross-model reduction to hemisphere `moist.moisture_budget` on the NH at `φ_EFE=0`
**machine-exact** (eddy 0.0, full ~4e-11: symmetric `q` zeroes the equatorial face-flux, collapsing the
stencil to moist.py's equator BC); `∫(P−E)=0` machine-exact (symmetric **and** displaced); symmetric ⟹ even
`P−E` peaking at the equator; two-cell Ψ structure. *Real-but-loose* — co-location-as-a-check (~2.6× margin);
~C–C warming response. *Plumbing/negatives* — co-location *at all* is by-construction (ascent **placed** at
the EFE, NOT derived); dipole direction by-construction; **the geometric-not-`q` peak negative**; `strength=0`
⟹ eddy-only bit-for-bit.

**Walls (carried).** Cell is **prescribed** (`HADLEY_STRENGTH` the wall; fully-emergent `Ψ∝H/GMS` = **rung
3+** — double-counts the mean transport in the EBM's `D`, needs a GMS closure); asymmetry **imposed** (Q-flux/
albedo); **subtropical desert stays mislocated** (hyper-peaked C–C `q`, unfixed); `P−E` not full `P` (no
honest zonal `E` keeps `P≥0`). **Demo+CI-guard** (`planet/demo_sphere_moist.py` →
`docs/figures/planet-sphere-moist.png`; `slow` `test_demo_reproduces_the_banked_headline`). Tests
`planet/tests/test_sphere_moist.py` (10 fast + 1 slow). `sphere_ebm.py` breadcrumb added; `moist.py`/`precip.py`
untouched. Extends [[planet-rung2x-itcz]], [[planet-rung2-hadley-fix]], [[moist-ebm-source]].
