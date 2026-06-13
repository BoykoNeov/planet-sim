---
name: planet-rung3-phaseB-outcropping
description: "Planet rung 3 Phase B (saturated baroclinic eddy flux = the headline payoff) SPIKE FINDING 2026-06-13: the free-surface two-layer SW engine OUTCROPS at saturation (h->0, PV detonates) — robust across all params/drag; control = Froude U_s/sqrt(g'H), overshoot drives eta/H~12*Fr. The saturated payoff does NOT live in this engine; it needs the pre-named two-layer QG (rigid-lid) solver (FFT 2x2 spectral PV inversion, ~200 lines, Held-Larichev). Phase B STILL the open bet, now sharpened + re-routed. Phase A (linear growth) stays banked."
metadata:
  node_type: memory
  type: project
  originSessionId: rung3-phaseB-session
---

Planet **rung 3 Phase B** (the **nonlinear saturated eddy heat flux = the headline payoff**: rung-1's
reduction-to-EBM finally non-vacuous via a strong, irreversible baroclinic `⟨v'·interface'⟩`) was
attempted **spike-first** 2026-06-13 (`outputs/rung3_baroclinic_phaseB_spike.py` +
`outputs/rung3_phaseB_debug.py`, both gitignored) and the spike returned a **clear negative result that
re-routes the build** — exactly the "open bet" the [[planet-rung3-scoped]] done-check flagged (Phase A
moved the payoff risk **zero**). **Phase A (the linear growth rate) stays fully banked**; this is a Phase-B
blocker, not a Phase-A regression.

**THE BLOCKER — the free-surface SW engine OUTCROPS at saturation.** At finite amplitude the
saturated/overshoot **interface displacement reaches the layer depth `H`** → the layer thickness `h→0` →
`PV=(f+ζ)/h` detonates (EKE→1e295). The control parameter is the **Froude ratio `Fr = U_s/√(g'H) ≈
η_sat/H`**. The killer is the **first-saturation OVERSHOOT**: empirically the peak displacement is
**`η/H ≈ 12·Fr`** — *far* above the ~4–5× RMS vortex-tail factor, because the first baroclinic life cycle
dumps the mean APE into eddies before equilibrating (EKE overshoots to ~5–19 then collapses). **Robust
across `g'∈{0.2,0.8,1,2}`, `H∈{400,500}`, `U_s∈{2,4}`, drag `r∈{0.5…3σ}`:** stronger drag only **delays**
the overshoot (`r=3σ` still `η/H=0.91`), stronger stratification gives a **bigger** overshoot (more stored
APE → `g'=2` overshot to EKE≈19 vs `g'=1`'s ≈9). Avoiding it needs `Fr≲0.04` ⟹ `U_s≲1` ⟹ e-fold `≳370 h`
— **full QG-regime cost in the explicit free-surface tool worst-suited to it.** The engine reproduces the
**analytic linear growth rate at every config tried (growth-err 2.2–2.6 %)**, so this is **not** leaving
the Phase-A-validated regime — it is the **finite-amplitude free-surface wall**.

**THE DISSIPATION (built in the spike, the advisor's load-bearing under-specification CONFIRMED).** The
plan named only hyperviscosity; that is **insufficient**. A **fixed-`G` background is an infinite APE
reservoir**, so hyperviscosity (a *small-scale* enstrophy sink) does nothing about the **inverse energy
cascade** condensing into the box-scale barotropic mode — **linear bottom Ekman drag `−r·u₂` on the lower
layer (Held–Larichev 1996) is required** to arrest it. The spike built **both** in a
`DissipativeLayered(LayeredShallowWater)` subclass: **hyperviscosity `−ν₄∇⁴` on momentum AND thickness +
bottom drag**, all **default-off** (`ν₄=r=0` → the Phase-A engine **bit-for-bit**, verified), **momentum
biharmonic + flat-bottom mass form keep per-layer mass machine-exact with dissipation ON** (verified).
Hyperviscosity has its **own CFL** `dt≲0.035·dx⁴/ν₄` (min'd with the gravity CFL in `max_dt`). These
dissipation operators are correct and reusable — they are the named Phase-B prerequisite, now built and
de-risked; only the **host model** is wrong.

**INTERPRETATION + RE-ROUTE (advisor, my lean A).** This is the **empirical reason two-layer turbulence is
done in QG** — Held–Larichev 1996 *is* two-layer QG; the free-surface explicit SW model has a **hard
thickness floor QG does not**, so **the saturated payoff does not live in this engine**. The plan
**pre-named** the rigid-lid/QG fork as the within-rung upgrade "*if the saturated runs make the penalty
bite*" → **they bit, via OUTCROPPING (not just the external-mode CFL).** Executing the plan's own
contingency, not a surprise. **`Fr = U_s/√(g'H)` is the banked control number; `η/H≈12·Fr` the overshoot
law.**

**P2 (the pre-outcrop finite-amplitude flux) REJECTED as the headline (advisor).** Measuring `κ` during
the growing/overshoot phase is **near-vacuous**: an unstable baroclinic mode fluxes buoyancy down-gradient
**by definition** (already implied by the Phase-A growth rate), and with **no sloshing to overcome
`irr~1` is trivial** — *not* the meaningful irreversibility contrast with rung-1's ~90 %-reversible
barotropic flux. Keep as a one-line note, do **not** bank it.

**THE NEXT BUILD (Phase B's real home) — two-layer QG (rigid-lid).** An **FFT-based 2×2 spectral PV
inversion** model (~200 lines, the **canonical Held–Larichev** doubly-periodic setup), **no outcropping**
(no thickness floor) and **no external-mode CFL** (the fast barotropic mode is filtered). Its **linear
anchor already exists** = the spike's **QG-Phillips cross-check** ([[planet-rung3-scoped]]: rigid-lid
`g→∞` σ_SW/σ_Phillips→1.004). **Honest cost (advisor):** a **new model OUTSIDE `engines/fluid`** (no
C-grid reuse, no bit-for-bit single-layer reduction) → **Phase A and Phase B would validate *different*
models** (the coherence price). And **QG makes the experiment POSSIBLE, it does NOT pre-guarantee** the
emergent `κ` comes out irreversible / down-gradient / well-scaled — **that is still the open bet, finally
*testable***. **USER PICKED A (bank + scope QG)** 2026-06-13 over build-QG-now / deep-layer-SW — so the **QG build is now
SCOPED in plan §10** (model = `q_k=∇²ψ_k+(−1)^k F_k(ψ_1−ψ_2)+βy`, 2×2 spectral PV inversion, β-plane for
Rhines arrest + critical shear `U_crit=β/F`, porting the spike's dissipation; triad/edges named) but **NOT
built this session**. Home TBD next session: `planet/baroclinic_qg.py` (single-consumer, lean) vs
`engines/spectral/` (reusable). Pre-registered discriminators carried forward: **irr fraction O(1)** (vs
rung-1 ~0.1) and **`κ_eff/(v'_rms·L_d)` O(0.1–10)** (vs rung-1 ~1e-3) — validate
**dimensionless** (idealized `κ_ML~v'·L_d` is intrinsically 15–60× below Earth's `κ₀=2.2e6`; a direct
compare manufactures a false failure). Sources: **Held & Larichev 1996** (two-layer QG turbulence),
Phillips 1954 / Eady 1949, Vallis 2017. See [[planet-rung3-scoped]], [[planet-rung1-two-way-coupler]],
[[planet-plan]]; plan §10.
