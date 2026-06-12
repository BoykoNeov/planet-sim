---
name: planet-rung25-mse-diffusion
description: "Planet rung 2.5 (MSE-diffusing moist EBM where T responds) BUILT 2026-06-12 (planet/moist_ebm.py); headline = emergent polar amplification (~1.5×) from moisture transport alone, framed as REDISTRIBUTION around a PINNED ⟨δT⟩=ΔA/B; D_eff=D_s(1+β(T)) frozen each Strang step like α; recalibrate D_s (the double-count wall); B fixed=rung-4 wall"
metadata:
  node_type: memory
  type: project
  originSessionId: rung25-session
---

Planet **rung 2.5** (the MSE-diffusing moist EBM where **T responds**) — **BUILT 2026-06-12**,
`planet/moist_ebm.py` + `tests/test_moist_ebm.py` (12 fast tests). Continues
[[planet-rung2-scoped]] (rung 2 Phase A was a *pure diagnostic*; 2.5's whole point is **T responds** →
an emergent *climate response* the diagnostic couldn't make). Full plan record
`docs/plans/planet-earth-system.md` §10 "Rung 2.5 — BUILT". Built **spike-first**
(`outputs/rung25_moist_ebm_spike.py`, gitignored) + advisor-pressure-tested — the **advisor caught a
load-bearing math error and supplied the attribution backbone** (record these, they're the corrections).
A **separate model ALONGSIDE rung-0, NOT a replacement** (dry EBM stays default; opt-in moist sibling,
like circ-precip / energy-constrained rate).

**The mechanism (the headline = emergent POLAR AMPLIFICATION).** A moist atmosphere diffuses moist
static energy `m=c_pT+L·q`, `q=RH·q_sat(T)`; in temp-equivalent units this is **temperature diffusion
with a moisture-amplified `D_eff(T)=D_s·(1+β(T))`**, `β=(L/c_p)·RH·dq_sat/dT` (`moisture_amplification`).
C–C `q_sat` is steep ⟹ **β large in the warm tropics, ≈0 at the cold pole** (Earth: `D_eff`~1.3 eq →
~0.35 pole). Warm the planet ⟹ tropical β grows fastest ⟹ tropics export more poleward ⟹ **poles warm
more** — moisture ALONE, no ice-albedo feedback, no change in `D_s`. Dry EBM (constant D) warms
**EXACTLY uniformly** = the clean null. Earth PA factor `pole/eq δT` ≈ **1.5** (RH 0.8); **direction
BANKED** (PA>1 robust), **magnitude LOOSE** (1.43→1.50 across RH 0.6→0.8; observed ~2–3× also needs
ice-albedo+lapse-rate feedbacks, out of scope).

**THE ADVISOR'S MATH CATCH (load-bearing — I had it wrong).** `D_eff` MUST sit **INSIDE** the divergence
`∂/∂x[(1−x²)·D_eff·∂T/∂x]`, NOT outside (`D_eff·∂/∂x[…]`): the cross-term `(1−x²)(∂D_eff/∂x)(∂T/∂x)` is
**same-order** as the main term and *is where part of the PA mechanism lives*; outside also breaks
conservation. **`ebm.py` already places the callable `D(x)` inside** (the rung-1 array-`D` path) → pass
`D_eff` as that callable, let the engine place it (spike-confirmed at machine precision: varying-`D_eff`
transport conserves ∫T dx + direct operator match rel-err ~2e-12).

**Design = 5th reuse of the diffusion spine; D_eff frozen each Strang step like α(T).** ONE nonlinear
relaxation (not nested solves): each substep **freezes `D_eff` at the current T** and rebuilds the
conservative operator — the **identical idiom the ice-albedo α(T) already uses** in `ebm.equilibrate`.
**Self-contained — does NOT edit `ebm.py`** (rung-0 hot path untouched; the ~1 duplicated radiation
helper is the correct price). **`face="harmonic"` PINNED** (bit-for-bit reduction needs the per-step
cells `(D_s/C)(1−x²)` to match the dry model's once-built harmonic cells; `face="exact"` would silently
break it).

**THE HEADLINE FRAMING (advisor) — REDISTRIBUTION around a PINNED mean, not added pole warmth.** Constant
albedo ⟹ `⟨δT⟩=ΔA/B` to machine precision for **any** D (diffusion conserves ∫T dx → transport can't move
the global mean). So moisture **redistributes a fixed `⟨δT⟩` poleward**. This is asserted **TIGHT**
(conservation), and PA is the *shape* of that redistribution. (Forcing = uniform `ΔA` = CO₂ proxy, NOT
`ΔS₀` — S₀ is equator-weighted by `(1+s₂P₂)`, imposing tropical structure that FIGHTS PA.)

**THE ATTRIBUTION NULL (advisor — the backbone of the claim).** The moist model differs from dry two
ways (recalibrated `D`-shape AND T-dependent `D_eff`) — which makes PA? **Freeze `D_eff` at its present
profile and warm → PA = 1.0 EXACTLY** (uniform `δT` solves the perturbation for *any* frozen `D(x)`) ⟹
the PA is **100% the `dD_eff/dT` feedback, 0% the `D`-shape** (test: PA=1.0000, spread~3e-10, via the
genuine array-`D` `EnergyBalanceModel`). Also a bug-catch (a leak would make it ≠1).

**The recalibration = the named WALL (the double-count).** Rung-0's `D=0.555` is an *effective*
diffusivity already absorbing latent transport; explicit MSE diffusion double-counts → `recalibrate_
sensible_D` re-derives a **smaller sensible `D_s≈0.30`** matching the dry present **equator-pole
contrast** (`⟨T⟩` auto-equal — energy balance fixes it from `A,B,ᾱ` independent of D — so contrast is the
natural single scalar). **THE TRADE (not a win):** a single scalar can't reproduce all of dry `T(x)` —
same mean+contrast but a **higher-moment SHAPE residual** (matched moist profile flatter interior,
curvature to the edges). **Target = a modeling CHOICE (named):** matching P₂ amplitude `T₂` instead moves
PA **<5%** (PA set by β-SHAPE, not scaling).

**Scope refinement: "rung 2.5 re-opens `(A,B,D)`" → only `D` re-opened.** `A`=forcing knob, **`B` FIXED**
— re-deriving B's water-vapour content is *local radiation = the rung-4 wall*, not opened. Fixed RH;
constant albedo for the clean Hwang–Frierson experiment. **Separate from the rung-2 `P−E` budget** — a
T-response model, doesn't touch `moist.py`; rung-2's `test_subtropical_evaporative_belt_is_not_reproduced`
did **NOT** flip (a different deliverable, the advisor flagged this regression expectation).

**Triad.** *tight* = exact analytic `dq_sat/dT` (vs finite-diff + ~7%/K) + `D_eff`-inside conservation +
pinned `⟨δT⟩=ΔA/B` + the frozen-`D_eff` PA=1 null; *real-but-loose (unlock)* = the PA (~1.5, direction
banked/magnitude loose); *plumbing* = RH=0 ∧ `D_s=0.555` ⟹ genuine `EnergyBalanceModel` rung-0
**bit-for-bit**; *named choices* = recalibration to present contrast (`D_s<0.555`) + target-invariance.
Cited sources unchanged from [[moist-ebm-source]] (Hwang & Frierson 2010 / Flannery 1984 /
Siler–Roe–Armour 2018). Full planet gate **261 passed, 1 skip**. No engine edit; `uses` unchanged.
**Deferred (still):** ITCZ/Hadley, rung-3 resolved storm-track pattern (the vertical), rung-4 radiation
(the `B` wall). A demo PA figure is an unbuilt nicety. See [[planet-rung2-scoped]],
[[planet-rung1-two-way-coupler]], [[ebm-radiation-source]], [[planet-plan]].
