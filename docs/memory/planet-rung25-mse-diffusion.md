---
name: planet-rung25-mse-diffusion
description: "Planet rung 2.5 (MSE-diffusing moist EBM where T responds) BUILT 2026-06-12 (planet/moist_ebm.py); headline = emergent polar amplification from moisture transport alone, framed as REDISTRIBUTION around a PINNED ⟨δT⟩=ΔA/B; D_eff=D_s(1+β(T)); recalibrate D_s (double-count wall); B fixed=rung-4 wall. dt-FREE RE-BANK 2026-06-14: the ~1.5/~1.4 headline was an O(Δt) operator-splitting artifact (relaxation @ n_tau=0.5); the dt-free Picard solve moist_steady_direct gives endpoint ~2.05 / band ~1.80"
metadata:
  node_type: memory
  type: project
  originSessionId: rung25-session
---

Planet **rung 2.5** (the MSE-diffusing moist EBM where **T responds**) — **BUILT 2026-06-12**
(dt-free re-bank 2026-06-14), `planet/moist_ebm.py` + `tests/test_moist_ebm.py` (14 fast tests). Continues
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
**EXACTLY uniformly** = the clean null. Earth PA factor `pole/eq δT` ≈ **2.05** (dt-free, RH 0.8);
**direction BANKED** (PA>1 robust), **magnitude LOOSE** (observed ~2–3× also needs ice-albedo+lapse-rate
feedbacks, out of scope; ~2.05 is the model's *converged* number, not a *validated* one). **NB the ~2.05
is the dt-FREE value — see the re-bank note below; the originally-banked ~1.5 was a splitting artifact.**

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
`ΔS₀` — S₀ is equator-weighted by `(1+s₂P₂)`, imposing tropical structure that FIGHTS PA.) **NAME THE
METRIC (advisor catch #2):** the headline (dt-free) ~2.05 is the **single-endpoint** ratio
`δT(pole)/δT(equator)` (most generous, polar cell on the harmonic-face bias); the **area-band** ratio
`mean(δT|≥60°)/mean(δT|≤30°)` ≈ **1.80** (less generous); both honest, both stored on
`PolarAmplification` (`pa_moist`/`pa_moist_band`). A `.converged` guard was also added (the one failure
mode the suite was blind to).

**THE dt-FREE RE-BANK (2026-06-14 — a splitting-error finding, user-flagged).** The headline was first
banked as endpoint **≈1.5** / band **≈1.4** — those were a **first-order operator-splitting ARTIFACT**.
The PA was read off the Strang relaxation at the default `n_tau=0.5`, whose **shape** carries an **O(Δt)
splitting bias** (backward-Euler transport split against the *exact* radiation half-step; the global
**mean** stays exact — `pa_dry≡1.0` exactly is the clean control), and that bias **suppresses** the
amplification. Recomputing dt-free lifts it to **endpoint ≈2.05 / band ≈1.80** (a Richardson dt→0 sweep
and a direct solve agree to 3 digits; the `D_s` recalibration drift is negligible — the fixed-`D_s` sweep
reproduces the drift). **The fix** = `moist_steady_direct`: a **Picard iteration on the frozen-`D_eff`
linear solve `(L_T−B·I)T=A−S(1−α)`**, the *nonlinear generalisation of `ebm.steady_linear`* (exactly as
rung-4 went Newton), **no time-stepping ⟹ no splitting bias**, ~20 iters / ~10 ms — now the headline path
(`polar_amplification` + `recalibrate_sensible_D` use it; the relaxation `moist_equilibrium` is kept for
the bit-for-bit rung-0 reduction + animations). **Direction was NEVER in doubt** — only the magnitude
moved (the way an under-converged number always under-reports). **Caveat (advisor): converged ≠
validated** — ~2.05 is *this model's* number; whether moisture-transport-alone "should" give ~1.5 or ~2.0
is a formulation property, NOT checked against the cited Hwang–Frierson / Flannery / Siler–Roe–Armour PA
values (a deferred discriminator). **Same finding corrected the O(Δt²)→O(Δt) mischaracterisation** in
`sphere_ebm`/`radiative_ebm`/`ebm` docstrings (the *shape* bias is first-order from backward-Euler
transport; the nonlinear-OLR rung-4 case left untouched = unmeasured). **sphere_ebm + rung-4 headlines are
SAFE** (steady_linear / Newton = dt-free already). Spike+findings: `outputs/rung25_splitting_dt_spike.py`,
`rung25_picard_prototype.py`, `rung25_splitting_dt_FINDINGS.md` (gitignored).

**THE ATTRIBUTION NULL (advisor — the backbone of the claim).** The moist model differs from dry two
ways (recalibrated `D`-shape AND T-dependent `D_eff`) — which makes PA? **Freeze `D_eff` at its present
profile and warm → PA = 1.0 EXACTLY** (uniform `δT` solves the perturbation for *any* frozen `D(x)`) ⟹
the PA is **100% the `dD_eff/dT` feedback, 0% the `D`-shape** (test: PA=1.0000, spread~3e-10, via the
genuine array-`D` `EnergyBalanceModel`). Also a bug-catch (a leak would make it ≠1).

**The recalibration = the named WALL (the double-count).** Rung-0's `D=0.555` is an *effective*
diffusivity already absorbing latent transport; explicit MSE diffusion double-counts → `recalibrate_
sensible_D` re-derives a **smaller sensible `D_s≈0.28`** matching the dry present **equator-pole
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
pinned `⟨δT⟩=ΔA/B` + the frozen-`D_eff` PA=1 null; *real-but-loose (unlock)* = the PA (dt-free ~2.05,
direction banked/magnitude loose); *plumbing* = RH=0 ∧ `D_s=0.555` ⟹ genuine `EnergyBalanceModel` rung-0
**bit-for-bit** (relaxation) and **machine-precision** (`moist_steady_direct`→`steady_linear`), + a test
pinning the splitting artifact (`n_tau=0.5` reads ~1.5, climbs to the dt-free value as `n_tau→0`); *named
choices* = recalibration to present contrast (`D_s<0.555`) + target-invariance. Cited sources unchanged
from [[moist-ebm-source]] (Hwang & Frierson 2010 / Flannery 1984 / Siler–Roe–Armour 2018). No engine edit;
`uses` unchanged.
**Deferred (still):** ITCZ/Hadley, rung-3 resolved storm-track pattern (the vertical), rung-4 radiation
(the `B` wall). A demo PA figure is an unbuilt nicety. See [[planet-rung2-scoped]],
[[planet-rung1-two-way-coupler]], [[ebm-radiation-source]], [[planet-plan]].
