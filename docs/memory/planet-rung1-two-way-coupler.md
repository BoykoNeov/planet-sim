---
name: planet-rung1-two-way-coupler
description: Planet rung 1 (two-way coupler) — ALL 3 STEPS BUILT, RUNG 1 COMPLETE 2026-06-11; step 1 tracer, step 2 emergent eddy flux (state-dependent κ_eff, magnitude/reduction honesty), step 3 circ-informed precip (storm-track band ← emergent jet_lat; seam+reduction+migration-mechanism banked, the "trade not accuracy" finding); next = rung 2 (moist dynamics, emergent precip)
metadata: 
  node_type: memory
  type: project
  originSessionId: 88f3e2c8-5355-42ce-bbb4-2fe776581b65
---

**Rung 1 of the GCM staircase (the two-way coupler) STARTED 2026-06-10.** Plan §5 splits it into
3 sub-deliverables: (1) advected tracer in `engines/fluid`, (2) close the heat budget back to the
EBM, (3) circulation-informed precip. Freeze-before-reuse... no longer applies ([[engines-living-contracts]]),
but the dependency order still makes the engine extension step 1. See [[planet-plan]], [[planet-phase4-coupler]].

**Step 1 BUILT (commit f7ed2e5):** `engines/fluid` advects its `tracer` slot — a scalar `θ` carried
in **flux form** (`∂(hθ)/∂t = −∇·(hθ u)`) through the same SSP-RK3, **strictly passive** (zero
back-reaction on `h,u,v` → dry dynamics bit-for-bit unchanged, the re-seal). `tracer_mass` (`∫hθ`,
machine-exact — the anchor) + `tracer_variance` (`∫½hθ²`, bounded) diagnostics; `test_tracer.py`.
**Honesty correction vs the plan table:** variance is **bounded** (enstrophy honesty class — round-off
smooth, cascades only under filamentation), **NOT "dt-convergent"** (measured: dt-halving sat at
round-off ~1e-14). Scheme is **not monotone** by default (no flux limiter → Gibbs over/undershoot on sharp fronts
= the named scope edge). **TVD limiter BUILT 2026-06-14** (the §12.2 within-rung slice;
`engines/fluid/shallowwater.py`, opt-in `tracer_limiter=` ∈ {minmod, vanleer, mc, superbee}, default-off
→ centered byte-for-bit; `test_tracer_limiter.py`, full-repo gate 432 fast): `θ_face=θ_up+½ψ(r)(θ_down−θ_up)`
with **ψ≡1 = the existing centered average** (unlimited = the ψ≡1 special case); anchor met = a uniform-flow
step develops **no new extrema** under all 4 limiters (rigorous 1-D TVD, h frozen by an exact steady state)
while `∫hθ` stays machine-exact. Honest edge (advisor): strict TVD is **1-D only** — 2-D dimension-split
limiting (Goodman–LeVeque) gives no maximum principle, so the 2-D test asserts only *reduced* overshoot vs
centered (limiter is dissipative ⇒ variance one-sided-decreases); Sweby `0≤ψ≤min(2,2r)` + raise on unknown
name (ψ≡2 is bounded yet non-monotone). WENO not needed — TVD meets the anchor at a fraction of the code.
**Advisor done-check caught a real sign bug** in the `U<0`/`V<0` smoothness ratio (`θ−θ_plus` not `θ_plus−θ`;
pinned by the linear-ramp `r=+1` criterion both directions) that the `+x`-only anchor missed — but its
signature is **direction-asymmetric OVER-DIFFUSION** (neg-flow peak 0.83 vs 0.96), **not** overshoot (a
step stays monotone under the bug too), so the regression is a **reflection-symmetry** test (peak retention
bit-identical for ±x/±y on the correct scheme) + the ±-parametrized monotone step. Lesson: a limiter sign
error can surface as asymmetric *accuracy* loss rather than overshoot (boundedness-preservation here was
observed, not guaranteed — at an extremum the buggy ψ>0 *can* inject one), so **test direction symmetry, not
just over/undershoot**.

**THE de-risking finding (step-0 probe, advisor-driven; `outputs/rung1_stability_probe.py`, gitignored):**
a passive tracer on a *steady zonal* jet (v̄≈0) transports ZERO meridional heat — real transport needs
eddies, and a dry single layer has NO baroclinic instability (rung 3). BUT the existing Phase-4 jet is
**barotropically unstable** (Rayleigh–Kuo met; a 0.5 m/s v-perturbation grows ~200× in meridional KE
over ~20 inertial periods then saturates). So **emergent eddies DO exist** → a passive tracer gets a
real `⟨v'θ'⟩` flux; step 2 does **not** need an imposed stationary-wave crutch. (Consistency: the
instability needs ~10 periods to bite, exactly why Phase-4's 10-period release test still saw a
persistent jet.) Caveat: eddies feed off the *windowed* jet's edge-curvature → transport **magnitude is
window/forcing-tuned**, echoing Phase-4's "flanking-easterly is window-shaped" edge.

**THE anchor reframe (advisor + the physics; governs step 2 — do not re-litigate):** the rung-1 anchor
is **reduction-to-EBM** (resolved flux → effective diffusivity `D_eff`; down-gradient limit recovers the
rung-0 diffusive EBM — sibling of the coupler's release test), **NOT** the plan-table "~5–6 PW" (an
eddy/baroclinic = rung-3 number; here magnitude is tuning, not validation). Rung 1 banks *the two-way
feedback loop closing*, transport magnitude named as a tuned scope edge.

**Step 2 — Phase A BUILT 2026-06-10 (the two-way feedback MACHINERY; `planet/transport.py`).** Advisor-driven
**A/B split**: (A) the anchored feedback spine — given *any* flux `⟨v'θ'⟩(φ)` — lands *independent* of the
eddy sim by driving it with a **synthetic** down-gradient flux (the Phase-4 synthetic-gradient playbook);
(B) the emergent eddy sim is Phase B (NOT built). Pieces: the **κ→D bridge** `D = C_atm·κ/a²` (physical/citable;
`C_atm = c_p·p_s/g ≈ 1.04e7`; rung-0 `D=0.555` ⟺ `κ≈2.2e6 m²/s` = observed midlat eddy-diffusivity order —
the advisor verified ocean-C drops at equilibrium so the transport carries the ATMOSPHERE's column capacity);
**band-bulk** down-gradient diffusivity `κ_bulk=−Σ(F·g)/Σ(g²)` over the **window-flat interior**; re-equilibrate
the EBM at uniform `D_eff` (headline) — `EnergyBalanceModel` now also takes a **callable `D(x)`** (the band-limited
diagnostic; scalar path bit-for-bit unchanged = refactor hygiene, NOT the anchor). The *design* anchor = **reduction-to-EBM**
(structural form-match: `⟨v'θ'⟩=−D_eff·∂θ̄/∂y` has the SAME FORM as the EBM transport term → two-way model with constant
flow-diagnosed `D_eff` *is* a rung-0 diffusive EBM with that `D`). **ADVISOR CAUGHT MY OVERCLAIM (fixed in a follow-up
commit): in Phase A the reduction reduces to rung-0 BY CONSTRUCTION** (`two_way_pass` literally re-runs the scalar-`D`
rung-0 EBM at `D_eff`) → it's PLUMBING not an independent anchor; and my original reduction test applied the κ→D bridge on
BOTH sides so a wrong `a²`/`C_atm` would CANCEL (not caught). **What Phase A actually validates TIGHT:** (i) the bridge
**PINNED ABSOLUTELY** (`C_atm≈1.037e7`, `κ₀≈2.17e6` — not just round-tripped); (ii) **right-signed response** (stronger
flux ⇒ flatter contrast — the EBM's genuine physical response, non-tautological); (iii) plumbing (κ recovered + `D_eff`
routed + up-gradient sign rejected). `tests/test_transport.py`; planet gate **151 passed, 1 skip** (no shared-engine edit
→ planet gate not full-repo). **The bridge is PHYSICAL; only the eddy-flux magnitude (Phase B) is forcing/window-tuned.**
**Phase-B geometry forward-flag (advisor):** the bridge `D=C_atm·κ/a²` is derived for UNIFORM κ on the sphere; Phase B's
latitude-varying `D_eff(φ)` diagnosed on a flat Cartesian β-plane then fed into the spherical `∂/∂x[(1−x²)∂/∂x]` operator
is where the geometry correspondence stops being free — that's where the genuinely tight reduction (flux-divergence = EBM
operator) belongs, NOT inherited by-construction.

**Forks settled (advisor):** Fork 1 = **release-primary** (NOT "parameter-light" — tuning moves to seed/window;
removed by measuring **post-saturation** via a **life-cycle integral** `D_eff=−∫F̄dt/∫θ̄_y dt`, the temporal half
of Fork 2 I'd left open). **"Loop closes" scoped to ONE feedback pass** (burst-averaged `D_eff` too noisy to claim a
converged fixed point — that's a `slow` demo if it converges cleanly).

**Step 2 — Phase B BUILT 2026-06-11 (the EMERGENT eddy flux; `planet/eddy_flux.py` + geometry in `transport.py`).**
`eddy_life_cycle`: spin up the jet **dry** → init `θ`=windowed-EBM profile → deterministic `cos(kx)` v-perturbation →
forcing-**OFF** release → band-bulk `κ_eff=−∫F̄dt/∫θ̄_y dt` over the window-flat interior. **De-risked in throwaway
spikes BEFORE the module (advisor-gated — `outputs/rung1_*spike*.py`, gitignored), and the spikes overturned 3 of my
priors:** (1) the instantaneous `⟨v'θ'⟩` is **largely REVERSIBLE** (oscillates sign every ~2 periods with the meander);
the life-cycle time-integral CANCELS it, leaving a small irreversible down-gradient residual (**irreversible fraction
~0.1**, NOT ~1e-3 — I conflated that with κ/κ_expected). (2) **PURE RELEASE is correct; tracer RELAXATION during release
BREAKS the climate ordering** (a τ_θ artifact: relaxed→flat>steep FAIL, pure-release→flat<steep PASS) — rejected, even
though the prompt said "θ relaxed to the EBM target" (init-at-release is equivalent since steady jet v̄≈0). (3) **WARMING
via CO₂/S₀ barely flattens the CHANNEL gradient (~3.6%)** — linear OLR shifts `T₀` uniformly, ice retreat is poleward of
the channel — so the non-circularity knob is **`s₂`/obliquity** (flattens the channel ~80% cleanly), NOT warming.

**What Phase B BANKS, by honesty class:** **(headline, DIRECTION) state-dependent diffusivity** — flat (`s₂=−0.32`, jet
~14) `κ_eff≈0.5–0.6×` steep (`s₂=−0.48`, jet ~20), `α` HELD FIXED → a real right-signed feedback (flatter→weaker jet→
weaker eddies→smaller κ), PASSES robustly across all windows. Mechanism (advisor): the gradient **cancels** in
`−∫F̄/∫θ̄_y` so `κ_eff≈v'·ℓ` (stirring rate) tracks climate ONLY through the jet → **must hold α fixed** (renormalizing to
fix jet speed = cosmetic by construction). **(MAGNITUDE — named, NOT banked)** `κ~10³ m²/s`, ~1000× below rung-0's
`2.2e6`; **resolution-CONVERGED** (nx80≈96) but **don't claim "barotropic is intrinsically 1000× weak"** (advisor:
suppressed by single coherent seed mixing more reversibly than broadband + band-bulk smearing the jet peak — can't
separate config from physics; window/forcing/config-tuned like Phase-4 jet *speed*; don't chase a bigger number,
unbankable). **(the tight reduction = a FINDING, NOT a match)** `reduction_to_ebm_operator` TESTS (not assumes) the
resolved `−∂F̄/∂y` shape vs smooth down-gradient diffusion built from the band-bulk SCALAR κ (not circular, normalized
shape only) → **partial corr ~0.6**, and near-VACUOUS (uniform-κ diffusion of the near-linear midlat gradient ≈ 0
divergence) → the tight reduction is **non-vacuous only at rung 3** (strong baroclinic flux). **(GEOMETRY — DELIVERED,
the genuinely tight part)** `transport.spherical_transport_tendency` = `(1/cosφ)∂/∂y[κ cosφ ∂θ/∂y]` (the EBM
`D·∂/∂x[(1−x²)∂T/∂x]` in β-plane `y`), **anchored on the P₂ eigenvalue** (both forms → `−6(κ/a²)P₂`, the analytic
Legendre check — NOT a self-comparison of 2 FD operators), `cos φ` metric gap **order-unity ~0.6** over the φ≈19–61°
channel → "not inherited for free"; β-plane tangent across 42° = named scope edge. **(seam) `close_loop`** routes the
emergent `D_eff` through the Phase-A bridge: CONVERGES, right sign (contrast 56→122°C steeper), degenerate
near-radiative-eq climate NOT banked. Tests: geometry **FAST** (`test_transport.py` — P₂ eigenvalue + order-unity
metric), eddy-sim **`slow`** (`test_eddy_flux.py`). No engine edit; `uses` unchanged.

**Step 3 BUILT 2026-06-11 (circulation-informed precip; `planet/circ_precip.py` + a `midlat_center_deg` param on
`precip.precip_pattern`). RUNG 1 COMPLETE.** Wires the precip **storm-track band centre** to the **emergent jet
latitude** instead of the hardcoded 50° (the §3 deep-end hook "rain where the flow puts the storm track, no moisture
physics"). `precip_pattern(lat, midlat_center_deg=MIDLAT_CENTER_DEG)` default = rung-0 **bit-for-bit BY CONSTRUCTION**
(the `two_way_pass` plumbing-reduction honesty); `circ_precip.circulation_informed_precip(state, jet)` feeds it
`jet.jet_lat`; `relocate()` returns the rung-0-vs-circ **trade**. **De-risked in 2 throwaway spikes FIRST**
(`outputs/rung1_circprecip*`, gitignored) and the spikes set the headline + KILLED an anchor — the advisor's expected
"jet_lat anchor + convergence co-location validates it" structure was OVERTURNED by the data (see below). **What is
BANKED:** (1) the **seam** (centre ← emergent circ latitude), (2) the **reduction** (jet at 50° ⇒ rung-0 exactly), (3)
the **migration MECHANISM** — band tracks a *dynamically-selected* latitude, shown via the coupler's synthetic-gradient
playbook, **anchored to `jet_lat` NOT `gradient_peak_lat`** (the gap opens only off-centre → a *flow* response not an
EBM-gradient diagnostic; this is the leg that makes it "circulation"-informed not "EBM-gradient"-informed). **The
rung-1 FINDING — a TRADE, NOT an accuracy gain (advisor reframe; don't claim better/worse):** rung-0's 50° is
**observation-calibrated**, the circ centre is the model's OWN jet at ~42° (~8° **equatorward**, production nx=96 jet_lat=42.4° — the Phase-4 channel's
known equatorward bias, it excludes the ice cliff) → relocation **trades observational calibration for internal
consistency**; so **rung-0 `precip.py` stays the DEFAULT** in the biome map/demos, circ-informed is **opt-in** (don't
regress the calibrated map). **Spike #1 (instant, no sim):** `gradient_peak_lat` PINNED ~43–46° across realistic
obliquity/CO₂/S₀ (moves only when a near-snowball cooling drags the ice cliff into the channel) → realistic migration is
**mechanism-only** (decisive only under a synthetic gradient); AND `jet_lat ≈ gradient_peak_lat` (gap ~2–3°: nx=96 jet 42.4° vs grad 45.1°) at
present-day. **Spike #2 (slow, the killed anchor):** the literal "rain where the flow **CONVERGES**" — centre on the
eddy heat-flux convergence `−∂F̄/∂y` — is **REJECTED**: the resolved convergence is **near-vacuous in the channel
interior + a window-taper EDGE DIPOLE**, NOT a physical storm-track convergence (don't build a band on a taper
artifact) → the SAME rung-3 boundary [[planet-rung1-two-way-coupler]] step 2 found (`eddy_flux` non-vacuous only under
a strong baroclinic flux). So rung 1 wires the **position** seam; the **shape/amplitude** refinement (wet-get-wetter)
stays deferred, **ITCZ/subtropics** stay prescribed (Hadley out of the midlat channel). **Scope edge:** large
equatorward displacement shallows the subtropical trough toward merging with the ITCZ (band-tracking asserted only
where the structure survives, centre ≳ 36°). Tests: fast reduction/migration/structure + **one** `slow` composition
(band follows the emergent jet on a synthetic gradient — coupler's jet-tracks-gradient proof NOT re-tested),
`test_circ_precip.py`. Full repo gate **214 passed, 1 skip**; no engine edit; `uses` unchanged. See
[[precip-parameterization-source]], [[planet-phase4-coupler]].

**Within-rung — WET-GET-WETTER, DRY-GET-DRIER BUILT 2026-06-14** (the §12.2 amplitude slice; `planet/moist.py`:
`wet_get_wetter_precip_field` + `_amplify_contrast` + `WetGetWetter`/`wet_get_wetter`; `test_moist.py` 7 fast +
1 slow guard; demo `demo_wet_get_wetter.py` → `docs/figures/planet-wet-get-wetter.png`; planet+engines gate 451
fast green). The thermodynamic contrast-sharpening rung-0's uniform `CC(T̄)` omits. **HOME CORRECTED by advisor:
`moist.py` NOT the §12-stated `circ_precip.py`** — §12 pointer was STALE (predated the rung-2 build); this is
literally the GENERALIZATION of the one-line `energy_constrained_precip_field` (scales mean+anomaly *together* at
one rate → this SPLITS them: `P=⟨P⟩·M(T̄)+(pattern−⟨P⟩)·W(T̄)`, mean `M`=energy ~2.5%/K, anomaly `W`=C–C ~7%/K).
Held & Soden 2006 "rich-get-richer" on the precip pattern: warming intensifies the ITCZ/storm-track while DRYING
the deserts (rung-0 uniform `CC` wrongly wettens them too). **Honesty (advisor) = better PRESCRIBED param, NOT
derived:** the split DIRECTION + the two rates are calibrated/cited; STRUCTURALLY EXACT = the mean-zero anomaly
split (`⟨pattern−mean⟩=0` machine-precision on the equal-area `x=sinφ` grid ⟹ ⟨P⟩ scales at energy rate =
PLUMBING by-construction, not a finding) + the reduction to BOTH fields when rates coincide (`M=W=CC`→rung-0;
`M=W=energy`→energy field). Opt-in/default-off (rung-0 `precip.py` untouched); **deliberately NOT fused** with the
storm-track position seam (moist.py's non-composition rule — trade×trade). Numbers (ΔT=6): ITCZ +92, desert@25°
−13 (vs rung-0 uniform +14). Edges: global-`T̄` (local-`q_sat(T(φ))`=richer upgrade), `P≥0` floor (deep warming→
total aridity), thermodynamic-only (dynamic circ-driven amplification = rung 3+). See [[precip-parameterization-source]].

**NEXT — rung 2 (moist dynamics):** now **SCOPED (not built), 2026-06-11** → see [[planet-rung2-scoped]].
