---
name: planet-rung3-scoped
description: "Planet rung 3 (vertical structure → baroclinic instability) — Phase A BUILT 2026-06-12 (linear growth rate). N-layer SIBLING engine LayeredShallowWater + TwoLayerStability anchor; single-layer untouched → nl=1 reduction BYTE-IDENTICAL. Periodicity resolved (advisor): basic state enters as constant BACKGROUND COEFFICIENTS, fields are perturbations, full nonlinear engine default-off background. Engine σ matches the 6×6 dispersion matrix within ~4% converging. Phase B (saturated irreversible flux = the headline payoff) STILL the open bet, needs unbuilt hyperviscosity."
metadata:
  node_type: memory
  type: project
  originSessionId: rung3-session
---

Planet **rung 3** (vertical structure → **baroclinic instability** = real storms, the §5 staircase's
biggest jump). **PHASE A BUILT 2026-06-12** (the linear growth rate); the scoping below remains the plan of
record for Phase B. Full plan `docs/plans/planet-earth-system.md` §10. Continues
[[planet-rung1-two-way-coupler]] (which named the payoff: the eddy-flux→EBM reduction is "non-vacuous only
at rung 3"). Advisor-pressure-tested (3 load-bearing refinements) + spike (`outputs/rung3_baroclinic_spike.py`,
gitignored) that answered all three risks BEFORE any build.

**PHASE A BUILT 2026-06-12 — `engines/fluid/layered.py` (`LayeredShallowWater` + `LayeredState` +
`ThermalWindBackground`) + `engines/fluid/stability.py` (`TwoLayerStability`), 22 tests
(`test_layered.py`/`test_stability.py`).** Headline: a full nonlinear N-layer vector-invariant C-grid
SSP-RK3 engine whose unstable two-layer mode **grows at the analytic linear rate within ~4%, converging
with resolution** (5.4→4.2→3.9% at nx 32→64→128, matching the spike's 5.2→4.1→3.8%). **THE ARCHITECTURE
FORK (advisor): a SIBLING engine, `ShallowWater`/`SWState` left UNTOUCHED — NOT an in-place restructure of
`SWState`.** Reason: the single-layer tests are all loose (`rel=2e-3`), so a shared-kernel refactor that
shifted the trajectory a few ULP would pass SILENTLY → no byte-exact guard; untouched = zero-risk AND makes
the `nl=1` reduction a *meaningful cross-engine* check, not a near-vacuous same-kernel one. Duplicate the
~25-line C-grid stencil vectorized over a **leading layer axis** (`np.roll` on `axis=-1/-2`); the **only**
inter-layer coupling is the **Montgomery pressure stack** (`M₀=g·η_top`, `M_k=M_{k-1}+g'_k·η_{k-1}`). The
**single-layer reduction is BYTE-IDENTICAL** (`np.array_equal`, max|Δh|=0.0 over 60 steps incl. topography)
— `M₁=g(h_b+h)` is IEEE-bit-commutative with the single-layer Bernoulli pressure, a length-1 leading axis is
the same float ops.

**THE PERIODICITY RESOLUTION (advisor, the one load-bearing blocker — was: nonlinear engine vs the
thermal-wind basic state that breaks y-periodicity).** A doubly-periodic domain **categorically cannot carry
a uniform meridional gradient as a field** (it seams at the wrap — physics, not an engine limit). So the
unstable basic state `(U_k, G_k)` enters as **constant BACKGROUND COEFFICIENTS** (Doppler `−U_k ∂_x` +
baroclinic `−G_k v'` in continuity), with the **prognostic fields being the PERTURBATIONS** — exactly the
spike's `TwoLayerLinear`. **NOT "linearized instead of nonlinear":** the engine is the FULL nonlinear
vector-invariant solver with an **optional, default-off** `background`; `background=None` → plain nonlinear
→ the bit-for-bit reduction; `background` set → the spike's linear terms PLUS the engine's O(ε²) eddy-eddy
nonlinearity (grows at the linear rate at ε~1e-3 = Phase A; saturates at large ε = Phase B). This "fixed
mean shear + prognostic eddies on a doubly-periodic plane" is *the* standard two-layer geostrophic-turbulence
setup (Held–Larichev/Salmon) = **the exact substrate Phase B needs** (the background feature is in-scope, not
a Phase-A hack). `thermal_wind(U)` derives `G` so the engine injects exactly the `G_k` the analytic rate
assumes (pinned: `bg.G == TwoLayerStability.basic_state_gradients`). f-PLANE only (β=0, matching the f-plane
operator; finite critical shear = the named β-capable extension).

**THE ANCHOR — `TwoLayerStability` (ported from the spike, first-principles validated as tests).** Roots the
**6×6 dispersion matrix** built from the linearized eqs. Self-validation pinned: zero-shear **neutral to
machine precision** (`max|Im ω|≈3e-20`); **recovers BOTH two-layer Poincaré dispersions** (external
`√(gH_tot)`, internal `√(g'H_e)`) — but only to **~0.2% (free-surface coupling, NOT exact** → assert
`rel=1e-2`; the spike *printed* these side-by-side, never asserted tight equality; the decisive exact check
is the neutrality); **short-wave cutoff** (`K²=2F`); **Eady coefficient 0.304 vs 0.310** — *load-bearing
gotcha: use the LAYER-depth radius `L_d=√(g'H₁)/f₀`, NOT the equivalent-depth `√(g'H_e)/f₀`* (the latter gives
0.215 — `Ld_int` on the class is the gravity-wave radius, a DIFFERENT quantity); **f-plane no-critical-shear**
(σ/U_s ≈ const). The ENGINE also reproduces both two-layer Poincaré dispersions to 0.06%/0.01% (the tight
Montgomery-coupling check) — these + the growth test are `@pytest.mark.slow` (the external-mode CFL = the
first compute wall; ~3–13s each).

**Two gotchas BANKED (advisor):** (1) conservation tests **only** on the `background=None` path — a background
extracts APE so perturbation energy *grows* (that IS the signal), and per-layer mass `∫h_k` is machine-exact
only there (the `−G_k v'` source breaks it when `∫v≠0`; the `l=0` mode has `∫v=0`). (2) **no passive tracer**
in the layered engine — the interface displacement IS the buoyancy now (the leap past the rung-1 *passive*
tracer). **Gate: full fast lane 295 passed; the standalone repo has NO `tools/`/`GATES` import-drift guard
(it stayed in the BigSim monorepo) → the "full-repo gate" here = the whole pytest suite.**

**PHASE B IS STILL THE OPEN BET (unchanged by Phase A — exactly the advisor's done-check hedge).** Phase A
de-risked BUILDABILITY (textbook-guaranteed) at high confidence; it moves the *quantitative payoff* risk
**zero**. The saturated, irreversible, down-gradient baroclinic eddy heat flux `⟨v'·interface'⟩` at a
**realistic magnitude → the rung-1 reduction-to-EBM finally non-vacuous** is untouched, and is exactly the
claim class downgraded@rung1 / overturned@rung2. Needs the **named-but-unbuilt hyperviscosity** (a turbulent
layered run cascades enstrophy to the grid — the inherited no-limiter behaviour, now load-bearing), long
post-saturation runs, + the external-mode cost over many cycles. The Phase-A RHS is structured to leave room
for that dissipation operator.

**The model (the fork, settled): a two-layer free-surface shallow-water model (Phillips 1954).** A single
layer **categorically cannot** be baroclinically unstable (no available potential energy, no vertical
shear). Two stacked SW layers of slightly different density, **coupled only through pressure**
(`P₁=g(h₁+h₂)`, `P₂=g(h₁+h₂)+g'h₂`). The table's "multi-level / 3-D" is realized **minimally as 2 layers**;
**N-layer = the within-rung upgrade**, not a new rung (no cleaner anchor above 2, just cost). **The
interface displacement IS the dynamical temperature/buoyancy** → heat transport becomes **intrinsic to the
dynamics** (the leap past rung-1's *passive* tracer; the CONTRACT already names "active buoyancy tracer
feeding h = a two-layer model"). It is the **first *structural* edit to `engines/fluid`** since Phase 3 (the
rung-1 tracer was an additive edit; this restructures `SWState` to a leading layer axis) → triggers the
**full-repo gate + import-drift guard**. And the **first compute wall**.

**THE ANCHOR — the two-layer SW equations' OWN linear stability, NOT QG Phillips (advisor's #2, load-
bearing).** QG Phillips as the primary anchor would import an approximation gap (ageostrophy, free-surface
coupling) into the *tight* leg. The rung-appropriate move (matching Phase 3's Poincaré/Rossby anchor + rung
2.5's exact `dq_sat/dT`): **linearize the two-layer SW eqs → root the 6×6 dispersion matrix** for the
growing mode (built from the equations, **not a recalled quartic**). σ = max Im(ω). **First-principles
validated in the spike:** zero-shear → **neutral to machine precision** (`max|Im ω|≈3e-20`) AND **recovers
both Poincaré dispersions** (external `ω²=f₀²+gH_tot k²`, internal `ω²=f₀²+g'H_e k²`). Gives σ(k) with a
**short-wave cutoff**, most-unstable `λ*≈6.8×` the layer deformation radius (`≈680 km` idealized). **The
baroclinic `G_k` terms are EXTERNALLY anchored** (zero-shear Poincaré leaves them untested → else they rest
only on the self-consistent Phillips convergence = two hand-derivations, advisor's circularity catch): max-
growth coefficient `σ_max≈0.304·U_s/L_d` matches the **EADY** model (independent derivation) `0.310` to
**2%**, and the critical-shear *formula* `U_s,crit=β·g'H/f₀²` matches the literature `β/k²_int`
(Pedlosky/Vallis, web-confirmed). **The solver is f-PLANE** (β not in the perturbation operator → correctly
**Eady-like: NO critical shear, unstable for all shear**); a *finite* critical shear needs β = a Phase-A
(β-capable QG/PV-gradient) item, NOT a spike claim (the advisor caught my overclaim here).

**QG Phillips = the LOOSE cross-check that CONFIRMS the tight leg.** Derived 2×2; in the **rigid-lid limit
(`g→∞`, external mode infinitely fast) σ_SW/σ_Phillips → 1.004 (<0.5%)** = mutual validation; at free-
surface params the gap is a clean **~4%** (SW slightly *more* unstable, the extra free-surface DOF). **The
spike caught + fixed exactly the recalled-coefficient trap the advisor flagged**: β=0 cutoff is `K²=2F`,
`F=f₀²/(g'H)` — **not `2F`** (√2 too short; that mismatch is what first surfaced the gap, then closed it).

**THE EXTERNAL-MODE CFL = the real cost (advisor's #1) — quantified + AFFORDABLE.** Free-surface two-layer
carries a fast **barotropic** gravity wave `√(gH_tot)` that sets the explicit RK3 step while the **slow
baroclinic mode** is all that matters → penalty `c_ext/c_int≈14×` at idealized params, but **4 e-folds =
~3–16 s wall** (nx 32→128) → **GO on the free-surface engine extension**. The **rigid-lid fork** (a
barotropic **elliptic solve** = a *structural* change to the explicit engine, "a different animal") is the
**named within-rung upgrade** if Phase-B's saturated, higher-res, many-life-cycle runs make the penalty
bite. **Idealized `(g,H,Δρ/ρ)`** chosen for modest `√(gH)` + resolvable internal `L_d` — **honest at rung
3** (validates *mechanism + growth rate*, not Earth jet speeds = the same config-tuned honesty banked at
rungs 1–2). **The route is sound:** a linearized two-layer C-grid SSP-RK3 solver (mirrors `engines/fluid`)
reproduces the analytic σ **within ~4%, monotonically converging** (5.2→4.1→3.8% at nx 32→64→128).

**CALIBRATE THE GO (advisor done-check — the hedge is load-bearing).** The spike de-risks **buildability**
(linear instability exists, engine reproduces it, CFL affordable) at HIGH confidence — but that's the
**textbook-guaranteed** part. **The Phase-B QUANTITATIVE PAYOFF is an OPEN BET** (irreversible flux at
realistic magnitude → non-vacuous reduction-to-EBM): untouched by a linear spike, and *exactly* the claim
class downgraded at rung 1 (named-not-banked) + overturned at rung 2. Phase A succeeding moves that risk
**ZERO**. "Spike-validated/GO" = *the route is sound to build*, NOT *the rung is de-risked* — the saturation
reality-check is owed before banking the payoff.

**THE A/B SPLIT (advisor's #3 — same machinery-vs-emergent split as rungs 1 & 2):**
- **Phase A = the LINEAR growth rate (tight; the current energy-conserving engine).** Extend `engines/fluid`
  to **N layers** (the CONTRACT's promised leading-axis seam): layers couple *only* through the
  Montgomery/Bernoulli pressure term; vector-invariant momentum generalizes per-layer; **single-layer
  `tracer=None` stays bit-for-bit** (the by-construction reduction held every rung). Anchor: perturbation on
  a thermal-wind-balanced basic state grows at the analytic two-layer SW rate within a few %.
- **Phase B = the NONLINEAR saturated eddy field → THE HEADLINE PAYOFF (loose magnitude; needs dissipation).**
  Rung 1 found the eddy-flux→EBM "tight reduction" **near-vacuous** (barotropic ~1000× too weak, ~90%
  reversible) and named it "non-vacuous **only at rung 3** (strong baroclinic flux)." Phase B integrates the
  unstable two-layer flow to saturation, measures the **emergent irreversible down-gradient baroclinic eddy
  heat flux** `⟨v'·(interface)'⟩` (now genuine), and re-runs the rung-1 reduction → expecting a **non-vacuous
  `D_eff` near a realistic magnitude**. Needs the **named-but-unbuilt hyperviscosity** (the engine has no
  limiter; turbulent runs cascade enstrophy to the grid — a CONTRACT non-goal now forced), a long post-
  saturation run, + the external-mode cost over many cycles. Phase A's extension must **leave room for the
  Phase-B dissipation operator**.

**Held–Suarez is DEFERRED, NOT the rung-3 anchor** — it's a **3-D sphere primitive-equation** benchmark
(Newtonian relaxation + Rayleigh friction → statistically-steady storm track) = **rung 5**; naming it here
over-reaches a 2-layer β-plane laptop model. Rung-3 anchor = **Phillips/SW growth rate (tight) + emergent
flux (Phase B)**. **Named deferrals/walls:** hyperviscosity (Phase B prereq); rigid-lid/external-mode fork;
**rigid channel walls in y** (the classic lifecycle geometry — the spike used the doubly-periodic `l=0` mode
to sidestep them; the saturated field may want a channel = the named BC extension); the sphere (rung 5).
**Sources to pin at build** (extending [[shallow-water-source]]): Phillips 1954 / Eady 1949 / Charney 1947;
two-layer SW formulation → Vallis 2017 *AOFD* / Cushman-Roisin & Beckers; Held–Suarez 1994 (named-deferred).
See [[planet-rung1-two-way-coupler]], [[planet-rung25-mse-diffusion]], [[planet-plan]].
