---
name: planet-rung3-qg-built
description: "Planet rung 3 Phase B BUILT 2026-06-13 (planet/baroclinic_qg.py, two-layer QG turbulence) and THE OPEN BET IS WON: the saturated baroclinic eddy thickness flux is down-gradient, irreversible (irr~1 vs rung-1 ~0.1), at order-unity dimensionless mixing efficiency kappa/(v'L_d)~0.7-1.3 (vs rung-1 ~1e-3) -> the rung-1 reduction-to-EBM is finally NON-VACUOUS. Advisor's load-bearing catch: down-gradient+irr~1 are GUARANTEED for any sustained baroclinic state; the bet is won only by SHOWING turbulence (inverse-cascade condensate), which the weak-drag case does. Banked DIMENSIONLESS+qualitative; dimensional kappa landing on Earth is coincidental."
metadata:
  node_type: memory
  type: project
  originSessionId: rung3-qg-session
---

Planet **rung 3 Phase B** (the saturated baroclinic eddy heat flux = the headline payoff) is **BUILT
2026-06-13 and THE OPEN BET IS WON**. The re-routed engine (`planet/baroclinic_qg.py`, ~430 lines)
produces a saturated, irreversible, **down-gradient** baroclinic eddy thickness flux at an **order-unity
dimensionless mixing efficiency** → the rung-1 reduction-to-diffusive-EBM is **finally non-vacuous** (the
claim class **downgraded@rung1** = barotropic flux ~1000× too weak + ~90 % reversible, **overturned@rung2**;
*won* at rung 3 exactly as the §5 staircase predicted). Continues [[planet-rung3-phaseB-outcropping]] (the
spike that re-routed free-surface SW → QG) and [[planet-rung3-scoped]] (Phase A, banked).

**THE ENGINE.** Doubly-periodic two-layer QG (Held & Larichev 1996): the PV anomaly `q_k` (a `QGState`,
leading layer axis) advected by its own flow; **ψ recovered by a 2×2 spectral PV inversion**
(`det=K²(K²+F₁+F₂)`, the `K=0` domain-mean gauge → `ψ=0`); a **pseudospectral 2/3-dealiased Jacobian** for
the eddy–eddy nonlinearity; **SSP-RK3** (mirrors `engines/fluid`); the mean shear `(U_k, ∂q̄_k/∂y)` as
**background coefficients** (eddies = the prognostic fields — the SW engine's `background` design); **β-plane**;
the spike's **hyperviscosity `−ν₄∇⁴q` + lower-layer bottom Ekman drag `−r∇²ψ_2`**, all **default-off**.
**Home scope-decision SETTLED: the lean single-consumer planet module, NOT `engines/spectral/`** (rule-of-three
unmet). **Why QG won where free-surface SW outcropped:** linearized in the thickness → **no layer-depth floor
(never outcrops**, pinned by a test); rigid lid **filters the external gravity wave** (no external-mode CFL →
cheap long runs). Both stiffnesses the SW engine paid are gone.

**THE TIGHT LEG (linear anchor FIRST — advisor's day-one sequencing).** Rooted 2×2 QG dispersion = the
**analytic Phillips closed form to 2e-15** (equal layers — *the same equations*; **advisor's load-bearing fix:
the bar is ~1e-9 NOT "a few %"** — the ~4 % I'd carried forward was the free-surface-SW-vs-Phillips gap, which
does NOT exist here; a loose tolerance would hide a partially-compensated PV-gradient sign bug). Plus `K²=2F`
short-wave cutoff, **zero-shear neutral** (machine ε), and **β's re-entry** (the thing the f-plane SW solver
*could not* test): a finite **Charney–Stern critical shear `U_crit=β/F`** (lower-layer `∂q̄₂/∂y=β−F·U_s`
reverses sign; sub-critical neutral, super-critical growing). **Cross-model bridge to Phase A:** the SW 6×6
`TwoLayerStability` in the rigid-lid limit `g→∞` → **σ_SW/σ_QG → 1 to <0.5 %** — *asserted* (the only tie
between the two models; no bit-for-bit reduction across the model boundary). The full **nonlinear** engine
reduces to this linear operator: a single growing eigenmode grows at σ to **0.1 %**.

**THE WIN — and what actually decides it (advisor, the load-bearing reframe).** Drag sweep `r∈{0.5,1,2}σ`:
`κ>0` (down-gradient), `irr=0.96–1.00` (vs rung-1 ~0.1), **`κ/(v'_rms·L_d)=0.71–1.27`** (vs rung-1 ~1e-3),
robust across drag. **BUT down-gradient + `irr≈1` are GUARANTEED for ANY sustained baroclinic state** (the
flux *is* the APE→EKE conversion powering the eddies; sign-pinned + spatially averaged ⟹ `irr≈1`
automatically) — **necessary NOT sufficient**, and they cannot distinguish developed turbulence from a
quasi-steady wave (exactly the **P2 "irr~1 trivial *without sloshing*"** rejection carried from the spike). So
the bet is won **only by SHOWING genuine turbulence**, which the **weak-drag (`r=0.5σ`) condensate candidate**
(`v'_rms≈16 ≫ U_s=4`) does on three independent diagnostics: (1) **EKE(t) irregular + drifting**
(`std/mean≈0.25`, not flat); (2) an **isotropic KE spectrum = inverse-cascade condensate** — 84 % of the
energy *below* the injection band, the peak migrated to the **box scale `0.33k*`**, broadband-continuous with a
clean dissipation tail (NOT spikes at `k*,2k*,3k*`); (3) a **PV snapshot of coherent vortices + rolled-up
filaments across scales** (not a wave train). The **dimensionless `κ/(v'L_d)~O(1)`** is the discriminating
number; the inverse cascade is what makes it turbulence.

**HONEST EDGES (advisor — do NOT overclaim).** (1) Banked claim is **DIMENSIONLESS + qualitative**; the
**dimensional `κ≈0.7–4×10⁶ m²/s`** lands in Earth's band (1–5e6) but that is **COINCIDENTAL + box/drag-dependent**
(5× across the sweep) — *not* an Earth-κ reproduction. (2) The condensate **dominates `v'_rms`** → read
`κ/(v'L_d)~1` as "mixing length ~ L_d-ish," not a precise efficiency. (3) **`κ₁=κ₂` is an estimator IDENTITY**
(`F₁−F₂=⟨½∂_x τ²⟩=0` periodic), *not* an independent cross-layer check (it was for the SW spike's separate
per-layer thicknesses) — so I dropped it as a validation leg. (4) **New model outside `engines/fluid`**
(pseudospectral) → Phase A and B validate *different* models, bridged only by the <0.5 % rigid-lid linear
cross-check. (5) **Homogeneous box → a domain-bulk κ**, not `κ(y)`; a meridional channel (the operator-shape
test) is the named BC extension. (6) **Resolution nx=96 ≈ 3.5 pts/L_d** marginal at the deformation scale
(condensate lives at well-resolved large scales; clean dissipation tail = no grid-scale pileup; the
dimensionless ratio is already robust across the drag sweep); **nx=128 firm-up (4.7 pts/L_d) HOLDS the
dimensionless ratio** (`κ/(v'L_d)=1.10`, `irr=1.00`) while the dimensional κ shifts ~45 % with res
(`1.44→2.10e6`) — exactly why only the dimensionless ratio is banked, never the magnitude.

**Tests:** `planet/tests/test_baroclinic_qg.py` — tight (Phillips 2e-15, cutoff, `U_crit=β/F`, zero-shear
neutral, SW↔QG rigid-lid <0.5 %) + plumbing (`q↔ψ` machine-exact, zero-shear decay, CFL guard) fast; the
nonlinear→linear reduction (0.1 %) + the saturated-flux discriminators (κ>0, irr>0.8, 0.1<ratio<10,
v'>U_s, never outcrops) slow. **Full fast-lane gate green** (~349 tests). The turbulence-vs-wave
characterization (spectrum + PV snapshot) lives in the gitignored spike `outputs/rung3_qg_turbulence_check.py`.
Sources: **Held & Larichev 1996** (two-layer QG turbulence + diffusivity scaling), Phillips 1954 / Eady 1949,
Vallis 2017 *AOFD*; extends [[shallow-water-source]]. Plan §10. **Rung 3 (vertical structure → baroclinic
instability) is now COMPLETE: Phase A (linear, banked) + Phase B (saturated flux, bet won).**
