---
name: planet-rung2-scoped
description: "Planet rung 2 (moist dynamics → emergent precip) — Phase A (the column moist-EBM diagnostic) BUILT 2026-06-11 (planet/moist.py); headline = energy-constrained ~2.5%/K rate (opt-in), trade = extratropical-ONLY P−E budget (the spike OVERTURNED the scoped \"subtropics good\" claim); rung 2.5 / rung 3 deferred"
metadata: 
  node_type: memory
  type: project
  originSessionId: 0369b58c-6c5c-4ddf-a80e-44cc29351f78
---

Planet **rung 2** (moist dynamics → emergent precipitation, the §5 staircase). Scoped 2026-06-11
then **Phase A BUILT same day** (`planet/moist.py` + `tests/test_moist.py`, 15 fast tests; full plan
record `docs/plans/planet-earth-system.md` §10 "Rung 2, Phase A — BUILT"). Continues
[[planet-rung1-two-way-coupler]] (rung 1 COMPLETE). Built **spike-first**
(`outputs/rung2_phaseA_moistebm_spike.py`, gitignored) + **advisor-pressure-tested, which reshaped the
deliverable and OVERTURNED one scoped claim** — record that, it's the load-bearing correction.

**The fork (settled at scoping, still holds): a COLUMN moist budget, NOT fluid-channel moisture.** A
throwaway spike (`outputs/rung2_moisture_convergence_spike.py`) showed `−∂⟨v'q'⟩/∂y` on the released
unstable jet is near-vacuous interior + window-edge artifact (interior/edge **0.32**, worse than θ's
0.50) → a *resolved* storm-track precip PATTERN is **rung 3** (the vertical), not rung 2.

**Phase A = the moist-EBM diagnostic, reuses the diffusion spine a 4th time, does NOT touch the T-eq
(Phase-1 stays green — asserted).** The advisor **split it in two** (don't fuse — a full emergent `P`
field forces an unphysical `E`):
- **(1) the RATE = the headline unlock, robust, OPT-IN.** `energy_constrained_factor(T̄)` replaces
  `precip.py`'s C–C **7%/K** with the **energy-constrained ~2.5%/K** (`L⟨P⟩=R_atm−SH`, `⟨P⟩₀≈100 cm/yr`)
  → closes `precip.py` scope-edge #3. Energy budget is **LINEAR** in T̄ → factor is **linear** (NOT a
  smaller C–C exponent — honest functional difference), `=1` at present ref (map unchanged bit-for-bit),
  floored ≥0. `precip.py` **untouched** (its 7%/K warmed tests stay green); rung-0 stays default.
- **(2) the emergent budget `P−E` = the trade, a DIAGNOSTIC not the default.** `P−E=(D/c_p)·∂/∂x[(1−x²)
  ∂q/∂x]`, `q=RH·q_sat(T)`. **The latent heat L CANCELS** (same eddies stir heat+moisture → latent
  transport = EBM operator on `L·q` with the **same D** via the rung-1 κ→D bridge; mass flux = energy
  conv /L → L drops) ⟹ **NO new `D_q`**. Conservative face-flux form → `∫(P−E)=0` **machine-exact** under
  the area-MEAN rule (NOT `np.trapezoid` — breaks telescoping; that was a real bug I hit). **Pure
  diagnostic.** **Why P−E not full-P (advisor catch):** eq export ~−2.4 m/yr vs `⟨P⟩~1` m/yr → no honest
  zonal `E` keeps `P≥0` (uniform→P(eq)<0; `E∝q_sat`→absurd ~6 m/yr eq evaporation) → report P−E, skip
  full-P.

**THE OVERTURNED CLAIM (the spike corrected memory — do not regress to the old scope):** the budget is
**extratropical-ONLY**, NOT the scoped "extratropics good + subtropical E>P". Deep equator is
**backwards** (diffusion EXPORTS from the moist equator; real ITCZ = up-gradient Hadley, deferred) **AND
the subtropical evaporative belt is NOT reproduced** — the steep eq–pole contrast hyper-peaks C–C `q` at
the equator, pushing the moisture-flux max equatorward (~20° zero-crossing) so subtropics come out
**`P>E`** (production nx=180: eq −267, subtrop 25–35° +83, midlat +101, polar +28 cm/yr). Only the
**extratropics (poleward ~40°) are right**. So the benchmark test asserts ONLY eq-export +
extratropical-convergence, and **one test pins the subtropical mislocation** as the honest limitation
(guards a silent "fix"). Same "trade not a win" the staircase keeps banking.

**The sub-grid WALL = the prescribed `R_atm` slope** (`R_ATM_SLOPE=2 W/m²/K`, cited **Held & Soden 2006 /
Allen & Ingram 2002**; the rate is a *cited-closure* result, **NOT derived**, and explicitly **not
`B_OLR`** — a different 2 W/m²/K quantity, the advisor's named trap). **Triad:** *tight* = `q_sat` exact
C–C (textbook 3.8/14.7 g/kg + ~7%/K log-slope) + the operator reproduces the **P₂ eigenvalue −6** (it
IS the EBM transport operator, the `transport.py` anchor); *real-but-loose (unlock)* = ~2–3%/K rate,
slower than C–C, doubling with the closure slope; *plumbing* = `∫(P−E)=0` + `q→0` reduction +
unity-at-ref; *benchmark (loose)* = the extratropical-only trade. Sources pinned (the `[[…-source]]`
discipline): C–C `q_sat` → Hartmann GPC / Bohren & Albrecht; rate → Held & Soden / Allen & Ingram;
diffusive moist EBM → Flannery 1984 / Hwang & Frierson 2010 / Siler–Roe–Armour 2018. See
[[moist-ebm-source]], [[precip-parameterization-source]].

**DEFERRED (named):** ITCZ/Hadley (mean circulation); resolved storm-track precip *pattern* = **rung 3**
(the vertical, spike-confirmed); the full **MSE-diffusing moist EBM where `T` responds** (emergent polar
amplification) = **rung 2.5** — re-opens Phase-1 `(A,B,D)` calibration (`D=0.555` already absorbs latent
transport). No engine edit; `uses` unchanged; full planet gate **179 passed, 1 skip**.
