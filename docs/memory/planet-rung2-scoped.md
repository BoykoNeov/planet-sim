---
name: planet-rung2-scoped
description: Planet rung 2 (moist dynamics → emergent precip) is SCOPED not built; fork settled (column moist budget, NOT fluid-channel — spike-confirmed Phase B = rung 3); core = moist-EBM diagnostic, headline unlock = energy-constrained ~2-3%/K rate
metadata:
  type: project
---

Planet **rung 2** (moist dynamics → emergent precipitation, the §5 staircase) is **SCOPED — design +
fork settled, NOT built** (2026-06-11; advisor-pressure-tested twice). Full record is in the plan
(`docs/plans/planet-earth-system.md` §10, the "Rung 2 — SCOPED" build-record paragraph); this is the
recall pointer. Continues [[planet-rung1-two-way-coupler]] (rung 1 COMPLETE).

**The fork — *where moisture lives* — is settled empirically: a COLUMN MOIST BUDGET, not fluid-channel
moisture transport.** A throwaway spike (`outputs/rung2_moisture_convergence_spike.py`, gitignored)
advected a steep C–C-shaped **moisture** tracer on the *same* released barotropically-unstable Phase-4
jet `eddy_flux` uses; its flux **convergence** `−∂⟨v'q'⟩/∂y` is **near-vacuous in the interior +
window-edge-artifact-dominated** (interior/edge RMS **0.32**, *worse* than the temperature tracer's 0.50),
and the mean balanced jet is near-nondivergent (eddy/mean 11–23×, `∇·(hu)=−∂h/∂t≈0`). So a *resolved
storm-track precip PATTERN* is **rung 3** (needs the vertical = ascent→condensation), **not rung 2** —
corroborates rung-1 step-3's `reduction_to_ebm_operator` finding. One realization suffices ("don't
gold-plate the spike", advisor).

**The bankable rung-2 core = a moist-EBM DIAGNOSTIC** (reuses the diffusion spine a **4th** time; does
**not** perturb the validated Phase-1 climate). `q(φ)=RH·q_sat(T)` (fixed RH over rung-0 `T`; `q_sat` from
C–C) + down-gradient latent transport whose diffusivity is **tied to rung-1's eddy `κ` via the existing
`transport.py` κ→D bridge** (no new `D_q`); `P = E − ∇·(moisture transport)`, a **pure diagnostic** (does
NOT enter the T-equation → Phase-1 triad stays green).

- **HEADLINE UNLOCK = the RATE, not the pattern:** replace `precip.py`'s prescribed C–C **7 %/K** global
  amplitude with the **energy-constrained ~2–3 %/K** rate — closes `precip.py`'s *own* scope-edge #3.
- **The first sub-grid closure (the rung-2 wall) = the atmospheric-energy closure** `R_atm(T)` with
  `L⟨P⟩≈R_atm−SH` — NOT the benign fixed-RH. ("Energy-limited evaporation" overstates it: a *constraint*
  on `⟨P⟩`, not a closed surface budget.)
- **Emergent pattern = a TRADE, named not a win:** extratropical `P−E` good, but the deep-tropical **ITCZ
  comes out backwards** (diffusion exports moisture from the equator; the real ITCZ is up-gradient Hadley
  = mean circulation, deferred). Rung-0 `precip.py` stays the default.
- **Triad re-classed:** *tight* = `q_sat` exact C–C function; *real-but-loose (the unlock)* = the 2–3 %/K
  rate; *plumbing* = `∫E=∫P` + the `L·q→0` reduction (by-construction, the `two_way_pass` honesty class).
- **Deferred:** the full moist EBM (diffuse MSE `m=c_pT+L·q` so `T` responds) = **rung 2.5** — it re-opens
  Phase-1 `(A,B,D)` calibration (rung-0 `D=0.555` already absorbs latent transport — MSE diffusion
  double-counts).
- **Sources to pin at build** (named now, NOT carried from memory): rate → Held & Soden 2006 / Allen &
  Ingram 2002 (extends [[precip-parameterization-source]]); diffusive-moist-EBM → Flannery 1984 / Hwang &
  Frierson 2010 / Siler–Roe–Armour 2018.

**Next action = build Phase A** (the moist-EBM diagnostic), spike-first per discipline; user said "before
any code" for the scoping pass, so building is the follow-on go.
