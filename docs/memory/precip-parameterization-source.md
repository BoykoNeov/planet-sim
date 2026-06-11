---
name: precip-parameterization-source
description: Planet P2 — cited zonal-mean precip-by-latitude structure + the Clausius-Clapeyron 7%/K (moisture) vs 2-3%/K (energy-constrained global) rates precip.py pins
metadata: 
  node_type: memory
  type: reference
  originSessionId: 65378612-4f52-44d4-8842-fe8158fc4aed
---

**Pinned source for Planet Phase 2's diagnostic precipitation parameterization** (`projects/planet/precip.py`),
the `[[…-source]]` discipline. Part of [[bigsim-program]]; see [[planet-plan]].

**The cited facts.**
- **Zonal-mean precip-by-latitude structure** (the *pattern*, set by the general circulation —
  Hadley/Ferrel cells): a wet equatorial **ITCZ**, dry **subtropics ~25–30°** (the descending Hadley
  branch = the world's great deserts), wet **midlatitude ~50°** storm tracks, dry **poles**. This is
  the **known observed structure**, named as a parameterization — NOT derived (no moisture variable →
  a real water cycle is moist thermo = the GCM rung-2 deferral). Units **cm/yr** (to match the
  Whittaker classifier [[whittaker-biome-source]] — the units trap was mm/yr).
- **Clausius–Clapeyron ~7 %/K** — saturation water-vapour content rises exponentially ≈ 7 %/K (the
  *moisture-capacity* rate). pinned as `CC_RATE_PER_K = 0.07`.
- **THE named gap (advisor honesty flag):** global-*mean* precipitation is **energy-constrained** to a
  slower **≈ 2–3 %/K** (the atmosphere can only radiate away so much latent heating), NOT 7 %/K. v1
  scales at the C–C 7 %/K as a moisture proxy and **does not claim** the energy-constrained global
  rate — that closure is the rung-2 moist energetics, named not modelled.

**The design (advisor-reviewed 2026-06-09).** `P(φ, T̄) = pattern(φ) · CC(T̄)`: a Gaussian-band
pattern (ITCZ amp 215 cm @ σ12°, midlat amp 75 cm @ 50°/σ15°, baseline 20 cm; subtropical dry trough
*emerges* between the bands) × a **global-mean-T̄** C–C amplitude `exp(0.07·(T̄−15))`, ref 15 °C
(≈ present global mean, CC=1 there).

**Deliberate refinement of advisor's "scale by local T":** scaling the present pattern by *local* T
against a fixed ref **over-amplifies** the warm equator (×2.6 at 28 °C) and breaks calibration — the
equatorial wet belt is the ITCZ (circulation), not local warmth. So **pattern held fixed (circulation),
only the global moisture amplitude responds to T̄** — the honest pattern/amplitude decomposition.
(Advisor later OWNED "local T" was imprecise — global-T̄ is correct; proper local-T = local-*anomaly*
scaling, of which global-T̄ is the simpler member.)

**TWO distinct deferred spatial responses (both named scope edges):** (1) band **migration** (ITCZ
shift / Hadley widening — centres move) = rung-1/2; (2) **pattern amplification** = the observed
**wet-get-wetter, dry-get-drier** contrast sharpening — but v1's *uniform* CC(T̄) multiplies the whole
pattern by one factor, so the **dry subtropics also get wetter** under warming (the opposite). Also
rung-1/2 (needs circulation-set moisture convergence). v1 intensifies bands uniformly in place.

*Validated tight:* the band structure (ITCZ-wet/subtropics-dry/midlat-wet/poles-dry — the subtropical
min is a *local* trough between the ITCZ and storm track, NOT the global floor which is the pole) +
**monotone** ⟨P⟩(T̄) (a *consistency* check, not a conservation law — prescribed field has no water
budget). *Loose:* band amplitudes/centres + the 7 %/K rate.

**Rung-2 update (2026-06-11):** the named **2–3 %/K energy-constrained** rate is now BUILT as an opt-in
in `planet/moist.py` (`energy_constrained_factor`, a **linear** factor) — closing this module's
scope-edge #3. `precip.py` itself is **unchanged** (the C–C 7 %/K stays the default); the energy rate is
opt-in. Sources for the rate pinned there: Held & Soden 2006 / Allen & Ingram 2002. See
[[moist-ebm-source]], [[planet-rung2-scoped]].
