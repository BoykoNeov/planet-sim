---
name: planet-rung5a-orographic
description: "Rung 5A (the first step off the zonal mean toward the north star) BUILT 2026-07-09: planet/orographic.py = Smith & Barstad (2004) linear orographic precip — a diagnostic 2-D rain shadow that wakes the inert elevation seam; a TRADE not the engine leaving the zonal mean; 5A.2 (sphere placement + jet source + serialization + demo) deferred"
metadata:
  node_type: memory
  type: project
  originSessionId: fdb095c9-cdad-43a3-ac05-2e4a17bc3da3
---

**Rung 5A = the first laptop-scale step toward the "north star"** (true 2-D longitudinal geography →
regional climate, orographic precip, rain shadows; plan §12.5 / §5). Literal rung 5 (a full idealized
GCM) is infeasible-tier, so it is climbed as a **spike-first sub-ladder of reduced models** — the same
decomposition the gas-giant sketch uses ([[gas-giant-feasibility]]). Continues [[planet-plan]];
follows the diagnostic-precip precedent [[planet-phase2]] one rung out.

**BUILT 2026-07-09** — `planet/orographic.py` + `planet/tests/test_orographic.py` (10 anchors green,
fast lane). **Smith & Barstad (2004) linear theory of orographic precipitation:** a *diagnostic* on a
**prescribed uniform** background wind over a **2-D terrain** → a wavenumber-space transfer function
(one FFT, laptop-trivial) → windward rain + a lee **rain shadow**. Numerics + constants pinned in
[[smith-barstad-orographic-source]].

**What it banks (the payoff):**
- **Wakes the dormant elevation seam.** Since v1 `planet_spec` carried an `elevation` layer tagged
  `inert=True` (carried/displayed/round-tripped, climate-inert; §9.3). This is the module the plan
  named as the one that finally makes elevation *do something*.
- **Tight anchor = convergence to the closed-form triangle-ridge solution** — but honestly scoped
  (advisor-caught): in that limit `(1 − i·m·H_w) → 1`, so `m` and its `sgn(σ)` branch DROP OUT — the
  triangle anchor pins only the *reduced* transfer function (`C_w` scaling, upslope `iσ`, fallout
  `τ_f`), NOT the branch or `H_w`. **The branch is guarded by exactly one test — the rain-shadow
  direction** (windward wetter than lee, peak upwind); empirically, flipping the branch sign reddens
  *only* that test (shadow flips to the windward side). The wind-reversal **mirror** (machine precision,
  1.8e-15) is a *reflection self-consistency* check — invariant under a branch flip, so NOT a branch
  test. Plus: upslope limit `C_w·max(0,U·∇h)` to the analytic Gaussian gradient; flat-ground null.

**The honesty flag (a TRADE, not a win — do not overclaim):** it makes the **precipitation** 2-D, NOT
the engine. The **temperature** climate underneath stays **zonal-mean** (the EBM). So Rung 5A does
*not* claim the engine has left the zonal mean — only the rain shadow has. Regional Cartesian patch,
uniform prescribed wind; cross-mountain flow is **prescribed, not emergent**.

**Advisor-caught traps (all now handled):** (1) the upslope limit is `H_w=τ_c=τ_f=0`, NOT `U→∞`
(large U drives precip to *zero* via the τ factors); (2) the `sgn(σ)` branch is the #1 bug — its ONLY
validator is the **rain-shadow direction** test (NOT the triangle anchor, where the branch drops out;
NOT wind-reversal, which is branch-flip-invariant); (3) zero the `σ=0` locus explicitly; (4) zero-pad
the FFT so lee drying doesn't wrap. Also: compare the **truncated** model to the exact solution
(untruncated has a −C anti-rain spike at the downwind kink).

**5A.2 — DEFERRED + named** (the integration + beauty pass): sphere placement of the patch; **where the
cross-mountain wind comes from on a zonal-jet globe** ([[planet-phase4-coupler]] gives a purely zonal
jet — the honesty caveat is that cross-mountain flow is prescribed); cm/yr↔mm/hr + lat×lon grid
integration with [[planet-phase2]]'s biome map; serialization ([[planet-interactive-map-design]]
schema); the demo/figure (a rain shadow behind a *real* mountain range on the globe). Provisional per
the living-staircase rule. Then Rung 5B/5C = the 2-D EBM engine step (needs seasonality for
continentality) / Charney–Eliassen stationary waves.
