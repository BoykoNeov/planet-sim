---
name: planet-phase1
description: "Planet Phase 1 — latitudinal energy-balance model (EBM) & the Snowball Earth bifurcation; diffusion-engine reuse as latitudinal heat transport, Jominy Strang-splitting, the operator-splitting-not-the-engine design-fork lesson, three interchangeable solver modes, finite-cap bistability / multiple climates at today's sun"
metadata: 
  node_type: memory
  type: project
  originSessionId: 4a878190-f97f-4c5b-84c9-1e0d35fd121f
---

**Planet (project #3, the capstone) Phase 1 BUILT 2026-06-09** — `projects/planet/`
(ebm + albedo + climate_reference + plots + demo_snowball, 26-test triad +1 skip,
whole-repo fast lane **362 green** at build): the **latitudinal energy-balance model
& the Snowball bifurcation**.

**Spine reuse #3** = the frozen diffusion engine as a *sphere's latitudinal heat
transport* (heat mode, array diffusivity `(D/C)(1−x²)`, Neumann0 both poles,
x=sinφ→uniform area weight→`total()`=global mean) + the **Jominy Strang-splitting
idiom reused** (radiation split AROUND the engine: linear −B·T = exact-exp half-step;
the **albedo threshold = the nonlinear step that CREATES the bistability**).

**DESIGN-FORK saga (the durable lesson):** the "5.5 °C North-two-mode error" was
**operator-SPLITTING steady-state error** (O(dt), an oversized steady-state timestep),
**NOT the engine** — the first-instinct harmonic-mean diagnosis was **FALSIFIED
empirically** (workarounds didn't move it; a tiny dt did); the ~0.1 °C residual IS the
real harmonic-mean polar bias (engine face-averaging where (1−x²)→0). **ENGINE NOT
MODIFIED.**

**THREE interchangeable modes (USER-directed "can we have all three", 2 orthogonal
knobs):** `face` = harmonic | **exact** (pre-distort cell-D via reciprocal recursion so
the engine's harmonic-mean reproduces the TRUE face coeff → removes the polar floor) ×
`method` = relax (Strang, nonlinear-capable, the snowball path) | **direct** (dt-free
linear tridiag solve, **CONSTANT-ALBEDO ONLY**, guarded-raises on feedback; its operator
**PINNED to the engine** by a test: `engine.step == solve(I−dt·L_T/C)` to 1e-14). North
reproduced: **exact+direct = 1e-4 °C** (the tight anchor), harmonic+direct ~0.1 °C (floor
NAMED), relax→direct proves the splitting converges.

**BISTABILITY / branch (advisor-caught):** present-day Earth = the **finite-cap branch**
(needs a **capped IC** warm-equator/frozen-pole) ice line **72.9°/14.7 °C**; the
dimming-sweep branch at present S₀ is the **warmer near-ice-free** one (82.9°/15.75 °C) →
**multiple stable climates at today's sun** (the figure ★ = the colder capped one, NOT on
the down-branch — caption fixed).

**Conservation:** net-TOA = the *discrete* energy balance **machine-exact ~1e-12**; the
0-D mean matches the *continuous* T̄ only to **O(1/n²)** (point-sampled insolation ≠ ∫P2=0
on a grid) — the "machine-exact mean" claim was CORRECTED.

**Banked:** present 14.7 °C / 72.9°, **freezes at −8.3 % dimming** (Voigt–Marotzke GCM
6–9 % ✓), **re-melts at +580 W/m²**, snowball −57.5 °C (`docs/figures/planet-snowball.png`).
**UNITS = SI/climlab** (W/m², °C, x=sinφ) unlike chip's per-module native. **climlab NOT
installed** → the live cross-check skips (the frozen reference table is the committed
benchmark). Source → [[ebm-radiation-source]]. Gate: planet row added = single-engine
`{diffusion}` (fluid joins at Phase 3 → the first multi-engine row). See [[planet-plan]];
later phases → [[planet-phase2]] / [[planet-phase3-engine]] / [[planet-phase4-coupler]].
