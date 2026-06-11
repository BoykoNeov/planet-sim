---
name: planet-phase3-engine
description: "Planet P3 BUILT — engines/fluid (2nd shared engine, shallow-water C-grid) + circulation.py; the design forks, what the scheme actually conserves, and the gate/import-drift consequences"
metadata: 
  node_type: memory
  type: project
  originSessionId: d28b5a93-2ea3-4982-9548-32d186e4dee6
---

**Planet Phase 3 BUILT 2026-06-09** — `engines/fluid` (the program's **second shared
engine**) + `projects/planet/circulation.py` + demo + banked figure. Part of
[[bigsim-program]] / [[planet-plan]]; sources [[shallow-water-source]].

**What it is.** A rotating shallow-water solver on a **doubly-periodic β-plane**, **Arakawa
C-grid**, **vector-invariant** form, **explicit SSP-RK3** — deliberately a *different solver
class* from the parabolic-implicit `engines/diffusion` (hyperbolic, CFL-limited). Built
standalone, **frozen behind `engines/fluid/CONTRACT.md`** before any coupling. State =
`SWState(h,u,v, tracer=None)` stacked plain 2-D arrays (the ADR-0001 boundary). Files:
`engines/fluid/{shallowwater.py, __init__.py, CONTRACT.md, tests/}` (24 tests, 4 slow),
`circulation.py` + `tests/test_circulation.py`, `demo_shallowwater.py` + test,
`plots.shallowwater_figure` → `docs/figures/planet-shallowwater.png`.

**Durable design calls (advisor-guided):**
- **Built nonlinear, NOT linear.** The advisor's gate: a *finite-amplitude* PV/enstrophy seal
  (a balanced vortex at Rossby ~0.5, where advection genuinely moves PV) — at small amplitude
  that leg is vacuous. The seal came out **crisp** (enstrophy drift ~1e-4 over 3 inertial
  periods, PV extrema growth <5e-3), so the engine stayed nonlinear. The fallback (downscope to
  a linear core, defer advection to rung-1) was NOT needed. Don't write FROZEN until that seal passes.
- **The symmetric scheme conserves ENERGY, not enstrophy** (one Sadourny-class scheme conserves
  one or the other; Arakawa–Lamb both, not built). Empirically: **mass machine-exact (~1e-13)**;
  **energy drift dt³-convergent** (5e-10→6e-11→8e-12 per dt halving → it's the semi-discrete
  invariant); **enstrophy drift flat ~1e-7** (spatial-discretization-limited, shrinks with Δx
  not dt). Contract states conservation **as measured**, not aspirational. Claim machine-precision
  ONLY for mass.
- **CFL guard fires at the true stability threshold (Courant≈1, `max_dt(safety=1.0)`), not the
  recommended 0.3× step** — else passing `max_dt()` false-trips on float rounding.
- **tracer slot = "seam, not machinery"** (the planetmap `vector_overlay` idiom): declared on
  `SWState`, `step` raises `NotImplementedError` naming rung 1; not advected in v1.

**Validation triad (all green):** gravity `√(gH)` & **Poincaré ω²=f₀²+gHk²** to ~1e-3 (rotation
check); **Rossby** westward+dispersive, loose **and converging to analytic with resolution**
(0.96→0.978→0.982 at N=48/96/144 — a named numerical-dispersion edge of the slow balanced mode,
NOT a bug; β-independent so not the β-plane approx); **geostrophic adjustment → Helmholtz over
L_R to ~1%**; balanced zonal jet steady (exact nonlinear steady state since v=0). circulation.py:
L_R(45°)≈960 km, 86% of an adjustment bump radiates away, Rossby c≈−4.4 m/s (analytic −4.7).

**Gate/infra consequences realized:** `tools/gate.py` — `planet` `uses` = **{engines/diffusion,
engines/fluid}**, the manifest's **first genuinely multi-engine row**. The **import-drift guard**
(deferred to engine #2 per ADR 0003) is **BUILT & live** in `tools/tests/test_gate.py` (ast-scans
each project's source `engines.*` imports ⊆ declared; runs inside the per-project gate). **No
pyproject change** (testpaths `engines`/`projects` already collect `engines/fluid/tests`
recursively). Fast lane: ~2.5s engine + cheap planetary-number tests; integration/demo/dispersion
= `slow`. Full planet gate **129 passed, 1 skip** (live-climlab) in 130s. **Steel untouched**
(another agent may be working steel concurrently).

**NEXT = Phase 4** (the one-way EBM→circulation coupler: a geostrophically-balanced jet emerges;
the map then registers a `vector_overlay` circulation layer, no renderer edit). Two-way = rung 1.
Engine extension seam: stacked fields → N layers (rung 3) / advected tracer (rung 1); rigid
channel walls + sphere are named-unbuilt.

**Phase-4 carry-forward (advisor forward-notes, none blocked P3):** (1) the domain is
**doubly-periodic**, but the EBM's equator→pole gradient is **monotonic/non-periodic** → force the
jet with a **periodic meridional profile** (the validated balanced sinusoidal zonal jet proves this
works), a deliberate choice, not a wall-bounded channel; rigid walls remain an additive BC extension
that won't un-freeze periodic behaviour. (2) `energy()` is **total** energy (incl. the resting ½gH²
PE background) → `|ΔE|/E` flatters the *dynamical* drift; if a Phase-4 figure claims dynamical-energy
conservation, subtract the resting state (the dt³ test is unaffected). (3) `solve(dt=None)` fixes dt
at the initial speed → a jet **spun up from rest** will trip the per-step CFL raise; pass explicit dt
or step in segments (now documented).
