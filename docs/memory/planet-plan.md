---
name: planet-plan
description: Planet (project
metadata: 
  node_type: memory
  type: project
  originSessionId: 003f02e6-29a0-4a8e-8c97-0976ac82b58b
---

**Planet / Earth-system capstone plan WRITTEN 2026-06-09** (`docs/plans/planet-earth-system.md`,
per-project plan #3, §10 template). This is the program's capstone and the **first project to build a
*second* shared engine**. Part of the [[bigsim-program]]. **STATUS 2026-06-09: Phases 1–2 built +
the deep-end interactive map + `planet.ipynb` + Phase 3 (`engines/fluid`, the 2nd shared engine) all
BUILT → see [[planet-phase3-engine]]; NEXT = Phase 4 (the one-way EBM→circulation coupler).**

**Four scope decisions LOCKED (user, via AskUserQuestion — do not re-litigate):**
1. **Biomes banked EARLY** — phase order is EBM → **biomes (the payoff)** → fluid engine → coupler,
   *not* the portfolio-arc order (which lists biomes last). Rationale: the dramatic end-to-end win
   (knobs → biome-band map) needs only Phase-1 EBM temp + a diagnostic precip field, **not** the
   fluid engine/coupler, so it banks before the risky new-engine work (invariant 2, graceful
   degradation as deliberate risk management).
2. **One-way coupler in v1** (EBM forces shallow-water → an emergent **geostrophically-balanced**
   jet; anchor = geostrophic balance of the coupled jet + jet-latitude benchmark). **Two-way**
   (advected tracer, close heat budget back to EBM) = **rung 1**, seamed not built (no clean anchor
   → "doesn't bank cleanly"). **Caught in review (advisor):** thermal-wind anchor (needs vertical
   shear → multi-layer = rung 3) and poleward-heat-transport / reduction-to-EBM anchors (need an
   advected thermodynamic variable the dry single layer lacks = rung 1) are NOT v1 Phase-4 legs —
   dry single-layer SW has neither vertical shear nor a `v·T` flux. Those legs live on the staircase.
3. **v1 = "rung 0"** with the full **GCM staircase DOCUMENTED** (§5) as the growth axis — the user
   wants the GCM tar pit *eventually*. The ceiling is written as a 7-rung ladder (0=v1 reduced model
   … 6=full GCM/ESM), each rung slotting at the ADR-0001 array seam, each with a weakening anchor;
   the climb is "one validated rung at a time," gated by the three walls (compute / context-coherence
   / validation). Rungs: 1=two-way coupler, 2=moist dynamics (emergent precip), 3=multi-layer
   baroclinic, 4=real radiation, 5=idealized GCM (sphere), 6=full ESM.
4. **Deep-end INTERACTIVE map** as the viz surface (`planetmap.py`, `[webviz]` extra) — beyond the
   floor + one teaching notebook. Planet is the *one* trio project that earns ADR-0002's selective
   deep-end ("planet maps"). NOT a Streamlit dual-surface like steel.

**Non-obvious technical calls baked in (advisor-flagged, 4 issues):**
- **Phase 1 EBM reuses the Jominy-2a Strang-splitting idiom, not just the engine.** The frozen
  diffusion engine carries the *transport* `D·(1−x²)∂T/∂x` (array diffusivity, x=sin φ), but the
  radiative source/sink `S(1−α(T)) − (A+BT)` is **state-dependent** and the engine's `S` is only
  `S(x,t)` not `S(u)` → composed *around* the engine by operator splitting, exactly as Jominy split
  its `−h(T−T_air)` lateral sink. The albedo threshold makes the local step a nonlinear relaxation =
  *what creates the bistability*. Reuses [[matcalc-mc-fe-database-source]]-era discipline; spine reuse #3.
- **Snowball hysteresis = a parameter-continuation sweep** (ramp solar constant up, then down, track
  the branch), not a single solve.
- **Precipitation is a PRESCRIBED kinematic parameterization, NOT a simulated water cycle** —
  single-layer dry shallow-water has no moisture; real precip = moist thermo = GCM tar pit. v1 uses
  physically-motivated latitude bands (ITCZ wet / ~30° dry deserts / midlat wet / poles dry) C–C-scaled
  by EBM temp. Named not derived, or Phase 2 over-claims. Emergent precip = rung 2.
- **`engines/fluid` is hyperbolic/explicit (C-grid, CFL)** — shares NO machinery with the
  parabolic-implicit `engines/diffusion`. Discriminating triad legs = **PV conservation + Rossby
  radius/geostrophic adjustment** (catch a wrong Coriolis), not just gravity-wave `√(gH)`. Built
  **extension-ready**: `state` = stacked fields (1→N layers = a contract extension) + a tracer slot
  (moisture for rungs 1–2). β-plane in v1; sphere = rung 5.

**Gate/infra consequences of engine #2** (see [[test-execution-policy]]): planet's `tools/gate.py`
manifest row is the **first multi-engine `uses` entry** `{diffusion, fluid}` (steel/chip were
single-engine) — the case the per-project gate was designed to validate; and the **import-drift
guard (deferred to "engine #2") is now built** with Phase 3. `climlab` = the pycalphad pattern
(optional `[climate]` extra + frozen reference table). **§9 diligence:** an observed
biome/topography/reanalysis benchmark dataset is the lone license-check item (CALPHAD-DB analogue).

**Phase-1 sources to pin at build** (`[[…-source]]` discipline): `[[ebm-radiation-source]]`
(Budyko 1969 / North 1975 / climlab `A,B,D,α,T_freeze`), later `[[whittaker-biome-source]]`,
`[[precip-parameterization-source]]`, `[[shallow-water-source]]`.
