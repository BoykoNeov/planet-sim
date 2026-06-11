---
name: planet-phase4-coupler
description: "Planet Phase 4 (the capstone-completing one-way EBM→shallow-water coupler) — design forks, the y-periodicity crux, the conservation reframe"
metadata: 
  node_type: memory
  type: project
  originSessionId: 2b3be86f-e0d4-4f62-8d8d-b2451bac82ae
---

**Planet Phase 4 BUILT 2026-06-09 — the capstone is complete (all 4 planet phases done).**
`projects/planet/coupler.py` (+`demo_coupler.py`, `plots.coupler_figure`, 9-test `test_coupler` +
2-test `test_demo_coupler`) **couples the two frozen shared engines**: the EBM's equilibrium
meridional temperature gradient forces the frozen shallow-water engine and a
**geostrophically-balanced westerly jet emerges**. Banked: jet **16.5 m/s @ ~42°**, core geostrophic
residual **~0.6%**, `docs/figures/planet-coupler.png` + interactive `planet-coupler-map.html`. Full
planet gate **140 passed, 1 skip**; no engine modified; planet `uses` unchanged `{diffusion,fluid}`.

**Design fork (advisor-decided): Option A — sustained thermal relaxation + weak Rayleigh drag, split
around the BARE frozen engine** (NOT Option B init-and-adjust, which = Phase-3's "adjustment in
isolation" the plan explicitly distinguishes from). The frozen `engines.fluid.step` has no forcing
term, so forcing is composed *around* it by **operator splitting = the THIRD reuse of the
EBM/Jominy `_radiation_half` idiom**: exact-exponential `h→h_target+(h−h_target)e^{−½dt/τ_relax}` +
drag `(u,v)·e^{−½dt/τ_drag}`, half/full-step/half. `τ_relax=3`, `τ_drag=15` inertial periods
(`τ_drag≫1/f` ⇒ near-geostrophic, advisor's note). `φ_ref=40°`, `n_LR=4.5`, `α=HEIGHT_PER_KELVIN=28`.

**THE Y-PERIODICITY CRUX (advisor's blocking blind-spot, the core lesson).** The engine is
**doubly-periodic** (walls in y = its named, *unbuilt* BC extension). A monotonic warm→cold height
target jumps by the full equator–pole contrast at the y-seam → a **spurious giant boundary jet** that
dominates the real one. Fix = a **windowed (Tukey), discretely zero-mean** height anomaly
(`height_target`): C¹-periodic (matched value+slope at edges) AND mass-neutral. The periodic channel
then *necessarily* exacts a **flanking easterly return** (∮u≈0, zero net zonal momentum). **DON'T
overclaim it** (advisor's final catch — I papered the easterly optics with an overclaim first): the
return's *existence* is geometric+physical and its E–W–E *sign* banding resembles the general
circulation, but this single-layer periodic channel does **NOT** reproduce observed westerly-dominant
magnitudes — in the banked run the **poleward easterly (−27) is actually the STRONGEST band** (> the
+16.5 westerly), and its concentration is **window-construction-dependent**, not observed. So: jet =
validated/benchmarked feature; the return = a **named scope edge**, NOT a trade-wind/polar-easterly
reconstruction. (Also: the 0.6% geostrophic residual is the **weakest** leg — near-self-fulfilling since
the flow locks onto a balanced target; **emergence + release are the load-bearing legs**.)
- **DEAD-END tried & rejected:** de-trending the anomaly (subtract the linear trend for periodicity)
  removes the *bulk* gradient → jet ∝ curvature only → **~2 m/s, far too weak**. The EBM profile is
  near-linear over the band, so its deviation-from-linear is tiny. Windowed-anomaly is right; de-trend
  is wrong. (Windowed-*gradient*-integrated also works but is less intuitive; the easterly came out
  comparable either way — it's geometry, not the method.)
- **Cliff exclusion:** the channel brackets the *smooth* midlatitude baroclinic zone and **excludes
  the ice-line albedo cliff** (~73°) — including it lets the cliff's sharp gradient dominate the
  forcing (an over-strong poleward easterly).

**CONSERVATION REFRAME (advisor-blessed build-time honesty call; the plan's wording was internally
inconsistent).** A **forced–dissipative** system does NOT conserve energy or PV — the forcing injects
APE, the drag removes KE, and *that balance is what selects the steady jet*. So the plan's
"(mass,PV,energy) preserved under the steady forcing" is reframed (same class as Phase-3's "energy OR
enstrophy, as measured") to what is true: **(a) mass machine-exact under forcing** (zero-mean target +
engine's mass invariant, ~1e-13), **(b) a RELEASE test** — switch forcing & drag OFF, run the bare
engine: mass/energy/enstrophy conserved (Phase-3 guarantees re-confirmed) **AND the jet persists**
(16.5→16.5) ⇒ a genuine balanced state, not forcing-propped. Plan §3 Phase-4 bullet + coupler docstring
updated; this memory is the third leg the advisor required. See [[planet-phase3-engine]].

**EMERGENCE PROOF (advisor's non-circularity caveat — don't hand-place the jet at channel centre).**
The jet sits at the **EBM gradient maximum (~45°)**, *poleward of the channel centre (40°)* → emergent,
not placed. Climate changes (warm/dim) barely move it (the gradient peak is robust), so the **decisive
test = a SYNTHETIC off-centre gradient** (`test_jet_latitude_tracks_the_climate_gradient_not_the_channel`):
a step at 30° → jet at 36°, at 50° → jet at 50° — the jet *follows the gradient*. Non-circularity split:
jet **latitude + geostrophic balance = amplitude-independent (validated)**; **height-per-K `α` → jet
speed = calibrated (tuning)** (`test_jet_speed_scales_with_amplitude_but_latitude_does_not`).

**Scope edge (dry single layer, one-way):** NO poleward heat transport / reduction-to-EBM (rung 1,
needs the `tracer`), NO thermal wind (rung 3, needs vertical shear). Two-way = rung 1, seamed not built.

**`vector_overlay` seam now PAINTED (the deferred Phase-4 machinery).** `planetmap.circulation_layer`
maps the jet's `u(φ)` onto the full globe (mirrored both hemispheres, 0 outside the band) as a stacked
`(2,n_lat,n_lon)` `[u,v]` layer; `_vector_overlay_trace` paints it as **Plotly cones** on the
sphere-tangent plane (eastward cones at midlats = the westerly jet); `build_view(jet=…)` registers it.
**Computed-then-viewed, NOT in the live-slider loop** (the jet integration is the first compute too
heavy for the rung-0 instant remap, §9.2). Round-trips through `planet_spec` (arbitrary array shape).
The render test (was "raises NotImplementedError naming Phase 4") flipped to "paints a Cone"; a fake
`SimpleNamespace` jet keeps it fast (no integration). See [[planet-interactive-map-design]].
