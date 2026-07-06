# Memory index

Per-project memories for **planet-sim** (planet simulator), extracted from the BigSim
program's shared memory on 2026-06-10 when the monorepo was split into standalone repos.
Cross-cutting *program* memories (bigsim-program, parallel-dev, test-execution-policy,
the notebook-discipline notes, the github-repo pointer) stayed in the BigSim archive store;
`[[links]]` to those will dangle here by design.

One line per entry — the hook, not the record. Open the topic file for the full detail.

## Working preferences (how the user wants me to work)
- [Always push commits](always-push-commits.md) — push every commit to the remote immediately; commit+push at the end of every work batch (straight to `main`, linear history)
- [Detailed test instructions](detailed-test-instructions.md) — when handing a manual test/verify to the user (browser play-through = the hand-off boundary), give numbered steps: what to open, how to start, what to look at, expected pass/fail, diagnostic to copy back
- [No pointer spike lines](no-pointer-spike-lines.md) — standing viz pref: `showspikes=False` on EVERY Plotly figure (spikes default ON for 3-D); tooltips/modebar kept
- [Viz prose novice→intermediate](viz-prose-novice-intermediate.md) — standing viz pref: richer novice→intermediate prose on every figure+notebook; appended glossaries/definitions/`<abbr>` hover now welcome (relabel-don't-append REMOVED 2026-06-12); figure sweep done

## Repo / infra
- [Memory vendored in repo](memory-vendored-in-repo.md) — auto-memory lives in the repo at `docs/memory/` (committed), junctioned to the hardcoded global Claude path; on a fresh clone recreate the junction or recall misses this history
- [Flatten repo-root gotcha](flatten-repo-root-gotcha.md) — post-flatten `parents[N]` off-by-one (fixed 2026-06-10); CI's `pip install -e` MASKS the bare-checkout failure → verify the suite without an editable install
- [Engines are living contracts](engines-living-contracts.md) — freeze-before-reuse DROPPED 2026-06-10 (ADR 0005): extend engines directly + test + changelog, no re-seal; kept ADR-0001 array boundary + validation triad
- [Teaching surfaces re-synced](teaching-surfaces-resynced.md) — notebook §8 "Up the staircase" + README rung log through rung 4 (2026-06-15); gotchas: `index.html` is GENERATED (byte-pinned by a slow test), `planet.ipynb` too big for Read (edit the JSON)

## Capstone (the 4 phases) + the interactive map
- [Planet plan](planet-plan.md) — project #3 capstone plan (2026-06-09; ALL 4 phases + map + notebook BUILT = CAPSTONE COMPLETE): the 4 locked scope decisions, the GCM staircase documented; → [[planet-phase3-engine]] [[planet-phase4-coupler]]
- [Planet Phase 1 EBM](planet-phase1.md) — latitudinal EBM + Snowball bifurcation (diffusion-spine reuse #3, Strang split); the "5.5 °C error" was operator-SPLITTING not the engine; finite-cap bistability; [[ebm-radiation-source]]
- [Planet Phase 2 biomes](planet-phase2.md) — climate→biome map: diagnostic precip + an original sloped-boundary Whittaker (T,P) partition (the Irvin precedent, not plotbiomes); C-C scales by global-T̄; warming migrates bands poleward
- [Planet Phase 3 engine](planet-phase3-engine.md) — `engines/fluid` (2nd shared engine: rotating shallow-water C-grid / SSP-RK3) BUILT 2026-06-09; conserves energy-not-enstrophy (mass exact only), CFL guard; → [[planet-phase4-coupler]]
- [Planet Phase 4 coupler](planet-phase4-coupler.md) — one-way EBM→shallow-water coupler BUILT 2026-06-09 (emergent jet 16.5 m/s @ 42°): windowed target, emergence via a synthetic gradient, `vector_overlay` now painted
- [Planet interactive-map design](planet-interactive-map-design.md) — deep-end map CONVERGED + ADR 0004; `planetmap.py`/`planet_spec.py` BUILT 2026-06-09; exoplanet+obliquity knobs + design-a-world bench + two-world diff BUILT 2026-06-10; layer registry + pin-the-schema round-trip
- [Interactive what-if](interactive-what-if.md) — no-install browser what-if (`docs/interactive/index.html`, precomputed grid) + `explain.py` shared prose engine + notebook enrichment; built 2026-06-12
- [Pedagogy: novice/intermediate](pedagogy-novice-intermediate.md) — `planet.ipynb` pedagogy: bucket C (design-a-world sandbox) BUILT 2026-06-10; buckets A (predict-then-check) + B (mission cards) DEFERRED with drafted prompts

## The rung staircase (rung 0 = the capstone; rungs 1→4 climb toward a GCM)
- [Planet rung 1 two-way coupler](planet-rung1-two-way-coupler.md) — rung 1 COMPLETE 2026-06-11 (`transport.py`/`eddy_flux.py`/`circ_precip.py`): emergent eddy κ→D bridge + reduction-to-EBM + circ-informed precip; κ~10³ config-tuned (named-not-banked), reduction near-vacuous → rung 3; → [[planet-rung2-scoped]]
- [Planet rung 2 scoped](planet-rung2-scoped.md) — rung 2 moist Phase A BUILT 2026-06-11 (`moist.py`): opt-in ~2.5%/K linear rate + a diagnostic P−E budget (extratropical-only; eq/subtropics wrong); wall = prescribed `R_ATM_SLOPE`; → [[moist-ebm-source]] [[planet-rung25-mse-diffusion]]
- [Planet rung 2.5 MSE-diffusion](planet-rung25-mse-diffusion.md) — rung 2.5 MSE-diffusing moist EBM BUILT 2026-06-12 (`moist_ebm.py`): emergent polar amplification ~2.05×/1.80× (dt-free) from moisture transport alone; `D_eff(T)` inside the divergence; direction banked / magnitude loose; wall = recalibrate `D_s`
- [Planet rung 2.x ITCZ](planet-rung2x-itcz.md) — full-sphere EBM + energetic ITCZ BUILT 2026-06-14 (`sphere_ebm.py`): ITCZ = the energy-flux-equator, migrates toward the warm hemisphere; the sensitivity is closed-form (not emergent); → [[planet-rung2-hadley-fix]]
- [Planet rung 2 Hadley fix](planet-rung2-hadley-fix.md) — deep-tropical backwards-P−E fixed 2026-06-14 (`moist.py`): opt-in prescribed Hadley cell flips the ITCZ sign; convergence = plumbing, amplitude ~C–C 6.6%/K = the nugget; does NOT relocate the desert
- [Planet rung 2.x emergent ITCZ rain](planet-rung2x-emergent-itcz.md) — full-sphere moisture budget (eddy + two-cell Hadley on the EFE), ∫(P−E)=0 exact → ITCZ rain RAINED not painted; intensity geometric not q; co-location as a falsifiable check
- [Planet rung 3 scoped](planet-rung3-scoped.md) — rung 3 baroclinic instability Phase A BUILT 2026-06-12 (`engines/fluid/layered.py`): two-layer free-surface SW linear growth; Phase-B payoff = the OPEN BET; anchors Phillips/Eady/Charney; plan §10
- [Planet rung 3 Phase B — outcropping finding](planet-rung3-phaseB-outcropping.md) — rung 3 Phase B SPIKE 2026-06-13 NEGATIVE: free-surface SW OUTCROPS at saturation (Froude control) → re-routes to two-layer QG (rigid-lid); dissipation BUILT default-off; → [[planet-rung3-scoped]]
- [Planet rung 3 QG built — bet WON](planet-rung3-qg-built.md) — rung 3 Phase B BUILT 2026-06-13 (`baroclinic_qg.py`) BET WON: saturated flux down-gradient + irreversible (irr 0.96–1.0) at κ/(v'L_d)~O(1) → reduction-to-EBM finally NON-VACUOUS; banked dimensionless+qualitative; rung 3 COMPLETE
- [Planet rung 4 radiation](planet-rung4-radiation.md) — rung 4 gray radiative transfer BUILT 2026-06-14 (`radiation.py`, sibling): emergent OLR retires the B-fixed wall; B=2 ≈ Planck−WV+LR (non-circular decomposition); per-latitude wire TROPICAL at Earth loading; spectral CO₂ log law; plan §10/§12.2

## Visualization showcase (the flow-globe ladder)
- [Planet viz animation rungs](planet-viz-animation-rungs.md) — eddy-flow viz rungs A/B/C BUILT 2026-06-11/12/13: A=matplotlib GIF, B=Plotly globe, C=`flow_globe.py` three.js particle globe; renderer-agnostic `FlowField` contract; honest-by-disclosure carve-out (ADR 0002); three.js vendored inline
- [Ocean currents viz rungs](ocean-currents-viz-rungs.md) — ALL FIVE O-rungs BUILT + BANKED 2026-07-06 (§9.6): O1 mask + O2 OSCAR producer + O3 beauty pass + O4 seasonal frames (browser-verified) + O5 `flow_field_from_qg` (the QG emergent producer; §9.4 rule-of-three re-affirmed HOLD); → [[planet-spinout-roadmap]]
- [Planet spin-out roadmap](planet-spinout-roadmap.md) — editable-ocean GPU spin-out (a Julia/ClimaOcean+Makie repo born FROM planet-sim); R1+R2 BUILT 2026-06-14 (`flow_serialize.py`); two seams (viz-out vs forcing-in), rungs R1–R3, spinoff S1–S5, ECCO the S1 anchor; plan §11

## Feasibility sketches (not built)
- [Gas-giant feasibility](gas-giant-feasibility.md) — sketch 2026-06-13 (NOT built): gas-giant atmosphere/vorticity on these engines = 3 tiers (β-plane jets ≈1 rung on `baroclinic_qg.py` → new sphere geometry → deep interior); full record `docs/explorations/gas-giant-atmosphere.md`

## Cited-source pins (the `[[…-source]]` discipline — numbers each module pins to literature)
- [EBM radiation source](ebm-radiation-source.md) — P1 cited climlab EBM defaults (A/B/D/Tf/albedos/S0/s2) + benchmark outputs (ice line ~70°, snowball 6–9%) that `ebm.py`/`climate_reference.py` pin
- [Precip parameterization source](precip-parameterization-source.md) — P2 cited zonal-mean precip-by-latitude + C-C 7%/K vs 2–3%/K energy-constrained (named gap); pattern×CC(global-T̄); rung-2 rate now built → [[moist-ebm-source]]
- [Whittaker biome source](whittaker-biome-source.md) — P2 cited Whittaker 1975 biome diagram (T °C, P cm/yr) + pinned (T,P) thresholds; the GRAPHICAL→original-rule-partition decision (Irvin precedent, not plotbiomes)
- [Shallow-water source](shallow-water-source.md) — P3 pinned Earth Ω/a/g + β-plane f/β, Poincaré/Rossby dispersion, L_R; H=1000 m → L_R(45°)≈960 km; Gill/Vallis/Sadourny/Arakawa-Lamb
- [Moist-EBM source](moist-ebm-source.md) — rung-2 cited sources `moist.py` pins: C-C q_sat + energy-constrained ~2–3%/K (Held&Soden2006) + the fixed-RH diffusive moist EBM (Flannery1984/Hwang&Frierson2010)
- [Obliquity insolation source](obliquity-insolation-source.md) — §9.1 (2026-06-10): cited daily-mean-insolation + mean-annual s₂(ε) DERIVED; recovers climlab s₂=−0.48 @23.44°; s₂(0)=−5/8 exact
- [Stellar-spectrum ice-albedo source](stellar-spectrum-ice-albedo-source.md) — §9.1: Joshi&Haberle2012/Shields host-SED→ice-albedo + Warren1982 snow optics; two-band ratio to solar (Sun recovers ai=0.62); size knob D∝1/size²
- [Ocean albedo+transport source](ocean-albedo-transport-source.md) — the 4th interactive knob (`ocean.py`, 2026-06-15): sea fraction → a0 (firm) + D (loose); Earth 0.71 = identity anchor; → [[interactive-what-if]]
