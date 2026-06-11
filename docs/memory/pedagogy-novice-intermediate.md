---
name: pedagogy-novice-intermediate
description: planet.ipynb pedagogy roadmap — bucket C (design-a-world sandbox) BUILT 2026-06-10; buckets A (predict-then-check) + B (mission cards) DEFERRED with their concrete content
metadata: 
  node_type: memory
  type: project
  originSessionId: 81e88905-468f-45b8-ad33-1bde80e06f27
---

After the expert-tier "Going deeper" collapsibles landed (commit f64ffeb tiered §1–§5
*upward*), the user asked for the **novice/intermediate + experimentation/design** direction.
The notebook already had a "three depths" scaffold (narrative → sliders → ▸ Going deeper); the
gap was that everything is a **guided tour** (read top-to-bottom, free-twiddle knobs). Three
buckets were proposed to add **active**, **goal-directed**, and **generative** modes. User: "I
like all three — mark A and B for the future, now work on C."

**Bucket C — "Design your own world" sandbox = BUILT 2026-06-10** (§7 of `planet.ipynb`, 5 cells
between §6 and the provenance table). Surfaces the previously-invisible `planet_spec` save/load as
a *creative workflow*: all six knobs at once (S₀, CO₂→A, D, obliquity, T_star, size) → live biome-map
preview → name + 💾 Save button → a shareable planet-spec file → `load_world()` re-runs it. Module
support added (both with real tests, full suite green):
- `planetmap.climate_params(...)` — extracted the knobs→`EBMParams` composition (obliquity `s2` +
  exoplanet `ai`/`D`) into the **single source of truth** `climate_view`, `planet_spec.build_spec`,
  and the notebook bench all ride on (behavior-preserving refactor; the canary
  `test_round_trip_identity_of_a_real_present_day_view` stayed green untouched).
- `planet_spec.build_spec(**knobs)` — the "knobs → saveable `PlanetSpec`" entry point. **The trap it
  fixes:** `climate_view` builds the composed params *and discards them*, so a naive
  `from_view(climate_view(...), EBMParams(S0,A,D))` stores the **un-perturbed** knobs for any
  non-default star/size/tilt → the exported world re-runs to the WRONG climate. `build_spec` solves
  and serializes the *same* composed params, so `to_params()` reconstructs the exact world. Stays
  headless (no Plotly pull) — the `test_importing_planet_spec_stays_headless` subprocess guard holds.
Design calls (advisor-endorsed): save fires **only on the button** (Run-All never writes a file → the
notebook smoke test + CI stay clean); saved worlds go in **`outputs/worlds/`** (gitignored); load
re-runs from `to_params()` (deterministic, teaches the API) and the prose is honest that the stored
*arrays* are the manifest the §6 globe / external tools consume, not what the matplotlib bench reloads.
See [[planet-interactive-map-design]].

**Bucket A — "predict, then check" (DEFERRED; ~free, markdown-only).** Before each §'s slider, one
prompt turning passive reading into a hypothesis-test loop. Concrete prompts to drop in:
- §1: "Before you brighten the sun — will the *pole* or the *equator* warm more? Drag S₀ and watch."
- §2: "You're at today's sun. Predict: if you dim then re-brighten to exactly the same S₀, do you get
  the same climate back? Now do it." (the hysteresis 'gotcha').
- §3: "Will warming move the *deserts* toward the equator or the poles? Crank −A and see."
- §4/§5: "Which way does the jet blow — east or west? Predict before the arrows render."

**Bucket B — goal-directed "mission" cards (DEFERRED; cheap, reuses existing knobs).** Short challenge
cards, each with a success criterion read straight off the figure:
- *Freeze the planet* — find the Snowball threshold by hand with S₀ (§2).
- *Make a desert world* — push −A / knobs until the temperate biome bands collapse (§3).
- *Find the habitable-zone edges* — the inner/outer S₀ where biomes appear / vanish.
- *Build an M-dwarf world* — use the star-spectrum knob; watch the ice line move (the §9.1 exoplanet
  work). Same framing for the size and obliquity knobs.
**Deferred follow-on to C — "two-world diff" (a stretch, NOT core C).** Load two saved specs and show
Earth vs. your-world side by side (the layer registry already supports multiple views; a small
side-by-side builder, or just two `load_world` panels). Was framed as a follow-on, correctly not built
in the first cut.

A and B are nearly free (markdown + framing over machinery that already works) and high-leverage for
novice engagement; C was the bigger "turn the tour into a lab" piece and is done. [[planet-plan]].
