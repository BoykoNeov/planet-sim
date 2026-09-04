---
name: pedagogy-novice-intermediate
description: planet.ipynb pedagogy — ALL THREE buckets now BUILT (C the design bench 2026-06-10; A predict-then-check + B mission cards 2026-09-04, after the drafted content failed against the code)
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

**Bucket A — "predict, then check" = BUILT 2026-09-04** (five markdown cells in `planet.ipynb`, one
immediately before each live slider in §§1–5). Shape: a **"🔮 Predict first"** question that asks for a
committed guess, then a `<details>` **"🔎 What you should see"** answer in the notebook's existing
collapsible idiom — which is *also* why the addition earns its place: the payoff survives on GitHub /
nbviewer, where none of the widgets run.

**The prompts drafted in 2026-06-10 did not survive contact with the code — check every one against its
widget before writing it down.**
Two were wrong (advisor-caught): the drafted §2 prompt ("dim then re-brighten — do you get the same climate
back?") describes a control the §2 widget *does not have* (it is one Sun × two starting climates), and the
drafted §4 prompt ("which way does the jet blow") belongs to §5's widget — §4's live cell is grid
refinement. Same staleness class as the `circ_precip.py` pointer and the `∝1/D` line in plan §12.2: the
drafts predate several builds. **Every answer was measured before being written**, not recalled:
- §1 — the pole warms ~2× the equator while ice retreats (+4 vs +9 °C, 1361→1400) and **the sign flips once
  ice-free** (+4.0 vs +2.6 °C, 1480→1520). "It was never about latitude, it was about ice"; cross-links to
  §8.2, which gets polar amplification from moisture with no ice at all.
- §2 — reframed to a question cell 10 does *not* already answer: **how many** Suns give both starting
  climates the same answer. Exactly one (the dim closure ≈1253); the frozen branch does not melt out until
  **≈1800 W/m², ~30 % brighter than today** — the asymmetry *is* the loop width.
- §3 — asks about **size only**. The first cut asked "which way *and* does it grow?", but §3's heading cell
  already says the bands migrate poleward, so half the question was spoiled before the reader reached it
  (advisor-caught on a second pass): **the spoiler check is against the surrounding prose, not only against the
  widget** — the same move §2's prompt needed. Answer: **smaller**, dry share 24 %→10 %, tundra leaves the
  legend; the shrinking is the counter-intuitive half and the hook into Mission 2.
- §4 — the sharpest of the five: refining the grid moves the **numerical error** (5.9 %→0.8 %) and **not the
  physics** (radiated fraction ~86 % at every resolution), while the *bump-size* slider moves the physics
  (95 %→53 % across 0.15→1.0 `L_R`). "A number that keeps moving as you refine has not measured anything."
- §5 — the jet **slows** as the Sun brightens (20→16 m/s: melting ice weakens the very gradient that drives
  it) and its latitude barely moves (44°→43°) — §1's prediction resurfacing one engine downstream.

**Bucket B — mission cards = BUILT 2026-09-04** (one markdown cell at the end of §7, before §8 — the first
point where every knob the missions need is on the page). Four cards in a fixed *where · allowed knobs ·
goal · done-when* shape. **All four success criteria were verified reachable on the actual slider ranges**,
and two of the drafted missions were not:
- *Make a desert world* — the draft ("push −A until the temperate bands collapse") is **backwards**. Warming
  *shrinks* the dry bands (a warmer atmosphere rains more: 24 %→10 %); the desert maximum is a **dim** Sun
  (~38 % at 1250 W/m² plus a few greenhouse notches), bounded by the Snowball cliff. The card is now built
  *around* that reversal, with §3's prediction box as the set-up.
- *Build an M-dwarf world* — the draft was **unreachable**: the red-star freeze cliff sits at ≈1150 W/m²,
  below the bench's 1250 minimum. Rebuilt as a one-notch demonstration at fixed `S₀ = 1250` — Sun-like
  (5772 K) → frozen at −46 °C; dragged to 2600 K → alive at +5 °C, ice line 59°.
- *Find the habitable-zone edges* — reframed as **one edge and an honest absence**: the cold edge brackets
  between the bench's two dimmest notches (1250 frozen / 1260 not), and the hot edge **is not in this
  model** (the OLR is linear in T; the runaway is §8.5's territory). "Find where the model stops modelling"
  turns the project's honesty gradient into an exercise.
- *Freeze the planet* — §2's slider, the ≈1253 cliff, unchanged from the draft.

**The third deferred item, the two-world diff, was ALREADY BUILT** — §7 carries both
`planet_spec.diff` ("Compare two worlds") and the `[webviz]` A · B · Δ comparison globes
(`planet_spec.delta_view` → `planetmap.render_comparison`). The plan §12.4 checkbox was stale; struck.

All three buckets are now built. A and B were indeed nearly free in *mechanism* (markdown over machinery
that already works) but **not** in content: the cheap half is the prose, the expensive half is measuring
every answer and re-deriving every prompt against the widget it sits in front of — half the drafted content
was wrong. Markdown-only means no output re-banking and no `docs/index.html` regeneration; still re-run
`test_planet_notebook.py` **solo** on the final state ([[teaching-surfaces-resynced]]). [[planet-plan]].
