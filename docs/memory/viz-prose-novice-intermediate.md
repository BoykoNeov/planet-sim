---
name: viz-prose-novice-intermediate
description: "User wants richer, novice→intermediate prose on EVERY visualization, now and future: plain-language labels (technical term parenthetical) + a compact caption defining the domain jargon, so a standalone figure stands on its own — not 'a glossary on every plot'"
metadata:
  type: feedback
---

**Standing preference (user, 2026-06-12): every visualization must carry richer, novice→intermediate
prose so a reader who doesn't know the domain jargon can still read it.** The user pointed at the Rung-B
eddy globe and listed the terms it failed to explain: *what is a "β-plane band", "planet-wide flow",
"instantaneous flux" (and what does "reversible" mean, to what %), "band streams", "κ residual",
"swirls range"* (the last is them reading `(the swirls rage)` off the legend — proof the cryptic
editorializing labels were the problem).

**The vehicle = relabel, don't append** (advisor; resolves the tension with the *same-day* "remove
overlapping text" request — see [[no-pointer-spike-lines]] for that one). Two moves, in order:
1. **Plain-language labels with the technical term parenthetical.** Rewrite cryptic legend/axis/subtitle
   text into what it *means*, keep the jargon in `(…)` for the intermediate reader. This makes the
   *existing* text legible instead of adding new elements. (Eddy globe: `throughput Σ∫|F̄|dt (the swirls
   rage)` → `back-and-forth stirring (throughput)`; `net … (the small residual)` → `net heat moved
   poleward (κ residual)`; subtitle de-jargoned but KEEPS the words the test pins.)
2. **ONE caption** (a handful of plain sentences) for what can't be inlined — the genuinely domain terms
   (β-plane, "reversible" + the %, κ). Seat it in a **deep bottom margin WELL below the play/slider
   controls** (`xref/yref="paper"`, `y` ≈ −0.40, big `b` margin) so it adds explanation without
   re-cluttering the plot. Drive the % off `irreversible_fraction` so caption + subtitle stay in sync.

**Two corrections the user made on the 2026-06-12 follow-up pass (apply going forward):**
- **Formulas & coefficients are WELCOME — they just must be EXPLAINED, not stripped.** The first pass
  over-corrected by deleting `Σ∫|F̄|dt` etc.; the user pushed back. KEEP the formula and gloss it in
  words right there (eddy globe caption: "throughput, Σ∫|F̄|dt — the flux size summed over time and
  latitude, ignoring direction"; "net, Σ|∫F̄dt| — what survives after the poleward/equatorward parts
  cancel"). Plain-language label in the legend, formula+gloss in the caption.
- **Placement & size:** the caption must clear the slider's "… periods" tick labels (they hang to
  ≈ `y=-0.24`) — park it lower (`y≈-0.40`) in plenty of empty space below the figure, not crammed just
  under the controls. And **make the caption text BIGGER** (eddy globe: 11 → **14 px**). Manual `<br>`
  wraps (Plotly doesn't auto-wrap annotations) keep each line ≲108 visible chars so nothing overflows;
  inline `<span style='color:…'>` can tint a term to match its curve.

**NOT "a glossary on every figure."** It's a per-figure judgment call, not a one-line sweep like the
spike toggle. **Skip hover as the definition vehicle** (poor discoverability + it would overload the
lat/lon tooltip deliberately kept). The figure must **stand alone** in the banked standalone HTML, but
**deep pedagogy stays in the notebook** (`planet.ipynb` markdown) — link, don't duplicate; see
[[pedagogy-novice-intermediate]].

**Why:** the banked HTML globes are viewed with no surrounding context; cryptic shorthand reads as
gatekeeping to a novice and as ambiguous to an intermediate.

**How to apply / status:** wired 2026-06-12 on the Rung-B eddy globe (`eddy_globe.eddy_globe_figure`)
as the **exemplar**; pinned by `test_eddy_globe.test_demo_eddy_globe_renders` (asserts a caption
annotation containing both "β-plane" and "reversible"). **The owed sweep is now DONE (2026-06-12):**
`planetmap.render` (banks BOTH `planet-map.html` biome + `planet-coupler-map.html` temp+circulation),
`planetmap.render_comparison` (the A·B·Δ triptych), and the Rung-A GIF `plots.eddy_life_animation`
(`planet-eddy-life.gif`) all carry the relabel+caption pass; each pinned by a caption assertion. The
**six static Phase PNGs** (snowball/biomes/shallowwater/coupler/exoplanet/obliquity) were deliberately
NOT swept — they're the notebook-embedded matplotlib floor (deep pedagogy stays in the notebook, already
richly labelled); revisit only if the user asks. Any NEW Plotly/matplotlib figure must carry the pass.

**Lessons from the sweep (apply going forward):**
- **A globe with no formula = relabel-only.** The biome/temperature globes carry no formula, so the
  "formulas welcome but explained" rule simply has nothing to restore — the caption is pure plain-prose
  (`planetmap._field_caption`). Don't manufacture a formula to gloss; the rule is "keep+explain existing
  formulas", not "add formulas everywhere".
- **Layer-adaptive caption ⇒ sync it on every redraw.** `planetmap.render` builds a caption keyed by the
  active layer; `interactive_map`'s live `update()` synced only `data`+`title`, so a slider layer-switch
  would leave the PREVIOUS field's caption stale under the new globe (untested reach — nothing flags it).
  The advisor caught this; fix is one line: `fig.layout.annotations = new.layout.annotations`.
- **Test-pin the LAYER-INVARIANT edge, not a field word.** Pin "latitude bands" (on the caption for every
  active layer) in the biome-specific test — a "biome" word only holds for the biome render. For the Δ
  triptych pin "differs from"+"latitude bands" (wording that covers BOTH the continuous diff and the
  biome changed-mask — the advisor's catch: don't write a temperature-only "warmer/cooler").
- **`_wrap_html(text, width)` helper** (planetmap): inserts `<br>` at word boundaries counting VISIBLE
  chars only (regex-strips `<b>`/`<span>` tags) so a bold run can straddle a break — replaces the
  hand-placed `<br>`s of the exemplar for the longer/variable-length globe captions.
- **matplotlib caption ≠ Plotly caption.** `fig.text` can't tint inline spans, so the LEGEND's curve
  colours carry the colour cue (not the caption). Use **mathtext `$\bar{F}$`** for the F-bar, NOT the
  combining-macron unicode (fonts render it poorly). Reserve the bottom band for the figtext with the
  constrained-layout engine: `fig.get_layout_engine().set(rect=(0, 0.17, 1, 0.83))` (left,bottom,w,h) —
  verified empirically by rendering a frame to PNG (the slow GIF test only proves it doesn't crash).
- **Drive the % off the figure's OWN `irreversible_fraction`** in every caption — the banked GIF lands on
  irr≈0.081→92/8, matching the globe (truthful-to-its-own-figure; they agree because both sims converge).

Related: [[planet-viz-animation-rungs]], [[planet-interactive-map-design]].
