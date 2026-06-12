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
2. **ONE compact caption** (2–3 plain sentences) for what can't be inlined — the genuinely domain terms
   (β-plane, "reversible" + the %, κ). Seat it in an **expanded bottom margin BELOW the play/slider
   controls** (`xref/yref="paper"`, `y<0`) so it adds explanation without re-cluttering the plot. Drive
   the % off `irreversible_fraction` so caption + subtitle stay in sync.

**NOT "a glossary on every figure."** It's a per-figure judgment call, not a one-line sweep like the
spike toggle. **Skip hover as the definition vehicle** (poor discoverability + it would overload the
lat/lon tooltip deliberately kept). The figure must **stand alone** in the banked standalone HTML, but
**deep pedagogy stays in the notebook** (`planet.ipynb` markdown) — link, don't duplicate; see
[[pedagogy-novice-intermediate]].

**Why:** the banked HTML globes are viewed with no surrounding context; cryptic shorthand reads as
gatekeeping to a novice and as ambiguous to an intermediate.

**How to apply / status:** wired 2026-06-12 on the Rung-B eddy globe (`eddy_globe.eddy_globe_figure`)
as the **exemplar**; pinned by `test_eddy_globe.test_demo_eddy_globe_renders` (asserts a caption
annotation containing both "β-plane" and "reversible"). **Follow-up sweep still owed** (not done this
batch): `planetmap.render`, `planetmap.render_comparison` (the A·B·Δ triptych), and the matplotlib
GIFs (`plots.py`) each need the same labels+caption pass. Any NEW Plotly/matplotlib figure must carry
it from the start. Related: [[planet-viz-animation-rungs]], [[planet-interactive-map-design]].
