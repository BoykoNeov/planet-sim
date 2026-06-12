---
name: no-pointer-spike-lines
description: "User wants NO hover spike/crosshair lines (the pointer-following projection lines) on any visualization — set Plotly showspikes=False (3-D scene axes default it ON) on every figure, now and going forward"
metadata:
  type: feedback
---

**Standing preference (user, 2026-06-12): remove the hover spike/crosshair lines that track the
pointer from EVERY visualization — now and in the future.** In Plotly these are the "spikes": thin
projection lines from the hovered point to the axis planes (3-D scenes) / the axis crosshair (2-D),
which **default ON for 3-D scene axes** (`scene.{x,y,z}axis.showspikes`). The user called them
"interlocking circles that localize with the pointer" and confirmed (single-select) that the spike /
crosshair lines were the target — *not* the modebar toolbar and *not* the hover tooltips (the lat/lon
readout on the eddy band was deliberately kept).

**Why:** the spikes are visual clutter that overlap the figure and chase the mouse; the user wants
clean, static-looking globes/panels (the 3-D scene stays rotatable — only the pointer overlays go).

**How to apply:** in every Plotly figure builder, set `showspikes=False` on the scene axes (and the
2-D cartesian axes). The repo's idiom is the shared `no_axis` dict — add `showspikes=False` beside
`visible=False`. Wired 2026-06-12 in `planetmap.render` + `planetmap.render_comparison` (both
`no_axis` dicts) and `eddy_globe.eddy_globe_figure` (the scene `no_axis` + the 2-D flux-panel
`update_xaxes`/`update_yaxes`); pinned by assertions in `test_planetmap.test_render_builds_the_biome_globe`
and `test_eddy_globe.test_demo_eddy_globe_renders`. Matplotlib figures (`plots.py`, the GIFs) have no
hover, so the rule is Plotly-only. Any NEW Plotly viz must carry it. Related: [[planet-viz-animation-rungs]],
[[planet-interactive-map-design]].
