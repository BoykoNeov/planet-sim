---
name: pending-interactive-hover-check
description: "OPEN: the browser play-through of the three what-if UI fixes (knob row / no-reflow / disk hover) is still unverified by the user"
metadata: 
  node_type: memory
  type: project
  originSessionId: e70e5026-9f6c-4a79-be0f-d41a9982163b
---

**OPEN TASK (deferred 2026-06-14, user said "do later"):** the manual browser
play-through of the three what-if fixes in `docs/interactive/index.html` is
still **unverified by the user**. The code is built, committed, and pushed
(`24b99fb`); only the human eyeball check remains.

What to verify when picked up (full numbered steps were given in session
`e56c4c87`):
1. **Knob row** — Sun / Greenhouse / Tilt sliders sit on one aligned row (the
   middle one used to drop); stays aligned while dragging CO2 (value-string
   width changes must not jog a slider).
2. **No reflow** — dragging any slider must not shift the stats/legend/footer.
   **The one decision I owe the user:** is the dark gap below the one-sentence
   explanation *too large* in the common single-knob case? I reserved
   worst-case height (`.oneline { min-height: 9.5rem }` ≈ 6 lines, all three
   knobs moved). If it looks too empty, dial the reservation back to ~4 lines
   (still stable for single/double-knob drags). See `planet/interactive.py`.
3. **Disk hover** — hovering the planet disk shows a tooltip naming the biome
   band + |lat|°, matching the band under the cursor; disappears off the globe.

Apply the [[detailed-test-instructions]] standing rule when re-issuing the
hand-off. Part of [[interactive-what-if]].
