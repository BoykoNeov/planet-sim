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
   **RESOLVED 2026-06-15 (`c9580ee`):** the gap-size decision is moot — there is
   no `min-height` to dial any more. The reserved-height approach was fragile
   both ways (too tall for the common single-knob case, too short for the worst
   three-knob/Snowball one — which still overflowed to 7–8 lines and shoved the
   legend) AND the `<details>` "Why" had no reserved space so expanding it always
   pushed the legend. Fix = **decouple, don't reserve**: stats+legend moved into a
   new `.vizcol` under the disk/curve; `.panel` is now pure prose with `<details>`
   last; `min-height` removed. A growing one-liner or an expanded "Why" now pushes
   only the footer. Still owed: a browser eyeball that the legend truly holds still
   while dragging and while toggling "Why".
3. **Disk hover** — hovering the planet disk shows a tooltip naming the biome
   band + |lat|°, matching the band under the cursor; disappears off the globe.

**Disk band-latitude BUG found + fixed 2026-06-14 (`7fc976b`, pushed):** the
user spotted boreal forest painted *poleward* of the ice-line ring (real cells:
S0=1265, CO2=+2..+4, tilt=35; and S0=1265, tilt=23, CO2=+0..+9 — dim/cold
worlds, low ice line). Root cause = the disk colored pixels with a **linear
index** `k=round(|lat|/latMax·(N−1))`, but the model grid is **equal-area
(uniform in sin φ)** so stored latitudes are non-uniform (~2° steps at equator →
~6° at pole, last band only ~76°). The linear map read a too-warm band and
dragged warm bands poleward, while the ice ring is drawn at *true* latitude →
forest on the ice cap. (The temperature curve plots at true lat = was correct,
so curve & disk disagreed — the tell.) Fix = shared `bandForLat(phi)`
nearest-latitude search used by BOTH the disk fill and the hover (the hover had
the **same** broken formula → it was mis-naming bands near the pole too, so this
also fixes verification item 3). `drawCurve` untouched. Scoped to the index map
(a residual ~1-cell band-edge-vs-ring gap is expected). Guarded by a string test
+ a render-path replay over the disk pixels for the flagged cells (no warm biome
poleward of the ring; boreal entirely equatorward). **Page regenerated** — so
the browser play-through now ALSO needs a re-eyeball of the disk bands vs the
ring on a dim/tilted cell (e.g. 1265/+2/35): boreal must stop equatorward of the
white ring, only tundra beyond it. NOTE: commit subject has a stray leading
"`@ `" (PowerShell here-string leaked into the Bash tool); body is clean;
cleaning it needs a force-push the user must authorize.

**Tilt-cooling "Why" prose added 2026-06-15 (`c9580ee`, same commit):** user asked
why, past a threshold, *more* tilt *cools* the planet. Mechanism (validated to the
digit): tilt only **redistributes** the fixed `S₀/4` — it never changes the total —
onto the model's intrinsically brighter high latitudes (the fixed `a₂` poleward-albedo
term), so once the ice cap is gone (the ice-albedo *warming* amplifier off) more tilt
absorbs slightly *less* → the closed form `Δ⟨T⟩ = −(S₀/20)·a₂·Δs₂/B` (e.g. ice-free
1465/+16: 0→45° tilt drops ⟨T⟩ 33.0→31.6 °C, matches `−1.349` exactly). With ice
present the sign is OPPOSITE (more tilt melts ice → warms). Baked into `explain.py`'s
paragraph via `_obliquity_cooling_note`, **gated ice-free + tilt-raised** (4 regime
tests in `test_explain.py`). Re-eyeball owed: drag tilt up on a warm/ice-free cell and
read the "Why" paragraph. Heads-up surfaced to user: the lookup cell `1265/+2/40°`
snowballs (⟨T⟩≈−44) while 35°/45° next to it don't — a bistable cell flipping basins at
the Snowball cliff, **not a new bug** (don't chase).

Apply the [[detailed-test-instructions]] standing rule when re-issuing the
hand-off. Part of [[interactive-what-if]].
