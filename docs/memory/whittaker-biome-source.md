---
name: whittaker-biome-source
description: "Planet P2 — cited Whittaker/Ricklefs biome-diagram source + the pinned (T,P) thresholds biomes.py uses; the GRAPHICAL→independent-partition (Irvin) decision"
metadata: 
  node_type: memory
  type: reference
  originSessionId: 65378612-4f52-44d4-8842-fe8158fc4aed
---

**Microchip-style pinned source for Planet Phase 2's biome classifier** (`projects/planet/biomes.py`),
the `[[…-source]]` discipline. Part of [[bigsim-program]]; see [[planet-plan]].

**The cited fact (a GRAPHICAL diagram, like [[irvin-sheet-resistance-source]]).** The **Whittaker
biome diagram** — biomes as a partition of the (mean annual temperature, mean annual precipitation)
plane. Primary: **Whittaker, R.H. (1975), *Communities and Ecosystems*** (2nd ed., Macmillan); the
widely-reproduced redraw is **Ricklefs, R.E. (2008), *The Economy of Nature*, Fig 5.5**. Units on the
diagram: **T in °C, P in cm/yr** (NOT mm/yr — the units trap; biomes.py + precip.py both use cm/yr).
A machine-readable digitization exists (`plotbiomes` R package, **MIT-licensed code**, 775 polygon
vertices digitized from Ricklefs Fig 5.5) — useful as a *cross-check reference*, but the figure
itself is copyrighted.

**The decision = the Irvin precedent (do NOT embed the polygons).** Per ARCHITECTURE §6 ("no verbatim
listings/figures") and exactly as Microchip P1a computed R_s independently rather than copying Irvin's
graphical curve: `biomes.py` is an **original, total, rule-based partition** of the (T,P) plane —
boundary thresholds **read off** the diagram and pinned as cited constants, NOT a digitized copy.
Total-by-construction (no nearest-polygon fallback → the "no unclassified gaps / area fractions sum to
1" consistency leg holds). Advisor-confirmed 2026-06-09.

**The pinned thresholds (loose / calibration-dependent — the non-circularity split).** Structure
mirrors the diagram's two physics:
- **Cold biomes = temperature-limited (vertical cuts):** tundra `T < −5 °C`; boreal forest
  `−5 ≤ T < 3 °C` (precip-independent — a named simplification).
- **Warm biomes = moisture-limited (DIAGONAL = precip thresholds LINEAR IN T;** a forest needs more
  rain when warmer): `p_arid = 25 + 0.8T`, `p_semiarid = 50 + 2.0T`, `p_humid = 150 + 2.0T` (cm/yr).
  Temperate band `3 ≤ T < 20`, tropical `T ≥ 20` — the warm T-boundary only **relabels** (temperate
  seasonal forest ↔ tropical savanna), no precip-threshold jump. 9 biomes total.

**Non-circularity (advisor's circularity guard).** Thresholds are *calibrated to* the diagram (loose);
the **probe points** the test asserts are **independent canonical facts** (rain forest = warm+wet,
tundra = cold, subtropical desert = warm+dry…), not points drawn inside lines the module invented.
*Validated tight:* partition totality/determinism + present-Earth band **ordering**. *Loose:* absolute
biome latitudes (depend on the calibrated precip param [[precip-parameterization-source]]).
**Scope edge:** Whittaker (annual T,P) not Köppen (needs seasonal precip the annual-mean v1 lacks).
