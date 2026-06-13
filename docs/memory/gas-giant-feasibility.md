---
name: gas-giant-feasibility
description: "Feasibility sketch (2026-06-13, NOT built) — simulating a gas-giant atmosphere/vorticity with these engines: three tiers (β-plane mechanism = ~1 rung on baroclinic_qg.py; sphere-correct globe = new geometry engine; deep interior = a steeper reach, NOT out of scope — Busse-annulus QG + rotating Rayleigh–Bénard are reduced entries), the isotropic-condensate≠jets correction, pyqg/Dedalus/EPIC/MagIC off-the-shelf. Full record docs/explorations/gas-giant-atmosphere.md"
metadata:
  type: project
---

A **scoping note, not code** (2026-06-13): can these engines simulate a gas-giant atmosphere + its
vorticity? Full record: `docs/explorations/gas-giant-atmosphere.md`.

**Verdict — three tiers.** (1) **β-plane mechanism** (banded jets + vortices, idealized) = feasible,
**~one rung on [[planet-rung3-qg-built]]** — two-layer QG β-plane *is* the Williams-1978 Jovian-jet
model (pyqg = its published twin); shallow water *is* the GRS-vortex model. (2) **Sphere-correct globe**
(global jet count, polar polygons, equatorial superrotation) = **NOT** with what we have — both engines
are doubly-periodic Cartesian β-planes; needs a real new spherical-geometry engine (Dedalus / EPIC). (3)
**Deep convective interior** = a **steeper reach, NOT out of scope** (user correction 2026-06-13 — aligns
with ARCHITECTURE.md §8 "deferrals, not foreclosures"): reduced **laptop-scale entry points** exist — the
**Busse annulus** (a *QG* model of deep convection; sloping ends = topographic β ⟹ same rotating-
turbulence / Rhines family as tier 1, so QG-adjacent to our machinery — the cleanest bridge) and
**rotating Rayleigh–Bénard** (Boussinesq, small box; a Dedalus problem). Only the *realistic* **anelastic
deep-shell + MHD dynamo** at planetary parameters is the frontier / HPC wall (MagIC / Rayleigh; Gastine /
Heimpel / Aurnou; Kaspi-2018 Juno deep-jet inversion).

**The load-bearing correction (advisor):** the rung-3 saturated **condensate** we banked (isotropic KE
pile-up at the box scale) is the regime where jets **FAILED** to form — *not* a step toward them. Jets
are the **anisotropic** alternative (β-arrest at the Rhines scale, *before* the box). So tier 1 is NOT
"turn up β" — it is three things: **scale separation** (`L_Rhines` between forcing and box scales,
zonostrophy window); **anisotropic diagnostics** (current `ke_spectrum` is azimuthally-binned ⟹
**jet-blind** — need zonal-mean ū(y) bands + PV-staircase); and **name the forcing** (current = fixed-
shear baroclinic APE = terrestrial picture; gas-giant-idiomatic = small-scale stochastic forcing +
large-scale drag, Scott–Polvani). Equatorial-superrotation sign = known-hard even for specialists.

Same "named, not banked" altitude as [[planet-spinout-roadmap]]; extends [[shallow-water-source]].
