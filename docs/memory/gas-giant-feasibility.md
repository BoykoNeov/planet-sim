---
name: gas-giant-feasibility
description: "Feasibility sketch (2026-06-13, NOT built) — simulating a gas-giant atmosphere/vorticity with these engines: three tiers (β-plane mechanism = ~1 rung on baroclinic_qg.py; sphere-correct globe = new geometry engine; deep interior = out of scope), the isotropic-condensate≠jets correction, pyqg/Dedalus/EPIC off-the-shelf. Full record docs/explorations/gas-giant-atmosphere.md"
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
**Deep convective interior** = **out of scope** (anelastic deep-shell convection, MagIC/Rayleigh class).

**The load-bearing correction (advisor):** the rung-3 saturated **condensate** we banked (isotropic KE
pile-up at the box scale) is the regime where jets **FAILED** to form — *not* a step toward them. Jets
are the **anisotropic** alternative (β-arrest at the Rhines scale, *before* the box). So tier 1 is NOT
"turn up β" — it is three things: **scale separation** (`L_Rhines` between forcing and box scales,
zonostrophy window); **anisotropic diagnostics** (current `ke_spectrum` is azimuthally-binned ⟹
**jet-blind** — need zonal-mean ū(y) bands + PV-staircase); and **name the forcing** (current = fixed-
shear baroclinic APE = terrestrial picture; gas-giant-idiomatic = small-scale stochastic forcing +
large-scale drag, Scott–Polvani). Equatorial-superrotation sign = known-hard even for specialists.

Same "named, not banked" altitude as [[planet-spinout-roadmap]]; extends [[shallow-water-source]].
