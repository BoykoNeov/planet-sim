---
name: shallow-water-source
description: "Planet P3 — pinned shallow-water constants/relations (Earth Ω/a/g, β-plane f, dispersion, L_R) + the C-grid scheme references that engines/fluid + circulation.py use"
metadata: 
  node_type: memory
  type: reference
  originSessionId: d28b5a93-2ea3-4982-9548-32d186e4dee6
---

**Planet Phase 3 reference facts** pinned at build for the shallow-water engine
(`engines/fluid`) and its planetary consumer (`projects/planet/circulation.py`). Standard
GFD textbook results — Gill 1982 *Atmosphere–Ocean Dynamics*; Vallis 2017 *Atmospheric &
Oceanic Fluid Dynamics*; the C-grid conserving schemes are Sadourny 1975 (JAS, energy *or*
enstrophy) and Arakawa–Lamb 1981 (both). Part of [[bigsim-program]] / [[planet-plan]].

**Earth constants** (universal): Ω = 7.292e-5 rad/s, a = 6.371e6 m, g = 9.81 m/s².
**β-plane Coriolis:** f = f₀ + β(y−y_ref); f₀ = 2Ω sin φ, β = 2Ω cos φ / a. At φ=45°:
f₀ ≈ 1.031e-4 1/s, β ≈ 1.619e-11 1/(m·s).

**Dispersion relations the engine reproduces** (linearized, mean depth H, c=√(gH)):
- gravity wave: ω = c·k; **Poincaré/inertia-gravity:** ω² = f₀² + gH·(k²+l²) (the rotation check).
- **Rossby (β-plane):** ω = −βk/(k²+l²+1/L_R²) — westward, dispersive.
- **deformation radius** L_R = √(gH)/f₀; **geostrophic adjustment** of a height anomaly →
  the Helmholtz-balanced state (1 − L_R²∇²)η_adj = η_init (Rossby's problem, exact at small amp).

**Equivalent depth** H = 1000 m is a **calibration (loose)** chosen so L_R(45°) ≈ 960 km
matches the cited extratropical deformation radius (~1000 km); c = 99 m/s follows. The
qualitative adjustment/Rossby behaviour is what's validated, not the exact H.

**C-grid scheme:** the engine uses a symmetric vector-invariant corner-PV Coriolis flux
(Sadourny-1975 lineage). A single such scheme conserves **energy OR potential enstrophy**,
not both — this realization conserves **energy** semi-discretely (Arakawa–Lamb conserves
both, not built). See [[planet-phase3-engine]] for the engineering detail.
