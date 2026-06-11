---
name: obliquity-insolation-source
description: "Planet §9.1: cited daily-mean-insolation formula + the mean-annual Legendre s₂(ε) derivation that obliquity.py uses (the last deferred exoplanet knob, now wired)"
metadata:
  node_type: memory
  type: reference
  originSessionId: planet-obliquity-knob
---

The **obliquity → annual-mean-insolation `s₂(ε)`** derivation BigSim's Planet §9.1 uses
(`projects/planet/obliquity.py`) — the **last deferred exoplanet knob, now wired**. Axial
tilt `ε` feeds the EBM's insolation P₂ coefficient `EBMParams.s2`; a pure parameter-feed
(no engine, no EBM machinery, no gate-manifest change — numpy-only, planet-local).

**Derived from first principles, NOT a memorized coefficient:** integrate the pinned
**daily-mean-insolation formula** over a circular-orbit year and project onto P₂. Applied
as the **ratio to the Earth value** (mirroring [[stellar-spectrum-ice-albedo-source]]'s
`stellar_albedo_factor`), so Earth's tilt recovers the climlab `s₂ = −0.48` from
[[ebm-radiation-source]] **bit-for-bit** (a clean perturbation).

**Source (pinned at build, 2026-06-10):** the daily-mean-insolation formula = **Hartmann,
*Global Physical Climatology* §2.7 / Berger 1978 / Rose's climlab notes**; the mean-annual
Legendre context = **Nadeau & McGehee, *Icarus* 2017 (arXiv:1510.04542) / North 1975**.
ε_Earth ≈ **23.44°**.

**Validated facts (the independent cross-checks):**
- **s₂(0) = −5/8 exactly** (flat orbit, Q̄ ∝ √(1−x²), no polar cutoffs) — the TIGHT analytic anchor.
- The numerical projection reproduces the known closed form **s₂(ε) = −(5/8)(1 − 1.5·sin²ε)**
  across the whole range (the closed form is the cross-check; the integration is the definition).
- Non-circular: geometry independently lands on climlab's **−0.48 at 23.44°** (s₂≈−0.477, <1%).
- Relaxed-climate benchmark: more tilt → ice line poleward (Earth 23.44° → ~71°, the climlab
  ~70° benchmark; ice-free by ~40°). **Sign reversal ≈ 54.7°** (poles warmer) — surfaced as a
  loose bracket, NOT a pinned crossing.

**Scope edge (the `a₂`-fixed treatment, named in the docstring):** only the **insolation s₂**
responds to tilt — the ice-free albedo's poleward `a₂` structure is held fixed (the standard
EBM treatment; North 1975 varies the insolation with obliquity, not the albedo's zenith-angle
dependence). Single-P₂-mode truncation degrades at high tilt (real annual-mean grows an s₄
term); annual-mean only (no seasonal extremes); eccentricity/precession = a separate deferred
axis. Wired into `climate_view` + the interactive-map slider (`OBLIQUITY_FAITHFUL_MAX = 45°`,
a UI/scope cap, not a pinned physics number). See [[planet-interactive-map-design]].
