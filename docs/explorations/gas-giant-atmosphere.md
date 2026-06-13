# Exploration — simulating a gas-giant atmosphere & its vorticity with these engines

*Feasibility sketch, 2026-06-13. Conversation findings, advisor-reviewed. **Nothing built** — this is
a scoping note, the same "named, not banked" altitude as [[planet-spinout-roadmap]].*

**The question.** How feasible is it to simulate the atmosphere and vorticity of a gas giant
(Jupiter/Saturn — banded zonal jets, the Great Red Spot, polar cyclone clusters) with the engines
this repo already owns, or with off-the-shelf engines?

**The verdict in one line.** The *idealized mechanism* (banded jets + long-lived vortices) is genuinely
feasible — about **one rung** on top of [[planet-rung3-qg-built]]. A recognizably-Jovian *globe* is
**not** (it needs spherical geometry we don't have). The *deep convective interior* is a **heavier,
steeper reach — not out of scope**: it has reduced, laptop-scale entry points (the Busse annulus,
rotating Rayleigh–Bénard), with only the *realistic* anelastic deep-shell-plus-dynamo regime as the
frontier / HPC end.

## Why the match is close

Gas-giant atmospheric dynamics, in the canonical *idealized* literature, is modeled with exactly the
two model families we already have:

- **`planet/baroclinic_qg.py` — two-layer QG on a β-plane** is essentially the Williams-1978 Jovian-
  turbulence model: an inverse energy cascade arrested by β at the **Rhines scale** → zonal jets. The
  published open-source twin is **pyqg**; our engine is a hand-built sibling. It already carries β,
  hyperviscosity, bottom Ekman drag, spectral PV inversion, and 2/3 dealiasing.
- **`engines/fluid` — rotating / reduced-gravity shallow water** is the Great-Red-Spot-as-a-vortex
  model (Marcus 1988; Dowling & Ingersoll 1989). A thin weather layer over a deep quiescent interior is
  precisely the 1.5-/two-layer idealization these engines make — a natural fit, not a stretch.
  See [[shallow-water-source]].

## Three feasibility tiers

### Tier 1 — mechanism on a β-plane (jets + vortices, idealized): feasible, ~one rung

Reuses `baroclinic_qg.py` wholesale. **But the load-bearing correction (advisor):** the saturated
condensate we already banked at rung 3 (KE peak at 0.33 k\*, v′≈16, isotropic pile-up at the box scale)
is the regime where jets **failed** to form — *not* a step toward them. It is an **isotropic** 2-D
condensate: energy ran upscale isotropically to the box and stopped. Zonal jets are the **anisotropic**
alternative — when β is strong enough the cascade is arrested at the Rhines scale and steered into
zonal bands *before* it reaches the box. So "we've shown the inverse cascade ⟹ jets are one knob away"
**overstates it**: the cascade landed in the wrong attractor for jets.

The honest tier-1 build is therefore three specific things, not "turn up β":

1. **Scale separation.** Size the domain and β so `L_Rhines` sits *between* the forcing/deformation
   scale and the box scale (the zonostrophy window). Too-weak β → the condensate we have; too-small box
   → steady waves; jets live in the middle (the zonostrophy index `Rβ` has to land in the jet window).
2. **Anisotropic diagnostics — the current ones are jet-blind.** `TwoLayerQG.ke_spectrum` is azimuthally
   binned (isotropic); it *literally cannot* distinguish a condensate from jets. The deliverable needs
   **zonal-mean ū(y)** showing alternating prograde/retrograde bands, plus a zonal-vs-meridional spectral
   split or a **PV-staircase** diagnostic. A real (modest) addition, not a free byproduct.
3. **Name the forcing.** Current driving is sustained baroclinic instability of a *fixed mean shear* (an
   APE reservoir — the terrestrial / Panetta picture). The gas-giant-idiomatic setup is **small-scale
   stochastic forcing + large-scale drag** (Scott–Polvani forced-dissipative turbulence, representing
   moist-convective storms). Both make jets but claim *different physics* — decide which before building,
   because the honest caption differs.

### Tier 2 — global, sphere-correct: NOT with what we have

Both engines are **doubly-periodic Cartesian β-planes**. The global jet count, the Juno polar-polygon
cyclone clusters, and the prograde equatorial jet all need **spherical geometry** — spherical-harmonic /
cubed-sphere shallow water or QG. That is a genuine new geometry engine (a real multi-step build, not a
tweak). And **equatorial-superrotation *sign* is a known-hard problem** even for the specialists.

### Tier 3 — deep convective interior: a heavier reach, but explorable (not a wall)

The Juno "jets ~3000 km deep" / Busse-column picture is rotating convection in a deep rotating shell —
a *different physics container* (3-D convection, not 2-D layered geostrophy). But the project's own ethos
is **"deferrals, not foreclosures"** (ARCHITECTURE.md §8), and this regime has genuine **reduced,
laptop-scale entry points**:

- **The Busse annulus** (Busse 1976 / 1986) — a *quasi-geostrophic* model of deep convection in a
  rotating annulus whose **sloping ends supply a topographic β**; convective Reynolds stresses drive
  zonal jets. This is **QG-adjacent** — the same rotating-turbulence / Rhines family as tier 1, just in a
  convective container — so it sits *closer* to our existing spectral-QG machinery than "a different
  class" implies. A real, explorable reduced model and the cleanest conceptual bridge: the *deep* jets
  are still Rhines/zonostrophic physics, in a 3-D rotating container.
- **Rotating Rayleigh–Bénard convection** (Boussinesq, rotating plane layer or box) — the canonical
  idealized deep-convection regime (Taylor columns, geostrophic turbulence). Runs in 2-D or a small box
  at modest Ekman & Rayleigh numbers; a tractable hand build or a Dedalus problem.

What stays **frontier / leaves the laptop** is the *realistic* end: full **anelastic (compressible)
convection in a deep spherical shell** at planetary parameters, plus the **MHD dynamo** coupling (the
Juno gravity-inferred deep-jet truncation by magnetic drag) — Gastine / Heimpel / Aurnou / Jones
territory, the **MagIC / Rayleigh** codes. So tier 3 is a *steep* rung with a heavy top, **not out of
scope**: the reduced entry is explorable now; only the realistic deep-shell + dynamo regime is the wall.

## Already-ready engines, if hand-rolling is not the point

| Tool | Tier | Note |
|---|---|---|
| **pyqg** | 1 | Drop-in two-layer QG spectral turbulence; the published version of `baroclinic_qg.py`. |
| **Dedalus** | 2 & 3-reduced | Python spectral PDE framework; shallow water / QG *on the sphere* (tier 2) **and** Boussinesq/anelastic rotating convection incl. spherical shells (the reduced tier-3 entry — rotating Rayleigh–Bénard / Busse-annulus-style). The most reusable single tool here. |
| **EPIC** (Dowling) | 1–2 | The specialist planetary isentropic-coordinate GCM, purpose-built for Jupiter/Saturn. |
| **MagIC / Rayleigh** | 3-realistic | Deep-shell **anelastic** convection + MHD dynamo at planetary parameters — the frontier / HPC end of tier 3. |

## Recommendation

If the goal is to **demonstrate gas-giant jets and vorticity as a mechanism**, that is a clean next rung
on top of `baroclinic_qg.py`: re-tune into the zonostrophic regime, pick and add the forcing, and build
the anisotropic jet diagnostics (ū(y) / Rhines / PV-staircase) that can actually *see* jets — staying at
this project's honest idealized altitude. If the goal is a recognizably-Jovian **globe** (global jet
count, GRS, superrotation), that is a sphere — a new engine — and at that point reaching for
**Dedalus** or **EPIC** likely beats hand-rolling spherical harmonics.

## Sources to pin if this is ever built

Williams 1978 (planetary geostrophic turbulence / Jovian jets); Rhines 1975 (the β-arrest scale); Cho &
Polvani 1996 and Scott & Polvani 2007 (forced-dissipative shallow-water turbulence on the sphere);
Marcus 1988 / Dowling & Ingersoll 1989 (GRS as a shallow-water vortex); Held & Larichev 1996 (the two-
layer QG basis we already cite); Galperin et al. (the zonostrophy index). For tier 3: Busse 1976/1986
(the rotating-annulus QG deep-convection model); Kaspi et al. 2018 (the Juno gravity inversion → ~3000-km
deep jets); Gastine / Heimpel / Aurnou (deep-shell anelastic convection). Extends [[shallow-water-source]].
