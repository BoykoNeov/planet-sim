---
name: stellar-spectrum-ice-albedo-source
description: "Planet §9.1 exoplanet knob — cited host-SED→ice-albedo feedback (Joshi & Haberle 2012, Shields 2013-14, Warren 1982 snow optics) + the two-band values exoplanet.py pins; the SIZE knob is a derivation, not a source"
metadata: 
  node_type: memory
  type: reference
  originSessionId: 1d49d915-1c78-4254-aa64-b817b4a2940a
---

Source pinned for **Planet §9.1** (`projects/planet/exoplanet.py`, the stellar-spectrum
exoplanet knob): **Joshi & Haberle 2012** (*Astrobiology* 12:3) and **Shields et al.
2013/2014** — around M-dwarfs the host's emission shifts into the near-IR, where snow/ice
are *dark* (Warren 1982 snow optics), so the **ice-albedo feedback weakens** and the
planet is **harder to snowball**. Citeable, modest effect.

**The two-band model `stellar_ice_albedo`:** bright-visible / dark-near-IR ice albedo
weighted by the host's blackbody spectrum — **a_vis = 0.80, a_nir = 0.35, λ_c = 0.7 µm**.
Used only as a **RATIO to the solar broadband value (≈0.57)**, scaling the climlab
`ai = 0.62` ([[ebm-radiation-source]]) so the **Sun recovers ai = 0.62 exactly**.
M-dwarf broadband comes out **≈0.4**, matching Joshi & Haberle. An `a_nir` floor keeps
ice albedo > ocean albedo for every star (never inverts).

**The size knob is a DERIVATION, not a new source:** `transport_for_size` gives
**D ∝ 1/size²** from the spherical Laplacian in x = sinφ; the D value reuses
[[ebm-radiation-source]]. The analytic 0-D mean is size-invariant (size enters only
transport); validation anchors on the constant-albedo two-mode solution (T₂ ∝ 1/(6D+B)
exact), NOT the relaxed ice-cap mean. See [[planet-interactive-map-design]] §9.1.
