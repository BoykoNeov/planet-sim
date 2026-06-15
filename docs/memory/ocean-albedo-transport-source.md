---
name: ocean-albedo-transport-source
description: "Cited sources + chosen magnitudes for the ocean-fraction knob's two channels (planetary albedo a0, transport D) in planet/ocean.py"
metadata:
  type: reference
---

The ocean-fraction knob (`planet/ocean.py`, the 4th interactive what-if axis, built
2026-06-15) maps sea fraction → two EBM params the model already accepts, each pinned
to a source per the project's `[[…-source]]` discipline. Earth's 0.71 is the ratio
anchor (both channels the identity there → bit-for-bit model recovery).

**Channel 1 — planetary albedo `a0` (the FIRM leg).** Donohoe & Battisti 2011,
*J. Climate* **24**, 4402 ("Atmospheric and Surface Contributions to Planetary Albedo"):
most of Earth's `a0≈0.30` is atmosphere (clouds/Rayleigh); the surface contributes only
`≈0.07`, because a surface-albedo change is attenuated on its two-way atmospheric trip by
~the square of the transmissivity. Surface albedos (Hartmann, *Global Physical
Climatology*): open ocean `≈0.06`, mean land `≈0.25` → surface swing `≈0.19`, attenuated
to a **planetary**-albedo swing `≈0.07` all-ocean→all-land. Used as `A0_LAND_MINUS_OCEAN
= 0.07` (the conservative end of the cited range). `a0(w)=ALBEDO_A0−(w−0.71)·0.07`. The
poleward `a2` term is held FIXED (same restraint obliquity shows with the insolation modes).

**Channel 2 — transport `D` (the LOOSE leg, flagged).** Trenberth & Caron 2001,
*J. Climate* **14**, 3433: the ocean carries a MINORITY of poleward heat transport,
peaking near a quarter of the total in the subtropics; the atmosphere dominates. The
EBM's `D` is a *bulk* coefficient (mostly atmospheric eddies), and removing ocean also
removes moisture that fuels part of the atmosphere's latent transport, so `D`'s true
sensitivity to ocean cover conflates effects and is NOT sharply constrained. Taken as one
modest monotonic coefficient `D_OCEAN_SENSITIVITY=0.35` (all-land `≈−25%` D, all-ocean
`≈+10%`), explicitly **order-of-magnitude not pinned** — the analogue of obliquity's loose
bracket. The rigorous transport story lives in [[planet-rung25-mse-diffusion]].

**Named ceilings (honest, surfaced in the explain.py prose + slider hint):** the ocean's
heat capacity / seasonal damping is INVISIBLE (equilibrium annual-mean climate), and the
rain *pattern* doesn't move with the water (precip = prescribed latitude pattern × global
C–C amplitude, no explicit ocean source) — a wetter world rains a little more everywhere,
not somewhere new. Part of [[interactive-what-if]].
