"""Ocean-fraction knob — how much of the surface is sea → its albedo and heat transport.

The fourth interactive climate knob (alongside the Sun, the greenhouse and the axial tilt). Like
:mod:`planet.obliquity` and :mod:`planet.exoplanet`, it adds **no physics**: it is a **parameter
derivation** that maps a single "what fraction of the planet is ocean" dial onto two parameters the
EBM already accepts — the ice-free planetary albedo :attr:`~planet.albedo.EBMParams.a0` and the
meridional heat-transport coefficient :attr:`~planet.albedo.EBMParams.D`. No engine, no EBM
machinery, no validated radiation constant is touched. Set the ocean fraction to Earth's and the
present-day model is recovered **exactly** (a clean perturbation, asserted in
:mod:`tests.test_ocean`).

The physics — two channels by which "more sea" changes the climate
------------------------------------------------------------------
**1. Albedo (the firm leg).** The sea is dark; land is bright. So a wetter planet reflects less of
the Sun's light and warms. The knob moves the *planetary* (top-of-atmosphere) albedo ``a0``, not the
surface albedo — and the two differ a lot. Most of Earth's ``a0 ≈ 0.30`` is the atmosphere (clouds,
Rayleigh scattering); the surface contributes only ``≈ 0.07`` of it, because a change in surface
brightness is attenuated on its two-way trip through the atmosphere by roughly the square of the
atmospheric transmissivity (Donohoe & Battisti 2011, *J. Climate* **24**, 4402 — [[ocean-albedo-transport-source]]).
Surface albedos: open ocean ``≈ 0.06``, mean land ``≈ 0.25`` (Hartmann, *Global Physical Climatology*),
a surface swing ``≈ 0.19``; attenuated to the top of the atmosphere that is an all-ocean→all-land
**planetary**-albedo swing of ``≈ 0.07`` — the conservative end of the cited surface-contribution
range, used here as :data:`A0_LAND_MINUS_OCEAN`. (The poleward ``a2`` albedo structure
:data:`~planet.ebm.ALBEDO_A2` is held fixed — only the ``a0`` offset moves, the same restraint the
obliquity knob shows with the insolation modes.)

**2. Heat transport (the loose leg, flagged).** Oceans carry heat poleward, so a wetter planet is
more equable (flatter equator-to-pole gradient). The knob nudges the EBM's diffusive ``D`` with ocean
fraction. This is the *soft* mapping: the EBM's ``D`` is a **bulk** meridional-transport coefficient
that on Earth is mostly atmospheric eddy transport, with the ocean carrying a minority share that
peaks near a quarter of the total in the subtropics (Trenberth & Caron 2001, *J. Climate* **14**,
3433). Removing the ocean also removes the moisture that fuels part of the atmosphere's *latent*
transport, so the true sensitivity of the bulk ``D`` to ocean cover conflates several effects and is
not sharply constrained. It is therefore taken as a single modest, monotonic coefficient
(:data:`D_OCEAN_SENSITIVITY`, an all-land world ``≈ −25 %`` in ``D``, an all-ocean world ``≈ +10 %``)
and **flagged as order-of-magnitude, not pinned** — the analogue of the obliquity knob's "loose
bracket". The richer, validated transport story lives in the moist-EBM rungs (:mod:`planet.moist_ebm`).

What the knob deliberately does **not** show (the honest ceiling)
-----------------------------------------------------------------
* **Thermal inertia / the seasons.** The ocean's single biggest real role — its huge heat capacity
  damping the seasonal swing and the day-night cycle — is **invisible here**, because the interactive
  surfaces show the *equilibrium* annual-mean climate, where heat capacity has dropped out. ``a0`` and
  ``D`` are the two channels that survive into the steady state.
* **The rain pattern.** Precipitation in :mod:`planet.precip` is a prescribed latitude pattern times a
  global Clausius–Clapeyron amplitude keyed to mean temperature; it has no explicit ocean source. So a
  wetter world rains a little *more everywhere* (via the warming) but does not move *where* it rains —
  an honest limitation surfaced in the explanation prose, not a modelled ocean–rain coupling.

Units — ocean fraction dimensionless on [0, 1]; ``a0`` dimensionless; ``D`` in W m⁻² K⁻¹.
"""
from __future__ import annotations

from dataclasses import replace

from .albedo import EBMParams
from .ebm import ALBEDO_A0, D_TRANSPORT

# --------------------------------------------------------------------------- #
# Pinned constants ([[ocean-albedo-transport-source]]). Earth's ocean fraction is the ratio anchor
# (the slider default + the value at which both channels are the identity); the two sensitivities are
# the cited mappings — the albedo swing firm, the transport sensitivity flagged loose (see docstring).
# --------------------------------------------------------------------------- #
OCEAN_FRACTION_EARTH = 0.71    # —  present-day Earth sea fraction (~71 %); the anchor + slider default
OCEAN_FRACTION_MIN = 0.0       # —  knob clamp: a bone-dry land world
OCEAN_FRACTION_MAX = 1.0       # —  knob clamp: a global ocean (water world)

# Channel 1 (firm): the all-ocean → all-land swing in the *planetary* albedo a0. Conservative end of
# the Donohoe & Battisti surface-contribution range (a ~0.19 surface swing attenuated to ~0.07 at TOA).
A0_LAND_MINUS_OCEAN = 0.07     # —  a0(land) − a0(ocean), the full 0→1 ocean-fraction range

# Channel 2 (loose, flagged): the fractional change in the bulk transport D per unit ocean fraction,
# applied about the Earth anchor. 0.35 ⇒ an all-land world has ~25 % less D, an all-ocean world ~10 %
# more — modest and monotonic; the magnitude is order-of-magnitude, not pinned (see docstring).
D_OCEAN_SENSITIVITY = 0.35     # —  d(D)/D per unit ocean fraction, about OCEAN_FRACTION_EARTH


def _clamp(ocean_fraction: float) -> float:
    return float(min(max(ocean_fraction, OCEAN_FRACTION_MIN), OCEAN_FRACTION_MAX))


def ocean_albedo_a0(ocean_fraction: float = OCEAN_FRACTION_EARTH,
                    a0_base: float = ALBEDO_A0) -> float:
    """Ice-free planetary albedo ``a0`` at this ocean fraction — channel 1 (the firm leg).

    A linear blend in ocean fraction anchored at Earth: ``a0(w) = a0_base − (w − w_Earth)·ΔA``,
    with ``ΔA = A0_LAND_MINUS_OCEAN``. Earth's fraction returns ``a0_base`` **exactly** (a clean
    perturbation); a wetter world (``w`` up) darkens (``a0`` down, warmer), a drier world brightens.
    """
    return a0_base - (_clamp(ocean_fraction) - OCEAN_FRACTION_EARTH) * A0_LAND_MINUS_OCEAN


def ocean_transport_D(ocean_fraction: float = OCEAN_FRACTION_EARTH,
                      D_base: float = D_TRANSPORT) -> float:
    """Meridional heat-transport coefficient ``D`` at this ocean fraction — channel 2 (the loose leg).

    A modest multiplicative tilt about the Earth anchor: ``D(w) = D_base·(1 + (w − w_Earth)·s)`` with
    ``s = D_OCEAN_SENSITIVITY``. Earth's fraction returns ``D_base`` **exactly**; a wetter world
    transports a little more heat poleward (larger ``D`` → flatter gradient). Magnitude is flagged
    order-of-magnitude (see the module docstring), so callers narrate it as the soft channel.
    """
    return D_base * (1.0 + (_clamp(ocean_fraction) - OCEAN_FRACTION_EARTH) * D_OCEAN_SENSITIVITY)


def ocean_params(ocean_fraction: float = OCEAN_FRACTION_EARTH,
                 base: EBMParams | None = None) -> EBMParams:
    """An :class:`~planet.albedo.EBMParams` for a planet of this ocean fraction.

    Returns ``base`` with **only** its ``a0`` (channel 1) and ``D`` (channel 2) replaced by the
    knob-derived values; every other parameter is untouched. The Earth default
    (``ocean_fraction=OCEAN_FRACTION_EARTH``) returns an ``EBMParams`` **equal to ``base``** — a clean
    perturbation that composes on top of whatever the caller already set (it commutes with the
    obliquity and exoplanet knobs, which replace ``s2`` and ``ai`` respectively). The bundle the
    interactive what-if (:mod:`planet.interactive`) consumes for its fourth axis.
    """
    if base is None:
        base = EBMParams()
    return replace(base,
                   a0=ocean_albedo_a0(ocean_fraction, base.a0),
                   D=ocean_transport_D(ocean_fraction, base.D))
