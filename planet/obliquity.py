"""Obliquity knob — axial tilt → the annual-mean-insolation gradient (Planet §9.1).

The §9.1 plan named **obliquity** as a core climate knob but left it *deferred* — a disabled
slider — until its ``s₂(obliquity)`` relation could be **pinned to a source**, the
``[[…-source]]`` discipline that gates every climate number in this project
([[obliquity-insolation-source]]). This module supplies it, and like the two exoplanet knobs
(:mod:`projects.planet.exoplanet`) it adds **no physics**: it is a **parameter derivation** that
computes the insolation P₂ coefficient ``s₂`` the EBM already accepts
(:attr:`~projects.planet.albedo.EBMParams.s2`). No engine, no EBM machinery, no validated constant
is touched — set the obliquity to Earth's and the present-day model is recovered **exactly** (a clean
perturbation, asserted in :mod:`tests.test_obliquity`).

The physics — axial tilt sets how evenly the year's sunlight is spread in latitude
------------------------------------------------------------------------------------
The EBM's insolation is the single-mode ``S(x) = (S₀/4)(1 + s₂·P₂(x))`` (:func:`~projects.planet.ebm.insolation`,
``x = sin φ``), with ``s₂ < 0`` concentrating sunlight at the equator. That ``s₂`` is **not** a free
constant — it is the second-Legendre coefficient of the **annual-mean insolation**, which is fixed by
the planet's **obliquity** ``ε`` (the tilt of the spin axis to the orbital plane). This module computes
it from first principles rather than pinning a memorized number:

1. **Daily-mean insolation** (the one pinned formula — Hartmann *Global Physical Climatology* §2.7;
   Berger 1978; Rose's climlab notes — [[obliquity-insolation-source]]):

       Q̄_day(φ, δ) ∝ H₀·sin φ·sin δ + cos φ·cos δ·sin H₀,      cos H₀ = −tan φ·tan δ

   with ``H₀`` the sunset hour angle, clamped to the polar limits — **polar day** (sun never sets,
   ``H₀ = π``) where ``−tan φ tan δ ≤ −1``, **polar night** (``H₀ = 0``) where ``≥ 1`` (the
   ``np.clip`` of the ``arccos`` argument does both). The leading ``S₀/π`` and the eccentricity
   distance factor are dropped — every use here is a *ratio*, so they cancel (which is also why the
   **circular-orbit** idealization is honest: eccentricity rescales the annual mean by a global
   ``1/√(1−e²)`` that divides out, and eccentricity/precession are a *separate* deferred Milankovitch
   axis, not this knob).
2. **Annual mean.** Over a circular-orbit year the solar longitude ``λ`` advances uniformly and the
   declination is ``sin δ = sin ε · sin λ``; averaging ``Q̄_day`` over ``λ ∈ [0, 2π)`` gives the
   annual-mean insolation ``Q̄(x; ε)`` (:func:`annual_mean_insolation`).
3. **Projection onto P₂.** With the equator-symmetric ``Q̄`` expanded in even Legendre modes
   ``Q̄ ∝ 1 + s₂P₂ + s₄P₄ + …``, the coefficient is ``s₂ = 5∫₀¹Q̄P₂dx / ∫₀¹Q̄dx``
   (:func:`insolation_p2_coefficient`) — the geometric/raw value, the analogue of
   :func:`~projects.planet.exoplanet.two_band_ice_albedo`.

As ``ε`` grows the summer sun reaches higher latitudes (even over the pole), so the annual sunlight
spreads out: ``s₂`` rises **toward zero**, and past a **critical obliquity ≈ 55°** it goes **positive**
— the poles then receive *more* annual-mean sunlight than the equator (the high-obliquity world, e.g.
Uranus). The exact analytic limits seal the construction: at ``ε = 0`` the orbit is flat, ``δ ≡ 0``,
``Q̄ ∝ √(1−x²)``, and the projection is **exactly ``s₂ = −5/8 = −0.625``** (no polar cutoffs anywhere);
at Earth's ``ε ≈ 23.44°`` the geometry independently lands on ``s₂ ≈ −0.48`` — the climlab/North-1975
fit the EBM already pins (:data:`~projects.planet.ebm.S2_INSOLATION`), a **non-circular cross-check**
(the two come from independent places and agree to <1%).

The knob — a ratio anchor, so Earth recovers the model exactly
--------------------------------------------------------------
The geometric ``s₂(ε)`` lands near but not exactly on the pinned climlab ``−0.48`` at Earth's tilt, so
the knob is applied as the **ratio to the Earth value** (:func:`obliquity_s2_factor`), exactly as the
stellar-albedo knob ratios to the solar value — the **Sun/Earth defaults recover the model bit-for-bit**
(``factor(ε_Earth) == 1`` by construction → :func:`insolation_s2` returns
:data:`~projects.planet.ebm.S2_INSOLATION` exactly). :func:`obliquity_params` composes it onto an
:class:`~projects.planet.albedo.EBMParams`, replacing only ``s2``.

Scope edge, named
-----------------
The knob feeds the EBM's **single P₂ mode only** — exactly consistent with the model's existing
one-mode ``insolation()``, so it is the honest extension, but the *real* annual-mean insolation grows a
significant ``s₄`` as ``ε`` rises (the latitude profile flattens and then humps at the poles), which a
P₂ truncation cannot carry. So the knob is most faithful near present tilt and **progressively
truncates at high obliquity** — the analogue of the EBM's "linear OLR accurate only near the present
climate". The sign reversal is real and surfaced, but asserted only as a **loose bracket** (negative at
45°, positive by 65°), not at a pinned crossing. Likewise **only the insolation ``s₂`` responds to
tilt — the ice-free albedo's poleward ``a₂`` structure (:data:`~projects.planet.ebm.ALBEDO_A2`) is held
fixed** (the standard EBM treatment: North 1975 varies the *insolation* with obliquity, not the albedo's
zenith-angle dependence). Combined with the **annual-mean** model (no seasonal extremes — the very thing
high obliquity makes dramatic), this is the named ceiling. Eccentricity and precession are a separate
deferred Milankovitch axis (above).

Units — SI/astronomical (W m⁻² implied, °C, x = sin φ dimensionless; obliquity in degrees)
------------------------------------------------------------------------------------------
Obliquity ``ε`` in **degrees** (Earth ≈ 23.44°); latitude coordinate ``x = sin φ`` dimensionless on
[0, 1]; ``s₂`` dimensionless. The insolation values are *relative* (every result is a ratio).
"""
from __future__ import annotations

import numpy as np

from .albedo import EBMParams
from .ebm import S2_INSOLATION, legendre_P2

# --------------------------------------------------------------------------- #
# Pinned constants ([[obliquity-insolation-source]]).
# The only pinned number is Earth's present obliquity; the s₂(ε) shape is COMPUTED (the
# daily-insolation formula is the pinned relation, validated by the exact −5/8 limit and the
# independent ≈−0.48 cross-check, NOT a memorized coefficient).
# --------------------------------------------------------------------------- #
OBLIQUITY_EARTH = 23.44        # °  — present-day Earth obliquity (the ratio anchor + the slider default)
OBLIQUITY_MIN = 0.0            # °  — knob clamp: an untilted world (sunlight pinned at the equator, s₂ = −5/8)
OBLIQUITY_MAX = 90.0           # °  — knob clamp: a pole-on world (the physics is defined to here)

_N_LAMBDA = 720                # year samples (0.5° of solar longitude) — the annual-mean integral
_N_PHI = 721                   # latitude samples (0.25°, equator→pole) — the P₂ projection integral


def annual_mean_insolation(phi: np.ndarray, obliquity_deg: float, n_lambda: int = _N_LAMBDA) -> np.ndarray:
    """Annual-mean insolation ``Q̄(φ; ε)`` (relative units) at latitudes ``phi`` (radians) for tilt ``ε``.

    Averages the daily-mean insolation ``H₀·sinφ·sinδ + cosφ·cosδ·sinH₀`` (the pinned formula, leading
    ``S₀/π`` dropped) over a circular-orbit year — ``sin δ = sin ε·sin λ`` with ``λ`` uniform on
    ``[0, 2π)``. The sunset hour angle ``H₀ = arccos(clip(−tanφ·tanδ, −1, 1))`` reduces to ``π`` in
    polar day and ``0`` in polar night through the clip, so the formula needs no branching. Returns one
    value per latitude (the year-average), the field :func:`insolation_p2_coefficient` projects onto P₂.
    """
    eps = np.radians(float(obliquity_deg))
    phi = np.asarray(phi, dtype=float)
    lam = np.linspace(0.0, 2.0 * np.pi, n_lambda, endpoint=False)
    sin_delta = np.sin(eps) * np.sin(lam)
    delta = np.arcsin(np.clip(sin_delta, -1.0, 1.0))
    # Broadcast latitude (rows) × solar longitude (cols); tan(π/2) is large-but-finite in float and
    # the arccos-argument clip turns it into the correct polar-day/-night limit.
    with np.errstate(invalid="ignore"):
        arg = -np.tan(phi)[:, None] * np.tan(delta)[None, :]
    H0 = np.arccos(np.clip(arg, -1.0, 1.0))
    daily = (H0 * np.sin(phi)[:, None] * sin_delta[None, :]
             + np.cos(phi)[:, None] * np.cos(delta)[None, :] * np.sin(H0))
    return daily.mean(axis=1)


def insolation_p2_coefficient(obliquity_deg: float, n_phi: int = _N_PHI) -> float:
    """The geometric/raw P₂ coefficient ``s₂(ε)`` of the annual-mean insolation (the un-anchored value).

    ``s₂ = 5∫₀¹Q̄·P₂dx / ∫₀¹Q̄dx`` with ``x = sin φ`` (equal-area), computed by integrating the
    annual-mean insolation (:func:`annual_mean_insolation`) over latitude. The integral is taken in
    ``φ`` (so ``dx = cos φ·dφ`` — a smooth weight that vanishes at the pole, avoiding the ``√`` cusp a
    direct ``x``-quadrature would hit at ``x = 1``). Exact analytic value at ``ε = 0``: **``−5/8``**
    (the tight test anchor); ``≈ −0.477`` at Earth's tilt (the climlab cross-check); **positive** above
    the ≈55° critical obliquity. This is the analogue of
    :func:`~projects.planet.exoplanet.two_band_ice_albedo` — the raw quantity the knob then ratios.
    """
    phi = np.linspace(0.0, 0.5 * np.pi, n_phi)
    Qbar = annual_mean_insolation(phi, obliquity_deg)
    weight = np.cos(phi)                                  # dx = cos φ dφ (the equal-area Jacobian)
    P2 = legendre_P2(np.sin(phi))
    num = np.trapezoid(Qbar * P2 * weight, phi)
    den = np.trapezoid(Qbar * weight, phi)
    return float(5.0 * num / den)


_S2_RAW_EARTH = insolation_p2_coefficient(OBLIQUITY_EARTH)   # the ratio denominator (the Earth anchor); ≈ −0.477


def obliquity_s2_factor(obliquity_deg: float) -> float:
    """The insolation modifier ``F(ε) = s₂(ε) / s₂(ε_Earth)`` — exactly ``1.0`` at Earth's obliquity.

    The dimensionless factor the obliquity knob scales the pinned ``s₂`` by. Built as a ratio so the
    geometry's small offset from the climlab fit cancels and **Earth's tilt recovers the model exactly**
    (``F(OBLIQUITY_EARTH) == 1`` by construction). ``F`` *rises* with obliquity (toward 0 and then
    negative — sunlight spreads poleward), and goes negative once ``s₂(ε)`` flips sign past the critical
    obliquity. Clamped to ``[OBLIQUITY_MIN, OBLIQUITY_MAX]``.
    """
    eps = float(np.clip(obliquity_deg, OBLIQUITY_MIN, OBLIQUITY_MAX))
    return insolation_p2_coefficient(eps) / _S2_RAW_EARTH


def insolation_s2(obliquity_deg: float = OBLIQUITY_EARTH, s2_base: float = S2_INSOLATION) -> float:
    """Effective insolation P₂ coefficient ``s₂`` at axial tilt ``ε`` (degrees) — the knob's output.

    ``s₂(ε) = s2_base · F(ε)``. The Earth anchor ``F(OBLIQUITY_EARTH) = 1`` makes
    ``insolation_s2(OBLIQUITY_EARTH)`` return ``s2_base`` (the climlab ``−0.48``) **exactly** — the
    clean-perturbation property. A *smaller* tilt → a *more negative* ``s₂`` (sun pinned at the equator,
    a steeper gradient, a colder pole); a *larger* tilt → ``s₂`` toward 0 and beyond (a flatter, then
    pole-warm planet). Feeds straight into :attr:`~projects.planet.albedo.EBMParams.s2`.
    """
    return s2_base * obliquity_s2_factor(obliquity_deg)


def obliquity_params(obliquity_deg: float = OBLIQUITY_EARTH, base: EBMParams | None = None) -> EBMParams:
    """An :class:`~projects.planet.albedo.EBMParams` for a planet of axial tilt ``ε`` (degrees).

    Returns ``base`` with **only** its insolation coefficient replaced by the knob-derived value
    (``s2 = insolation_s2(ε, base.s2)``); every other parameter is untouched. The Earth default
    (``obliquity_deg=OBLIQUITY_EARTH``) returns an ``EBMParams`` **equal to ``base``** — a clean
    perturbation, composing on top of whatever the caller already set (it commutes with the exoplanet
    knobs, which replace ``ai``/``D`` only). The bundle the demo (:mod:`projects.planet.demo_obliquity`)
    and the interactive map (:mod:`projects.planet.planetmap`) consume.
    """
    from dataclasses import replace
    if base is None:
        base = EBMParams()
    return replace(base, s2=insolation_s2(obliquity_deg, base.s2))
