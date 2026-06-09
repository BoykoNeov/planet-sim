"""Exoplanet knobs — stellar spectrum → ice albedo, and planet size → transport (Planet §9.1).

The §9.1 *exoplanet sandbox*: two knobs that turn the present-day Earth model into an
**other-world** climate. Both were *named* in the plan but **built only now that each one's
source is pinned** — the ``[[…-source]]`` discipline that gates every climate number in this
project. Crucially, neither adds physics: each is a **parameter derivation** that computes an
effective ice albedo ``a_ice`` and an effective meridional transport ``D`` which the existing
EBM already accepts (:class:`~projects.planet.albedo.EBMParams` ``ai`` / ``D``). No engine, no
EBM machinery, and no validated constant is touched — turn both knobs to their solar/Earth
values and the present-day model is recovered **exactly** (the knobs are clean perturbations,
asserted in :mod:`tests.test_exoplanet`).

Knob 1 — stellar spectrum / type → ice-albedo modifier ([[stellar-spectrum-ice-albedo-source]])
-----------------------------------------------------------------------------------------------
Snow and ice are **bright in the visible but dark in the near-IR** (Warren 1982 snow optics). A
cooler, redder star (an M-dwarf at ~3000 K vs the Sun at 5772 K) emits a far larger fraction of
its light in the near-IR, where ice barely reflects — so the **broadband ice albedo is lower**,
the **ice-albedo feedback weakens**, and the planet is **harder to snowball** (Joshi & Haberle
2012; Shields et al. 2013/2014 — a real, *modest* effect, not a runaway). This is captured by a
**two-band** ice albedo — a bright visible band ``a_vis`` and a dark near-IR band ``a_nir`` split
at a crossover wavelength — weighted by the fraction of the stellar blackbody flux in each band
(:func:`blackbody_visible_fraction`). The band values are **pinned empirical** numbers; the knob
itself is applied as the **ratio to the solar value** :func:`stellar_albedo_factor` so the Sun
recovers the climlab ``ai = 0.62`` *exactly* (the ratio is robust to the crossover choice — only
the calibrated band albedos enter, and only through a ratio). Only the bright **ice/snow** albedo
is modified; the dark ice-free ocean/land albedo is left unchanged — the dominant spectral effect
is on ice, and modelling more would over-claim (the named scope edge).

Knob 2 — planet size → meridional transport (a *derivation*; the ``D`` value is [[ebm-radiation-source]])
--------------------------------------------------------------------------------------------------------
In the EBM the transport term is ``D·∂/∂x[(1−x²)∂T/∂x]`` with ``x = sin φ``. That operator is the
spherical Laplacian written in the area coordinate, and each meridional derivative is
``∂/∂(a·φ) = (1/a)·∂/∂(·)`` — so the operator carries ``1/a²`` and the EBM coefficient is
``D = κ/a²`` for a physical (radius-independent) eddy diffusivity ``κ``. Hence a planet of radius
``a`` (in Earth radii) has ``D(a) = D_Earth / a²``: a **bigger planet transports heat less
effectively per unit area → a steeper equator-to-pole gradient**. This is a **derivation**, not a
published scaling law — North 1975/1981 source the *value* of ``D``, not its ``a``-dependence;
that comes from the geometry. The **0-D global mean is identically size-invariant** (``D`` enters
only the mean-preserving transport, never :func:`~projects.planet.ebm.equilibrium_temperature_0d`),
so size sharpens the gradient without *directly* shifting the mean — the clean analytic anchor.
With the ice feedback, the mean *does* drift, but only through the albedo response to the sharpened
gradient (a colder pole grows the ice cap) — a feedback-mediated shift, named as such.

**Rung-0 idealization, named.** ``κ`` is held radius-independent. Size's richer effects route
through **rotation** (the Rossby radius / β-plane — the real circulation lever), which lives in the
Phase-3 fluid engine, not the EBM (plan §5/§9.1) — a different rung, not faked here.

Units — SI/astrophysical (W m⁻², K, °C, Earth radii dimensionless)
------------------------------------------------------------------
Stellar effective temperature ``T_star`` in **K**; planet ``size`` in **Earth radii** (Earth = 1);
albedos dimensionless; ``D`` in **W m⁻² K⁻¹** (the EBM convention). Wavelengths in **µm**.
"""
from __future__ import annotations

from dataclasses import replace

import numpy as np

from .albedo import EBMParams
from .ebm import ALBEDO_A0, ALBEDO_ICE, D_TRANSPORT

# --------------------------------------------------------------------------- #
# Knob 1 — the two-band snow/ice spectral albedo ([[stellar-spectrum-ice-albedo-source]]).
# Pinned EMPIRICAL band values (Warren 1982 snow optics; Joshi & Haberle 2012; Shields 2013/2014).
# Used only through a RATIO to the solar value (the calibrated absolutes cancel — see module docstring).
# --------------------------------------------------------------------------- #
ICE_ALBEDO_VIS = 0.80          # —  snow/ice albedo shortward of the crossover (the bright visible band)
ICE_ALBEDO_NIR = 0.35          # —  snow/ice albedo longward of the crossover (the dark near-IR band)
ICE_BAND_CROSSOVER_UM = 0.7    # µm — the visible/near-IR crossover where snow albedo drops (Warren 1982)

T_SUN = 5772.0                 # K — the Sun's effective temperature (IAU 2015 nominal) — the solar anchor
STAR_TEFF_MIN = 2600.0         # K — knob clamp: a late-M dwarf (cooler → unpinned spectral regime)
STAR_TEFF_MAX = 10000.0        # K — knob clamp: an early-A star

# Main-sequence effective temperatures (Pecaut & Mamajek 2013) — the demo's labelled stellar types.
STELLAR_TYPES = {
    "M5V": 3050.0, "M3V": 3400.0, "M0V": 3850.0, "K5V": 4400.0,
    "G2V (Sun)": T_SUN, "F5V": 6500.0, "A5V": 8000.0,
}

_C2_UM_K = 1.438776877e4        # µm·K — second radiation constant c₂ = hc/k (the Planck-law exponent scale)
_PLANCK_LAMBDA_UM = (0.05, 60.0)   # µm — integration span (well below 0.1 µm to the far-IR tail)
_PLANCK_N = 12000                  # integration nodes — the fraction is stable to <1e-4 here (and ratio'd)


def _planck_lambda(lam_um: np.ndarray, T_star: float) -> np.ndarray:
    """Planck spectral radiance ``∝ B_λ`` (arbitrary units) vs wavelength in µm at temperature ``T_star`` (K).

    ``B_λ ∝ λ⁻⁵ / (exp(c₂/λT) − 1)``. Only relative values matter (every use is a flux *fraction*),
    so the leading constants are dropped. ``np.expm1`` keeps the denominator accurate in the
    Rayleigh–Jeans tail (small exponent) without overflow in the Wien tail (large exponent, float64).
    """
    lam = np.asarray(lam_um, dtype=float)
    return 1.0 / (lam ** 5 * np.expm1(_C2_UM_K / (lam * T_star)))


def blackbody_visible_fraction(T_star: float, crossover_um: float = ICE_BAND_CROSSOVER_UM) -> float:
    """Fraction of a blackbody's radiant power emitted shortward of ``crossover_um`` (the bright-ice band).

    The share of a star's flux that lands in the band where snow/ice is bright. It **increases
    monotonically with ``T_star``** (a hotter, bluer star emits more short-wavelength light) and
    depends only on the product ``crossover·T_star``. Computed by integrating Planck's law over
    wavelength, with the crossover inserted as an exact grid node so the split is clean.
    """
    lam = np.linspace(*_PLANCK_LAMBDA_UM, _PLANCK_N)
    lam = np.unique(np.concatenate([lam, [float(crossover_um)]]))   # crossover as an exact node
    B = _planck_lambda(lam, float(T_star))
    total = np.trapezoid(B, lam)
    below = lam <= crossover_um
    visible = np.trapezoid(B[below], lam[below])
    return float(visible / total)


def two_band_ice_albedo(T_star: float) -> float:
    """Blackbody-weighted broadband snow/ice albedo at stellar temperature ``T_star`` (the two-band model).

    ``a_ice(T★) = a_vis·f_vis(T★) + a_nir·(1 − f_vis(T★))`` with ``f_vis`` the visible-band flux
    fraction (:func:`blackbody_visible_fraction`). For the Sun this is ≈ 0.57 (the canonical broadband
    snow albedo); for an M-dwarf (~3000 K) it falls to ≈ 0.39 — matching Joshi & Haberle 2012. The
    knob uses the **ratio** of this to the solar value (:func:`stellar_albedo_factor`), so this
    absolute is a calibration that cancels.
    """
    f_vis = blackbody_visible_fraction(T_star)
    return ICE_ALBEDO_VIS * f_vis + ICE_ALBEDO_NIR * (1.0 - f_vis)


_ICE_ALBEDO_SUN = two_band_ice_albedo(T_SUN)   # the ratio denominator (the solar anchor); ≈ 0.57


def stellar_albedo_factor(T_star: float) -> float:
    """The ice-albedo modifier ``R(T★) = a_ice(T★) / a_ice(T_Sun)`` — exactly ``1.0`` at the Sun, ``<1`` redder.

    The dimensionless factor the stellar-spectrum knob scales the ice albedo by. Built as a ratio so
    the calibrated band absolutes cancel and the **Sun recovers the model exactly** (``R(T_SUN) == 1``
    by construction). Clamped to the pinned stellar range ``[STAR_TEFF_MIN, STAR_TEFF_MAX]``.
    """
    T = float(np.clip(T_star, STAR_TEFF_MIN, STAR_TEFF_MAX))
    return two_band_ice_albedo(T) / _ICE_ALBEDO_SUN


def stellar_ice_albedo(T_star: float = T_SUN, ai_base: float = ALBEDO_ICE) -> float:
    """Effective ice/snow albedo under a star of effective temperature ``T_star`` (K) — knob 1's output.

    ``a_ice(T★) = ai_base · R(T★)``. The solar anchor ``R(T_SUN) = 1`` makes ``stellar_ice_albedo(T_SUN)``
    return ``ai_base`` (the climlab ``0.62``) **exactly**. A cooler/redder star lowers ``a_ice`` →
    the ice-albedo feedback weakens → the planet is **harder to snowball**. The ``a_nir`` floor keeps
    ``a_ice`` bounded *above* the ice-free ocean/land albedo (:data:`~projects.planet.ebm.ALBEDO_A0`)
    for every stellar type, so the ice/ocean contrast **weakens but never inverts** — *modest*, not a
    sign flip (asserted in the triad). Only this bright-ice albedo is modified (the scope edge).
    """
    return ai_base * stellar_albedo_factor(T_star)


# --------------------------------------------------------------------------- #
# Knob 2 — planet size → meridional transport (a derivation; the D value is [[ebm-radiation-source]]).
# --------------------------------------------------------------------------- #
SIZE_MIN = 0.3                 # Earth radii — knob clamp (a small rocky world)
SIZE_MAX = 3.0                 # Earth radii — knob clamp (a super-Earth)


def size_transport_factor(size: float) -> float:
    """The transport modifier ``1/size²`` (``size`` in Earth radii) — exactly ``1.0`` at Earth, ``<1`` bigger.

    The dimensionless factor knob 2 scales ``D`` by, from the ``1/a²`` the spherical Laplacian carries
    in the ``x = sin φ`` coordinate (see the module docstring). Clamped to ``[SIZE_MIN, SIZE_MAX]``.
    """
    s = float(np.clip(size, SIZE_MIN, SIZE_MAX))
    return 1.0 / (s * s)


def transport_for_size(size: float = 1.0, D_base: float = D_TRANSPORT) -> float:
    """Effective meridional transport ``D(size) = D_base / size²`` for a planet of ``size`` Earth radii — knob 2.

    A larger planet has a smaller ``D`` → a **steeper equator-to-pole gradient** (the two-mode
    amplitude ``|T₂| ∝ 1/(6D + B)`` grows). ``transport_for_size(1.0) == D_base`` **exactly** (Earth
    recovers the model). Because ``D`` is absent from
    :func:`~projects.planet.ebm.equilibrium_temperature_0d`, the **0-D global mean is identically
    size-invariant** — the clean analytic anchor; the relaxed ice-cap mean drifts only through the
    feedback (named in the module docstring).
    """
    return D_base * size_transport_factor(size)


# --------------------------------------------------------------------------- #
# The composition point — both knobs into one EBMParams (the loose-coupling currency).
# --------------------------------------------------------------------------- #
def exoplanet_params(T_star: float = T_SUN, size: float = 1.0, base: EBMParams | None = None) -> EBMParams:
    """An :class:`~projects.planet.albedo.EBMParams` for a world around a ``T_star`` star at ``size`` Earth radii.

    Returns ``base`` with **only** its ice albedo and transport replaced by the knob-derived values
    (``ai = stellar_ice_albedo(T_star, base.ai)``, ``D = transport_for_size(size, base.D)``); every
    other parameter (S₀, A, B, the ice-free albedo, …) is untouched. The solar/Earth defaults
    (``T_star=T_SUN``, ``size=1``) return an ``EBMParams`` **equal to ``base``** — the knobs are clean
    perturbations that compose on top of whatever the caller already set. This is the bundle the demo
    (:mod:`projects.planet.demo_exoplanet`) and the interactive map (:mod:`projects.planet.planetmap`)
    consume.
    """
    if base is None:
        base = EBMParams()
    return replace(base,
                   ai=stellar_ice_albedo(T_star, base.ai),
                   D=transport_for_size(size, base.D))
