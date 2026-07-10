"""Lee-side moisture depletion — the rain shadow that drops *below* baseline (Rung 5A.3, plan §12.5).

Rung 5A.2 (:mod:`planet.orographic_scene`) placed the Smith & Barstad rain shadow on the sphere, but
its combination with the zonal-mean climate is **enhancement-only**: the windward slope is lifted above
the baseline and the lee is left *at* the baseline (``P_total = P_zonal + P_orographic``, with
``P_orographic ≥ 0``). That is honest to Smith & Barstad — a condensation-from-forced-ascent diagnostic
that never removes the rained-out water from the air column. But it does **not** reproduce the real
Columbia-Basin desert *behind* the Cascades, whose dryness is not merely "no orographic bonus" but a lee
baseline that has fallen **below** the zonal mean: the windward rainout has drained the passing air, so
the column arriving in the lee carries less moisture than the zonal-mean climate assumes.

This module adds that one missing piece — a **column moisture budget along the wind** — as an *opt-in*
refinement (default off keeps 5A.2's enhancement-only combination exactly; cf. the opt-in moist rate of
:mod:`planet.moist`). It is the smallest closure that turns the dry-shadow into a real rain-shadow
*desert*.

The reduced model — a 1-D along-wind flux budget (the L→∞ limit of the consistent budget)
-----------------------------------------------------------------------------------------
The wind on the patch is the (purely zonal) mid-latitude westerly, so the along-wind axis is longitude
(``x``) and each latitude row is an independent streamline. Track the vertically-integrated moisture
flux ``F = U·W`` (kg m⁻¹ s⁻¹), where ``W`` is the column precipitable water (kg m⁻² ≡ mm). The full
streamline budget with sources is::

    d(U·W)/dx = S − P_total

In the zonal-mean *equilibrium* (no mountain) evaporation balances precipitation, ``S₀ = P_base``. Take
that same source into the mountain case (``S = P_base``, evaporation keeps resupplying at the baseline
rate — the one modelling assumption) and write ``P_total = g·P_base + P_orographic``. The budget for the
dimensionless **depletion factor** ``g = W/W₀`` becomes::

    U·W₀ · dg/dx = P_base·(1 − g) − P_orographic
                   └── refill ──┘   └ depletion ┘

Two consequences that are *forced*, not chosen (advisor-caught):

* The **refill** term ``P_base·(1 − g)`` and the **depletion** term ``−P_orographic`` come from the
  *same* premise — you cannot keep one and drop the other. Their competition has a length scale, the
  **evaporative refill length** ``L = U·W₀ / P_base`` (:func:`refill_length_m`): the distance over which
  evaporation refills a drained column.
* On Earthlike numbers (``U ≈ 15 m/s``, ``W₀ ≈ 3 cm``, ``P_base ≈ 100 cm/yr``) ``L`` is **thousands of
  km** — far larger than the few-hundred-km regional patch. Over the patch the refill term is therefore
  negligible, and the budget collapses to its ``L → ∞`` limit::

    g(x) = 1 − 1/(U·W₀) · ∫  P_orographic(x') dx'          (cumulative, integrated DOWNWIND)

  which is what :func:`depletion_factor` computes. This is **not** "we chose not to model refill"; it is
  the large-``L`` limit of the consistent budget, and ``L`` is a *derived* quantity we assert is ≫ the
  patch (:func:`refill_length_m`; the test pins ``L ≫ patch``). The honest scope note is the *converse*:
  because refill is dropped, the modelled desert does **not** relax back to the zonal mean within the
  patch — a real one does, over ``~L``, far downwind of the window shown.

The drying ratio — the cited calibration target (the loose tier)
----------------------------------------------------------------
The lee depletion is set by one loose-magnitude knob, the incoming column water :data:`PWV_IN_MM`
(``W₀``). Rather than pin ``W₀`` arbitrarily, it is calibrated through the **drying ratio**

    ``DR ≡ ∫ P_orographic dx / (U·W₀) = 1 − g_lee``          (:func:`drying_ratio`)

the fraction of the passing moisture flux that rains out crossing the range — a standard, *observed*
orographic quantity (Smith et al. 2003/2005; Kirshbaum & Smith 2008 report ``DR ≈ 0.3–0.5`` across the
Southern Alps of New Zealand and the coastal ranges of western North America). ``W₀`` is tuned so the
demo ``DR`` lands in that band; the absolute lee dryness stays in the **loose** validation tier, exactly
as the S&B amplitudes and :data:`~planet.orographic_scene.OROGRAPHIC_HOURS_PER_YEAR` do.

Validation triad (plan §3) — what is tight vs loose
---------------------------------------------------
* **Tight (exact).** *Conservation*: the water removed from the flux equals the orographic water rained,
  ``U·W₀·(1 − g_out) = ∫ P_orographic dx`` per streamline — exact by construction of the cumulative sum
  (so ``DR = 1 − g_lee`` identically). *Reduction*: with zero orographic precip ``g ≡ 1``, recovering
  5A.2's enhancement-only combination bit-for-bit (the opt-in is a strict superset).
* **Tight (structural).** ``g`` is **monotone non-increasing downwind** (a flux that only loses water);
  the depletion lands in the **lee, not the windward side** — the integration-direction guard that is
  the exact analogue of Rung 5A's ``sgn(σ)`` branch (integrate the wrong way and the desert appears
  *upwind* of the range).
* **Directional (the payoff).** With depletion on, the lee total drops **below** the zonal-mean baseline
  (``min P_total < P_base``) — a real rain-shadow *desert*, the thing enhancement-only structurally
  cannot produce. Falsifiable: turn depletion off and the minimum returns to the baseline.
* **Loose (magnitude).** The absolute lee dryness, set by :data:`PWV_IN_MM` via the drying ratio.

Honest scope (named, not fixed)
-------------------------------
* **Per-streamline** — valid *only* because the demo wind is purely zonal (``v ≈ 0``); a meridional wind
  component would advect moisture across rows and couple them. Named, not handled.
* **The orographic bonus is not itself depleted** — ``P_orographic`` keeps drawing on Smith & Barstad's
  own reference saturation density while it drains the *baseline*'s column. A fully consistent budget
  would deplete both; simplest-first depletes the baseline and names the caveat.
* **No on-patch refill** — see above: the ``L → ∞`` limit, honest because ``L ≫ patch``, at the cost of a
  desert that does not relax back within the window.

See [[planet-rung5a-orographic]]; the engine is :mod:`planet.orographic`, the 5A.2 scene
:mod:`planet.orographic_scene`. Constants pinned in [[smith-barstad-orographic-source]].
"""
from __future__ import annotations

import numpy as np

from planet.orographic import SECONDS_PER_HOUR

# The one loose-magnitude calibration knob (§3 / advisor-caught): the incoming column precipitable water
# W₀ (kg/m² ≡ mm). Tuned so the demo drying ratio DR = ∫P_oro dx/(U·W₀) lands in the cited observed band
# (~0.3–0.5; Smith et al. 2003/2005, Kirshbaum & Smith 2008). ~30 mm ≈ 3 cm is a typical mid-latitude
# column water vapour; at the Cascades demo it gives DR ≈ 0.46. A calibration band, not a pinned constant.
PWV_IN_MM = 30.0                     # kg/m² ≡ mm — incoming column precipitable water W₀
_MM_PER_CM = 10.0                    # mm per cm
SECONDS_PER_YEAR = 3.156e7           # s/yr — for the refill-length P_base conversion (cm/yr → kg/m²/s)


def depletion_factor(orographic_mm_hr: np.ndarray, dx_m: float, u_m_s: float, *,
                     pwv_in_mm: float = PWV_IN_MM) -> np.ndarray:
    """The along-wind moisture **depletion factor** ``g(y, x) ∈ (0, 1]`` — the fraction of the column
    moisture that remains after the upwind orographic rainout.

    The ``L → ∞`` limit of the consistent streamline budget (module docstring)::

        g(x) = 1 − 1/(U·W₀) · ∫  P_orographic(x') dx'          (cumulative, DOWNWIND)

    computed as a discrete cumulative sum of the **instantaneous** Smith & Barstad rate along the wind.

    Parameters
    ----------
    orographic_mm_hr : ``(n_lat, n_lon)`` orographic precipitation in **mm/hr** — the *instantaneous* S&B
        rate (:func:`planet.orographic.orographic_precip`), **before** any cm/yr annualisation. The
        budget must run in physical (instantaneous) units; the ``OROGRAPHIC_HOURS_PER_YEAR`` accumulation
        knob must **not** enter here (advisor-caught unit trap).
    dx_m : the along-wind (longitude) grid spacing in metres (:func:`planet.orographic_scene.patch_spacings`).
    u_m_s : the along-wind wind component (m/s). Its **sign sets the downwind direction**: ``u > 0``
        (a westerly, +x) integrates west→east; ``u < 0`` (an easterly) integrates east→west. Getting this
        wrong puts the desert *upwind* of the range — the analogue of Rung 5A's ``sgn(σ)`` branch bug.
    pwv_in_mm : the incoming column precipitable water ``W₀`` (kg/m² ≡ mm) — the calibration knob
        (:data:`PWV_IN_MM`).

    Returns ``g`` clamped to ``[0, 1]``. ``U = 0`` (a mountain off the westerly band) or ``W₀ ≤ 0`` gives
    ``g ≡ 1`` (no advection → no depletion). By construction ``g`` is monotone non-increasing downwind.
    """
    p = np.clip(np.asarray(orographic_mm_hr, dtype=float), 0.0, None)
    U = abs(float(u_m_s))
    if U == 0.0 or pwv_in_mm <= 0.0:
        return np.ones_like(p)

    p_si = p / SECONDS_PER_HOUR                          # mm/hr → mm/s ≡ kg/(m²·s)
    # Cumulative water removed from the flux, integrated DOWNWIND along x (axis=1), in kg/(m·s).
    if float(u_m_s) >= 0.0:                              # westerly: downwind = +x (increasing lon index)
        removed = np.cumsum(p_si, axis=1) * dx_m
    else:                                               # easterly: downwind = −x
        removed = np.cumsum(p_si[:, ::-1], axis=1)[:, ::-1] * dx_m
    g = 1.0 - removed / (U * pwv_in_mm)
    return np.clip(g, 0.0, 1.0)


def drying_ratio(orographic_mm_hr: np.ndarray, dx_m: float, u_m_s: float, *,
                 pwv_in_mm: float = PWV_IN_MM) -> float:
    """The **drying ratio** ``DR = ∫ P_orographic dx / (U·W₀)`` — the cited calibration target (~0.3–0.5).

    The fraction of the passing moisture flux that rains out crossing the range, taken as the maximum
    over streamlines (the wettest row, i.e. the ridge crest line). By construction ``DR = 1 − g_lee`` for
    that streamline. Reported **unclamped**, so a value ``> 1`` flags a column that fully dries out
    (``g`` saturating at 0) — a signal the calibration ``W₀`` is too small. See :data:`PWV_IN_MM`.
    """
    p = np.clip(np.asarray(orographic_mm_hr, dtype=float), 0.0, None)
    U = abs(float(u_m_s))
    if U == 0.0 or pwv_in_mm <= 0.0:
        return 0.0
    removed = (p / SECONDS_PER_HOUR).sum(axis=1) * dx_m         # kg/(m·s) rained per streamline
    return float(removed.max() / (U * pwv_in_mm))


def refill_length_m(u_m_s: float, baseline_cm_yr: float, *, pwv_in_mm: float = PWV_IN_MM) -> float:
    """The **evaporative refill length** ``L = U·W₀ / P_base`` (m) — the scale over which evaporation
    refills a drained column.

    This is the length that justifies dropping the refill term: the on-patch budget is the ``L → ∞``
    limit (:func:`depletion_factor`), honest precisely when ``L`` is much larger than the patch. On
    Earthlike numbers (``U ≈ 15 m/s``, ``W₀ ≈ 3 cm``, ``P_base ≈ 100 cm/yr``) ``L`` is ~thousands of km,
    ≫ a few-hundred-km patch. ``baseline_cm_yr`` is the zonal-mean precipitation (cm/yr); ``P_base ≤ 0``
    returns ``∞`` (no sink → no refill demand).
    """
    p_base_si = float(baseline_cm_yr) * _MM_PER_CM / SECONDS_PER_YEAR   # cm/yr → mm/yr ≡ kg/m²/yr → /s
    if p_base_si <= 0.0:
        return float("inf")
    return abs(float(u_m_s)) * pwv_in_mm / p_base_si
