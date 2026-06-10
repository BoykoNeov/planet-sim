"""Planet §9.1 obliquity-knob triad — axial tilt → the annual-mean-insolation P₂ coefficient s₂.

Mirrors ``test_exoplanet``: the knob is a **parameter derivation**, so the triad is analytic +
structural, not a new conservation law (the EBM it feeds is already sealed by ``test_ebm`` /
``test_albedo``). What is asserted **tight**: the exact ``ε=0 → −5/8`` limit, the closed-form P₂
cross-check across the whole range, the clean-perturbation property, and monotonicity. What is
**loose** (cited / calibration-adjacent): the ``≈−0.48`` climlab cross-check at Earth's tilt and the
relaxed ice-line response. Every test here is cheap (an insolation integral, or one short
equilibration).
"""
import numpy as np
import pytest

from planet import obliquity as ob
from planet.albedo import EBMParams, present_day_climate
from planet.ebm import S2_INSOLATION


# --------------------------------------------------------------------------- #
# The geometric s₂(ε): exact analytic limits + the closed-form cross-check
# --------------------------------------------------------------------------- #
def test_zero_obliquity_gives_exactly_minus_five_eighths():
    # The tight analytic anchor: at ε=0 the orbit is flat (δ≡0), Q̄ ∝ √(1−x²) with no polar-day/night
    # cutoffs anywhere, and the P₂ projection is EXACTLY −5/8. Seals the geometry + normalization in one
    # assertion — the project's signature tight-limit leg.
    assert ob.insolation_p2_coefficient(0.0) == pytest.approx(-0.625, abs=1e-4)


def test_geometric_s2_matches_the_closed_form_across_the_range():
    # The numerical projection of the annual-mean insolation reproduces the known closed-form P₂
    # coefficient s₂(ε) = −(5/8)(1 − 3/2·sin²ε) to the integration's accuracy — an independent,
    # MULTI-point analytic check. The integration is the definition (derived from the pinned
    # daily-insolation formula); the closed form is only the cross-check it is verified against.
    for eps in (0.0, 15.0, 23.44, 45.0, 60.0, 90.0):
        closed = -(5.0 / 8.0) * (1.0 - 1.5 * np.sin(np.radians(eps)) ** 2)
        assert ob.insolation_p2_coefficient(eps) == pytest.approx(closed, abs=5e-4)


def test_earth_obliquity_matches_the_climlab_s2_cross_check():
    # The non-circular benchmark: the geometry independently lands on the climlab/North-1975 fit
    # (S2_INSOLATION = −0.48) at Earth's tilt, to <1% — two independent sources agreeing (the geometry
    # knows nothing of climlab's fitted value).
    assert ob.insolation_p2_coefficient(ob.OBLIQUITY_EARTH) == pytest.approx(S2_INSOLATION, abs=0.005)


def test_s2_rises_monotonically_toward_zero_with_tilt():
    # More tilt → the year's sunlight spreads poleward → s₂ rises (less negative). Monotone over the
    # whole physical range, from exactly −5/8 at no tilt to positive past the critical obliquity.
    s2 = [ob.insolation_p2_coefficient(e) for e in (0, 10, 23.44, 35, 45, 60, 75, 90)]
    assert np.all(np.diff(s2) > 0.0)
    assert s2[0] == pytest.approx(-0.625, abs=1e-4)


def test_gradient_reverses_at_high_obliquity_loose_bracket():
    # The teaching gem (and the worst-truncation regime): past a critical tilt the poles receive MORE
    # annual sun than the equator → s₂ flips sign. Asserted as a LOOSE bracket (negative at 45°, positive
    # by 65°), NOT a pinned 54.7° crossing — the honest claim is the sign change, not its exact location.
    assert ob.insolation_p2_coefficient(45.0) < 0.0
    assert ob.insolation_p2_coefficient(65.0) > 0.0


# --------------------------------------------------------------------------- #
# The knob: clean perturbation, composition, clamping
# --------------------------------------------------------------------------- #
def test_earth_obliquity_recovers_the_model_exactly():
    # The clean-perturbation property: at Earth's tilt the ratio is exactly 1, so the climlab s₂ is
    # recovered bit-for-bit (the present-day map/figure cannot move when the knob sits at its default).
    assert ob.obliquity_s2_factor(ob.OBLIQUITY_EARTH) == 1.0
    assert ob.insolation_s2(ob.OBLIQUITY_EARTH) == S2_INSOLATION


def test_smaller_tilt_steepens_larger_tilt_flattens():
    # The knob's direction: a smaller axial tilt pins sun at the equator (a more negative s₂, steeper
    # gradient); a larger tilt spreads it poleward (a less negative s₂, flatter planet).
    assert ob.insolation_s2(10.0) < S2_INSOLATION < ob.insolation_s2(40.0)


def test_obliquity_params_apply_only_s2():
    base = EBMParams()
    assert ob.obliquity_params(ob.OBLIQUITY_EARTH, base) == base          # clean perturbation: no drift
    p = ob.obliquity_params(40.0)
    assert p.s2 > S2_INSOLATION                                           # higher tilt → less negative s₂
    assert p.S0 == base.S0 and p.A == base.A and p.D == base.D and p.ai == base.ai   # nothing else moved


def test_knob_clamps_outside_the_physical_range():
    # The clamp guards the arccos/tan (no nan below 0° or above 90°): out-of-range values pin to the ends.
    assert ob.obliquity_s2_factor(-5.0) == ob.obliquity_s2_factor(ob.OBLIQUITY_MIN)
    assert ob.obliquity_s2_factor(120.0) == ob.obliquity_s2_factor(ob.OBLIQUITY_MAX)


# --------------------------------------------------------------------------- #
# The relaxed-climate benchmark (loose, mirrors the exoplanet size-knob ice-line test)
# --------------------------------------------------------------------------- #
def test_higher_tilt_pushes_the_ice_line_poleward():
    # Loose qualitative benchmark on the relaxed present-day climate: more tilt → flatter insolation →
    # warmer pole → the ice cap retreats poleward (Earth's 23.44° lands at ~70°, the climlab benchmark;
    # by 40° the cap is gone). The downstream proof the s₂ knob actually moves the climate.
    icelines = [present_day_climate(EBMParams(s2=ob.insolation_s2(e), n_cells=90), n_tau=0.05).ice_line_lat
                for e in (10.0, ob.OBLIQUITY_EARTH, 40.0)]
    assert icelines[0] < icelines[1] < icelines[2]
