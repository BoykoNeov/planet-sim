"""Planet §9.1 exoplanet-knob triad — stellar spectrum → ice albedo, planet size → transport.

The two knobs are **parameter derivations**, so the triad is structural + analytic, not a new
conservation law (the EBM machinery they feed is already sealed by ``test_ebm`` / ``test_albedo``).
What is asserted **tight**: the clean-perturbation property (solar/Earth recovers the model exactly),
the exact transport/two-mode scaling, and the bounded "modest, never inverts" property. What is
**loose** (calibration-dependent, cited): the absolute M-dwarf ice albedo (≈0.4, Joshi & Haberle 2012)
and the snowball-resistance threshold. The expensive end-to-end demo is ``slow`` in
``test_demo_exoplanet``; here every test is one or two equilibrations at most.
"""
import numpy as np
import pytest

from projects.planet import exoplanet as ex
from projects.planet.albedo import EBMParams
from projects.planet.ebm import (
    ALBEDO_A0, ALBEDO_ICE, B_OLR, D_TRANSPORT, S0_EARTH, S2_INSOLATION,
    equilibrium_temperature_0d, two_mode_solution,
)


# --------------------------------------------------------------------------- #
# Knob 1 — stellar spectrum → ice albedo
# --------------------------------------------------------------------------- #
def test_solar_star_recovers_the_model_exactly():
    # The clean-perturbation property: at the Sun the ratio is exactly 1, so the climlab ai is recovered
    # bit-for-bit (the map/figure cannot move when the knob sits at its default).
    assert ex.stellar_albedo_factor(ex.T_SUN) == 1.0
    assert ex.stellar_ice_albedo(ex.T_SUN) == ALBEDO_ICE


def test_visible_fraction_is_a_monotone_fraction_of_temperature():
    fracs = [ex.blackbody_visible_fraction(T) for T in (2800, 3500, 4500, 5772, 7000, 9000)]
    assert all(0.0 < f < 1.0 for f in fracs)
    assert np.all(np.diff(fracs) > 0.0)          # hotter, bluer star → more visible-band flux


def test_cooler_star_lowers_the_ice_albedo_monotonically():
    ais = [ex.stellar_ice_albedo(T) for T in (3000, 3800, 4400, ex.T_SUN, 6500, 8000)]
    assert np.all(np.diff(ais) > 0.0)            # redder (cooler) → lower ice albedo
    assert ex.stellar_ice_albedo(3000) < ALBEDO_ICE < ex.stellar_ice_albedo(8000)


def test_mdwarf_broadband_ice_albedo_matches_the_literature_band():
    # Loose benchmark (Joshi & Haberle 2012): an M-dwarf's broadband snow/ice albedo ≈ 0.4 (vs ≈0.57
    # solar). The two-band model lands there from independently-pinned visible/near-IR band values.
    assert 0.35 < ex.two_band_ice_albedo(3000.0) < 0.45
    assert ex.two_band_ice_albedo(ex.T_SUN) == pytest.approx(0.57, abs=0.02)


def test_ice_albedo_stays_above_the_ice_free_albedo_for_every_star():
    # The bounded "modest, never inverts" property (the a_nir floor): for ANY stellar type the ice
    # albedo stays above the ice-free ocean/land albedo, so the ice/ocean contrast weakens but never
    # flips sign — the effect is "harder to snowball", not "ice darker than ocean".
    for T in np.linspace(ex.STAR_TEFF_MIN, ex.STAR_TEFF_MAX, 25):
        assert ex.stellar_ice_albedo(T) > ALBEDO_A0


def test_redder_star_resists_snowball_where_the_sun_freezes():
    # The cheap, decisive proxy (advisor): dim the present (temperate, finite-cap) planet to S₀=1200.
    # The Sun-like planet snowballs (ice to the equator); the SAME planet around an M-dwarf — weaker
    # ice feedback — keeps a temperate cap. One equilibration each (the finite-cap start, comfortably
    # off the cold-start bifurcation edge); the full hysteresis-shift is the slow demo.
    from projects.planet.albedo import present_day_climate
    S0 = 1200.0
    sun_st = present_day_climate(EBMParams(S0=S0, ai=ALBEDO_ICE, n_cells=90), n_tau=0.05)
    mdwarf_st = present_day_climate(EBMParams(S0=S0, ai=ex.stellar_ice_albedo(3000.0), n_cells=90), n_tau=0.05)
    assert sun_st.ice_line_lat == 0.0                       # the Sun-like planet has frozen over
    assert mdwarf_st.ice_line_lat > 30.0                    # the M-dwarf planet stays temperate
    assert mdwarf_st.global_mean_T > sun_st.global_mean_T + 20.0


# --------------------------------------------------------------------------- #
# Knob 2 — planet size → transport
# --------------------------------------------------------------------------- #
def test_earth_size_recovers_the_transport_exactly():
    assert ex.transport_for_size(1.0) == D_TRANSPORT
    assert ex.size_transport_factor(1.0) == 1.0


def test_transport_scales_as_inverse_size_squared():
    for size in (0.5, 0.8, 1.5, 2.0):
        assert ex.transport_for_size(size) == pytest.approx(D_TRANSPORT / size**2, rel=1e-15)


def test_zero_d_mean_is_exactly_size_invariant_but_the_gradient_sharpens():
    # The analytic anchor (constant-albedo, no feedback): the two-mode mean component T0 is the
    # D-free 0-D equilibrium, so it does not move with size; the gradient amplitude T2 ∝ 1/(6D+B)
    # scales EXACTLY. Anchoring on two_mode (not the relaxed ice-cap mean) keeps this leg exact.
    T0 = equilibrium_temperature_0d()
    invariant = []                                          # T2·(6D+B) must be the same constant for all sizes
    gradients = []
    for size in (0.5, 1.0, 2.0):
        D = ex.transport_for_size(size)
        T2 = float(two_mode_solution(1.0, D=D)) - T0        # P₂(1)=1 ⇒ two_mode(1) = T0 + T2
        invariant.append(T2 * (6.0 * D + B_OLR))
        gradients.append(abs(T2))
    assert np.allclose(invariant, invariant[0], rtol=1e-12)               # the exact T2 ∝ 1/(6D+B) scaling
    assert gradients[0] < gradients[1] < gradients[2]                      # bigger planet → steeper gradient


def test_bigger_planet_pushes_the_ice_line_equatorward():
    # Loose qualitative benchmark on the relaxed present-day climate: a bigger planet's weaker transport
    # leaves a colder pole, so the ice cap reaches further toward the equator.
    from projects.planet.albedo import present_day_climate
    icelines = [present_day_climate(EBMParams(D=ex.transport_for_size(s), n_cells=90), n_tau=0.05).ice_line_lat
                for s in (0.5, 1.0, 2.0)]
    assert icelines[0] > icelines[1] > icelines[2]


# --------------------------------------------------------------------------- #
# Composition — both knobs into one EBMParams
# --------------------------------------------------------------------------- #
def test_solar_earth_params_equal_the_base_exactly():
    base = EBMParams()
    assert ex.exoplanet_params(ex.T_SUN, 1.0, base) == base          # clean perturbation: no drift at default


def test_exoplanet_params_apply_both_knobs():
    p = ex.exoplanet_params(T_star=3000.0, size=2.0)
    assert p.ai < ALBEDO_ICE                                          # cooler star lowered the ice albedo
    assert p.D == pytest.approx(D_TRANSPORT / 4.0, rel=1e-15)         # 2 Earth radii → D/4
    assert p.S0 == S0_EARTH and p.s2 == S2_INSOLATION and p.A == EBMParams().A   # nothing else moved
