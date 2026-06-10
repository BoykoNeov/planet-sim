"""Planet Phase-2 validation: the diagnostic precipitation parameterization (plan §3, the precip leg).

The precip half of the biome-map triad (planet-earth-system.md). What is asserted **tight** is the
*qualitative band structure* (the ITCZ-wet / subtropics-dry / midlat-wet / poles-dry partition of
latitude) and the *monotone* global-water response to warming; what is **loose** is the *calibrated
band amplitudes/centres* and the 7 %/K rate (cited [[precip-parameterization-source]]). A prescribed
field has no water-mass conservation law, so the "conservation" leg is a **consistency** check
(``⟨P⟩`` moves monotonically with T̄ as designed) — named as such, not dressed up as a law.
"""
import numpy as np
import pytest

from planet import precip
from planet.albedo import present_day_climate


def _area_grid(n: int = 400):
    """Latitudes (degrees) sampled uniformly in x = sin φ — the equal-area coordinate (mean = area mean)."""
    x = (np.arange(n) + 0.5) / n
    return np.degrees(np.arcsin(x))


# --------------------------------------------------------------------------- #
# Pinned parameters ([[precip-parameterization-source]])
# --------------------------------------------------------------------------- #
def test_pinned_constants_are_cited_values():
    assert precip.CC_RATE_PER_K == pytest.approx(0.07)          # Clausius–Clapeyron ~7 %/K
    assert precip.PRECIP_REF_TEMP_C == pytest.approx(15.0)      # ≈ present-day global mean
    assert precip.ITCZ_CENTER_DEG == pytest.approx(0.0)         # ITCZ on the equator
    assert precip.MIDLAT_CENTER_DEG == pytest.approx(50.0)      # storm tracks near 50°


# --------------------------------------------------------------------------- #
# Structural (tight): the latitudinal band structure
# --------------------------------------------------------------------------- #
def test_pattern_band_structure_itcz_subtropics_midlat_poles():
    p_eq = float(precip.precip_pattern(0.0))                   # ITCZ
    p_sub = float(precip.precip_pattern(25.0))                 # subtropics
    p_mid = float(precip.precip_pattern(50.0))                 # midlatitude storm track
    p_pole = float(precip.precip_pattern(90.0))               # pole
    assert p_eq > p_mid > p_sub                                # equator wettest; subtropics a local min
    assert p_sub < p_mid and p_sub < p_eq                      # the dry subtropics sit BETWEEN two wet belts
    assert p_pole < p_mid                                      # poles dry
    assert p_eq > 150.0 and p_sub < 60.0                       # loose magnitude bands (cm/yr)


def test_subtropical_minimum_is_a_local_trough_between_itcz_and_storm_track():
    # the desert belt is a LOCAL trough between the equatorial ITCZ and the midlatitude storm track
    # (the poles are drier in absolute terms — the global floor — but that is not the subtropical dip).
    phi = np.linspace(0.0, 90.0, 181)
    p = precip.precip_pattern(phi)
    belt = (phi > 8.0) & (phi < precip.MIDLAT_CENTER_DEG)       # equator → midlat peak
    phi_min = phi[belt][int(np.argmin(p[belt]))]
    assert 20.0 < phi_min < 38.0                                # the trough sits in the subtropics
    # and it really is a local minimum: drier than both its equatorward and poleward neighbours
    assert float(precip.precip_pattern(phi_min)) < float(precip.precip_pattern(phi_min - 10.0))
    assert float(precip.precip_pattern(phi_min)) < float(precip.precip_pattern(phi_min + 15.0))


def test_precipitation_is_nonnegative_everywhere():
    phi = _area_grid()
    for Tbar in (-30.0, 0.0, 15.0, 40.0):
        assert np.all(precip.precipitation(phi, Tbar) >= 0.0)


# --------------------------------------------------------------------------- #
# Clausius–Clapeyron amplitude
# --------------------------------------------------------------------------- #
def test_cc_factor_unity_at_reference_and_monotone():
    assert precip.clausius_clapeyron_factor(precip.PRECIP_REF_TEMP_C) == pytest.approx(1.0)
    warmer = precip.clausius_clapeyron_factor(20.0)
    cooler = precip.clausius_clapeyron_factor(10.0)
    assert cooler < 1.0 < warmer                               # warmer atmosphere holds more moisture
    # the exact exp form at the cited 7 %/K rate
    assert warmer == pytest.approx(np.exp(0.07 * (20.0 - 15.0)))


# --------------------------------------------------------------------------- #
# Conservation = a *consistency* check (honestly weaker): the global-water budget
# --------------------------------------------------------------------------- #
def test_global_mean_precip_increases_monotonically_with_warming():
    # NOT a conservation law (prescribed field, no water budget) — a consistency check: the C–C-scaled
    # area-mean precip moves monotonically with T̄ exactly as designed. x = sin φ is the equal-area grid.
    phi = _area_grid()
    Tbars = np.linspace(-5.0, 35.0, 9)
    means = np.array([float(np.mean(precip.precipitation(phi, T))) for T in Tbars])
    assert np.all(np.diff(means) > 0.0)                        # strictly increasing with warming
    # and it scales as the pattern mean times the C–C factor (the pattern shape is T̄-independent)
    base = float(np.mean(precip.precip_pattern(phi)))
    assert means[-1] == pytest.approx(base * precip.clausius_clapeyron_factor(Tbars[-1]))


# --------------------------------------------------------------------------- #
# The state convenience wrapper
# --------------------------------------------------------------------------- #
def test_precip_field_matches_precipitation_on_the_state_grid():
    st = present_day_climate(n_tau=0.02)
    got = precip.precip_field(st)
    want = precip.precipitation(st.latitude_deg(), st.global_mean_T)
    assert np.allclose(got, want)
    assert got.shape == st.T.shape                             # paired with T on the EBM grid
