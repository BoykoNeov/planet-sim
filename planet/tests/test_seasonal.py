"""Triad for the seasonal zonal EBM + continentality (:mod:`planet.seasonal`, rung 5B.1).

*Tight* = (a) the **0-D slab** analytic amplitude ``F₁/√(B²+ω²C²)`` and phase lag ``arctan(ωC/B)``
reproduced by both solvers with transport off; (b) the **reduction** — the spectral ``n=0`` harmonic
equals the annual-mean :class:`~planet.sphere_ebm.SphereEBM` steady solve to machine precision, and
``⟨T_L⟩ = ⟨T_O⟩`` (the annual means of land and ocean are identical — continentality lives *only* in the
seasonal amplitude); (c) **marcher → spectral** at first order in ``dt`` (the anti-damping cross-check —
the backward-Euler transport substep converges to the exact frequency-domain solution, it does not quietly
flatten the swing); (d) **hemispheric antisymmetry** ``T(x,t) = T(−x, t+½yr)``. *Conservation* = the
global-and-annual-mean net TOA ≈ 0. *Loose (calibrated)* = continentality — land seasonal range ≫ ocean,
ocean lags more — banked in direction, the magnitude riding the calibrated heat capacities.
"""
import numpy as np
import pytest

from planet import ebm, seasonal as sea
from planet.sphere_ebm import SphereEBM


# --------------------------------------------------------------------------- #
# TIGHT (a) — the 0-D slab analytic response (the mechanism of continentality).
# --------------------------------------------------------------------------- #
def test_slab_analytic_amplitude_and_lag_both_solvers():
    # Transport OFF (D=0), a pure sinusoidal forcing F0 + F1 cos ωt at every latitude → the exact slab
    # response. Both the spectral solve and the time-marcher must reproduce the closed-form amplitude and
    # phase lag — this pins the time integration and the C-dependence that CREATES continentality.
    m = sea.SeasonalEBM(D=0.0, n_cells=20, n_steps=720, land_fraction=1.0)
    F0, F1 = 300.0, 120.0
    forcing = F0 + F1 * np.cos(2.0 * np.pi * m.day_seconds / sea.SECONDS_PER_YEAR)
    absorbed = np.broadcast_to(forcing, (m.n_cells, m.n_steps)).copy()
    amp_an, lag_an = sea.slab_amplitude_lag(m.C_land, F1)

    sp = m.spectral(absorbed=absorbed)
    mr = m.march(absorbed=absorbed, tol=1e-8, max_years=40)
    lag_sp = sea._time_of_max(sp.T_land[0], m.days) - sea._time_of_max(absorbed[0], m.days)
    lag_mr = sea._time_of_max(mr.T_land[0], m.days) - sea._time_of_max(absorbed[0], m.days)

    assert sp.amplitude("land")[0] == pytest.approx(amp_an, abs=1e-3)
    assert mr.amplitude("land")[0] == pytest.approx(amp_an, abs=5e-3)
    assert lag_sp == pytest.approx(lag_an, abs=1e-2)
    assert lag_mr == pytest.approx(lag_an, abs=0.2)           # marcher lag limited by the 0.5-day dt


def test_slab_larger_C_smaller_amplitude_larger_lag():
    # The continentality mechanism as a monotone law: a bigger heat capacity (ocean) shrinks the amplitude
    # and grows the lag toward the quarter-year quadrature limit.
    ampL, lagL = sea.slab_amplitude_lag(sea.land_heat_capacity(), 120.0)
    ampO, lagO = sea.slab_amplitude_lag(sea.ocean_heat_capacity(), 120.0)
    assert ampL > ampO
    assert lagO > lagL
    assert lagO < 0.25 * sea.SECONDS_PER_YEAR / 86400.0       # below the π/2 (91-day) quadrature ceiling


# --------------------------------------------------------------------------- #
# TIGHT (b) — the reduction to the annual-mean parent (the n=0 harmonic).
# --------------------------------------------------------------------------- #
def test_spectral_dc_harmonic_is_the_annual_mean_ebm():
    # The spectral n=0 harmonic IS the annual-mean EBM forced by ⟨S⟩(1−α). Build the reference from the
    # YEAR-AVERAGE of the actual seasonal insolation (not the P2-truncated ebm.insolation) so both sides
    # see identical forcing — then the reduction is machine-tight, not ~1e-2.
    m = sea.SeasonalEBM(n_cells=180, n_steps=360)
    sp = m.spectral()
    Sbar = m.absorbed_series().mean(axis=1)                   # ⟨S(x,t)⟩(1−α)
    ref = SphereEBM(A=m.A, B=m.B, D=m.D, n_cells=m.n_cells).steady_linear(lambda x, T: Sbar)
    assert np.max(np.abs(sp.annual_mean("mean") - ref.T)) < 1e-9


def test_annual_means_of_land_and_ocean_are_identical():
    # THE core insight: at the annual mean C cancels, so land and ocean at the same latitude reach the
    # SAME temperature — continentality is entirely a seasonal-amplitude phenomenon, zero in the mean.
    m = sea.SeasonalEBM(n_cells=120, n_steps=360)
    sp = m.spectral()
    assert np.max(np.abs(sp.annual_mean("land") - sp.annual_mean("ocean"))) < 1e-9


# --------------------------------------------------------------------------- #
# TIGHT (c) — marcher → spectral at first order (the anti-damping cross-check).
# --------------------------------------------------------------------------- #
def test_marcher_converges_to_spectral_first_order_in_dt():
    # The tight anchors above are all blind to time-accuracy in the transport substep (slab turns it off;
    # reduction/conservation time-average). This is the check that proves the backward-Euler transport is
    # NOT damping the seasonal swing: as dt halves the marcher→spectral gap halves (clean 1st order → a
    # convergent discretization of the exact solution, not a systematic amplitude deficit).
    errs = []
    for ns in (180, 360, 720):
        m = sea.SeasonalEBM(n_cells=120, n_steps=ns)
        sp, mr = m.spectral(), m.march(tol=1e-7, max_years=120)
        assert mr.converged
        errs.append(max(np.max(np.abs(mr.T_land - sp.T_land)),
                        np.max(np.abs(mr.T_ocean - sp.T_ocean))))
    ratio1 = errs[0] / errs[1]
    ratio2 = errs[1] / errs[2]
    assert 1.7 < ratio1 < 2.3 and 1.7 < ratio2 < 2.3          # halving dt halves the error (O(dt))
    # and the PAYOFF (midlatitude amplitude) is essentially exact regardless of dt
    m = sea.SeasonalEBM(n_cells=120, n_steps=360)
    sp, mr = m.spectral(), m.march(tol=1e-7, max_years=120)
    i = m.nearest_index(45)
    assert mr.amplitude("mean")[i] == pytest.approx(sp.amplitude("mean")[i], abs=5e-3)


# --------------------------------------------------------------------------- #
# TIGHT (d) — hemispheric antisymmetry (the seasons are half a year out of phase).
# --------------------------------------------------------------------------- #
def test_hemispheric_antisymmetry():
    # Symmetric geometry (uniform land fraction) ⟹ the limit cycle satisfies T(x,t) = T(−x, t+½yr):
    # NH summer IS SH summer half a year later. A structural check on the seasonal forcing + solver.
    m = sea.SeasonalEBM(n_cells=180, n_steps=360)
    sp = m.spectral()
    half = m.n_steps // 2
    flipped_shifted = np.roll(sp.T_mean[::-1, :], half, axis=1)
    assert np.max(np.abs(sp.T_mean - flipped_shifted)) < 1e-9


# --------------------------------------------------------------------------- #
# CONSERVATION — annual+global energy balance.
# --------------------------------------------------------------------------- #
def test_global_annual_energy_balance():
    # Over one converged year the global-and-annual-mean net TOA ⟨⟨S(1−α) − A − B T̄⟩⟩ ≈ 0 (transport
    # conserves column energy each step: uniform f_L ⟹ constant C_a ⟹ the engine's ∫T̄dx invariant).
    m = sea.SeasonalEBM(n_cells=180, n_steps=360)
    sp = m.spectral()
    net = m.absorbed_series() - m.A - m.B * sp.T_mean
    assert abs(float(net.mean())) < 1e-9


# --------------------------------------------------------------------------- #
# LOOSE (calibrated) — continentality: the payoff, banked in direction.
# --------------------------------------------------------------------------- #
def test_continentality_land_swings_more_and_ocean_lags_more():
    # The headline: at a midlatitude the small-C land tile swings far more than the large-C ocean tile,
    # and the ocean lags the sun more. Direction is banked; the magnitude bands are loose (calibrated C).
    m = sea.SeasonalEBM(n_cells=180, n_steps=360)
    sp = m.spectral()
    i = m.nearest_index(45)
    ampL, ampO = sp.amplitude("land")[i], sp.amplitude("ocean")[i]
    lagL, lagO = m.phase_lag_days(sp.T_land)[i], m.phase_lag_days(sp.T_ocean)[i]
    assert ampL > 4.0 * ampO                                  # land swings several× the ocean
    assert lagO > lagL                                        # the ocean lags the sun more
    assert 15.0 < ampL < 45.0                                 # strong-continental amplitude (loose band)
    assert 1.0 < ampO < 6.0                                   # a maritime seasonal range (loose band)
    assert 40.0 < lagO < 95.0                                 # ~1.5–3 month ocean lag (below the ¼-yr ceiling)


def test_land_fraction_end_members_bracket_the_mixed_planet():
    # f_L = 0 → an aquaplanet (ocean's small amplitude); f_L = 1 → a land world (large amplitude); a mixed
    # planet's zonal-mean amplitude sits between — the continentality knob, reducing to each end-member.
    def amp45(fL):
        m = sea.SeasonalEBM(n_cells=120, n_steps=360, land_fraction=fL)
        return m.spectral().amplitude("mean")[m.nearest_index(45)]
    a_ocean, a_mixed, a_land = amp45(0.0), amp45(0.3), amp45(1.0)
    assert a_ocean < a_mixed < a_land


# --------------------------------------------------------------------------- #
# Forcing sanity — the seasonal insolation reduces to the pinned global mean.
# --------------------------------------------------------------------------- #
def test_seasonal_insolation_global_annual_mean_is_S0_over_4():
    # The absolute seasonal insolation (S0/π · daily kernel) must average — over the globe and the year —
    # to the textbook S0/4 (the disk/sphere factor), the anchor that the reused kernel is normalized right.
    m = sea.SeasonalEBM(n_cells=360, n_steps=365)
    S = m.insolation_series()                                 # [n_x, n_steps], W/m²
    global_annual = np.trapezoid(S.mean(axis=1), m.x) / 2.0   # ∫dx/2 area mean, time mean
    assert global_annual == pytest.approx(ebm.S0_EARTH / 4.0, rel=2e-3)


def test_zero_obliquity_kills_the_seasonal_cycle():
    # An untilted planet (ε=0) has δ≡0 all year → no seasons: the amplitude is ~0 everywhere.
    m = sea.SeasonalEBM(n_cells=90, n_steps=180, obliquity_deg=0.0)
    sp = m.spectral()
    assert np.max(sp.amplitude("land")) < 1e-9
