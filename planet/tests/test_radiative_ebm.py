"""Triad for the rung-4-completion gray-OLR EBM (:mod:`planet.radiative_ebm`).

Classed for honesty (plan §10 / the module docstring):

* **tight** — a **linear** ``olr_fn`` makes the steady residual affine, so the Newton solve reproduces
  rung-0's direct solve (:meth:`planet.ebm.EnergyBalanceModel.steady_linear`) **bit-for-bit**; at
  convergence the residual is ~machine, so **net-TOA = 0** (``⟨S(1−α)⟩ = ⟨OLR⟩``); and ``column_olr`` is
  the radiation forward model itself (the °C↔K seam pinned).
* **the discriminator (tight, pure column)** — the local slope ``B_loc(Ts)`` is **smallest at the warm
  equator** (the water-vapour feedback beating Planck), which is the fact that fixes the headline *sign*.
* **the unlock (real but loose)** — **tropical** amplification (``δT(pole)/δT(equator) < 1``), the mirror of
  rung-2.5's polar amplification; and the emergent present climate **flatter and warmer-mean** than rung-0.
  Crucially the global-mean warming is **not** pinned at ``ΔA/B`` (Jensen — the moist-EBM identity does not
  transfer): ``⟨δT⟩ > ΔA/B_tan``.
* **named operating point** — the **climlab-matched loading** (water-vapour fraction giving global-mean
  ``B = 2``), the principled, sub-runaway point the wire runs at.
"""
import numpy as np
import pytest

from planet import radiation as rad
from planet import radiative_ebm as r
from planet.albedo import EBMParams
from planet.ebm import D_TRANSPORT, EnergyBalanceModel
from planet.moist_ebm import constant_albedo_absorbed


@pytest.fixture(scope="module")
def params():
    return EBMParams()


@pytest.fixture(scope="module")
def column():
    # The climlab-matched loading (global-mean B = 2), reused across tests.
    return r.climlab_matched_column()


@pytest.fixture(scope="module")
def ta(column, params):
    # The headline experiment once (gray + present-tangent null, warmed by ΔA = 10), reused.
    return r.tropical_amplification(column, params, forcing=10.0)


# --------------------------------------------------------------------------- #
# TIGHT — reduction to rung-0, conservation, the radiation seam.
# --------------------------------------------------------------------------- #
def test_linear_olr_reduces_to_rung0_direct_solve_bit_for_bit(params):
    # A linear olr_fn ≡ rung-0's OLR, so the Newton solve must reproduce EnergyBalanceModel.steady_linear
    # to machine precision — the whole gray departure is then attributable to the OLR curvature alone.
    absorbed = constant_albedo_absorbed(params)
    st_lin = r.RadiativeEBM(r.linear_olr_fn(), D=D_TRANSPORT, n_cells=params.n_cells).equilibrium(absorbed)
    st_dir = EnergyBalanceModel(D=D_TRANSPORT, n_cells=params.n_cells, face="harmonic").steady_linear(absorbed)
    assert np.allclose(st_lin.T, st_dir.T, atol=1e-9, rtol=0.0)
    assert abs(st_lin.global_mean_T - st_dir.global_mean_T) < 1e-9


def test_net_toa_machine_exact_at_equilibrium(ta):
    # At Newton convergence the residual is ~machine and transport conserves ∫T dx, so the global energy
    # balance ⟨absorbed⟩ = ⟨OLR(T)⟩ holds to that tolerance — for both the present and the warmed climate.
    assert ta.gray_present.converged and ta.gray_warm.converged
    assert abs(ta.gray_present.net_toa) < 1e-8
    assert abs(ta.gray_warm.net_toa) < 1e-8


def test_column_olr_is_the_radiation_forward_model_with_the_degC_seam(column):
    # column_olr is a loop over radiation.outgoing_longwave (single source of truth); this pins the °C→K
    # conversion at the seam (the EBM works in °C, the gray column in kelvin).
    Ts_C = np.array([28.0, 14.85, -25.0])
    looped = r.column_olr(column, Ts_C)
    direct = np.array([column.outgoing_longwave(float(t) + 273.15) for t in Ts_C])
    assert np.allclose(looped, direct, rtol=0.0, atol=0.0)


# --------------------------------------------------------------------------- #
# THE DISCRIMINATOR — B_loc smallest at the warm equator (the sign-fixing fact).
# --------------------------------------------------------------------------- #
def test_local_slope_smallest_at_the_warm_equator(column):
    # The advisor's discriminator: whichever latitude has the smallest B_loc warms most. The water-vapour
    # feedback makes B_loc small where it is warm, so a 30 °C equator has a smaller slope than the 15 °C
    # global mean — which then has a smaller slope than nothing… the point is the warm end is the minimum.
    B_eq = float(r.local_radiative_slope(column, [30.0])[0])
    B_mean = float(r.local_radiative_slope(column, [14.85])[0])
    assert B_eq < B_mean                      # equator damps less than the mean → tropical amplification
    assert B_eq > 0.0                         # but still positive (sub-runaway at this loading)


def test_water_vapour_flips_the_planck_ordering(column):
    # The mechanism, cleanly: the *bare Planck* slope 4σTs³ RISES with temperature, so on its own the warm
    # equator would damp MORE than the cold pole (→ polar amplification). The actual B_loc does the OPPOSITE
    # (smaller at the equator) — the water-vapour feedback reverses the Planck ordering, and that flip is the
    # whole reason the amplification comes out tropical.
    sigma = rad.STEFAN_BOLTZMANN
    eq_K, pole_K = 28.0 + 273.15, -25.0 + 273.15
    planck_eq, planck_pole = 4 * sigma * eq_K ** 3, 4 * sigma * pole_K ** 3
    B_eq = float(r.local_radiative_slope(column, [28.0])[0])
    B_pole = float(r.local_radiative_slope(column, [-25.0])[0])
    assert planck_eq > planck_pole            # Planck alone: warm damps more → would be polar
    assert B_eq < B_pole                       # actual: warm damps less → tropical (the WV flip)


def test_local_slope_minimised_at_the_warmest_latitude(column):
    # Sweeping equator→pole, the minimum local damping is at the WARMEST end — so that is where warming
    # concentrates. (B_loc is non-monotone at the cold end, but the global minimum is the warm equator.)
    Ts_C = np.linspace(32.0, -30.0, 40)        # warm equator first → cold pole
    B = r.local_radiative_slope(column, Ts_C)
    assert int(np.argmin(B)) == 0              # the warmest cell has the smallest slope
    assert B.max() - B.min() > 0.5             # and the slope genuinely varies across the planet (not B≈2)


# --------------------------------------------------------------------------- #
# NAMED OPERATING POINT — the climlab-matched loading.
# --------------------------------------------------------------------------- #
def test_climlab_matched_loading_gives_global_slope_two(column):
    # By construction the present-day (288 K) single-column slope is 2, the value rung-0 prescribes; the
    # loading lands near 0.35 (below rung-4's default 0.5, i.e. less greenhouse water vapour).
    assert abs(float(r.local_radiative_slope(column, [rad.PRESENT_SURFACE_T - 273.15])[0]) - 2.0) < 1e-3
    assert 0.30 < column.wv_fraction < 0.40


# --------------------------------------------------------------------------- #
# THE UNLOCK (real but loose) — tropical amplification + the un-pinned mean.
# --------------------------------------------------------------------------- #
def test_warming_is_tropically_amplified(ta):
    # The headline AT THE EARTH-CALIBRATED LOADING: the gray model concentrates warming in the tropics
    # (amp < 1, both endpoint and band) — the mirror of rung-2.5's polar amplification. Both sign and
    # magnitude ride the loading (test_amplification_sign_rides_the_water_vapour_loading); here it is tropical.
    assert ta.amp_gray < 0.95
    assert ta.amp_gray_band < 0.95
    assert 0.5 < ta.amp_gray < 0.9            # loose band around the ~0.68 endpoint ratio
    assert np.all(ta.delta_T_gray > 0)        # everywhere warms (it is a redistribution of *more* warming)


def test_amplification_sign_rides_the_water_vapour_loading(params):
    # The sign is NOT loading-independent (the single-loading suite's blind spot): it is a water-vapour vs
    # Planck competition. A DRY, Planck-dominated planet damps most at the warm equator → POLAR (amp > 1);
    # Earth's calibrated loading is well into the tropical regime (amp < 1). Both sign and magnitude ride
    # the loading (the wall) — unlike rung-2.5's polar direction, which is robust to its RH wall.
    import planet.radiation as rad
    dry = rad.calibrate_column(wv_fraction=0.10)              # Planck-dominated
    earth = r.climlab_matched_column()                        # ~0.35, water-vapour-dominated
    amp_dry = r.tropical_amplification(dry, params, forcing=10.0).amp_gray
    amp_earth = r.tropical_amplification(earth, params, forcing=10.0).amp_gray
    assert amp_dry > 1.0                                      # dry planet: polar
    assert amp_earth < 1.0                                    # Earth loading: tropical
    assert amp_dry > amp_earth                                # monotone: more water vapour → more tropical


def test_warmed_state_is_stable_but_near_the_runaway_edge(ta, column):
    # Newton finds a ROOT whether or not it is a stable equilibrium, so confirm the warmed climate has
    # positive local damping everywhere (stable) — but note the warmed equator (~39 °C) sits at a small
    # B_loc (~0.5), so part of the amplification magnitude is runaway-proximity (the named hot edge).
    B_warm = r.local_radiative_slope(column, ta.gray_warm.T)
    assert np.all(B_warm > 0.0)                               # stable: positive damping everywhere
    assert B_warm[0] < 0.7                                    # but the warm equator is close to runaway


def test_uniform_B_null_warms_uniformly(ta):
    # The clean baseline: the present-tangent uniform-B model warms exactly uniformly (transport of a
    # uniform field is zero), so amp_null = 1 — isolating the nonlinearity as the source of the signal.
    assert abs(ta.amp_null - 1.0) < 1e-6
    assert abs(ta.amp_null_band - 1.0) < 1e-6


def test_global_mean_warming_is_not_pinned_jensen(ta):
    # Unlike rung-2.5 the mean is NOT pinned at ΔA/B: OLR is concave (Jensen) and the WV feedback amplifies
    # the mean response, so the measured ⟨δT⟩ exceeds the naive present-tangent ΔA/B_tan.
    assert ta.mean_delta_T_gray > ta.dA_over_B_tan
    assert ta.mean_delta_T_gray > ta.mean_delta_T_null   # the gray mean warms more than the uniform-B null


def test_present_climate_is_a_jensen_warm_shift_of_rung0(ta, params):
    # The emergent present-day structure at the SAME D, vs rung-0's dt-free reference (steady_linear, NOT the
    # relaxation default, which carries an O(Δt²) contrast error). The contrast is essentially unchanged (the
    # loading-matched mean slope ≈ 2), but the mean is lifted by Jensen (concave OLR) — a near-uniform warm
    # shift. Direction banked, magnitude loose.
    absorbed = constant_albedo_absorbed(params)
    rung0 = EnergyBalanceModel(D=params.D, n_cells=params.n_cells, face="harmonic").steady_linear(absorbed)
    gray_contrast = ta.gray_present.T[0] - ta.gray_present.T[-1]
    rung0_contrast = rung0.T[0] - rung0.T[-1]
    assert ta.gray_present.global_mean_T > rung0.global_mean_T + 1.0   # warmer mean (Jensen), ~2 °C
    assert abs(gray_contrast - rung0_contrast) < 3.0                   # contrast essentially unchanged
