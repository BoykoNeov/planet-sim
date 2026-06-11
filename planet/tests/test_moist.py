"""Triad for the rung-2 column moist-EBM diagnostic (:mod:`planet.moist`).

Re-classed for honesty (plan §10 / the module docstring): *tight* = the exact Clausius–Clapeyron
``q_sat`` + the conservative moisture operator reproducing the P₂ Legendre eigenvalue (it *is* the EBM
transport operator); *real-but-loose (the unlock)* = the energy-constrained ~2–3 %/K rate; *plumbing*
= machine-exact ``∫(P−E)=0`` + the ``q→0`` reduction + the energy factor's unity-at-reference; *bench-
mark (loose)* = the **named extratropical-only trade** — the equator exports (ITCZ backwards) and the
extratropics converge (``P>E``), but the subtropical evaporative belt is **not** reproduced (so no test
asserts subtropical ``E>P``; one test pins that limitation as the honest finding).
"""
from dataclasses import replace

import numpy as np
import pytest

from planet import moist, precip, transport
from planet.albedo import EBMParams, present_day_climate
from planet.ebm import legendre_P2


@pytest.fixture(scope="module")
def climate():
    return present_day_climate(EBMParams())


# --------------------------------------------------------------------------- #
# TIGHT — q_sat is the exact Clausius–Clapeyron function; the operator is the EBM transport operator.
# --------------------------------------------------------------------------- #
def test_qsat_matches_the_analytic_clausius_clapeyron_formula():
    T = np.array([-20.0, 0.0, 15.0, 30.0])
    e_sat = moist.E_SAT_0 * np.exp((moist.L_VAPOR / moist.R_VAPOR)
                                   * (1.0 / moist.T0_KELVIN - 1.0 / (T + moist.T0_KELVIN)))
    want = moist.EPSILON * e_sat / (moist.P_SURFACE - (1.0 - moist.EPSILON) * e_sat)
    assert np.allclose(moist.saturation_specific_humidity(T), want)
    # exact function, not a fit: monotone increasing with temperature
    assert np.all(np.diff(moist.saturation_specific_humidity(T)) > 0.0)


def test_qsat_reproduces_known_textbook_values():
    # Loose benchmark on the absolute function (Hartmann GPC): ~3.8 g/kg at 0 °C, ~14.7 at 20 °C.
    assert moist.saturation_specific_humidity(0.0) * 1e3 == pytest.approx(3.8, abs=0.3)
    assert moist.saturation_specific_humidity(20.0) * 1e3 == pytest.approx(14.7, abs=0.6)


def test_surface_pressure_matches_transport_no_silent_drift():
    # moist.P_SURFACE (q_sat's surface pressure) and transport.P_SURFACE (the κ→D bridge's C_atm) are
    # the *same* global-mean surface pressure declared in both modules — pin them as a real drift guard
    # (the inline "== transport.P_SURFACE" comment was the only guard; a constant edit must trip a test).
    assert moist.P_SURFACE == transport.P_SURFACE


def test_qsat_local_rate_is_the_clausius_clapeyron_seven_percent():
    # d ln q_sat/dT ≈ L_v/(R_v T²) — the moisture-capacity ~7 %/K precip.py cites (distinct from the
    # energy-constrained *precip* rate). Check the numerical log-slope against the C–C prediction.
    T = 10.0
    h = 1e-3
    rate = (np.log(moist.saturation_specific_humidity(T + h))
            - np.log(moist.saturation_specific_humidity(T - h))) / (2 * h)
    cc = moist.L_VAPOR / (moist.R_VAPOR * (T + moist.T0_KELVIN) ** 2)
    assert rate == pytest.approx(cc, rel=0.02)
    assert 0.06 < rate < 0.075                                   # ~6–7.5 %/K


def test_moisture_operator_reproduces_the_P2_legendre_eigenvalue(climate):
    # ∂/∂x[(1−x²) ∂P₂/∂x] = −6 P₂ (Legendre): the same anchor planet.transport uses for the channel
    # geometry — proves the conservative flux divergence *is* the EBM spherical transport operator.
    x = climate.x
    P2 = legendre_P2(x)
    op = moist._spherical_flux_divergence(P2, x)
    interior = slice(3, -3)                                      # boundary cells are O(Δx)-biased
    # absolute form (not op/P2): robust through the P₂ zero-crossing at φ≈35°, where a ratio would
    # divide by ~0. Truncation is ~2e-5 at n=180; atol=1e-3 is a comfortable margin below the −6 signal.
    assert np.allclose(op[interior], -6.0 * P2[interior], atol=1e-3)


# --------------------------------------------------------------------------- #
# REAL-BUT-LOOSE (the unlock) — the energy-constrained rate, slower than C–C, set by the closure.
# --------------------------------------------------------------------------- #
def test_energy_constrained_rate_is_2_to_3_percent_and_slower_than_cc(climate):
    rate = moist.energy_constrained_rate()
    assert 0.015 < rate < 0.035                                 # ~2–3 %/K (Held & Soden 2006)
    # strictly slower than the column-water-vapour C–C ~7 %/K (the named gap precip.py flags)
    cc = moist.L_VAPOR / (moist.R_VAPOR * (climate.global_mean_T + moist.T0_KELVIN) ** 2)
    assert rate < 0.5 * cc
    assert rate < precip.CC_RATE_PER_K


def test_energy_rate_is_set_by_the_cited_closure_slope_not_derived():
    # The rate is the *cited closure* R_ATM_SLOPE / (L⟨P⟩₀): doubling the prescribed slope doubles it
    # (a free closure choice, not a first-principles constant). R_ATM_SLOPE is the atmospheric-column
    # cooling sensitivity (Held & Soden); it is *coincidentally* 2 W m⁻² K⁻¹ like B_OLR but a different
    # quantity (the module derives the slope from neither B nor the EBM — guarded by the doubling).
    base = moist.energy_constrained_rate(slope=2.0)
    assert moist.energy_constrained_rate(slope=4.0) == pytest.approx(2.0 * base)
    assert moist.R_ATM_SLOPE == pytest.approx(2.0)


def test_energy_constrained_factor_is_linear_not_exponential():
    ref = precip.PRECIP_REF_TEMP_C
    f1 = moist.energy_constrained_factor(ref + 1.0) - 1.0
    f2 = moist.energy_constrained_factor(ref + 2.0) - 1.0
    assert f2 == pytest.approx(2.0 * f1)                        # linear in ΔT (energy budget is linear)
    # the C–C factor (exp) is strictly convex, so it is NOT linear — the honest functional difference
    c1 = precip.clausius_clapeyron_factor(ref + 1.0) - 1.0
    c2 = precip.clausius_clapeyron_factor(ref + 2.0) - 1.0
    assert c2 > 2.0 * c1


# --------------------------------------------------------------------------- #
# PLUMBING (by-construction) — ∫(P−E)=0 machine-exact, the q→0 reduction, unity-at-reference.
# --------------------------------------------------------------------------- #
def test_moisture_convergence_integrates_to_zero_machine_exact(climate):
    pme = moist.moisture_convergence(climate)
    assert abs(float(np.mean(pme))) < 1e-9                       # ∫(P−E) dx = 0 (conservative flux form)


def test_reduction_to_dry_limit_q_to_zero(climate):
    # q → 0 (RH → 0): a vanishing moisture layer ⟹ P − E → 0 (the rung-0 dry limit, by construction).
    assert np.allclose(moist.moisture_convergence(climate, RH=0.0), 0.0)


def test_energy_factor_unity_at_reference_and_field_reduces_to_rung0():
    ref = precip.PRECIP_REF_TEMP_C
    assert moist.energy_constrained_factor(ref) == pytest.approx(1.0)
    # at a reference-temperature climate both amplitudes are 1, so the opt-in field is the rung-0
    # precip field bit-for-bit (the by-construction reduction; replace() sets the frozen dataclass).
    st = replace(present_day_climate(EBMParams()), global_mean_T=ref)
    assert np.allclose(moist.energy_constrained_precip_field(st), precip.precip_field(st))


def test_energy_factor_floored_nonnegative_in_deep_cooling():
    # a deep-Snowball cooling must not drive the global amplitude negative (physical floor)
    assert moist.energy_constrained_factor(-100.0) == 0.0


# --------------------------------------------------------------------------- #
# TRANSPORT D — the consistent path for a non-default-D world (the §9.1 size knob / two_way_pass).
# ClimateState carries no D, so the diagnostic defaults to rung-0 D; the optional param threads the
# climate's own D so moisture diffuses with the coefficient its temperature did.
# --------------------------------------------------------------------------- #
def test_moisture_convergence_default_D_is_rung0_bit_for_bit(climate):
    # the optional D defaults to rung-0 D_TRANSPORT — the no-arg call is unchanged (the consistent-path
    # param must not move the default output; the whole Phase-1 reduction depends on it).
    assert np.array_equal(moist.moisture_convergence(climate),
                          moist.moisture_convergence(climate, D=moist.D_TRANSPORT))


def test_moisture_convergence_scales_linearly_with_transport_D(climate):
    # P − E = (D/c_p)·∂/∂x[(1−x²)∂q/∂x] is linear in D: doubling D doubles the convergence everywhere.
    # This is the consistent path for a non-default-D climate (the §9.1 size knob D∝1/size²,
    # two_way_pass's D_eff) — pass that D so moisture diffuses with the coefficient T did.
    base = moist.moisture_convergence(climate, D=moist.D_TRANSPORT)
    doubled = moist.moisture_convergence(climate, D=2.0 * moist.D_TRANSPORT)
    assert np.allclose(doubled, 2.0 * base)
    assert np.allclose(moist.moisture_budget(climate, D=2.0 * moist.D_TRANSPORT).p_minus_e, doubled)


def test_moisture_convergence_rejects_callable_D(climate):
    # scalar D only — the array-D(x) EBM the rung-1 feedback can drive is a non-goal for the column
    # diagnostic; a callable must raise clearly, not silently broadcast into the (D/c_p) scalar factor.
    with pytest.raises(TypeError):
        moist.moisture_convergence(climate, D=lambda x: moist.D_TRANSPORT * np.ones_like(x))


# --------------------------------------------------------------------------- #
# BENCHMARK (loose) — the named extratropical-only trade.
# --------------------------------------------------------------------------- #
def test_equator_exports_moisture_the_itcz_is_backwards(climate):
    # Down-gradient diffusion EXPORTS moisture from the moist equator → P < E there (the real ITCZ is
    # up-gradient Hadley convergence, deferred). The named trade, not a win.
    b = moist.moisture_budget(climate)
    assert b.equatorial_export < 0.0
    assert b.p_minus_e[0] < 0.0


def test_extratropics_converge_precip_exceeds_evaporation(climate):
    # Poleward of ~40° the column gains moisture (P > E) — the extratropical budget Phase A gets right.
    b = moist.moisture_budget(climate)
    assert b.extratropical_convergence > 0.0
    phi = b.phi
    assert b.p_minus_e[(phi >= 40.0) & (phi <= 60.0)].mean() > 0.0


def test_subtropical_evaporative_belt_is_not_reproduced(climate):
    # HONEST LIMITATION (pinned): the steep equator–pole contrast hyper-peaks C–C q at the equator, so
    # the moisture-flux maximum sits equatorward of the canonical subtropics and the subtropics come out
    # as P > E — NOT the observed evaporative E > P. This is why no test asserts subtropical E > P.
    # NOTE: this pins a KNOWN-WRONG sign as current behaviour. When rung 2.5/3 reproduces the evaporative
    # belt (correct E > P here), this test SHOULD fail — update it *deliberately* then; a red here after
    # that work is the guard doing its job, not a regression.
    b = moist.moisture_budget(climate)
    phi = b.phi
    subtropics = (phi >= 25.0) & (phi <= 35.0)
    assert b.p_minus_e[subtropics].mean() > 0.0                 # model says P>E (wrong sign vs observed)


def test_moisture_budget_is_a_pure_diagnostic(climate):
    # The diagnostic must not perturb the rung-0 climate (Phase-1 triad stays green).
    T_before = climate.T.copy()
    b = moist.moisture_budget(climate)
    assert np.array_equal(climate.T, T_before)
    # the recorded global mean is the energy-constrained ⟨P⟩, and the rate is the unlock
    assert b.mean_precip == pytest.approx(
        moist.GLOBAL_PRECIP_REF_CMYR * moist.energy_constrained_factor(climate.global_mean_T))
    assert b.energy_rate == pytest.approx(moist.energy_constrained_rate())
