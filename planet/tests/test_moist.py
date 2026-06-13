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
    # HONEST LIMITATION (pinned) of the EDDY-ONLY DEFAULT: the steep equator–pole contrast hyper-peaks C–C
    # q at the equator, so the moisture-flux maximum sits equatorward of the canonical subtropics and the
    # subtropics come out as P > E — NOT the observed evaporative E > P. This is why no test asserts
    # subtropical E > P for the default path.
    # NOTE: the mislocation PERSISTS past the Hadley fix. The opt-in Hadley path
    # (moisture_budget(..., hadley=True)) flips the ITCZ sign and DOES make a dry belt — but equatorward of
    # 25–35° (~10–15°), because the hyper-peaked fixed-RH q pulls the moisture flux equatorward. So 25–35°
    # stays P>E on BOTH paths (test_hadley_creates_an_equatorward_dry_belt_but_does_not_relocate_the_desert).
    # Relocating the desert to the canonical subtropics needs a realistic q (moist dynamics / vertical) =
    # rung 3+. This test pins the limitation for the default; it stays green by design.
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


# --------------------------------------------------------------------------- #
# THE HADLEY MOISTURE-CONVERGENCE FIX (opt-in) — the deferred deep-tropical mean-circulation term.
# Honesty classification (the advisor's framing): convergence-at-ITCZ / divergence-in-subtropics is
# BY-CONSTRUCTION (plumbing, not a win — HADLEY_STRENGTH is the prescribed wall); the genuinely emergent,
# non-vacuous claim is the AMPLITUDE — q(T) from the EBM ⟹ the ITCZ convergence intensifies at the ~C–C
# rate ("rich-get-richer" P−E, faster than the energy-constrained global mean).
# --------------------------------------------------------------------------- #
def test_hadley_streamfunction_is_a_tropics_confined_cell():
    # PLUMBING (by-construction): ψ(x) = 0 at the equatorial ascent AND the subtropical-edge descent,
    # positive between, exactly 0 in the extratropics (the prescribed single cell).
    x = np.linspace(0.0, 1.0, 200)
    psi = moist.hadley_streamfunction(x)
    x_edge = np.sin(np.radians(moist.HADLEY_EDGE_DEG))
    assert psi[0] == pytest.approx(0.0)                         # ascent at the equator
    assert np.all(psi[x >= x_edge] == 0.0)                      # nothing poleward of the cell edge
    assert np.all(psi[(x > 0.02) & (x < x_edge - 0.02)] > 0.0)  # a positive overturning in between
    assert psi.max() == pytest.approx(1.0, abs=1e-3)            # normalized peak


def test_hadley_convergence_conserves_water_machine_exact(climate):
    # PLUMBING: the conservative face form ⟹ ∫(P−E)_Hadley dx = 0 (the ITCZ convergence is exactly paid
    # for by subtropical divergence — a budget, not a painted band).
    had = moist.hadley_moisture_convergence(climate)
    assert abs(float(np.mean(had))) < 1e-9


def test_hadley_convergence_is_confined_to_the_tropics(climate):
    # The mean cell vanishes poleward of its edge ⟹ the extratropical eddy budget is UNTOUCHED.
    had = moist.hadley_moisture_convergence(climate)
    phi = climate.latitude_deg()
    assert np.all(had[phi >= moist.HADLEY_EDGE_DEG + 1.0] == 0.0)


def test_hadley_reduces_to_eddy_only_at_zero_strength(climate):
    # REDUCTION (by-construction): strength = 0 ⟹ the Hadley term is identically zero, and the opt-in
    # budget is the eddy-only default bit-for-bit (the independent-diff discipline).
    assert np.array_equal(moist.hadley_moisture_convergence(climate, strength=0.0),
                          np.zeros_like(climate.x))
    off = moist.moisture_budget(climate)
    on0 = moist.moisture_budget(climate, hadley=True, strength=0.0)
    assert np.array_equal(off.p_minus_e, on0.p_minus_e)


def test_hadley_flips_the_deep_tropical_sign_to_convergence(climate):
    # THE FIX (the headline): the eddy-only default EXPORTS at the equator (P<E, backwards); adding the
    # mean Hadley convergence flips it to CONVERGENCE (P>E) at the ITCZ. The SIGN FLIP is the deliverable;
    # the magnitude is calibrated (observed order only — HADLEY_STRENGTH is the prescribed wall).
    eddy = moist.moisture_budget(climate)
    full = moist.moisture_budget(climate, hadley=True)
    assert eddy.equatorial_export < 0.0                        # was backwards (the named trade)
    assert full.equatorial_export > 0.0                        # now converges (the fix)
    assert full.hadley is True
    # observed ITCZ order (~1–2 m/yr = 100–200 cm/yr); a loose band, NOT a tuned match
    assert 50.0 < full.equatorial_export < 350.0


def test_hadley_leaves_the_extratropical_budget_unchanged(climate):
    # The fix is local to the tropics: the extratropical convergence (the leg the eddy diffusion gets
    # right) is identical with and without the Hadley term.
    eddy = moist.moisture_budget(climate)
    full = moist.moisture_budget(climate, hadley=True)
    assert full.extratropical_convergence == pytest.approx(eddy.extratropical_convergence)


def test_hadley_creates_an_equatorward_dry_belt_but_does_not_relocate_the_desert(climate):
    # THE TRADE (named, honest): the descending branch produces a dry belt (E>P) — but EQUATORWARD of the
    # canonical 25–35° subtropics (~10–15°), because the hyper-peaked fixed-RH C–C q pulls the moisture
    # flux ψ·q equatorward (the SAME mislocation the eddy budget has). So the Hadley path flips the ITCZ
    # SIGN robustly but does NOT relocate the desert: the canonical subtropics stay P>E on BOTH paths.
    full = moist.moisture_budget(climate, hadley=True)
    phi = full.phi
    dry = (phi >= 8.0) & (phi <= 20.0)
    assert full.p_minus_e[dry].mean() < 0.0                    # an off-equator dry belt (E>P) emerges
    assert full.subtropical_balance > 0.0                      # but 25–35° stays P>E (desert NOT relocated)
    # the dry belt sits equatorward of the canonical subtropics (the persistent hyper-peaked-q limitation)
    desert_lat = phi[np.argmin(np.where(phi < 35.0, full.p_minus_e, np.inf))]
    assert desert_lat < 25.0


def test_hadley_itcz_convergence_intensifies_at_the_cc_rate(climate):
    # REAL-BUT-LOOSE — the non-vacuous EMERGENT nugget (advisor): the cell strength is FIXED, so any
    # warming response is carried entirely by q(T) from the EBM. The ITCZ convergence intensifies at the
    # ~Clausius–Clapeyron moisture rate (~7 %/K) — the "rich-get-richer" P−E scaling — and is FASTER than
    # the energy-constrained global-mean rate (~2.5 %/K). This is emergent, not prescribed.
    from dataclasses import replace
    dT = 4.0
    warm = replace(climate, T=climate.T + dT, global_mean_T=climate.global_mean_T + dT)
    eq_now = moist.hadley_moisture_convergence(climate)[0]
    eq_warm = moist.hadley_moisture_convergence(warm)[0]        # SAME strength — only q(T) moved
    rate = (eq_warm / eq_now - 1.0) / dT
    assert 0.04 < rate < 0.09                                   # C–C moisture order (~7 %/K)
    assert rate > moist.energy_constrained_rate()              # faster than the energy-constrained mean


def test_hadley_convergence_is_a_pure_diagnostic(climate):
    # Like the eddy budget: the mean-circulation term must not perturb the rung-0 climate.
    T_before = climate.T.copy()
    moist.hadley_moisture_convergence(climate)
    moist.moisture_budget(climate, hadley=True)
    assert np.array_equal(climate.T, T_before)


@pytest.mark.slow
def test_demo_reproduces_the_hadley_fix_headline():
    # Guards the committed figure (planet-hadley-moisture.png): a fresh clone reproduces the headline, not
    # just reads it — the deep-tropical sign FLIP (export → convergence), the conserving budget, and the
    # EMERGENT ~C–C-order intensification (faster than the energy-constrained global mean).
    from planet import demo_hadley_moisture as demo
    r = demo.compute()
    assert r.eq_eddy < 0.0                                      # eddy-only: backwards ITCZ (export)
    assert r.eq_full > 0.0                                      # + Hadley: the fix (convergence)
    assert 50.0 < r.eq_full < 350.0                             # observed ITCZ order (~1–2 m/yr), loose
    assert r.dry_belt_min < 0.0 and r.dry_belt_lat < 25.0       # a dry belt emerges, equatorward of 25–35°
    assert r.subtropics_full > 0.0                              # canonical 25–35° NOT relocated (the trade)
    assert abs(r.net_full) < 1e-9                               # water conserved (a budget, not a band)
    assert 0.04 < r.itcz_rate < 0.09                            # emergent ~C–C rate ("rich-get-richer")
    assert r.itcz_rate > r.energy_rate                          # faster than the global-mean rate
