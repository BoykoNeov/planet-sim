"""Triad for the full-sphere MSE-diffusing EBM (:mod:`planet.sphere_moist_ebm`, rung-2.x refinement).

The headline is a **negative**: moisture does NOT tighten the ITCZ-migration sensitivity (it moves the dry
``−6.3 deg/PW`` only ~10 %, to ``−5.7``), because the sensitivity is a **radiation** quantity
``δ/AHT = −1/(2π a² NEI(0))`` that the moisture-amplified transport ``D_eff`` cancels out of.

*Tight* = bit-for-bit reduction to the dry :meth:`SphereEBM.steady_linear` at ``RH = 0`` + EFE ``= 0`` for a
symmetric climate + the ``NEI``-form sensitivity matching the **measured** moist migration. *Real-but-loose
(the unlock — a negative)* = moist ``≈ −5.7``, the dry number barely moved and still ~2× observed; the whole
difference is the ``−B·ΔT_eq`` equatorial-cooling shift of ``NEI(0)``. *Plumbing* = the recalibration matches
the dry contrast; an imposed Q-flux shifts the EFE toward the warmed hemisphere.
"""
import numpy as np
import pytest

from planet import ebm
from planet.sphere_ebm import SphereEBM, itcz_sensitivity_from_nei
from planet import sphere_moist_ebm as sm


ABSORBED = sm.constant_albedo_absorbed()


# --------------------------------------------------------------------------- #
# TIGHT — the reduction, the symmetry, and the NEI identity.
# --------------------------------------------------------------------------- #
def test_rh0_reduces_to_dry_sphere_bit_for_bit():
    # RH = 0 (β ≡ 0) AND D_s = D_TRANSPORT ⟹ D_eff is the dry constant ⟹ the Picard solve collapses to the
    # dry direct solve: the moist model reproduces SphereEBM.steady_linear to machine precision (plumbing).
    dry = SphereEBM(D=ebm.D_TRANSPORT, n_cells=180).steady_linear(ABSORBED)
    moist0 = sm.SphereMoistEBM(D_s=ebm.D_TRANSPORT, RH=0.0, n_cells=180).steady(ABSORBED)
    assert np.max(np.abs(dry.T - moist0.T)) < 1e-11


def test_symmetric_moist_climate_has_efe_at_equator():
    # A symmetric moist climate ⟹ the moist transport H is odd ⟹ EFE = 0 and cross-equatorial AHT = 0.
    c = sm.SphereMoistEBM(D_s=0.28, RH=0.8, n_cells=360).steady(ABSORBED)
    assert c.phi_efe == pytest.approx(0.0, abs=1e-6)
    assert c.aht_eq == pytest.approx(0.0, abs=1e-9)


def test_moist_steady_rejects_ice_feedback():
    # The direct Picard solve is constant-albedo only (like the dry steady_linear / moist_steady_direct).
    m = sm.SphereMoistEBM(D_s=0.28, RH=0.8, n_cells=180)
    with pytest.raises(ValueError, match="state-independent"):
        m.steady(lambda x, T: ebm.insolation(x) * (1.0 - (0.3 + 0.1 * (T > 0))))


def test_moist_nei_identity_matches_measured_migration():
    # THE TIGHT ANCHOR: the moist sensitivity is STILL δ/AHT = −1/(2π a² NEI(0)) — the NEI form matches the
    # engine's MEASURED moist Q-sweep migration, and the product D_eff(0)·T̄ₓₓ(0) is pinned to −NEI(0).
    m = sm.SphereMoistEBM(D_s=0.28, RH=0.8, n_cells=360)
    slope, slope_nei = m.itcz_sensitivity(ABSORBED)
    assert slope_nei == pytest.approx(slope, rel=0.01)          # NEI form == measured moist migration
    c = m.steady(ABSORBED)
    nei0 = m.net_radiative_input_equator(ABSORBED, c.T)
    Teq = float(np.interp(0.0, c.x, c.T))
    Deff0 = sm.effective_diffusivity(Teq, m.D_s, m.RH)
    assert Deff0 * m.equatorial_curvature(c.T) == pytest.approx(-nei0, rel=0.02)  # the identity (fit-limited)


# --------------------------------------------------------------------------- #
# REAL-BUT-LOOSE (the unlock — a NEGATIVE result banked at the identity's altitude).
# --------------------------------------------------------------------------- #
def test_moisture_does_not_tighten_the_sensitivity():
    # THE HEADLINE NEGATIVE: the moist MSE upgrade moves the dry −6.3 by only ~10% (to ≈ −5.7), still ~2×
    # observed. Moisture is a transport intervention; the sensitivity is radiation — so it barely responds.
    r = sm.moist_vs_dry_sensitivity(n_cells=360)
    assert r.slope_dry == pytest.approx(-6.3, abs=0.3)
    assert r.slope_moist == pytest.approx(-5.7, abs=0.3)
    assert abs(r.slope_moist) > abs(r.slope_dry) * 0.85        # moved < 15% — did NOT halve toward observed
    assert abs(r.slope_moist) / 3.0 > 1.6                      # still well above observed ~3


def test_the_whole_moist_effect_is_the_equatorial_cooling_of_NEI():
    # WHY it barely moves (the identity, machine-tight): the ONLY thing the MSE upgrade changes about the
    # sensitivity is NEI(0), and it changes it by exactly −B·ΔT_eq (the moist equator runs ~1.7 K cooler) —
    # NOT through D_eff. D_eff(0) rises and T̄ₓₓ(0) flattens, but their product (= −NEI) is pinned.
    r = sm.moist_vs_dry_sensitivity(n_cells=360)
    nei_diff = r.nei_moist - r.nei_dry
    predicted = -ebm.B_OLR * (r.Teq_moist - r.Teq_dry)         # −B·ΔT_eq, the radiation shift
    assert nei_diff == pytest.approx(predicted, rel=1e-6)      # the entire difference IS the cooling
    assert r.Teq_moist < r.Teq_dry                             # the moist equator is cooler
    # and both sensitivities are exactly their NEI forms:
    assert r.slope_dry == pytest.approx(itcz_sensitivity_from_nei(r.nei_dry), rel=0.01)
    assert r.slope_moist == pytest.approx(itcz_sensitivity_from_nei(r.nei_moist), rel=0.01)


def test_moist_sensitivity_saturates_across_RH():
    # The ~10% is a CEILING, not a slope: raising RH from 0.4 to 0.95 barely moves the sensitivity (the
    # D_eff/curvature cancellation is structural, not a coincidence at one RH).
    slopes = []
    for RH in (0.4, 0.95):
        D_s = sm.recalibrate_sensible_D_sphere(RH=RH, n_cells=180)
        _, slope_nei = sm.SphereMoistEBM(D_s=D_s, RH=RH, n_cells=180).itcz_sensitivity(ABSORBED)
        slopes.append(slope_nei)
    assert abs(slopes[0] - slopes[1]) < 0.3                    # RH 0.4 → 0.95 moves it < 0.3 deg/PW


# --------------------------------------------------------------------------- #
# PLUMBING — the recalibration and the imposed-Q migration direction.
# --------------------------------------------------------------------------- #
def test_recalibration_matches_the_dry_contrast():
    # D_s is re-derived (the double-count discipline) so the moist present climate reproduces the dry sphere
    # equator-pole contrast — and lands on ~0.28 (< the dry 0.555), the hemisphere value.
    D_s = sm.recalibrate_sensible_D_sphere(RH=0.8, n_cells=180)
    assert D_s == pytest.approx(0.28, abs=0.03)
    absorbed = ABSORBED
    dry = SphereEBM(D=ebm.D_TRANSPORT, n_cells=180).steady_linear(absorbed)
    moist = sm.SphereMoistEBM(D_s=D_s, RH=0.8, n_cells=180).steady(absorbed)
    assert sm.equator_pole_contrast(moist.x, moist.T) == pytest.approx(
        sm.equator_pole_contrast(dry.x, dry.T), abs=0.05)


def test_imposed_qflux_shifts_moist_itcz_toward_warm_hemisphere():
    # An imposed cross-equatorial Q-flux (heats the NH) moves the moist EFE/ITCZ into the NH (φ>0) — the
    # direction is by-construction, exactly as in the dry sphere.
    m = sm.SphereMoistEBM(D_s=0.28, RH=0.8, n_cells=360)
    c0 = m.steady(ABSORBED)
    cq = m.steady(ABSORBED, Q=2.0 * m.x)
    assert c0.phi_efe == pytest.approx(0.0, abs=1e-6)
    assert cq.phi_efe > 0.3


@pytest.mark.slow
def test_demo_reproduces_the_radiation_limit_headline():
    # Guards the committed figure (planet-itcz-radiation-limit.png): a fresh clone reproduces the headline —
    # transport (dry→moist) barely moves the sensitivity, the isothermal ceiling is the floor, and observed
    # needs an NEI(0) ABOVE that ceiling (unreachable by any transport).
    from planet import demo_itcz_radiation_limit as demo
    r = demo.compute()
    assert r.slope_dry == pytest.approx(-6.3, abs=0.3)
    assert r.slope_moist == pytest.approx(-5.7, abs=0.3)
    assert r.slope_ceiling == pytest.approx(-3.9, abs=0.2)          # the D→∞ floor
    assert abs(r.slope_ceiling) > 3.0                              # floor sits above observed
    assert r.nei_observed > r.nei_ceiling                         # observed needs NEI beyond the ceiling
    # the D-sweep monotonically slides NEI up toward (never past) the ceiling:
    assert np.all(np.diff(r.d_sweep_nei) > 0)
    assert r.d_sweep_nei[-1] < r.nei_ceiling + 1e-6
