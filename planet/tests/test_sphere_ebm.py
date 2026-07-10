"""Triad for the full-sphere EBM + the energetic ITCZ (:mod:`planet.sphere_ebm`, rung 2.x).

*Tight* = the North two-mode on the full sphere (splitting-free direct solve) + the **reduction** to the
hemisphere ``ebm.py`` climate under symmetric forcing + the closed form ``δ/AHT = 1/(2π a² D T̄ₓₓ(0))``
reproduced by the engine + EFE ``= 0`` for symmetric forcing. *Real-but-loose (the unlock, lower altitude)*
= the ITCZ sensitivity is of the observed order (~3 deg/PW, a factor ~1.5–2 high) and is a pure function of
the **equatorial radiative surplus** ``NEI(0)`` (the identity ``D·T̄ₓₓ(0) = −NEI(0)`` ⟹
``δ/AHT = −1/(2π a² NEI(0))``; equivalently ``∝ (6 + B/D)``, **not** the naive ``∝ 1/D``) — banked as a
closed-form consequence of radiation, **not** an emergent prediction, and radiatively floored at ``−3.9``
*above* observed so no transport ``D`` reaches it (the forcing-independence is a linear-operator identity,
pinned by a test). *Plumbing* = the precip wiring reduces to rung-0 at
``φ_EFE = 0`` and the imposed Q-flux shifts the ITCZ toward the warm hemisphere.
"""
import numpy as np
import pytest

from planet import ebm, precip, sphere_ebm as se
from planet.albedo import absorbed_shortwave, present_day_climate, EBMParams


def _absorbed_const(x, T):
    """Constant-albedo (no ice) absorbed shortwave — the linear, splitting-free reference forcing."""
    return ebm.insolation(x) * (1.0 - ebm.ALBEDO_A0)


# --------------------------------------------------------------------------- #
# TIGHT — the new engine's analytic anchors.
# --------------------------------------------------------------------------- #
def test_reduces_to_hemisphere_ebm_under_symmetric_forcing():
    # The protected Phase-1 climate, re-validated as a CROSS-MODEL reduction (ebm.py untouched):
    # the full-sphere ice climate's northern hemisphere == the hemisphere ebm.py climate (same n_tau).
    m = se.SphereEBM(n_cells=360)
    T_init = 30.0 - 60.0 * np.abs(m.x)                       # warm equator → iced poles (matches ebm.py IC)
    cice = m.equilibrate(lambda x, T: absorbed_shortwave(x, T), T_init, n_tau=0.5)
    hemi = present_day_climate(EBMParams())                  # hemisphere n=180, same default n_tau
    assert np.max(np.abs(cice.T[m.x > 0] - hemi.T)) < 1e-6
    assert cice.global_mean_T == pytest.approx(hemi.global_mean_T, abs=1e-6)
    assert np.max(np.abs(cice.T - cice.T[::-1])) < 1e-9     # the full-sphere solution is even


def test_north_two_mode_on_full_sphere_constant_albedo():
    # The tight analytic anchor of the new engine, via the splitting-free direct solve. The error is
    # limited by the SAME harmonic-face polar floor ebm.py documents (~0.16 °C), not clean 2nd order — it
    # improves only weakly with resolution (face="exact" would remove it; out of scope for this diagnostic).
    errs = []
    for n in (180, 360, 720):
        m = se.SphereEBM(n_cells=n)
        c = m.steady_linear(_absorbed_const)
        errs.append(np.max(np.abs(c.T - ebm.two_mode_solution(m.x, albedo=ebm.ALBEDO_A0))))
    assert errs[-1] < 0.2                                    # the harmonic-face polar floor
    assert errs[-1] <= errs[0]                               # does not worsen with resolution


def test_symmetric_climate_has_efe_at_equator():
    # A hemispherically symmetric climate ⟹ H is odd ⟹ EFE = 0 and cross-equatorial transport = 0.
    m = se.SphereEBM(n_cells=360)
    c = m.steady_linear(_absorbed_const)
    assert c.phi_efe == pytest.approx(0.0, abs=1e-6)
    assert c.aht_eq == pytest.approx(0.0, abs=1e-9)


def test_closed_form_matches_engine_sensitivity():
    # The engine's measured deg/PW equals the closed form δ/AHT = 1/(2π a² D T̄ₓₓ(0)) — the tight identity.
    m = se.SphereEBM(n_cells=360)
    slope, closed = m.itcz_sensitivity(_absorbed_const, linear=True)
    assert slope == pytest.approx(closed, rel=0.02)


def test_global_mean_matches_0d_anchor():
    # The relaxed/solved mean obeys the 0-D global energy balance (transport is mean-preserving).
    m = se.SphereEBM(n_cells=360)
    c = m.steady_linear(_absorbed_const)
    assert c.global_mean_T == pytest.approx(ebm.equilibrium_temperature_0d(albedo=ebm.ALBEDO_A0), abs=0.05)


def test_global_energy_balance_closes_to_machine_precision():
    # The conservation plumbing leg (named in the plan): the no-flux transport only redistributes, so at
    # steady state ⟨S(1−α)⟩ = A + B⟨T⟩ to machine precision (the discrete net-TOA is ~0).
    m = se.SphereEBM(n_cells=360)
    c = m.steady_linear(_absorbed_const)
    net_toa = m.global_mean(_absorbed_const(m.x, c.T)) - m.A - m.B * c.global_mean_T
    assert abs(net_toa) < 1e-9


def test_steady_linear_rejects_ice_feedback():
    # The direct solve is constant-albedo only; the ice nonlinearity must go through equilibrate.
    m = se.SphereEBM(n_cells=180)
    with pytest.raises(ValueError, match="state-independent"):
        m.steady_linear(lambda x, T: absorbed_shortwave(x, T))


def test_nei_identity_pins_the_sensitivity_denominator():
    # THE TIGHT ANCHOR (upgraded loose→tight): at the symmetric steady state the equatorial energy balance
    # pins D·T̄ₓₓ(0) = −NEI(0), so δ/AHT = −1/(2π a² NEI(0)) is IDENTICALLY the curvature closed form — the
    # sensitivity is a RADIATION quantity (D cancelled), not a transport one. And the NEI form is the TIGHTER
    # reading: it matches the engine's MEASURED migration slope, free of the curvature-fit error.
    m = se.SphereEBM(n_cells=360)
    c = m.steady_linear(_absorbed_const)                        # the symmetric base state
    nei0 = m.net_radiative_input_equator(_absorbed_const, c.T)
    # the identity D·T̄ₓₓ(0) = −NEI(0), machine-tight (the dry base is exactly P₂ ⟹ the polyfit is exact):
    assert m.D * m.equatorial_curvature(c.T) == pytest.approx(-nei0, rel=1e-3)
    # the two closed forms are the same number, and both equal the measured migration:
    slope, closed_curv = m.itcz_sensitivity(_absorbed_const, linear=True)
    assert se.itcz_sensitivity_from_nei(nei0) == pytest.approx(slope, rel=5e-3)
    assert se.itcz_sensitivity_from_nei(nei0) == pytest.approx(closed_curv, rel=0.02)


def test_sensitivity_is_radiatively_floored_above_observed():
    # THE HANDOFF (why re-deriving D cannot tighten it): as D→∞ the equator cools to T̄ (isothermal), NEI(0)
    # rises to its CEILING S(0)(1−α)−A−B·T̄, and the sensitivity bottoms out at ≈ −3.9 deg/PW — ABOVE observed
    # −3. So no transport reaches observed; the lever is a stronger equatorial radiative surplus (rung 4).
    m = se.SphereEBM(n_cells=360)
    absorbed = _absorbed_const(m.x, np.zeros_like(m.x))
    Tbar = (m.global_mean(absorbed) - m.A) / m.B                # the isothermal-limit temperature
    absorbed0 = float(np.interp(0.0, m.x, absorbed))
    nei_ceiling = absorbed0 - m.A - m.B * Tbar                  # NEI(0) at the isothermal (D→∞) limit
    floor = se.itcz_sensitivity_from_nei(nei_ceiling)
    assert floor == pytest.approx(-3.9, abs=0.2)               # the radiative floor
    # a huge D nearly reaches it (confirms the floor is the D→∞ limit, not an asymptote we never approach):
    big = se.SphereEBM(n_cells=360, D=100.0 * ebm.D_TRANSPORT)
    slope_big, _ = big.itcz_sensitivity(_absorbed_const, linear=True)
    assert slope_big == pytest.approx(floor, abs=0.1)
    # observed −3 needs NEI(0) ABOVE the isothermal ceiling ⟹ unreachable by ANY transport:
    nei_for_observed = -(180.0 / np.pi) * se.PW / (se.AREA_FACTOR * -3.0)
    assert nei_for_observed > nei_ceiling                      # ~75 W/m² needed vs ~57 ceiling
    assert abs(floor) > 3.0                                    # the floor magnitude sits above observed


# --------------------------------------------------------------------------- #
# REAL-BUT-LOOSE (the unlock) — banked at the lower altitude: a closed-form consequence of D.
# --------------------------------------------------------------------------- #
def test_itcz_sensitivity_order_and_sign():
    # The no-ice sensitivity ≈ −6.3 deg/PW: same ORDER as the observed ~3 deg/PW but a factor ~1.5–2 high
    # (NOT a precise match — the lower-altitude claim), and negative (ITCZ toward the warm hemisphere).
    m = se.SphereEBM(n_cells=360)
    slope, _ = m.itcz_sensitivity(_absorbed_const, linear=True)
    assert slope < 0.0                                       # toward the warm hemisphere
    assert slope == pytest.approx(-6.3, abs=0.4)             # the splitting-free no-ice value
    assert 2.0 < abs(slope) / 3.0 < 2.5                      # same order as observed, factor ~2 high


def test_sensitivity_depends_on_calibrated_D():
    # EXPOSING the D-dependence: the sensitivity is a property of the transport D (through the mean-state
    # curvature), NOT of ITCZ dynamics. It tracks the closed form at every D, weakens as D grows (more
    # transport ⟹ flatter base state), and follows the analytic two-mode law deg/PW ∝ (6 + B/D) — the
    # curvature is ∝ 1/(6D+B), so it is a *pure function of D*, not "∝ 1/D" (the curvature moves with D too).
    slopes = []
    for fac in (0.8, 1.0, 1.2):
        m = se.SphereEBM(n_cells=360, D=ebm.D_TRANSPORT * fac)
        slope, closed = m.itcz_sensitivity(_absorbed_const, linear=True)
        assert slope == pytest.approx(closed, rel=0.02)             # the closed form holds at every D
        slopes.append(slope)
    assert abs(slopes[0]) > abs(slopes[1]) > abs(slopes[2])         # |deg/PW| weakens as D increases
    law = [-(6.0 + ebm.B_OLR / (ebm.D_TRANSPORT * f)) for f in (0.8, 1.0, 1.2)]
    ratios = [s / l for s, l in zip(slopes, law)]
    assert max(ratios) - min(ratios) < 0.02 * abs(np.mean(ratios))  # deg/PW ∝ (6 + B/D) exactly


def test_sensitivity_is_forcing_independent_a_linear_operator_identity():
    # Q-flux and an antisymmetric ALBEDO give the SAME deg/PW — NOT robustness but a linear-operator
    # tautology (both make an equatorial gradient anomaly through the same operator). Pinned so the honest
    # framing is not silently lost: this is why the number is a property of D, not of ITCZ physics.
    m = se.SphereEBM(n_cells=360)
    slope_q, _ = m.itcz_sensitivity(_absorbed_const, linear=True)
    # albedo-asymmetry forcing α = a0 + asym·x (odd); still state-independent ⟹ steady_linear applies.
    ahts, phis = [], []
    for asym in (0.0, 0.02, 0.04):
        def absorbed_asym(x, T, a=asym):
            return ebm.insolation(x) * (1.0 - (ebm.ALBEDO_A0 + a * x))
        c = m.steady_linear(absorbed_asym)
        phis.append(c.phi_efe); ahts.append(c.aht_eq)
    slope_alb = float(np.polyfit(ahts, phis, 1)[0])
    assert slope_alb == pytest.approx(slope_q, rel=0.03)


@pytest.mark.slow
def test_ice_climate_sensitivity_converged():
    # The present-day (ice) sensitivity from a CONVERGED profile (small n_tau — the splitting-error gotcha):
    # ≈ −4.9 deg/PW (a steeper equatorial curvature than no-ice), still ~2× the observed ~3.
    m = se.SphereEBM(n_cells=360)
    slope, closed = m.itcz_sensitivity(lambda x, T: absorbed_shortwave(x, T), linear=False, n_tau=0.02)
    assert slope == pytest.approx(closed, rel=0.05)          # closed form still holds on the ice base
    assert slope == pytest.approx(-4.9, abs=0.6)
    assert slope < 0.0


# --------------------------------------------------------------------------- #
# PLUMBING — the precip wiring + the imposed-asymmetry migration.
# --------------------------------------------------------------------------- #
def test_imposed_qflux_shifts_itcz_toward_warm_hemisphere():
    # An imposed cross-equatorial Q-flux Q(x)=q·x (heats the NH) moves the EFE/ITCZ into the NH (φ>0).
    m = se.SphereEBM(n_cells=360)
    c0 = m.steady_linear(_absorbed_const)
    cq = m.steady_linear(_absorbed_const, Q=2.0 * m.x)
    assert c0.phi_efe == pytest.approx(0.0, abs=1e-6)
    assert cq.phi_efe > 0.3                                  # shifted into the (warmed) northern hemisphere


def test_precip_wiring_reduces_to_rung0_at_equator():
    # itcz_informed_precip with φ_EFE = 0 recovers the rung-0 precip field BIT-FOR-BIT (the reduction).
    m = se.SphereEBM(n_cells=180)
    c0 = m.steady_linear(_absorbed_const)
    assert c0.phi_efe == pytest.approx(0.0, abs=1e-6)
    lat = c0.latitude_deg()
    rung0 = precip.precipitation(lat, c0.global_mean_T)      # itcz_center_deg defaults to 0
    # force exact φ=0 to isolate the wiring reduction from a ~1e-7 EFE residual
    wired0 = precip.precipitation(lat, c0.global_mean_T, itcz_center_deg=0.0)
    assert np.array_equal(rung0, wired0)


def test_itcz_informed_precip_moves_the_rain_belt():
    # With the ITCZ shifted off the equator, the precip maximum follows it (the "moist precip pattern" wire).
    lat = np.linspace(-90.0, 90.0, 361)
    p_sym = precip.precipitation(lat, 15.0, itcz_center_deg=0.0)
    p_shift = precip.precipitation(lat, 15.0, itcz_center_deg=6.0)
    assert lat[np.argmax(p_sym)] == pytest.approx(0.0, abs=0.5)
    assert lat[np.argmax(p_shift)] == pytest.approx(6.0, abs=0.7)


@pytest.mark.slow
def test_demo_reproduces_the_banked_headline():
    # Guards the committed figure (planet-sphere-itcz.png): a fresh clone reproduces the headline, not just
    # reads it. The migration is engine==closed-form, of the observed ORDER but a factor ~1.5–2 high, the
    # right sign, and the ITCZ band relocates off the equator under the imposed Q-flux.
    from planet import demo_sphere_itcz as demo
    r = demo.compute()
    assert r.slope == pytest.approx(r.slope_closed, rel=0.02)   # closed-form consequence of D
    assert r.slope < 0.0                                        # toward the warm hemisphere
    assert 1.5 < abs(r.slope) / abs(demo.OBSERVED_DEG_PER_PW) < 2.5   # observed order, factor ~2
    assert r.phi_efe_asym > 0.5                                 # ITCZ displaced into the warmed hemisphere
    assert r.lat[np.argmax(r.precip_shift)] == pytest.approx(r.phi_efe_asym, abs=1.0)
