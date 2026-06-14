"""Triad for the rung-4 gray radiative-transfer column (:mod:`planet.radiation`).

* **tight (the analytic anchor)** — the numerical two-stream solver reproduces the *derived* gray
  radiative-equilibrium profile, skin and **ground** temperatures (the surface-discontinuity
  coefficient) to ~2nd order in layer thickness, with ``OLR = σTₑ⁴`` machine-exact; and the
  no-feedback slope sits near the ``4σTₑ³`` Planck touchstone.
* **real-but-loose (the unlock)** — the emergent slope decomposition ``B ≈ Planck − water-vapour``:
  the dry slope ~3.4 drops *through* climlab's prescribed 2 to ~1.3 once water vapour is on; direction
  banked, magnitude loose (rides on the water-vapour optical-depth loading, the wall).
* **plumbing / reduction** — the calibration hits the present operating point ``OLR(288)=239`` by
  construction, the emergent ``OLR(Ts)`` is locally affine (rung-0's ``A + B·T`` is its tangent), and
  that operating point is consistent with climlab's ``(A, B) = (210, 2)``.
* **named edges** — the CO₂ forcing is **saturating, not logarithmic** (a gray band), pinned as the
  honest limitation that motivates the within-rung band upgrade.
"""
import numpy as np
import pytest

from planet import radiation as rad
from planet.ebm import A_OLR, B_OLR


@pytest.fixture(scope="module")
def column():
    return rad.calibrate_column()


# --------------------------------------------------------------------------- #
# TIGHT — the numerical solver reproduces the derived analytic gray-RE solution.
# --------------------------------------------------------------------------- #
def test_numerical_solver_converges_to_the_analytic_profile():
    """rel-err of the two-stream RE solver vs the analytic profile falls ~2nd order with resolution."""
    Te, tau_s = 255.0, 4.0
    errors = []
    for n in (10, 40, 160):
        tau_mid, T_num, _, _ = rad.solve_gray_equilibrium(tau_s, Te, n)
        T_ana = rad.gray_equilibrium_temperature(tau_mid, Te)
        errors.append(np.max(np.abs(T_num - T_ana) / T_ana))
    assert errors[0] > errors[1] > errors[2]                 # monotone convergence
    assert errors[1] / errors[2] > 8.0                       # ~2nd order (4× resolution → ~16× drop)
    assert errors[2] < 2e-5


def test_olr_is_machine_exact_sigma_Te4():
    """In radiative equilibrium the outgoing longwave equals σTₑ⁴ to machine precision (conservation)."""
    Te, tau_s = 255.0, 4.0
    _, _, _, olr = rad.solve_gray_equilibrium(tau_s, Te, 80)
    assert olr == pytest.approx(rad.STEFAN_BOLTZMANN * Te ** 4, rel=1e-9)


def test_solver_recovers_the_ground_discontinuity_coefficient():
    """The numerical ground temperature converges to the derived ½Tₑ⁴(2+τ_s) surface jump.

    The ground sits warmer than the air just above it; a wrong surface coefficient would converge to
    a constant offset instead of zero. This is the recalled-coefficient guard.
    """
    Te, tau_s = 255.0, 4.0
    errs = [abs(rad.solve_gray_equilibrium(tau_s, Te, n)[2] - rad.ground_temperature(tau_s, Te))
            for n in (40, 160)]
    assert errs[0] > errs[1]                                 # converging to the derived coefficient
    assert errs[1] < 5e-3


def test_skin_temperature_is_the_optically_thin_top():
    """The analytic profile at τ=0 is the skin temperature Tₑ/2^¼; the solver's top layer approaches it.

    The top *layer* sits at τ=Δτ/2 (just optically below the skin), so it is slightly warmer than the
    τ=0 skin value and converges down to it as the resolution grows.
    """
    Te = 255.0
    skin = rad.skin_temperature(Te)
    assert rad.gray_equilibrium_temperature(0.0, Te) == pytest.approx(skin)
    top_coarse = rad.solve_gray_equilibrium(4.0, Te, 80)[1][0]
    top_fine = rad.solve_gray_equilibrium(4.0, Te, 320)[1][0]
    assert top_coarse > skin and top_fine > skin                 # the top layer is below the skin in τ
    assert (top_fine - skin) < (top_coarse - skin)               # approaching the skin with resolution


def test_no_feedback_slope_sits_near_the_planck_touchstone(column):
    """The clear-sky no-WV slope is near 4σTₑ³ (the emission-level Planck slope), well above 2.

    The advisor's bug guard: a radiative core that is right lands ~3.2–3.8 here; ~1.5 would mean the
    emission level is mis-placed.
    """
    planck = 4.0 * rad.STEFAN_BOLTZMANN * rad.emission_temperature() ** 3
    b_dry = column.feedback_slope(water_vapour=False)
    assert abs(b_dry - planck) < 0.20 * planck
    assert b_dry > B_OLR + 1.0                                # clearly above climlab's prescribed 2


# --------------------------------------------------------------------------- #
# UNLOCK — the emergent B ≈ Planck − water-vapour decomposition.
# --------------------------------------------------------------------------- #
def test_water_vapour_is_a_positive_feedback_that_lowers_B(column):
    """Turning water vapour on (τ rises with Ts via C–C) reduces the OLR slope — a positive feedback."""
    b_dry, b_moist, wv_feedback = column.feedback_decomposition()
    assert b_moist < b_dry
    assert wv_feedback > 0.5                                  # a substantial feedback


def test_the_decomposition_brackets_climlabs_B(column):
    """climlab's prescribed B=2 sits between the dry Planck slope and the moist slope — the headline.

    B_with_wv < 2 < B_no_wv: the rung-0 number is the gray emission-level Planck slope minus the
    water-vapour feedback. Magnitude loose (the exact moist landing rides on the WV loading).
    """
    b_dry, b_moist, _ = column.feedback_decomposition()
    assert b_moist < B_OLR < b_dry
    assert 0.5 < b_moist < 2.5


def test_the_decomposition_is_order_validated_against_soden_held(column):
    """The emergent slopes land at Soden & Held (2006) feedback orders — order-validated, not tuned.

    no-WV ≈ Planck |λ₀| (3.2); the WV feedback ≈ λ_wv (1.8); and the gap from the gray net to climlab's
    prescribed 2 ≈ the lapse-rate feedback |λ_LR| (0.84) the fixed lapse rate omits — so
    climlab's B ≈ Planck − WV + lapse-rate, every term independently pinned (non-circular).
    """
    b_dry, b_moist, wv_feedback = column.feedback_decomposition()
    assert abs(b_dry - rad.SH_PLANCK) < 0.6                  # no-WV slope ≈ the Planck feedback
    assert abs(wv_feedback - rad.SH_WATER_VAPOUR) < 0.6      # WV feedback ≈ Soden–Held λ_wv (clear-sky order)
    assert b_moist + rad.SH_LAPSE_RATE == pytest.approx(B_OLR, abs=0.5)   # net + omitted λ_LR ≈ climlab's 2


def test_co2_forcing_is_positive_and_grows_with_co2(column):
    """Adding CO₂ optical depth reduces OLR → a positive forcing that grows with the CO₂ amount."""
    forcings = [column.co2_forcing(co2_factor=f) for f in (2.0, 4.0, 8.0)]
    assert forcings[0] > 0.0
    assert forcings[0] < forcings[1] < forcings[2]


def test_co2_forcing_saturates_it_is_not_logarithmic(column):
    """The per-doubling forcing DECREASES at high CO₂ — a gray (saturating) band, not the log law.

    A logarithmic law (Myhre) gives constant forcing per doubling; linear-in-τ would *grow* with each
    doubling. The gray band does neither — it saturates (concave OLR(τ)) — the named within-rung edge.
    """
    olr = {m: column.outgoing_longwave(rad.PRESENT_SURFACE_T, co2_factor=m, water_vapour=False)
           for m in (2, 4, 8, 16)}
    d_2_4 = olr[2] - olr[4]
    d_4_8 = olr[4] - olr[8]
    d_8_16 = olr[8] - olr[16]
    assert d_2_4 > d_4_8 > d_8_16                             # diminishing returns (saturation)
    assert d_8_16 < 0.6 * d_2_4                               # clearly non-constant → not logarithmic


# --------------------------------------------------------------------------- #
# THE SPECTRAL-BAND LOG LAW (within-rung upgrade) — exponential wings make the CO₂ forcing
# LOGARITHMIC where the gray band saturates. Gray's saturation test above is left untouched.
# --------------------------------------------------------------------------- #
@pytest.fixture(scope="module")
def band():
    return rad.SpectralCO2Band()


def test_planck_flux_integrates_to_sigma_T4():
    """π∫B_ν dν = σT⁴ — the spectral Planck source sums (over the whole spectrum) to the gray σT⁴.

    This is what lets the per-bin spectral emission reduce to the gray column: summing πB_ν·Δν over a
    grid spanning the spectrum reconstructs the Stefan–Boltzmann flux the gray kernel uses.
    """
    n = np.linspace(1.0, 4.0e5, 400_000)                     # 0..4000 cm⁻¹ in m⁻¹
    for T in (220.0, 255.0, 288.0):
        integ = np.trapezoid(rad.planck_flux_per_wavenumber(n, T), n)
        assert integ == pytest.approx(rad.STEFAN_BOLTZMANN * T ** 4, rel=1e-3)


def test_band_kernel_reduces_to_the_gray_olr_kernel(column):
    """REDUCTION (the independent anchor): the band kernel with a σT⁴ source == the gray _olr_from.

    :func:`radiation._transmission_emission` is written independently of
    :meth:`GrayRadiationColumn._olr_from`; feeding it the gray whole-spectrum ``σT⁴`` source over the
    same column reproduces ``_olr_from`` to machine precision (the residual is the float
    multiplication-order ULP — exactly what two independent implementations agreeing looks like).
    Collapsing the spectral resolution recovers gray — the cross-check that ties the band machinery to
    the validated gray column.
    """
    Ts = rad.PRESENT_SURFACE_T
    p, T = column._profile(Ts)
    for co2 in (1.0, 2.0, 8.0):
        tau = column._optical_depth(Ts, co2, water_vapour=False)
        surface = rad.STEFAN_BOLTZMANN * float(T[-1]) ** 4
        layer = rad.STEFAN_BOLTZMANN * (0.5 * (T[:-1] + T[1:])) ** 4
        band_kernel = rad._transmission_emission(tau, surface, layer)
        assert band_kernel == pytest.approx(column._olr_from(T, tau), rel=1e-12)


def test_spectral_co2_forcing_is_logarithmic(band):
    """THE UNLOCK — the mirror of gray's saturation test: per-doubling ΔF is CONSTANT (the log law).

    Exponential band wings make the τ=1 emission level spread by a constant spectral width per CO₂
    doubling, so ΔF is constant per doubling (the Myhre law) where the gray band's ΔF *decreases*.
    Tested in the flat-middle 0.5×–8× range (the named edges sit far outside it).
    """
    per = band.forcing_per_doubling((0.5, 1, 2, 4, 8))
    assert per.max() / per.min() < 1.05                      # constant per doubling (≈ 0.5% spread)
    # and far flatter than the gray band over the same range
    col = rad.calibrate_column()
    olr = {m: col.outgoing_longwave(rad.PRESENT_SURFACE_T, co2_factor=m, water_vapour=False)
           for m in (0.5, 1, 2, 4, 8)}
    gray_per = np.array([olr[0.5] - olr[1], olr[1] - olr[2], olr[2] - olr[4], olr[4] - olr[8]])
    assert gray_per.max() / gray_per.min() > 1.5             # gray clearly saturates over the same span


def test_spectral_forcing_lands_in_the_myhre_band(band):
    """The constant per-doubling forcing is order-validated against Myhre (~3.7 W m⁻²), not tuned.

    Loose by design (advisor): the magnitude rides the band parameters (wing scale, band-centre τ,
    half-width) — calibrated to order, the wall. The *functional form* (logarithmic) is the win.
    """
    per = band.forcing_per_doubling((0.5, 1, 2, 4, 8))
    assert np.all((per > 2.0) & (per < 6.0))                 # the Myhre band, not gray's 20–53
    assert 2.0 < per.mean() < 6.0
    # same order as Myhre's 5.35·ln2 ≈ 3.71 (within a factor of ~1.6, the loose magnitude)
    assert 0.6 < per.mean() / rad.MYHRE_PER_DOUBLING < 1.6


def test_uniform_wings_saturate_the_wing_is_the_ingredient():
    """REDUCTION/null: flatten the wings (uniform k) and the forcing SATURATES like gray.

    Isolates the exponential wing as the whole ingredient — same column, same band, only the wing
    shape removed. A moderate band-centre τ so the band transitions through τ=1 across the sweep.
    """
    flat = rad.SpectralCO2Band(uniform=True, band_centre_tau=8.0)
    per = flat.forcing_per_doubling((2, 4, 8, 16, 32))
    assert per[-1] < 0.6 * per[0]                            # saturates (gray's d_8_16 < 0.6·d_2_4 style)
    # the exponential-wing band over the same span stays ~constant
    wings = rad.SpectralCO2Band().forcing_per_doubling((2, 4, 8, 16, 32))
    assert wings[-1] > 0.9 * wings[0]


def test_spectral_log_law_is_range_limited_at_low_co2():
    """NAMED EDGE: below the band-centre-saturation threshold the forcing is linear/√, not log.

    With a weakly-absorbing band centre (small τ), low CO₂ leaves the whole band optically thin, so
    ΔF GROWS per doubling (linear-in-amount) rather than holding constant — the log law's low-CO₂ edge.
    """
    weak = rad.SpectralCO2Band(band_centre_tau=4.0)
    per = weak.forcing_per_doubling((0.0625, 0.125, 0.25, 0.5, 1, 2))
    assert per[-1] > 2.0 * per[0]                            # growing, not constant → not yet logarithmic


def test_log_law_coefficient_predicts_the_right_order(band):
    """DERIVATION (consistency): the τ=1 wing formula 2l·π[B(Ts)−B(T_strat)] predicts ΔF to order.

    The cold-to-space limit; the column's finite-layer emission realizes ~20–30% more, so this is a
    derivation/consistency leg (it assumes the same exponential wing), not an independent anchor.
    """
    analytic_per_doubling = band.log_law_coefficient() * np.log(2.0)
    column_per_doubling = band.forcing_per_doubling((1, 2)).item()
    # analytic is the sharp-limit floor; the column sits above it but within ~30%
    assert analytic_per_doubling < column_per_doubling
    assert abs(analytic_per_doubling - column_per_doubling) / column_per_doubling < 0.4


def test_band_path_reduces_to_gray_end_to_end():
    """REDUCTION (end-to-end): the full ``band_olr`` path — not just the kernel — collapses to gray.

    A full-spectrum band (uniform ``k`` = the column's CO₂ optical depth, the ``(p/p_s)`` shape, on a
    pure-CO₂ ``wv_fraction=0`` column so its τ is exactly that shape) drives the actual per-bin
    assembly; summed over the whole Planck spectrum it reproduces :meth:`GrayRadiationColumn.
    outgoing_longwave` to the Planck-grid truncation (~0.1%). Exercises the per-bin Planck weighting
    that the machine-precision *kernel* test and the Myhre-band magnitude test only constrain indirectly.
    """
    col = rad.GrayRadiationColumn(total_tau=2.0, wv_fraction=0.0)     # τ is a pure (p/p_s) CO₂ shape
    full = rad.SpectralCO2Band(column=col, band_centre_cm=2000.0, half_width_cm=1999.0,
                               band_centre_tau=col.total_tau, uniform=True, n_bins=2000)
    for co2 in (1.0, 2.0):
        gray = col.outgoing_longwave(rad.PRESENT_SURFACE_T, co2_factor=co2, water_vapour=False)
        band_full = full.band_olr(rad.PRESENT_SURFACE_T, co2_factor=co2)
        assert band_full == pytest.approx(gray, rel=3e-3)            # ~Planck-grid truncation


def test_spectral_forcing_is_resolution_converged():
    """The per-doubling forcing is stable across band resolution (the 'converged ≠ validated' reflex)."""
    vals = [rad.SpectralCO2Band(n_bins=n).forcing_per_doubling((1, 2)).item() for n in (150, 300, 600)]
    assert abs(vals[1] - vals[0]) / vals[1] < 0.01
    assert abs(vals[2] - vals[1]) / vals[1] < 0.01


# --------------------------------------------------------------------------- #
# PLUMBING / REDUCTION — the emergent OLR reduces to rung-0's A + B·T near present.
# --------------------------------------------------------------------------- #
def test_calibration_hits_the_present_operating_point(column):
    """OLR(288 K) = 239 W m⁻² by construction (the 33 K greenhouse), feedback-independent at present."""
    assert column.outgoing_longwave(rad.PRESENT_SURFACE_T) == pytest.approx(rad.PRESENT_OLR, abs=0.5)
    # water vapour is unity at the reference Ts → calibration is the same with/without it
    on = column.outgoing_longwave(rad.PRESENT_SURFACE_T, water_vapour=True)
    off = column.outgoing_longwave(rad.PRESENT_SURFACE_T, water_vapour=False)
    assert on == pytest.approx(off, abs=1e-9)


def test_emergent_olr_is_locally_affine_near_present(column):
    """Near present (±3 K) the emergent OLR(Ts) is well-fit by a line — A + B·T is a valid local reduction.

    Over a *wide* range the water-vapour feedback makes OLR(Ts) visibly curved (the curvature is the
    changing slope = the feedback), so the linear reduction is explicitly a *local* tangent — the named
    "linearization breaks far from present" edge.
    """
    Ts = rad.PRESENT_SURFACE_T + np.linspace(-3.0, 3.0, 7)
    olr = np.array([column.outgoing_longwave(t) for t in Ts])
    slope, intercept = np.polyfit(Ts, olr, 1)
    residual = np.max(np.abs(olr - (slope * Ts + intercept)))
    assert residual < 0.5                                    # tangent line fits to < 0.5 W m⁻² locally


def test_reduction_is_consistent_with_climlab_constants(column):
    """The forced operating point is consistent with climlab's (A, B): A = OLR(present) − B·T̄.

    The tangent shares the present operating point with rung-0; feeding climlab's B=2 through that
    point recovers climlab's A≈210 — confirming the A/B linkage the decomposition rests on.
    """
    A_recovered, B_recovered = column.linearized_olr()
    T_bar_c = rad.PRESENT_SURFACE_T - 273.15
    # the tangent passes through the present operating point
    assert A_recovered + B_recovered * T_bar_c == pytest.approx(column.outgoing_longwave(rad.PRESENT_SURFACE_T), abs=1e-6)
    # the same operating point with climlab's prescribed B recovers climlab's A
    A_at_climlab_B = rad.PRESENT_OLR - B_OLR * T_bar_c
    assert A_at_climlab_B == pytest.approx(A_OLR, abs=4.0)


def test_calibrate_column_is_deterministic_and_physical():
    """The calibration returns a positive optical depth and is reproducible."""
    c1 = rad.calibrate_column()
    c2 = rad.calibrate_column()
    assert c1.total_tau == c2.total_tau
    assert 0.5 < c1.total_tau < 10.0


# --------------------------------------------------------------------------- #
# The banked-figure guard (slow) — a fresh clone reproduces the headline, not just reads it.
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_demo_reproduces_the_radiation_headline():
    # Guards docs/figures/planet-radiation.png: the emission temperature, the present operating point
    # matched by construction, the B = Planck − water-vapour decomposition (climlab's 2 recovered at a
    # plausible water-vapour loading), and the saturating (non-logarithmic) CO₂ forcing.
    from planet import demo_radiation as demo
    r = demo.compute()
    assert r.Te == pytest.approx(255.0, abs=1.0)                 # emission temperature ~255 K
    assert r.olr_curve[np.argmin(np.abs(r.ts_curve - rad.PRESENT_SURFACE_T))] == pytest.approx(rad.PRESENT_OLR, abs=1.0)
    assert abs(r.b_dry - r.planck_slope) < 0.20 * r.planck_slope  # no-WV slope near the Planck touchstone
    assert r.b_moist < B_OLR < r.b_dry                           # the decomposition brackets climlab's B
    assert r.wv_feedback > 0.5                                   # a substantial water-vapour feedback
    assert 0.2 < r.wv_fraction_at_climlab < 0.5                  # climlab's 2 recovered at a plausible loading
    # the CO₂ forcing saturates: the last doubling adds less than the second (concave, not logarithmic)
    assert r.co2_per_doubling[-1] < r.co2_per_doubling[1]


@pytest.mark.slow
def test_demo_reproduces_the_spectral_log_law():
    # Guards docs/figures/planet-spectral-band.png: gray's per-doubling forcing SATURATES (decreasing)
    # while the spectral band's is CONSTANT (logarithmic) and lands in the Myhre band — and the
    # cumulative spectral forcing tracks Myhre's 5.35·ln(C/C₀) far better than gray does.
    from planet import demo_spectral_band as demo
    r = demo.compute()
    assert r.gray_per_doubling[-1] < r.gray_per_doubling[1]       # gray saturates (mirror of the gray demo)
    assert r.band_per_doubling_mid.max() / r.band_per_doubling_mid.min() < 1.05  # spectral constant
    assert 2.0 < r.band_per_doubling_mean < 6.0                   # the Myhre band, not gray's 20–53
    # the spectral cumulative forcing hugs Myhre; gray departs from it badly at high CO₂
    band_err = np.max(np.abs(r.band_forcing - r.myhre_forcing))
    gray_err = np.max(np.abs(r.gray_forcing - r.myhre_forcing))
    assert band_err < 0.25 * gray_err
