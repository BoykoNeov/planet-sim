"""Rung 5A.4 validation: the terrain cools its own air (:mod:`planet.elevation_temperature`, plan §12.5).

The triad, said plainly (plan §3):

* **Tight (exact).** The constant-``Γ`` path *is* the closed form ``T − Γz``; flat terrain returns the
  input temperature untouched; and — the structural anchor — the **moist integrator reduces to that
  closed form to machine precision** when handed a constant ``Γ`` callable, so the emergent path and the
  prescribed path are the *same code*, not two implementations that can drift apart.
* **Tight (convergent).** The moist march is Heun in height → **second order in the step** (halving the
  step quarters the error), and its **small-``z`` limit** is ``T₀ − Γ_m(T₀, p₀)·z`` with an O(z²)
  residual. The lapse rate itself is rung 4's own object (identity-checked), so the two rungs cannot drift.
* **Conservation → a consistency leg (honestly weaker, named).** Cooling a surface without re-solving the
  EBM has no energy law to close (the :mod:`planet.biomes` precedent). The substitute is the exact
  identity ``⟨T_sea − T⟩ = Γ·⟨z⟩`` over the patch, plus the biome partition still tiling it.
* **Loose (benchmark) — and it is the leg that picked the default.** The **freezing level**: against the
  observed deep-tropical ~4.5–5 km (Harris, Bowman & Shin 2000), the pinned 6.5 K/km constant lands
  4.38 km — just *below* the band — while the emergent moist adiabat overshoots it by ~45 %, because a
  *saturated parcel* adiabat is not the *environmental* lapse rate of an unsaturated mean column. The
  decisive part is the **ordering**, not the constant landing in the band, and the assertions below say
  exactly that. Hence ``moist=False`` stays the default and the emergent path ships as a diagnostic —
  the negative is asserted here, not hidden.
"""
import numpy as np
import pytest

from planet import elevation_temperature as et, radiation

TROPICAL_SURFACE_C = 28.2          # °C — the model's deep-tropical zonal-mean surface temperature
MIDLAT_SURFACE_C = 6.6             # °C — the model's 47° zonal-mean surface temperature (the demo range)
POLAR_SURFACE_C = -14.6            # °C — the model's high-latitude zonal-mean surface temperature
RIDGE_M = 2500.0                   # m  — the demo ridge crest


# --------------------------------------------------------------------------- #
# Tight (exact): the closed form, the flat-terrain null, and the reduction of the integrator to it
# --------------------------------------------------------------------------- #
def test_constant_path_is_the_closed_form_and_reuses_rung4s_constant():
    z = np.array([0.0, 500.0, 2500.0])
    assert np.allclose(et.surface_temperature(10.0, z), 10.0 - radiation.LAPSE_RATE * z, atol=0.0)
    # the default Γ is rung 4's pinned convective-adjustment rate, not a constant invented here
    assert et.surface_temperature(0.0, 1000.0) == pytest.approx(-radiation.LAPSE_RATE * 1000.0)


def test_flat_terrain_returns_the_input_temperature_exactly():
    T = np.array([-20.0, 0.0, 28.0])
    assert np.array_equal(et.terrain_temperature(T, 0.0), T)                     # constant path
    assert np.allclose(et.terrain_temperature(T, 0.0, moist=True), T, atol=1e-12)  # moist path


def test_integrator_reduces_to_the_closed_form_under_a_constant_gamma():
    # THE structural anchor: hand the emergent marcher a constant Γ and it *is* the prescribed path.
    z = np.array([[0.0, 800.0], [2500.0, 4000.0]])
    marched = et.column_temperature(12.0, z, lapse_rate_fn=et.constant_lapse_rate_fn(), n_steps=16)
    assert np.allclose(marched, et.surface_temperature(12.0, z), atol=1e-12)
    # and it holds for a non-default rate too (the reduction is in the code path, not in the number)
    marched9 = et.column_temperature(12.0, z, lapse_rate_fn=et.constant_lapse_rate_fn(9.8e-3), n_steps=16)
    assert np.allclose(marched9, et.surface_temperature(12.0, z, lapse_rate=9.8e-3), atol=1e-12)


def test_moist_lapse_rate_is_rung4s_own_object():
    # the two rungs share the function itself, so radiation.py cannot drift away from this module
    assert et.MOIST_LAPSE_RATE_FN is radiation.moist_adiabatic_lapse_rate
    # its cold limit is the dry adiabat g/c_p (no vapour left to release) — the physics sanity check
    dry = radiation.GRAVITY / radiation.C_P_AIR
    assert et.MOIST_LAPSE_RATE_FN(200.0, radiation.P_SURFACE) == pytest.approx(dry, rel=0.02)


# --------------------------------------------------------------------------- #
# Tight (convergent): Heun is second order in the height step, and the small-z limit is Γ_m(T₀, p₀)·z
# --------------------------------------------------------------------------- #
def test_moist_march_is_second_order_in_the_height_step():
    reference = float(et.column_temperature(MIDLAT_SURFACE_C, 3000.0, n_steps=4096))
    errs = [abs(float(et.column_temperature(MIDLAT_SURFACE_C, 3000.0, n_steps=n)) - reference)
            for n in (8, 16, 32, 64)]
    ratios = [a / b for a, b in zip(errs[:-1], errs[1:])]
    assert all(3.5 < r < 4.5 for r in ratios), ratios                      # halving the step quarters it
    assert errs[-1] < 1e-3                                                 # the shipped DEFAULT_STEPS=64


def test_small_height_limit_is_the_surface_moist_lapse_rate():
    gamma0 = float(radiation.moist_adiabatic_lapse_rate(MIDLAT_SURFACE_C + et.KELVIN_0C, radiation.P_SURFACE))
    resid = []
    for z in (400.0, 200.0, 100.0):
        marched = float(et.column_temperature(MIDLAT_SURFACE_C, z, n_steps=256))
        resid.append(abs(marched - (MIDLAT_SURFACE_C - gamma0 * z)))
    ratios = [a / b for a, b in zip(resid[:-1], resid[1:])]
    assert all(3.5 < r < 4.5 for r in ratios), ratios                      # O(z²) departure from linear


def test_temperature_falls_monotonically_with_height_on_both_paths():
    z = np.linspace(0.0, 6000.0, 41)
    for T in (et.surface_temperature(20.0, z), et.column_temperature(20.0, z)):
        assert np.all(np.diff(T) < 0.0)


# --------------------------------------------------------------------------- #
# Conservation → the consistency leg: the exact patch-mean identity (there is no energy law here)
# --------------------------------------------------------------------------- #
def test_patch_mean_cooling_is_exactly_gamma_times_mean_elevation():
    rng = np.random.default_rng(5)
    z = rng.uniform(0.0, 3000.0, size=(23, 17))
    T_sea = np.full(z.shape, 9.0)
    cooling = T_sea - et.surface_temperature(T_sea, z)
    assert float(cooling.mean()) == pytest.approx(radiation.LAPSE_RATE * float(z.mean()), rel=1e-12)


def test_elevation_out_of_range_is_refused_not_clamped():
    with pytest.raises(ValueError, match="MAX_ELEVATION_M"):
        et.surface_temperature(10.0, et.MAX_ELEVATION_M + 1.0)
    with pytest.raises(ValueError, match="must be >= 0"):
        et.terrain_temperature(10.0, -1.0)


# --------------------------------------------------------------------------- #
# Loose (calibrated): the emergent rate CONFIRMS the constant at mid-latitudes and diverges away from it
# --------------------------------------------------------------------------- #
def test_emergent_rate_confirms_the_pinned_constant_at_midlatitudes():
    # the bet was that making Γ emergent would retire the 6.5 K/km pin; at the demo range it reproduces it
    gamma_mid = float(et.effective_lapse_rate(MIDLAT_SURFACE_C, RIDGE_M)) * 1000.0
    assert 6.0 < gamma_mid < 6.6                                            # ≈ 6.31 — within ~3 % of 6.5
    assert abs(gamma_mid - radiation.LAPSE_RATE * 1000.0) < 0.5


def test_emergent_rate_diverges_strongly_away_from_midlatitudes():
    gamma_trop = float(et.effective_lapse_rate(TROPICAL_SURFACE_C, RIDGE_M)) * 1000.0
    gamma_pole = float(et.effective_lapse_rate(POLAR_SURFACE_C, RIDGE_M)) * 1000.0
    assert gamma_trop < 4.5                                                 # ≈ 3.7 — latent heat flattens it
    assert gamma_pole > 8.0                                                 # ≈ 8.8 — toward the dry adiabat
    assert gamma_pole / gamma_trop > 2.0                                    # the contrast this model predicts


def test_effective_rate_is_nan_where_there_is_no_column():
    assert np.isnan(float(et.effective_lapse_rate(10.0, 0.0)))


# --------------------------------------------------------------------------- #
# Loose (benchmark): the freezing level — the leg on which the constant beat the emergent adiabat
# --------------------------------------------------------------------------- #
def test_freezing_level_constant_path_is_the_closed_form_and_floors_below_zero():
    assert float(et.freezing_level(13.0)) == pytest.approx(13.0 / radiation.LAPSE_RATE)
    assert float(et.freezing_level(-3.0)) == 0.0                            # already freezing at the surface
    assert float(et.freezing_level(np.array(-3.0), moist=True)) == 0.0


def test_tropical_freezing_level_the_constant_is_close_and_the_moist_adiabat_is_far():
    """The negative that decided the default: the emergent path is *worse* against the one observation.

    The observed deep-tropical 0 °C isotherm is the highest on the planet, cited at ~5 km (Harris, Bowman
    & Shin 2000, a 20-yr NCEP + TRMM climatology — [[smith-barstad-orographic-source]]). The pinned
    6.5 K/km constant puts it at ~4.4 km — just *below* the band, close but not inside it; the saturated
    moist adiabat, effective ~3.7 K/km there, puts it above 7 km. A *parcel* adiabat is not the
    *environmental* lapse rate of an unsaturated mean column, so this is a scope failure of the
    idealisation, not a bug. What is asserted is the **ordering** plus loose bands — deliberately not
    "the constant lands in the band", which the numbers do not support.
    """
    lo, hi = et.OBSERVED_TROPICAL_FREEZING_LEVEL_M
    fixed = float(et.freezing_level(TROPICAL_SURFACE_C))
    moist = float(et.freezing_level(np.array(TROPICAL_SURFACE_C), moist=True))
    assert 0.9 * lo < fixed < hi                                            # ≈ 4.4 km — close, just under `lo`
    assert moist > 1.2 * hi                                                 # ≈ 7.1 km — ~45 % above the band
    assert moist - hi > lo - fixed                                          # THE verdict: one close, one far


def test_freezing_level_is_monotone_in_surface_temperature():
    T = np.array([-5.0, 0.0, 5.0, 15.0, 28.0])
    for z in (et.freezing_level(T), et.freezing_level(T, moist=True)):
        assert np.all(np.diff(z) >= 0.0)
