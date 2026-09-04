"""Rung 5A.2 validation: placing the rain shadow on the sphere + into the biome map (plan §12.5).

The 5A.2 integration triad (:mod:`planet.orographic_scene`), on top of the Rung-5A engine anchors
(:mod:`planet.tests.test_orographic`). What is **tight** is *exact* (the mm/hr ↔ cm/yr conversion
round-trips; the tangent-plane metric is closed-form) and *structural* (a **localized hill** places its
shadow in **longitude** with ~no latitude asymmetry — the pad-safe placement anchor that pins
``lon → x``; a mountain **outside the jet band** gives a strictly-zero anomaly; the enhancement-only
combination leaves the **lee at the baseline** and lifts the **windward** above it; the mountain
**changes the biome map** — the payoff). What is **loose** is the absolute cm/yr magnitude (the
:data:`~planet.orographic_scene.OROGRAPHIC_HOURS_PER_YEAR` calibration + the cited S&B constants).

**Rung 5A.4 (2026-09-04) added here too:** the opt-in ``lapse=True`` wiring of
:mod:`planet.elevation_temperature` — the terrain finally cools its *own* air, so the crest reads alpine.
Its own triad lives in :mod:`planet.tests.test_elevation_temperature`; what is pinned *here* is the
**integration**: the default-off and flat-terrain reductions (bit-for-bit 5A.2), the exact ``Γ·z``
correction the classifier saw, the alpine payoff, its robustness to which lapse rate, and the **named
degeneracy** that the Whittaker cold bands ignore precipitation — so on the crest the rain shadow and the
lapse rate act on the same cells but only one of them can change the answer.

Honest scope note (advisor-caught): the "zonal ridge casts ≈ no shadow" idealisation does *not* hold on
the finite zero-padded engine — a lon-uniform ridge is localised by the pad into a responding block. So
the pad-safe orientation anchor is a **compact hill's latitude-symmetry**, not a zonal ridge's
amplitude (see :func:`planet.orographic_scene.gaussian_hill`).
"""
import types

import numpy as np
import pytest

from planet import biomes, demo_biomes, elevation_temperature as et, orographic as og, orographic_scene as sc
from planet import planet_spec as ps
from planet.albedo import EBMParams

# Coarse climate keeps the fast lane snappy; the scene structure is resolution-robust (cf. test_planetmap).
COARSE = dict(n_cells=40, n_tau=0.25)


@pytest.fixture(scope="module")
def result():
    """A fast, coarse Phase-2 zonal-mean climate — the baseline the orographic bonus is added onto."""
    return demo_biomes.compute(EBMParams(n_cells=COARSE["n_cells"]), n_tau=COARSE["n_tau"])


def _cascades_grid(lat_span=6.0, lon_span=6.0, n_lat=41, n_lon=161):
    lat_c, lon_c = sc.DEMO_RANGES["cascades"]
    lat, lon = sc.regional_grid(lat_c, lon_c, lat_span, lon_span, n_lat, n_lon)
    return lat_c, lon_c, lat, lon


# --------------------------------------------------------------------------- #
# Tight (exact): the unit conversion round-trips, and the tangent-plane metric is closed-form
# --------------------------------------------------------------------------- #
def test_mm_hr_cm_yr_round_trips_exactly():
    p = np.array([0.0, 0.5, 4.32, 12.0])
    assert np.allclose(sc.cm_yr_to_mm_hr(sc.mm_hr_to_cm_yr(p)), p, atol=1e-12)
    # and the forward factor is the effective-duration annualisation (not the naive 8766 h/yr)
    assert sc.mm_hr_to_cm_yr(1.0) == pytest.approx(sc.OROGRAPHIC_HOURS_PER_YEAR / 10.0)
    assert sc.OROGRAPHIC_HOURS_PER_YEAR < 8766.0                       # a fraction of the year, not all of it


def test_patch_spacings_are_the_tangent_plane_metric():
    # dy = R·Δφ (latitude-independent); dx = R·cosφ·Δλ (shrinks toward the poles)
    dy, dx = sc.patch_spacings(0.0, 1.0, 1.0)
    assert dy == pytest.approx(sc.R_EARTH_M * np.radians(1.0))
    assert dx == pytest.approx(sc.R_EARTH_M * np.radians(1.0))         # cos 0 = 1 at the equator
    dy60, dx60 = sc.patch_spacings(60.0, 1.0, 1.0)
    assert dy60 == pytest.approx(dy)                                   # dy is latitude-independent
    assert dx60 == pytest.approx(dx * np.cos(np.radians(60.0)))        # dx halves at 60° (cos 60 = 0.5)


# --------------------------------------------------------------------------- #
# Tight (structural): placement — a localized hill shadows in LON, symmetric in LAT (pins lon → x)
# --------------------------------------------------------------------------- #
def test_localized_hill_places_the_shadow_in_lon_not_lat(result):
    lat_c, lon_c, lat, lon = _cascades_grid(lat_span=8.0, lon_span=8.0, n_lat=81, n_lon=161)
    hill = sc.gaussian_hill(lat, lon, lat_c, lon_c, amplitude_m=2500.0, half_width_deg=0.5)
    scene = sc.build_scene(result, lat, lon, hill, speed=15.0, direction_deg=270.0, lat_ref_deg=lat_c)
    P = scene.orographic_precip_cm
    ci_lat, ci_lon = int(np.argmin(np.abs(lat - lat_c))), int(np.argmin(np.abs(lon - lon_c)))

    # LON (along-wind): a strong windward/lee asymmetry, peak displaced UPWIND (to smaller lon)
    row = P[ci_lat, :]
    w_lon, l_lon = row[:ci_lon].sum(), row[ci_lon:].sum()
    assert (w_lon - l_lon) / (w_lon + l_lon) > 0.2                     # materially windward-weighted in lon
    assert lon[int(np.argmax(row))] < lon_c                            # peak sits upwind of the crest

    # LAT (cross-wind): ~symmetric, no north/south displacement — the discriminator that lon is the wind axis
    col = P[:, ci_lon]
    s_lat, n_lat = col[:ci_lat].sum(), col[ci_lat:].sum()
    assert abs(s_lat - n_lat) / (s_lat + n_lat) < 0.15                 # nearly lat-symmetric
    assert lat[int(np.argmax(col))] == pytest.approx(lat_c, abs=0.2)   # no lat displacement of the peak


# --------------------------------------------------------------------------- #
# Tight (structural): the rain shadow + enhancement-only combination
# --------------------------------------------------------------------------- #
def test_meridional_ridge_casts_a_lon_shadow(result):
    lat_c, lon_c, lat, lon = _cascades_grid()
    elev = sc.meridional_ridge(lat, lon, lon_center=lon_c, amplitude_m=2500.0)
    scene = sc.build_scene(result, lat, lon, elev, speed=15.0, direction_deg=270.0, lat_ref_deg=lat_c)
    mid = scene.orographic_precip_cm.shape[0] // 2
    row = scene.orographic_precip_cm[mid, :]
    crest = int(np.argmax(elev[mid, :]))
    windward, lee = row[:crest].sum(), row[crest:].sum()
    assert lee < 0.5 * windward                                       # a genuine shadow (lee materially drier)
    assert lon[int(np.argmax(row))] < lon[crest]                      # peak upwind of the crest (correct branch)


def test_enhancement_only_shadow_is_dry_windward_lifts_never_below_baseline(result):
    lat_c, lon_c, lat, lon = _cascades_grid()
    elev = sc.meridional_ridge(lat, lon, lon_center=lon_c, amplitude_m=2500.0)
    scene = sc.build_scene(result, lat, lon, elev, speed=15.0, direction_deg=270.0, lat_ref_deg=lat_c)
    mid = scene.precip_cm.shape[0] // 2
    crest = int(np.argmax(elev[mid, :]))
    base = scene.baseline_precip_cm[mid, 0]
    peak_bonus = scene.orographic_precip_cm[mid, :].max()
    # windward is lifted ABOVE the zonal baseline (the enhancement)
    assert scene.precip_cm[mid, :crest].max() > 1.5 * base
    # the IMMEDIATE lee (just behind the crest) is the dry rain shadow → total ≈ the zonal baseline there.
    # (Note: the *full* S&B model has weak downstream lee-wave rain bands further out — a real feature, not
    # the reduced-limit clean decay; so we anchor the dryness at the shadow, not the far edge.)
    imm_lee = scene.orographic_precip_cm[mid, crest + 2:crest + 12]
    assert imm_lee.mean() < 0.15 * peak_bonus                          # the shadow is dry vs the windward crest
    assert np.allclose(scene.precip_cm[mid, crest + 2:crest + 12], base, rtol=0.1)
    # enhancement-only: the total never dips BELOW the zonal baseline (no depletion is modelled)
    assert np.all(scene.precip_cm >= scene.baseline_precip_cm - 1e-9)


def _far_band(result, lon_span, n_lon, Hw):
    """The downstream secondary band's ``(peak_cm, downwind_deg)`` for a meridional ridge, or ``(0, nan)``.

    Looks *well downwind* of the immediate lee (> 0.5° past the crest) for a residual rain band.
    """
    lat_c, lon_c = sc.DEMO_RANGES["cascades"]
    lat, lon = sc.regional_grid(lat_c, lon_c, 6.0, lon_span, 41, n_lon)
    elev = sc.meridional_ridge(lat, lon, lon_center=lon_c, amplitude_m=2500.0)
    scene = sc.build_scene(result, lat, lon, elev, speed=15.0, direction_deg=270.0, lat_ref_deg=lat_c, Hw=Hw)
    mid = scene.orographic_precip_cm.shape[0] // 2
    row = scene.orographic_precip_cm[mid, :]
    crest = int(np.argmax(elev[mid, :]))
    downwind_deg = lon - lon[crest]
    far = downwind_deg > 0.5
    if far.any() and row[far].max() > 1.0:
        j = int(np.where(far)[0][int(np.argmax(row[far]))])
        return float(row[far].max()), float(downwind_deg[j])
    return 0.0, float("nan")


def test_downstream_band_is_propagating_mode_phase_not_numerical(result):
    # The full model shows a weak secondary rain band well downwind of the shadow. This is the
    # discriminator that it is the *propagating-mode phase* (the 1 − i m H_w factor), real to the linear
    # model — NOT FFT wrap-around or a pad artifact (which the docstrings claim). Two decisive checks:
    peak_full, dist_full = _far_band(result, lon_span=6.0, n_lon=161, Hw=sc.og.HW_M)
    peak_hw0, _ = _far_band(result, lon_span=6.0, n_lon=161, Hw=0.0)
    _, dist_wide = _far_band(result, lon_span=12.0, n_lon=321, Hw=sc.og.HW_M)  # 2× domain, SAME dx

    assert peak_full > 1.0                                             # the band exists in the full model
    assert peak_hw0 == 0.0                                            # (1) it VANISHES at H_w = 0 → phase-driven
    assert dist_full == pytest.approx(dist_wide, abs=0.2)             # (2) holds its physical distance under
    #                                                                     domain doubling → not wrap / not pad


# --------------------------------------------------------------------------- #
# Tight (structural): the wind is sourced from the zonal jet — zero outside the westerly band
# --------------------------------------------------------------------------- #
def _fake_jet():
    """A minimal stand-in for a CoupledJet — a westerly channel band peaking at 40°, zero outside it."""
    return types.SimpleNamespace(phi=np.array([20.0, 40.0, 60.0]), u_profile=np.array([0.0, 15.0, 0.0]))


def test_wind_from_jet_is_westerly_in_band_zero_outside():
    jet = _fake_jet()
    speed_in, dir_in = sc.wind_from_jet(jet, 40.0)                    # jet core
    assert speed_in == pytest.approx(15.0)
    assert dir_in == pytest.approx(og.DIRECTION_WESTERLY_DEG)         # 270 = westerly
    # a southern-hemisphere range reads the same jet by |lat| symmetry
    assert sc.wind_from_jet(jet, -40.0)[0] == pytest.approx(15.0)
    # outside the band → no cross-mountain wind
    assert sc.wind_from_jet(jet, 5.0)[0] == 0.0
    assert sc.wind_from_jet(jet, 80.0)[0] == 0.0


def test_mountain_outside_jet_band_gives_zero_anomaly(result):
    jet = _fake_jet()
    lat, lon = sc.regional_grid(5.0, 100.0, 6.0, 6.0, 41, 161)        # equatorial → outside the westerlies
    elev = sc.meridional_ridge(lat, lon, amplitude_m=2500.0)
    scene = sc.build_scene(result, lat, lon, elev, jet=jet, lat_ref_deg=5.0)
    assert scene.wind_speed == 0.0
    assert np.max(np.abs(scene.orographic_precip_cm)) == 0.0          # no wind → no orographic anomaly
    assert np.array_equal(scene.precip_cm, scene.baseline_precip_cm)  # total is exactly the zonal baseline


# --------------------------------------------------------------------------- #
# The payoff: the mountain changes the biome map; and the magnitude is sane (loose band)
# --------------------------------------------------------------------------- #
def test_biome_shift_is_the_payoff(result):
    lat_c, lon_c, lat, lon = _cascades_grid()
    elev = sc.meridional_ridge(lat, lon, lon_center=lon_c, amplitude_m=2500.0)
    scene = sc.build_scene(result, lat, lon, elev, speed=15.0, direction_deg=270.0, lat_ref_deg=lat_c)
    assert scene.biome_changed_fraction > 0.05                        # the mountain moves the biome map
    # the change is on the WINDWARD side (the enhancement) — the wetter windward biome differs from baseline
    mid = scene.biome_codes.shape[0] // 2
    crest = int(np.argmax(elev[mid, :]))
    changed = scene.biome_codes[mid, :crest] != scene.baseline_biome_codes[mid, :crest]
    assert changed.any()


def test_orographic_magnitude_is_sane(result):
    # the loose-magnitude anchor: the OROGRAPHIC_HOURS_PER_YEAR calibration lands the windward bonus in a
    # few-hundred cm/yr band — NOT the thousands a naive 8766-h/yr annualisation would give (which would
    # swamp the Whittaker classifier). A guard on the calibration, not a pinned value.
    lat_c, lon_c, lat, lon = _cascades_grid()
    elev = sc.meridional_ridge(lat, lon, lon_center=lon_c, amplitude_m=2500.0)
    scene = sc.build_scene(result, lat, lon, elev, speed=15.0, direction_deg=270.0, lat_ref_deg=lat_c)
    peak = scene.orographic_precip_cm.max()
    assert 50.0 < peak < 1500.0                                       # a few hundred cm/yr — a sane surplus


# --------------------------------------------------------------------------- #
# Rung 5A.4 — the terrain cools its own air (opt-in; the module triad is test_elevation_temperature)
# --------------------------------------------------------------------------- #
def _ridge_scene(result, **kw):
    lat_c, lon_c, lat, lon = _cascades_grid()
    elev = sc.meridional_ridge(lat, lon, lon_center=lon_c, amplitude_m=2500.0)
    return elev, sc.build_scene(result, lat, lon, elev, speed=15.0, direction_deg=270.0,
                                lat_ref_deg=lat_c, **kw)


def test_lapse_is_default_off_and_the_5a2_scene_is_bit_for_bit_unchanged(result):
    # the default-off reduction (ARCHITECTURE §4.2): with lapse off nothing about 5A.2/5A.3 moves
    _, scene = _ridge_scene(result)
    assert np.array_equal(scene.temperature_C, scene.sea_level_temperature_C)
    assert np.all(scene.elevation_cooling_K == 0.0)
    assert scene.alpine_fraction == 0.0
    assert np.array_equal(scene.biome_codes, scene.sea_level_biome_codes)


def test_flat_terrain_with_lapse_on_reproduces_the_lapse_off_scene(result):
    # the other exact null: turning the correction ON over zero terrain changes nothing, bit-for-bit
    lat_c, lon_c, lat, lon = _cascades_grid()
    flat = np.zeros((lat.size, lon.size))
    off = sc.build_scene(result, lat, lon, flat, speed=15.0, direction_deg=270.0, lat_ref_deg=lat_c)
    on = sc.build_scene(result, lat, lon, flat, speed=15.0, direction_deg=270.0, lat_ref_deg=lat_c,
                        lapse=True)
    assert np.array_equal(on.temperature_C, off.temperature_C)
    assert np.array_equal(on.biome_codes, off.biome_codes)


def test_patch_cooling_is_exactly_gamma_times_the_elevation(result):
    # the consistency leg on the scene: the correction the classifier saw IS Γ·z, cell by cell
    elev, scene = _ridge_scene(result, lapse=True)
    assert np.allclose(scene.elevation_cooling_K, et.LAPSE_RATE * elev, atol=1e-12)
    assert float(scene.elevation_cooling_K.mean()) == pytest.approx(et.LAPSE_RATE * float(elev.mean()))


def test_the_crest_turns_alpine_the_5a4_payoff(result):
    # THE payoff: a 2500 m crest under a temperate-forest latitude reads as a COLD-limited biome
    elev, scene = _ridge_scene(result, lapse=True)
    mid = scene.biome_codes.shape[0] // 2
    crest = int(np.argmax(elev[mid, :]))
    assert scene.biome_codes[mid, crest] in (biomes.Biome.TUNDRA, biomes.Biome.BOREAL_FOREST)
    _, warm = _ridge_scene(result)                                    # the same rain, no cooling
    assert warm.biome_codes[mid, crest] not in (biomes.Biome.TUNDRA, biomes.Biome.BOREAL_FOREST)
    assert scene.alpine_fraction > 0.2                                # the cooling alone moves a lot of map
    # and the partition still tiles the patch (the biomes.py consistency leg, carried through)
    assert sum(biomes.biome_area_fractions(scene.biome_codes).values()) == pytest.approx(1.0)


def test_the_payoff_is_robust_to_which_lapse_rate(result):
    # the emergent moist adiabat ≈ the pinned constant at this latitude, so the story does not ride on it
    _, fixed = _ridge_scene(result, lapse=True)
    _, moist = _ridge_scene(result, lapse=True, moist_lapse=True)
    assert abs(fixed.alpine_fraction - moist.alpine_fraction) < 0.05
    assert abs(float(fixed.elevation_cooling_K.max()) - float(moist.elevation_cooling_K.max())) < 2.0


def test_cold_crest_ignores_the_orographic_rain_the_named_degeneracy(result):
    """The classifier's cold bands are precip-independent, so on the crest 5A and 5A.4 degenerate.

    :mod:`planet.biomes` states it plainly ("a very wet sub-zero climate is still called boreal here").
    Where the terrain has cooled a cell below the cold-limited threshold, the Whittaker rule stops
    reading precipitation at all — so the rain shadow and the lapse rate act on the same cells but only
    one of them can change the answer there. Named in :mod:`planet.elevation_temperature`, pinned here.
    """
    elev, scene = _ridge_scene(result, lapse=True)
    cold = scene.temperature_C < biomes.BOREAL_MAX_C
    assert cold.any()
    with_rain = biomes.classify_field(scene.temperature_C[cold], scene.precip_cm[cold])
    without_rain = biomes.classify_field(scene.temperature_C[cold], scene.baseline_precip_cm[cold])
    assert np.array_equal(with_rain, without_rain)


# --------------------------------------------------------------------------- #
# Serialization: the regional scene rides the grid-agnostic planet_spec schema (round-trip identity)
# --------------------------------------------------------------------------- #
def test_scene_view_round_trips_through_planet_spec(result, tmp_path):
    lat_c, lon_c, lat, lon = _cascades_grid()
    elev = sc.meridional_ridge(lat, lon, lon_center=lon_c, amplitude_m=2500.0)
    scene = sc.build_scene(result, lat, lon, elev, speed=15.0, direction_deg=270.0, lat_ref_deg=lat_c)
    view = sc.scene_to_view(scene)
    spec = ps.from_view(view, result.params)
    ps.save(spec, tmp_path / "scene")
    assert ps.load(tmp_path / "scene") == spec                       # ADR 0004 #4 round-trip identity


def test_scene_elevation_is_live_not_inert(result):
    # on the regional patch the elevation finally DRIVES the rain — so it is NOT the inert globe seam
    lat_c, lon_c, lat, lon = _cascades_grid()
    elev = sc.meridional_ridge(lat, lon, lon_center=lon_c, amplitude_m=2500.0)
    scene = sc.build_scene(result, lat, lon, elev, speed=15.0, direction_deg=270.0, lat_ref_deg=lat_c)
    view = sc.scene_to_view(scene)
    assert view.layer("elevation").inert is False
    assert "orographic precipitation" in [ly.name for ly in view.layers]


# --------------------------------------------------------------------------- #
# Slow: the full emergent-jet wiring (drives the shallow-water spin-up) + the demo builds end-to-end
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_demo_scene_from_real_jet_builds_and_shadows():
    scene = sc.demo_scene("cascades", n_lat=31, n_lon=121, use_jet=True)
    assert scene.wind_speed > 5.0                                     # the emergent westerly reached the patch
    assert scene.wind_direction_deg == pytest.approx(og.DIRECTION_WESTERLY_DEG)
    assert scene.orographic_precip_cm.max() > 0.0                     # it rained a shadow
    assert scene.biome_changed_fraction > 0.0                        # and moved the biome map
