"""Rung 5A.3 validation: lee moisture depletion — the rain shadow that drops *below* baseline (§12.5).

The 5A.3 moisture-budget triad (:mod:`planet.orographic_depletion`), on top of the 5A.2 scene
(:mod:`planet.tests.test_orographic_scene`). What is **tight** is *exact* (conservation — the water
removed from the along-wind flux equals the orographic water rained, so ``DR = 1 − g_lee`` identically;
and *reduction* — zero orographic precip gives ``g ≡ 1``, recovering 5A.2's enhancement-only combination
bit-for-bit) and *structural* (``g`` monotone non-increasing downwind; the depletion lands in the
**lee, not the windward side** — the integration-direction guard, the analogue of Rung 5A's ``sgn(σ)``
branch). The **directional payoff** is the lee total falling **below** the zonal-mean baseline — a real
rain-shadow desert. What is **loose** is the absolute lee dryness, set by the incoming column water
:data:`~planet.orographic_depletion.PWV_IN_MM` and calibrated through the cited drying ratio.
"""
import numpy as np
import pytest

from planet import demo_biomes, orographic as og, orographic_depletion as ogd, orographic_scene as sc
from planet.albedo import EBMParams

COARSE = dict(n_cells=40, n_tau=0.25)


@pytest.fixture(scope="module")
def result():
    """A fast, coarse Phase-2 zonal-mean climate — the baseline the depletion drains."""
    return demo_biomes.compute(EBMParams(n_cells=COARSE["n_cells"]), n_tau=COARSE["n_tau"])


def _cascades_scene(result, *, deplete, **kw):
    lat_c, lon_c = sc.DEMO_RANGES["cascades"]
    lat, lon = sc.regional_grid(lat_c, lon_c, 6.0, 6.0, 41, 161)
    elev = sc.meridional_ridge(lat, lon, lon_center=lon_c, amplitude_m=2500.0)
    scene = sc.build_scene(result, lat, lon, elev, speed=15.0, direction_deg=270.0,
                           lat_ref_deg=lat_c, deplete=deplete, **kw)
    return scene, elev


# --------------------------------------------------------------------------- #
# Tight (exact): conservation, and reduction to the 5A.2 enhancement-only limit
# --------------------------------------------------------------------------- #
def test_conservation_water_removed_equals_water_rained():
    # Per streamline: U·W₀·(1 − g_out) == ∫ P_oro dx  (exact by the cumulative sum) → DR = 1 − g_lee.
    rng = np.random.default_rng(0)
    p = np.abs(rng.normal(size=(5, 40)))                  # a positive orographic rate (mm/hr)
    dx, U, W = 3000.0, 12.0, 25.0
    g = ogd.depletion_factor(p, dx, U, pwv_in_mm=W)
    removed = (p / og.SECONDS_PER_HOUR).sum(axis=1) * dx  # kg/(m·s) rained per row
    assert np.allclose(U * W * (1.0 - g[:, -1]), removed, rtol=1e-12)   # downwind end holds all the loss
    # DR (max over rows) is exactly 1 − g_lee for the wettest streamline
    DR = ogd.drying_ratio(p, dx, U, pwv_in_mm=W)
    assert DR == pytest.approx(1.0 - g.min(), rel=1e-12)


def test_zero_orography_gives_unit_g_and_recovers_5A2(result):
    # No orographic precip → g ≡ 1 (exactly) → the enhancement-only 5A.2 combination, bit-for-bit.
    flat = np.zeros((5, 40))
    assert np.array_equal(ogd.depletion_factor(flat, 3000.0, 12.0), np.ones_like(flat))
    # The default (deplete=False) is the strict enhancement-only limit: g ≡ 1, so the total never dips
    # below baseline and there is no lee desert — the opt-in is a strict superset of 5A.2.
    enh, _ = _cascades_scene(result, deplete=False)
    assert np.array_equal(enh.depletion_factor, np.ones_like(enh.precip_cm))
    assert enh.lee_desert_fraction == 0.0
    # And a huge incoming column (W₀ → ∞) drives the depletion term to zero → matches enhancement-only.
    dep_weak, _ = _cascades_scene(result, deplete=True, pwv_in_mm=1e9)
    assert np.allclose(dep_weak.precip_cm, enh.precip_cm, atol=1e-3)


# --------------------------------------------------------------------------- #
# Tight (structural): g monotone non-increasing downwind; depletion in the LEE not the windward side
# --------------------------------------------------------------------------- #
def test_g_is_monotone_non_increasing_downwind():
    rng = np.random.default_rng(1)
    p = np.abs(rng.normal(size=(4, 50)))
    g = ogd.depletion_factor(p, 3000.0, 12.0)             # westerly: downwind = +x
    assert np.all(np.diff(g, axis=1) <= 1e-12)            # never increases toward the lee


def test_depletion_lands_in_the_lee_not_the_windward_side():
    # A single ridge of rain at the centre. The drained column (g < 1) must sit DOWNWIND of it.
    p = np.zeros((1, 41)); p[0, 20] = 5.0                 # all the orographic rain at column 20
    g_west = ogd.depletion_factor(p, 3000.0, +12.0)[0]    # westerly: downwind = east (larger index)
    assert g_west[10] == pytest.approx(1.0)               # upwind (west) undrained
    assert g_west[30] < 0.999                             # lee (east) drained
    # reverse the wind → the desert must flip to the west (the integration-direction guard = new sgn(σ))
    g_east = ogd.depletion_factor(p, 3000.0, -12.0)[0]
    assert g_east[30] == pytest.approx(1.0)               # now the east is upwind → undrained
    assert g_east[10] < 0.999                             # and the west is the lee → drained


def test_off_band_wind_gives_no_depletion():
    p = np.abs(np.random.default_rng(2).normal(size=(3, 30)))
    assert np.array_equal(ogd.depletion_factor(p, 3000.0, 0.0), np.ones_like(p))   # U = 0 → g ≡ 1


# --------------------------------------------------------------------------- #
# The derived refill length ≫ patch — what makes the no-refill (L → ∞) limit honest
# --------------------------------------------------------------------------- #
def test_refill_length_is_much_larger_than_the_patch():
    # Earthlike numbers: L = U·W₀/P_base ~ thousands of km, ≫ a few-hundred-km patch → refill negligible.
    L = ogd.refill_length_m(15.0, 90.0)                   # U=15 m/s, P_base=90 cm/yr, W₀ default
    assert L > 5_000e3                                    # kilometres → many thousands of km
    patch_width_m = sc.patch_spacings(47.0, 0.0, 6.0 / 160)[1] * 160   # dx·(n_lon−1)
    assert L > 20.0 * patch_width_m                       # an order of magnitude clear of the patch
    assert ogd.refill_length_m(15.0, 0.0) == float("inf")  # no baseline sink → no refill demand


# --------------------------------------------------------------------------- #
# Loose (magnitude): the drying ratio lands in the cited observed band (~0.3–0.5)
# --------------------------------------------------------------------------- #
def test_drying_ratio_is_in_the_cited_band(result):
    lat_c, lon_c = sc.DEMO_RANGES["cascades"]
    lat, lon = sc.regional_grid(lat_c, lon_c, 6.0, 6.0, 41, 161)
    elev = sc.meridional_ridge(lat, lon, lon_center=lon_c, amplitude_m=2500.0)
    dy, dx = sc.patch_spacings(lat_c, lat[1] - lat[0], lon[1] - lon[0])
    p_mm_hr = og.orographic_precip(elev, dx, dy, speed=15.0, direction_deg=270.0, latitude_deg=lat_c)
    u, _ = og.wind_components(15.0, 270.0)
    DR = ogd.drying_ratio(p_mm_hr, dx, u)
    assert 0.3 < DR < 0.5                                 # Smith et al. 2003/2005; Kirshbaum & Smith 2008


# --------------------------------------------------------------------------- #
# The payoff: the depleted lee falls BELOW the baseline — a real rain-shadow desert
# --------------------------------------------------------------------------- #
def test_depletion_drives_the_lee_below_baseline(result):
    enh, elev = _cascades_scene(result, deplete=False)
    dep, _ = _cascades_scene(result, deplete=True)
    mid = dep.precip_cm.shape[0] // 2
    crest = int(np.argmax(elev[mid, :]))
    base = dep.baseline_precip_cm[mid, 0]

    # enhancement-only leaves the lee AT the baseline; depletion pulls it BELOW.
    assert enh.precip_cm[mid, crest:].min() == pytest.approx(base, abs=1.0)
    assert dep.precip_cm[mid, crest:].min() < 0.85 * base
    # the drain is confined to the lee: the far-upwind (west) edge column is undrained (g ≈ 1),
    # the far-downwind (lee) edge is drained (g < 1) — the depletion sits downwind, not upwind.
    assert dep.depletion_factor[mid, 0] == pytest.approx(1.0, abs=1e-3)
    assert dep.depletion_factor[mid, -1] < 0.9
    # and the payoff metrics move: a real lee desert exists only with depletion on.
    assert enh.lee_desert_fraction == 0.0
    assert dep.lee_desert_fraction > 0.1
    assert dep.biome_changed_fraction >= enh.biome_changed_fraction


def test_depleted_scene_still_round_trips_and_serializes(result, tmp_path):
    from planet import planet_spec as ps
    dep, _ = _cascades_scene(result, deplete=True)
    view = sc.scene_to_view(dep)
    spec = ps.from_view(view, result.params)
    ps.save(spec, tmp_path / "dep")
    assert ps.load(tmp_path / "dep") == spec             # ADR 0004 #4 round-trip identity holds with 5A.3
