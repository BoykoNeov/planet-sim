"""Triad for the seasonal ice-albedo feedback **on the 2-D map** + albedo maps (:mod:`planet.seasonal_map`, rung 5B.3).

*Tight (bit-identical / machine).* (a) **Warm-limit reduction** — an unreachable freezing threshold never
trips, so the ice march is bit-identical to the fixed-albedo 5B.2 march (and a per-latitude albedo passed
as a *map* is bit-identical to the same albedo passed per latitude). (b) **Reduction to 5B.1+** — an
all-land / all-ocean planet under the ice feedback collapses onto the 5B.1+ tile marcher
(``land_fraction`` 1 / 0). (c) **Zonal invariance** survives the feedback for a zonally-uniform mask.
(d) With a **fixed albedo map** the model is linear: at ``D = 0`` every cell's annual mean is its own
radiative equilibrium (machine), and with transport the zonal mean of the annual-mean map equals the 1-D
parent driven by the zonal-mean co-albedo. (e) **Hemispheric antisymmetry** under ice with a symmetric mask.
*Conservation.* Global-and-annual net TOA with the *realized* co-albedo ≈ 0.
*Loose (the payoff).* The seasonal-ice **map**: a continental interior freezes for a good part of the year
where the open ocean at the same latitude stays open, and — the new nugget — the annual mean is **no longer
blind to the mask**: the winter snow rectifies the seasonal cycle, so the continent ends colder in the annual
mean than the ocean at its latitude (zero for the fixed-albedo 5B.2 map, by its own theorem).
"""
import numpy as np
import pytest

from planet import seasonal as sea
from planet import seasonal_map as sm
from planet.ebm import ALBEDO_ICE, T_FREEZE
from planet.seasonal import SeasonalEBM
from planet.sphere_ebm import SphereEBM


def _lon(n):
    return (np.arange(n) + 0.5) * 2.0 * np.pi / n


# --------------------------------------------------------------------------- #
# TIGHT (a) — warm-limit reduction + the map form of a per-latitude albedo.
# --------------------------------------------------------------------------- #
def test_never_freezing_ice_march_is_bit_identical_to_fixed_march():
    x = SeasonalEBM(n_cells=36).x
    mask = sm.box_mask(x, _lon(12), (10.0, 60.0), (30.0, 150.0))
    m = sm.SeasonalMapEBM(land_mask=mask, n_cells=36, n_lon=12, n_steps=120)
    never_ice = lambda x, T: sea.ice_coalbedo(x, T, T_freeze=-1.0e3)
    fixed = m.march(tol=1e-9, max_years=40)
    icy = m.march(coalbedo_fn=never_ice, tol=1e-9, max_years=40)
    assert np.array_equal(icy.T, fixed.T)                    # bit-identical, same seed, same steps


def test_albedo_map_equals_per_latitude_albedo_bit_for_bit():
    # A [n_x, n_lon] map that is zonally uniform must reproduce the per-latitude path exactly — the map
    # path multiplies the same S(x,t) by the same co-albedo values.
    x = SeasonalEBM(n_cells=30).x
    mask = sm.box_mask(x, _lon(8), (0.0, 50.0), (0.0, 180.0))
    m = sm.SeasonalMapEBM(land_mask=mask, n_cells=30, n_lon=8, n_steps=90)
    alb_lat = 0.25 + 0.1 * x ** 2
    per_lat = m.march(albedo=alb_lat, tol=1e-9, max_years=40)
    as_map = m.march(albedo=np.broadcast_to(alb_lat[:, None], (30, 8)).copy(), tol=1e-9, max_years=40)
    assert np.array_equal(per_lat.T, as_map.T)


def test_masked_ice_coalbedo_with_offset_free_map_matches_the_tile_ice_coalbedo():
    x = SeasonalEBM(n_cells=20).x
    mask = sm.box_mask(x, _lon(6), (0.0, 40.0), (0.0, 120.0))
    fn = sm.masked_ice_coalbedo(sm.ice_free_albedo_map(x, mask))
    T = np.linspace(-30.0, 30.0, 20)[:, None] + np.zeros((1, 6))
    assert np.array_equal(fn(x[:, None], T), sea.ice_coalbedo(x[:, None], T))


# --------------------------------------------------------------------------- #
# TIGHT (b) — reduction to 5B.1+ (the tile marcher under ice).
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("land, field", [(True, "T_land"), (False, "T_ocean")])
def test_uniform_planet_under_ice_reduces_to_5b1_ice_marcher(land, field):
    n_cells, n_steps = 48, 120
    mask = sm.uniform_mask(n_cells, 4, land=land)
    cmap = sm.SeasonalMapEBM(land_mask=mask, n_cells=n_cells, n_lon=4, n_steps=n_steps).march(
        coalbedo_fn=sea.ice_coalbedo, T_init=15.0, tol=1e-8, max_years=200)
    ref = SeasonalEBM(land_fraction=1.0 if land else 0.0, n_cells=n_cells, n_steps=n_steps).march(
        coalbedo_fn=sea.ice_coalbedo, T_init=15.0, tol=1e-8, max_years=200)
    assert cmap.converged and ref.converged
    assert np.max(np.abs(cmap.T[:, 0, :] - getattr(ref, field))) < 1e-6
    assert (cmap.T[:, 0, :] < T_FREEZE).any()               # the reduction is exercised WITH ice present


# --------------------------------------------------------------------------- #
# TIGHT (c) — zonal invariance survives the feedback.
# --------------------------------------------------------------------------- #
def test_zonal_invariance_under_ice_for_zonal_mask():
    x = SeasonalEBM(n_cells=40).x
    mask = sm.zonal_band_mask(x, n_lon=10, bands_deg=[(20.0, 70.0), (-70.0, -20.0)])
    m = sm.SeasonalMapEBM(land_mask=mask, n_cells=40, n_lon=10, n_steps=120)
    c = m.march(coalbedo_fn=sea.ice_coalbedo, T_init=15.0, tol=1e-7, max_years=150)
    assert np.max(np.abs(c.T - c.T[:, 0:1, :])) < 1e-11


# --------------------------------------------------------------------------- #
# TIGHT (d) — the fixed albedo MAP is linear: two exact annual-mean anchors.
# --------------------------------------------------------------------------- #
def test_albedo_map_with_zero_transport_gives_per_cell_radiative_equilibrium_annual_mean():
    # D = 0 ⟹ every cell is an independent slab; with a fixed (land-brighter) albedo map its annual mean is
    # exactly ⟨S⟩(1−α_cell) − A over B — machine-tight, so the map path's forcing is pinned cell by cell.
    x = SeasonalEBM(n_cells=24).x
    mask = sm.box_mask(x, _lon(8), (-30.0, 60.0), (0.0, 180.0))
    m = sm.SeasonalMapEBM(land_mask=mask, D=0.0, n_cells=24, n_lon=8, n_steps=180)
    alb = sm.ice_free_albedo_map(x, mask, land_offset=0.08)
    c = m.march(albedo=alb, tol=1e-9, max_years=200)
    assert c.converged
    Sbar = m.zonal.insolation_series().mean(axis=1)
    expect = (Sbar[:, None] * (1.0 - alb) - m.A) / m.B
    assert np.max(np.abs(c.annual_mean() - expect)) < 2e-6


@pytest.mark.slow
def test_albedo_map_zonal_mean_of_annual_mean_equals_parent_with_zonal_mean_coalbedo():
    # With transport on, the (linear) annual-mean field solves the 2-D annual-mean EBM; its ZONAL mean solves
    # the 1-D parent driven by the ZONAL-MEAN co-albedo (the zonal operator averages to zero). And the mask is
    # now VISIBLE in the annual mean: the brighter continent is colder than the ocean at its latitude.
    n_cells = 48
    x = SeasonalEBM(n_cells=n_cells).x
    mask = sm.box_mask(x, _lon(24), (10.0, 70.0), (30.0, 150.0))
    m = sm.SeasonalMapEBM(land_mask=mask, n_cells=n_cells, n_lon=24, n_steps=240)
    alb = sm.ice_free_albedo_map(x, mask, land_offset=0.06)
    c = m.march(albedo=alb, tol=1e-6, max_years=80)
    Sbar = m.zonal.insolation_series().mean(axis=1)
    forcing = Sbar * (1.0 - alb).mean(axis=1)
    ref = SphereEBM(A=m.A, B=m.B, D=m.D, n_cells=n_cells).steady_linear(lambda xx, T: forcing)
    assert np.max(np.abs(c.annual_mean().mean(axis=1) - ref.T)) < 0.3
    i = m.nearest_index(40.0)
    anom = c.zonal_anomaly()
    assert anom[i, m.nearest_lon_index(90.0)] < -0.5          # land interior colder in the annual mean
    assert anom[i, m.nearest_lon_index(270.0)] > 0.0          # open ocean warmer


# --------------------------------------------------------------------------- #
# TIGHT (e) — hemispheric antisymmetry under ice with a symmetric mask.
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_hemispheric_antisymmetry_under_ice_with_symmetric_mask():
    x = SeasonalEBM(n_cells=48).x
    mask = sm.zonal_band_mask(x, n_lon=6, bands_deg=[(25.0, 65.0), (-65.0, -25.0)])
    m = sm.SeasonalMapEBM(land_mask=mask, n_cells=48, n_lon=6, n_steps=180)
    c = m.march(coalbedo_fn=sea.ice_coalbedo, T_init=15.0, tol=1e-7, max_years=200)
    assert c.converged and (c.T < T_FREEZE).any()
    half = m.n_steps // 2
    assert np.max(np.abs(c.T - np.roll(c.T[::-1, :, :], half, axis=2))) < 1e-3


# --------------------------------------------------------------------------- #
# CONSERVATION — global + annual energy balance with the realized co-albedo.
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_global_annual_net_toa_under_ice():
    x = SeasonalEBM(n_cells=40).x
    mask = sm.box_mask(x, _lon(16), (10.0, 70.0), (30.0, 150.0))
    m = sm.SeasonalMapEBM(land_mask=mask, n_cells=40, n_lon=16, n_steps=240)
    c = m.march(coalbedo_fn=sea.ice_coalbedo, T_init=15.0, tol=1e-7, max_years=200)
    assert c.converged
    S = m.zonal.insolation_series()[:, None, :]
    absorbed = S * sea.ice_coalbedo(m.x[:, None, None], c.T)
    net = absorbed - m.A - m.B * c.T                         # equal-area cells ⟹ cell mean = area mean
    assert abs(float(net.mean())) < 5e-3                     # the O(dt) diagnostic-sampling residual


# --------------------------------------------------------------------------- #
# LOOSE (the payoff) — the seasonal-ice map + the annual-mean rectification.
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_ice_map_land_freezes_where_ocean_stays_open_and_the_mean_sees_the_mask():
    n_cells = 45
    x = SeasonalEBM(n_cells=n_cells).x
    mask = sm.box_mask(x, _lon(36), (15.0, 75.0), (40.0, 160.0))    # one broad NH continent
    m = sm.SeasonalMapEBM(land_mask=mask, n_cells=n_cells, n_lon=36, n_steps=240)
    c = m.march(coalbedo_fn=sea.ice_coalbedo, T_init=15.0, tol=1e-6, max_years=200)
    assert c.converged
    i = m.nearest_index(55.0)
    interior, ocean = m.nearest_lon_index(100.0), m.nearest_lon_index(280.0)
    frac = c.ice_fraction()
    assert frac[i, interior] > 0.25                           # the interior freezes a good part of the year
    assert frac[i, ocean] < 0.05                              # the open ocean at 55° stays open
    # winter snow rectifies the cycle: the continent is COLDER in the annual mean than the ocean at 55°
    anom = c.zonal_anomaly()
    assert anom[i, interior] < anom[i, ocean] - 0.8
    # …an effect the fixed-albedo (5B.2) map cannot have, by its own theorem
    flat = m.march(tol=1e-6, max_years=80)
    assert abs(flat.zonal_anomaly()[i, interior] - flat.zonal_anomaly()[i, ocean]) < 0.5
    # land ice is seasonal (melts every summer) — no perennial land ice at this latitude
    assert not (c.T[i, interior] < T_FREEZE).all()


# --------------------------------------------------------------------------- #
# API guards.
# --------------------------------------------------------------------------- #
def test_api_guards():
    x = SeasonalEBM(n_cells=12).x
    m = sm.SeasonalMapEBM(n_cells=12, n_lon=4, n_steps=30)
    with pytest.raises(ValueError, match="exclusive"):
        m.march(coalbedo_fn=sea.ice_coalbedo, albedo=0.3)
    with pytest.raises(ValueError, match="albedo map must be shape"):
        m.march(albedo=np.zeros((12, 5)))
    with pytest.raises(ValueError, match="T_init field must be shape"):
        m.march(T_init=np.zeros((5, 4)))
    assert m.coalbedo_map(0.3).shape == (12, 4)
    assert sm.ice_free_albedo_map(x, sm.uniform_mask(12, 4), land_offset=0.1).shape == (12, 4)
