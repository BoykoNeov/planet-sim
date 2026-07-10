"""Triad for the 2-D (lat×lon) seasonal EBM + the continentality map (:mod:`planet.seasonal_map`, rung 5B.2).

*Tight* = (a) the **per-cell 0-D slab** (transport off): every grid point relaxes as an independent slab, so
a land cell and an ocean cell reproduce the analytic amplitude ``F₁/√(B²+ω²C²)`` for their own ``C``; (b)
**zonal invariance** — a zonally-uniform mask keeps the field exactly flat in longitude (the zonal sweep is
the identity), and all-land / all-ocean collapse **bit-for-bit** to the 5B.1 single-field marcher (which is
itself spectral-validated); (c) the **cyclic zonal solver** returns circulant eigenvectors ``cos mλ`` exactly
(machine — the periodic assembly + Sherman–Morrison, the one new numerical object); (d) the **annual-mean
reduction** — the NMS83 headline — ``⟨T⟩`` is zonally flat *for any mask* and equals the 1-D parent
:class:`~planet.sphere_ebm.SphereEBM` steady solve (the land/sea contrast lives entirely in the seasonal
amplitude), to the marcher's convergence; (e) **hemispheric antisymmetry** with a symmetric mask.
*Conservation* = the global-and-annual net TOA ≈ 0. *Loose (calibrated)* = the continentality **map** — a
continental interior swings far more than its coast, which swings more than open ocean, at the same latitude.
"""
import numpy as np
import pytest

from planet import seasonal_map as sm
from planet.seasonal import SeasonalEBM, slab_amplitude_lag
from planet.sphere_ebm import SphereEBM


# --------------------------------------------------------------------------- #
# TIGHT (a) — the per-cell 0-D slab (the mechanism, one cell at a time).
# --------------------------------------------------------------------------- #
def test_zero_transport_per_cell_slab():
    # D=0 ⟹ every (φ, λ) cell is an independent slab on its own C. Under a pure sinusoidal forcing a land
    # cell must reproduce the analytic amplitude for C_land and an ocean cell that for C_ocean — the
    # continentality mechanism, pinned per cell (and proving the marcher's per-cell radiation is exact).
    x = SeasonalEBM(n_cells=20).x
    lon = (np.arange(8) + 0.5) * 2.0 * np.pi / 8
    mask = sm.box_mask(x, lon, (20.0, 60.0), (0.0, 180.0))       # land in a lon half, ocean in the other
    m = sm.SeasonalMapEBM(land_mask=mask, D=0.0, n_cells=20, n_lon=8, n_steps=360)
    F0, F1 = 300.0, 120.0
    forcing = F0 + F1 * np.cos(2.0 * np.pi * np.arange(360) * m.dt / (365.25 * 86400.0))
    absorbed = np.broadcast_to(forcing, (20, 360)).copy()
    c = m.march(absorbed=absorbed, tol=1e-8, max_years=40)
    amp = c.amplitude()
    i = m.nearest_index(40)
    land_col = int(np.where(m.land_mask[i])[0][0])
    ocean_col = int(np.where(~m.land_mask[i])[0][0])
    ampL, _ = slab_amplitude_lag(m.C_land, F1)
    ampO, _ = slab_amplitude_lag(m.C_ocean, F1)
    assert amp[i, land_col] == pytest.approx(ampL, abs=5e-3)
    assert amp[i, ocean_col] == pytest.approx(ampO, abs=5e-3)


# --------------------------------------------------------------------------- #
# TIGHT (b) — zonal invariance + reduction to the 5B.1 single-field marcher.
# --------------------------------------------------------------------------- #
def test_zonal_invariance_preserved_for_zonally_uniform_mask():
    # A zonally-uniform mask has no longitude structure, so a longitude-flat field stays flat for all time:
    # the zonal sweep is the identity on it. Machine-tight — proves the periodic solve doesn't leak
    # spurious zonal structure.
    x = SeasonalEBM(n_cells=48).x
    mask = sm.zonal_band_mask(x, n_lon=12, bands_deg=[(20.0, 70.0), (-70.0, -20.0)])
    m = sm.SeasonalMapEBM(land_mask=mask, n_cells=48, n_lon=12, n_steps=180)
    c = m.march(tol=1e-7, max_years=80)
    assert np.max(np.abs(c.T - c.T[:, 0:1, :])) < 1e-11


@pytest.mark.parametrize("land, field", [(True, "T_land"), (False, "T_ocean")])
def test_uniform_planet_reduces_to_5b1_marcher(land, field):
    # THE bridge to rung 5B.1: an all-land (or all-ocean) 2-D planet is the 5B.1 single-field marcher at
    # land_fraction 1 (or 0), bit-for-bit — same grid, operator, insolation, dt (the 2-D model wraps the
    # 5B.1 model as its source of truth). Since 5B.1's marcher is validated against its exact spectral
    # solve, this transitively inherits the anti-damping guarantee for the meridional + time integration.
    n_cells, n_steps = 60, 180
    mask = sm.uniform_mask(n_cells, 6, land=land)
    cmap = sm.SeasonalMapEBM(land_mask=mask, n_cells=n_cells, n_lon=6, n_steps=n_steps).march(
        tol=1e-8, max_years=120)
    ref = SeasonalEBM(land_fraction=1.0 if land else 0.0, n_cells=n_cells, n_steps=n_steps).march(
        tol=1e-8, max_years=120)
    assert cmap.converged and ref.converged
    assert np.max(np.abs(cmap.T[:, 0, :] - getattr(ref, field))) < 1e-6


# --------------------------------------------------------------------------- #
# TIGHT (c) — the periodic zonal solver reproduces circulant eigenmodes exactly.
# --------------------------------------------------------------------------- #
def test_thomas_columns_matches_solve_banded_with_varying_diagonal():
    # The batched meridional solver, pinned DIRECTLY on the case its transitive anchors don't reach: a
    # non-uniform diagonal (the varying C the mask creates down a column). Compare to scipy's banded LU on
    # the same tridiagonal, machine-tight — so the hand-rolled Thomas sweep is proven, not just inferred.
    from scipy.linalg import solve_banded
    rng = np.random.default_rng(0)
    n, m = 25, 4
    lower = np.concatenate([[0.0], rng.uniform(-2, -0.5, n - 1)])   # sub (lower[0] unused)
    upper = np.concatenate([rng.uniform(-2, -0.5, n - 1), [0.0]])   # sup (upper[-1] unused)
    diag = rng.uniform(6.0, 10.0, (n, m))                          # distinct, dominant diagonal per column
    rhs = rng.standard_normal((n, m))
    got = sm._thomas_columns(lower, diag, upper, rhs)
    for k in range(m):
        ab = np.zeros((3, n))
        ab[0, 1:] = upper[:-1]; ab[1, :] = diag[:, k]; ab[2, :-1] = lower[1:]
        assert np.max(np.abs(got[:, k] - solve_banded((1, 1), ab, rhs[:, k]))) < 1e-10


def test_cyclic_solver_reproduces_circulant_eigenmodes():
    # cos(mλ) is an exact eigenvector of the constant-coefficient periodic Laplacian, so the cyclic solve
    # (diag(C) − Δt L_λ) x = C·cos(mλ) must return cos(mλ)/(C + 2a(1−cos mΔλ))·C — machine-tight. This is
    # the sole check on the Sherman–Morrison periodic solve (the one new numerical object). m=0 (the
    # constant, the operator's null space) also confirms a flat field is returned unchanged.
    n_lon = 32
    dlon = 2.0 * np.pi / n_lon
    lon = (np.arange(n_lon) + 0.5) * dlon
    C = np.full((1, n_lon), 2.2e8)
    a = 1.7e6
    for mwave in (0, 1, 3, 7):
        field = np.cos(mwave * lon)[None, :]
        x = sm._cyclic_thomas_rows(np.array([-a]), C + 2.0 * a, C * field)
        denom = C + 2.0 * a * (1.0 - np.cos(mwave * dlon))
        assert np.max(np.abs(x - C * field / denom)) < 1e-9


# --------------------------------------------------------------------------- #
# TIGHT (d) — the annual-mean reduction (the NMS83 headline).
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_annual_mean_is_zonally_flat_and_equals_parent():
    # NMS83 headline: average over the year and C cancels, so ⟨T⟩ solves the annual-mean EBM with
    # longitude-INDEPENDENT forcing ⟨S⟩(1−α) — hence zonally uniform and equal to the 1-D parent
    # SphereEBM, FOR ANY MASK. The land/sea contrast is entirely in the seasonal amplitude; the annual
    # mean is blind to the mask. (Convergence-limited, not machine: the marcher carries an O(dt) splitting
    # residual in the shape — a spectral 2-D solve would make it machine-tight; that is deferred.)
    n_cells = 60
    x = SeasonalEBM(n_cells=n_cells).x
    lon = (np.arange(36) + 0.5) * 2.0 * np.pi / 36
    mask = sm.box_mask(x, lon, (15.0, 75.0), (40.0, 160.0))     # a real continent (genuine λ structure)
    m = sm.SeasonalMapEBM(land_mask=mask, n_cells=n_cells, n_lon=36, n_steps=300)
    c = m.march(tol=1e-5, max_years=60)
    amean = c.annual_mean()                                     # [n_x, n_lon]
    zonal_spread = float(np.max(amean.max(axis=1) - amean.min(axis=1)))
    assert zonal_spread < 0.5                                   # ~zonally flat (vs the ~37 K seasonal range)
    Sbar = m.absorbed_series().mean(axis=1)
    ref = SphereEBM(A=m.A, B=m.B, D=m.D, n_cells=n_cells).steady_linear(lambda xx, T: Sbar)
    assert np.max(np.abs(amean.mean(axis=1) - ref.T)) < 0.3     # equals the annual-mean parent


# --------------------------------------------------------------------------- #
# TIGHT (e) — hemispheric antisymmetry under a symmetric mask.
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_hemispheric_antisymmetry_with_symmetric_mask():
    # A mask symmetric about the equator + the antisymmetric seasonal forcing ⟹ T(x,λ,t) = T(−x,λ,t+½yr):
    # NH summer is SH summer half a year later. A structural check the split preserves.
    x = SeasonalEBM(n_cells=60).x
    mask = sm.zonal_band_mask(x, n_lon=8, bands_deg=[(25.0, 65.0), (-65.0, -25.0)])
    m = sm.SeasonalMapEBM(land_mask=mask, n_cells=60, n_lon=8, n_steps=240)
    c = m.march(tol=1e-7, max_years=80)
    half = m.n_steps // 2
    flipped_shifted = np.roll(c.T[::-1, :, :], half, axis=2)
    assert np.max(np.abs(c.T - flipped_shifted)) < 1e-3


# --------------------------------------------------------------------------- #
# CONSERVATION — global + annual energy balance.
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_global_annual_energy_balance():
    # Over one converged year the global-and-annual-mean net TOA ⟨⟨S(1−α) − A − B T⟩⟩ ≈ 0: each backward-
    # Euler sweep conserves ∫ C T dA (no-flux / periodic operators, uniform area element a²dxdλ).
    x = SeasonalEBM(n_cells=48).x
    lon = (np.arange(24) + 0.5) * 2.0 * np.pi / 24
    m = sm.SeasonalMapEBM(land_mask=sm.box_mask(x, lon, (10.0, 70.0), (30.0, 150.0)),
                          n_cells=48, n_lon=24, n_steps=240)
    c = m.march(tol=1e-6, max_years=90)
    net = m.absorbed_series()[:, None, :] - m.A - m.B * c.T     # equal-area cells ⟹ cell-mean = area-mean
    assert abs(float(net.mean())) < 1e-5


# --------------------------------------------------------------------------- #
# LOOSE (calibrated) — the continentality MAP (the payoff).
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_continentality_map_interior_beats_coast_beats_ocean():
    # The headline of the 2-D model: at a fixed midlatitude, a continental INTERIOR point swings far more
    # than a COASTAL point (which the adjacent ocean moderates by diffusion) which swings more than open
    # OCEAN. The continentality now varies WITHIN a latitude — a map, not a single number. Direction
    # banked; the magnitudes ride the calibrated heat capacities (loose).
    n_cells = 60
    x = SeasonalEBM(n_cells=n_cells).x
    lon = (np.arange(48) + 0.5) * 2.0 * np.pi / 48
    mask = sm.box_mask(x, lon, (15.0, 75.0), (40.0, 160.0))     # a broad continent, lon 40..160
    m = sm.SeasonalMapEBM(land_mask=mask, n_cells=n_cells, n_lon=48, n_steps=300)
    c = m.march(tol=1e-5, max_years=60)
    rng = c.seasonal_range()
    i = m.nearest_index(45)
    interior = m.nearest_lon_index(100)                        # continent center
    coast = m.nearest_lon_index(45)                            # just inside the west edge
    ocean = m.nearest_lon_index(300)                           # open ocean
    assert m.land_mask[i, interior] and m.land_mask[i, coast] and not m.land_mask[i, ocean]
    assert rng[i, interior] > rng[i, coast] > rng[i, ocean]    # interior ≫ coast ≫ ocean
    assert rng[i, interior] > 2.5 * rng[i, ocean]              # a strong continental/maritime contrast
    # and the annual means barely differ across the same points (continentality is all amplitude)
    amean = c.annual_mean()
    assert np.max(np.abs(amean[i, [interior, coast, ocean]] - amean[i].mean())) < 0.5


# --------------------------------------------------------------------------- #
# Mask builders — sanity.
# --------------------------------------------------------------------------- #
def test_mask_builders_shapes_and_zonal_uniformity():
    x = SeasonalEBM(n_cells=40).x
    lon = (np.arange(16) + 0.5) * 2.0 * np.pi / 16
    assert sm.uniform_mask(40, 16, land=True).all()
    assert not sm.uniform_mask(40, 16, land=False).any()
    band = sm.zonal_band_mask(x, 16, [(20.0, 60.0)])
    assert band.shape == (40, 16)
    assert np.all(band == band[:, 0:1])                        # zonally uniform
    box = sm.box_mask(x, lon, (10.0, 50.0), (90.0, 180.0))
    assert box.shape == (40, 16) and box.any() and not box.all()
    earth = sm.earthlike_mask(x, lon)
    assert earth.shape == (40, 16) and 0.1 < earth.mean() < 0.6
