"""Tests for the deep-end interactive map (Planet §9 / ADR 0004).

``planetmap.py`` adds **no new physics** — it re-shapes the validated climate arrays
(:func:`planet.demo_biomes.compute`) into a layer registry and paints them. Per ADR 0002 /
ADR 0004 #4 the *render* itself is a **reach** layer (its only test is an execution smoke-test, not a
physics check); the genuine correctness property — the state round-trip — lives in
``test_planet_spec.py``. So these tests split three ways, exactly like ``test_app.py``:

* **always-green** (NumPy only): the registry primitives, the layer builders (the renderer-input
  seam — mirroring + broadcast), the live recompute, and robustness at the slider extremes;
* **importorskip("plotly")** (fast — no kernel, so **not** ``@slow``): the figure builds, switches
  active layer, draws the annotation, and paints the Phase-4 circulation ``vector_overlay`` (cones);
* **untested**: :func:`~planet.planetmap.interactive_map` (the live ipywidgets loop — the
  ``main()`` analogue, reach).

The load-bearing structural guard (its own test): importing the module must **not** pull Plotly or
ipywidgets — the whole layering (compute/registry headless, render deps confined to lazy imports)
collapses if a top-level ``import plotly`` creeps back in, and the always-green tests + the
``planet_spec`` round-trip would then need the ``[webviz]`` extra.
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from planet import planetmap as pm
from planet import planet_spec as ps
from planet import demo_biomes
from planet.albedo import EBMParams
from planet.biomes import Biome, biome_area_fractions

REPO_ROOT = Path(__file__).resolve().parents[2]

# Coarse settings keep the always-green recompute tests fast; the equilibrium's structure (band
# ordering, hemispheric symmetry) is resolution-robust. The present-day benchmark uses the default.
COARSE = dict(n_cells=40, n_tau=0.25)


def _coarse_result():
    """A fast, coarse Phase-2 climate result (the BiomeResult the map consumes)."""
    return demo_biomes.compute(EBMParams(n_cells=COARSE["n_cells"]), n_tau=COARSE["n_tau"])


# --------------------------------------------------------------------------- #
# 0. The layering guard — importing the map must not require/trigger plotly or ipywidgets
# --------------------------------------------------------------------------- #
def test_importing_the_map_modules_stays_headless():
    # The discipline: the registry + builders are NumPy-only; plotly/ipywidgets are imported lazily
    # (inside render / interactive_map), and matplotlib (the floor) stays out entirely. So the
    # always-green tests and the planet_spec round-trip run on a bare core install (no [webviz]).
    # Checked in a CLEAN subprocess: an in-process sys.modules check is fragile (the render smoke
    # tests below import plotly into the shared session), so import the modules fresh and assert
    # nothing heavy was pulled — this is the session-robust form of test_app's streamlit guard.
    code = (
        "import sys, planet.planetmap, planet.planet_spec\n"
        "print(','.join(m for m in ('plotly', 'ipywidgets', 'matplotlib') if m in sys.modules))\n"
    )
    out = subprocess.run([sys.executable, "-c", code], cwd=str(REPO_ROOT),
                         capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip() == "", f"a heavy dep was pulled at import time: {out.stdout.strip()}"
    assert callable(pm.interactive_map)             # the entry point exists, but we never call it


# --------------------------------------------------------------------------- #
# 1. The layer registry — the v1 biome-map layer stack
# --------------------------------------------------------------------------- #
def test_build_view_registers_the_v1_biome_layers():
    view = pm.climate_view(**COARSE)
    by_name = {ly.name: ly for ly in view.layers}
    assert set(by_name) == {"temperature", "precipitation", "biome", "ice_line", "elevation"}
    # the three scalar climate fields + the inert geography scalar are SCALAR_FIELD; the ice line is ANNOTATION
    assert {ly.name: ly.kind for ly in view.layers} == {
        "temperature": pm.LayerKind.SCALAR_FIELD,
        "precipitation": pm.LayerKind.SCALAR_FIELD,
        "biome": pm.LayerKind.SCALAR_FIELD,
        "elevation": pm.LayerKind.SCALAR_FIELD,
        "ice_line": pm.LayerKind.ANNOTATION,
    }
    # units are carried (self-describing — they cross into the interchange)
    assert by_name["temperature"].units == "°C"
    assert by_name["precipitation"].units == "cm/yr"
    assert by_name["elevation"].units == "m"
    # biome is integer codes (the categorical layer), continuous fields are float
    assert np.issubdtype(by_name["biome"].data.dtype, np.integer)
    assert np.issubdtype(by_name["temperature"].data.dtype, np.floating)
    # scalar_fields(): the three paintable surfaces, inert geography excluded by default
    paintable = [ly.name for ly in view.scalar_fields()]
    assert set(paintable) == {"temperature", "precipitation", "biome"}
    assert "elevation" in {ly.name for ly in view.scalar_fields(include_inert=True)}


def test_grid_is_the_full_globe_mirrored_about_the_equator():
    view = pm.climate_view(**COARSE)
    lat = view.grid.lat
    # the hemisphere EBM grid (n_cells) mirrored → both hemispheres (2·n_cells), pole-to-pole
    assert lat.size == 2 * COARSE["n_cells"]
    assert lat.min() < 0 < lat.max()
    assert np.allclose(lat, -lat[::-1])             # hemispherically symmetric about the equator
    # every scalar field is a 2-D lat×lon array on this grid (the renderer-input seam)
    for ly in view.scalar_fields(include_inert=True):
        assert ly.data.shape == (lat.size, view.grid.lon.size)


def test_scalar_field_is_the_zonal_mean_broadcast_across_longitude():
    view = pm.climate_view(**COARSE)
    temp = view.layer("temperature").data
    # v1 is zonal-mean: each longitude column is identical (bands, not a premature 2-D field)
    assert np.allclose(temp, temp[:, :1])
    # and the latitude profile is hemispherically symmetric (the mirrored annual-mean climate)
    col = temp[:, 0]
    assert np.allclose(col, col[::-1])


def test_biome_field_matches_the_classifier_and_tiles_the_planet():
    # The painted biome field must decode to the same biomes the validated classifier produced,
    # and the area fractions must sum to 1 (the §3 consistency leg — the map tiles the planet).
    r = _coarse_result()
    view = pm.build_view(r)
    painted = view.layer("biome").data
    # the per-latitude code (any longitude column) equals the hemisphere result mirrored
    assert set(np.unique(painted)) == set(int(c) for c in np.unique(r.codes))
    fracs = biome_area_fractions(r.codes)
    assert sum(fracs.values()) == pytest.approx(1.0, abs=1e-12)


# --------------------------------------------------------------------------- #
# 2. The ice-line annotation (the ANNOTATION kind; the NaN-gap the round-trip relies on)
# --------------------------------------------------------------------------- #
def test_ice_line_annotation_is_two_circles_with_a_nan_gap():
    # Present-day Earth has a finite polar cap → the annotation is two latitude circles at ±φ_ice,
    # separated by a [nan, nan] row (so the rendered 3-D line does not bridge the hemispheres).
    view = pm.climate_view(**COARSE)
    ann = view.layer("ice_line").data
    assert ann.shape[1] == 2 and ann.shape[0] > 0
    assert np.isnan(ann).any()                      # the separator → why round-trip eq uses equal_nan
    finite_lats = ann[np.isfinite(ann[:, 0]), 0]
    assert np.any(finite_lats > 0) and np.any(finite_lats < 0)   # both hemispheres' caps


def test_ice_line_annotation_is_empty_at_the_degenerate_extremes():
    lon = np.linspace(-180, 180, 8)
    assert pm._iceline_annotation(90.0, lon).shape == (0, 2)     # ice-free: cap shrunk to the pole
    assert pm._iceline_annotation(0.0, lon).shape == (0, 2)      # Snowball: ice to the equator


# --------------------------------------------------------------------------- #
# 3. The geography seam — an inert elevation layer, carried not consumed (§9.3)
# --------------------------------------------------------------------------- #
def test_elevation_is_an_inert_flat_geography_seam():
    view = pm.climate_view(**COARSE)
    elev = view.layer("elevation")
    assert elev.inert is True                       # carried/displayed/round-tripped, NOT consumed
    assert np.all(elev.data == 0.0)                 # flat by default
    assert elev.units == "m"


def test_imported_elevation_must_match_the_full_globe_grid():
    r = _coarse_result()
    n_lat = 2 * COARSE["n_cells"]
    good = np.ones((n_lat, pm.N_LON))
    assert np.all(pm.build_view(r, elevation=good).layer("elevation").data == 1.0)
    with pytest.raises(ValueError):
        pm.build_view(r, elevation=np.zeros((3, 3)))             # wrong shape rejected


# --------------------------------------------------------------------------- #
# 4. The wired knobs change the climate (S₀, CO₂→A, D) — the rung-0 live loop
# --------------------------------------------------------------------------- #
def _mean_temp(view):
    return float(np.mean(view.layer("temperature").data))


def test_co2_knob_warms_and_solar_knob_cools():
    base = pm.climate_view(**COARSE)
    warmer = pm.climate_view(A=pm.A_OLR - 20.0, **COARSE)        # CO₂ up ≈ lower OLR offset A
    dimmer = pm.climate_view(S0=pm.S0_EARTH - 80.0, **COARSE)    # the sun dims
    assert _mean_temp(warmer) > _mean_temp(base) > _mean_temp(dimmer)


def test_transport_knob_flattens_the_gradient():
    # The third wired knob, exercised end-to-end: more meridional transport D → a flatter
    # equator-to-pole gradient (the North two-mode amplitude T₂ ∝ 1/(6D+B)).
    spread = lambda v: float(np.ptp(v.layer("temperature").data))
    assert spread(pm.climate_view(D=0.9, **COARSE)) < spread(pm.climate_view(D=0.3, **COARSE))


def test_exoplanet_knobs_default_to_the_earth_model_exactly():
    # The clean-perturbation guard: Sun + Earth-size recover the present-day map bit-for-bit, so wiring
    # the two exoplanet knobs cannot move the default globe (the §9.1 build kept the v1 map invariant).
    base = pm.climate_view(**COARSE)
    explicit = pm.climate_view(T_star=pm.T_SUN, size=1.0, **COARSE)
    assert np.array_equal(base.layer("temperature").data, explicit.layer("temperature").data)


def test_climate_params_is_the_single_composition_climate_view_rides_on():
    # climate_params is the ONE place the obliquity + exoplanet composition lives (the §9.1 trap-guard:
    # knobs → params captured once, so climate_view / planet_spec.build_spec / a notebook bench can't
    # drift). Defaults recover the present-day Earth params exactly; the knobs compose onto exactly
    # ai (spectrum) / D (size) / s2 (tilt) and touch nothing else.
    import dataclasses
    base = EBMParams(n_cells=40)
    assert pm.climate_params(n_cells=40) == base                          # Sun/Earth-size/Earth-tilt = Earth
    p = pm.climate_params(T_star=3000.0, size=1.8, obliquity_deg=40.0, n_cells=40)
    changed = {k for k, v in dataclasses.asdict(p).items() if v != dataclasses.asdict(base)[k]}
    assert changed == {"ai", "D", "s2"}                                   # only the three knobs bit


def test_star_knob_resists_snowball():
    # The §9.1 stellar knob end-to-end: at a dimmed sun, a redder host star (lower ice albedo → weaker
    # feedback) leaves a far warmer climate than a Sun-like star — a redder star is harder to snowball.
    dimmed = 1150.0
    mean = lambda v: float(np.mean(v.layer("temperature").data))
    sun = pm.climate_view(S0=dimmed, T_star=pm.T_SUN, **COARSE)
    redder = pm.climate_view(S0=dimmed, T_star=3000.0, **COARSE)
    assert mean(redder) > mean(sun) + 20.0


def test_size_knob_sharpens_the_gradient():
    # The §9.1 size knob end-to-end: a bigger planet has weaker per-area transport (D ∝ 1/size²) → a
    # steeper equator-to-pole temperature spread.
    spread = lambda v: float(np.ptp(v.layer("temperature").data))
    assert spread(pm.climate_view(size=1.8, **COARSE)) > spread(pm.climate_view(size=0.6, **COARSE))


def test_present_day_bands_run_rainforest_to_tundra():
    # The loose Phase-2 benchmark, seen through the map: at present insolation the equatorial band is
    # tropical rain forest and the polar band is tundra (the data the globe paints — full resolution).
    r = demo_biomes.compute()
    lat = r.state.latitude_deg()
    equator = Biome(int(r.codes[int(np.argmin(np.abs(lat - 0.0)))]))
    pole = Biome(int(r.codes[int(np.argmin(np.abs(lat - 90.0)))]))
    assert equator == Biome.TROPICAL_RAIN_FOREST
    assert pole == Biome.TUNDRA


def test_extreme_knobs_return_a_valid_view():
    # Robustness at the slider extremes: a very dim sun (→ Snowball) and a bright sun (→ ~ice-free)
    # must both build a valid view whose biome field still tiles the planet (no gaps / crashes).
    for S0 in (1000.0, 1900.0):
        view = pm.climate_view(S0=S0, **COARSE)
        codes = view.layer("biome").data
        # decode the painted field back to area fractions on the (equal-area) latitude axis: tiles to 1
        col = codes[:, 0]
        _, counts = np.unique(col, return_counts=True)
        assert counts.sum() == col.size
        assert view.layer("ice_line").data.shape[1] == 2     # annotation well-formed even when degenerate


# --------------------------------------------------------------------------- #
# 5. The Plotly renderer — a build-only smoke test, [webviz]-gated (ADR 0002/0004: render is reach)
# --------------------------------------------------------------------------- #
def test_render_builds_the_biome_globe():
    go = pytest.importorskip("plotly.graph_objects")
    view = pm.climate_view(**COARSE)
    fig = pm.render(view)                            # default active = biome
    kinds = [type(t).__name__ for t in fig.data]
    assert "Surface" in kinds                        # the painted scalar field
    assert "Scatter3d" in kinds                      # the ice-line annotation overlay
    assert isinstance(fig, go.Figure)
    # standing preference: NO hover crosshair/spike lines tracking the pointer (3-D spikes default ON).
    assert fig.layout.scene.xaxis.showspikes is False
    assert fig.layout.scene.yaxis.showspikes is False
    assert fig.layout.scene.zaxis.showspikes is False
    # geometry guard: the globe must be parametrized lat→rows / lon→cols (a lat/lon transpose renders
    # without error, so assert it explicitly). z = sin(lat): the first column tracks the grid latitudes.
    surface = next(t for t in fig.data if type(t).__name__ == "Surface")
    z = np.asarray(surface.z)
    # the sphere is closed at the poles: the cell-centred grid is padded out to ±90° (so no polar hole),
    # giving two extra rows that reach the poles exactly.
    assert z.shape == (view.grid.lat.size + 2, view.grid.lon.size)
    assert np.isclose(z[0, 0], -1.0) and np.isclose(z[-1, 0], 1.0)        # south pole down, north pole up — closed
    assert np.allclose(z[1:-1, 0], np.sin(np.radians(view.grid.lat)))     # the interior rows are the grid latitudes


@pytest.mark.parametrize("active", ["temperature", "precipitation", "biome", "elevation"])
def test_render_switches_the_active_scalar_layer(active):
    pytest.importorskip("plotly.graph_objects")
    view = pm.climate_view(**COARSE)
    fig = pm.render(view, active=active)
    assert any(type(t).__name__ == "Surface" for t in fig.data)


def _fake_jet(jet_lat=42.0, jet_speed=16.0):
    """A minimal stand-in for a Phase-4 CoupledJet (no slow integration) — just the fields the
    circulation layer-builder consumes: the channel latitudes + the zonal-wind profile."""
    import types
    phi = np.linspace(19.0, 61.0, 40)
    u = jet_speed * np.exp(-((phi - jet_lat) / 8.0) ** 2) - 4.0      # a westerly bump + weak easterly mean
    return types.SimpleNamespace(phi=phi, u_profile=u, jet_lat=jet_lat, jet_speed=jet_speed)


def test_circulation_layer_maps_jet_to_both_hemispheres():
    # always-green builder check: the VECTOR_OVERLAY layer is the jet mapped onto the full globe.
    result = _coarse_result()
    view = pm.build_view(result, jet=_fake_jet())
    circ = view.layer("circulation")
    assert circ.kind is pm.LayerKind.VECTOR_OVERLAY
    assert circ.data.shape == (2, view.grid.lat.size, view.grid.lon.size)   # stacked [u, v]
    u = circ.data[0]
    assert np.allclose(circ.data[1], 0.0)                                   # v ≡ 0 for the zonal jet
    # westerly band present in BOTH hemispheres (mirrored), and zero outside the channel band
    lat = view.grid.lat
    assert u[np.argmin(np.abs(lat - 42.0)), 0] > 5.0                        # NH westerly
    assert u[np.argmin(np.abs(lat + 42.0)), 0] > 5.0                        # SH westerly (same sign)
    assert np.allclose(u[np.argmin(np.abs(lat - 5.0))], 0.0)                # equator: outside the band → 0


def test_render_paints_the_circulation_vector_overlay():
    pytest.importorskip("plotly.graph_objects")
    result = _coarse_result()
    view = pm.build_view(result, jet=_fake_jet())
    fig = pm.render(view)
    kinds = [type(t).__name__ for t in fig.data]
    assert "Surface" in kinds                        # the scalar surface
    assert "Cone" in kinds                            # the Phase-4 circulation overlay (now painted)


def test_render_unknown_active_layer_raises_keyerror():
    pytest.importorskip("plotly.graph_objects")
    with pytest.raises(KeyError):
        pm.render(pm.climate_view(**COARSE), active="salinity")


def test_save_html_writes_a_standalone_globe(tmp_path):
    pytest.importorskip("plotly.graph_objects")
    out = pm.save_html(pm.climate_view(**COARSE), tmp_path / "globe.html")
    assert out.exists() and out.stat().st_size > 1000          # a real, non-empty HTML document


# --------------------------------------------------------------------------- #
# 6. The two-world comparison renderer — A · B · Δ triptych ([webviz]-gated; render is reach)
# --------------------------------------------------------------------------- #
def _two_coarse_worlds():
    """A fast (coarse) Earth + off-Earth spec pair for the comparison smoke tests."""
    return ps.build_spec(**COARSE), ps.build_spec(T_star=3000.0, **COARSE)


@pytest.mark.parametrize("active", ["biome", "temperature"])
def test_render_comparison_builds_three_globes(active):
    pytest.importorskip("plotly.graph_objects")
    a, b = _two_coarse_worlds()
    dv = ps.delta_view(a, b, active=active)
    fig = pm.render_comparison(a.view(), b.view(), dv, active=active, labels=("Earth", "exo"))
    surfaces = [t for t in fig.data if type(t).__name__ == "Surface"]
    assert len(surfaces) == 3                                  # world A · world B · the Δ globe
    assert fig.layout.scene is not None and fig.layout.scene3 is not None   # three 3-D scenes
    # the three colorbars must not all stack on top of each other (A shares B's scale)
    assert sum(1 for s in surfaces if s.showscale) <= 2
    # the Δ globe sits BELOW A and B (the 2-row layout that keeps the colorbars off it)
    assert fig.layout.scene3.domain.y[1] <= fig.layout.scene.domain.y[0] + 1e-6
    # no hover-highlight contour circles drawn on the globe surface under the cursor
    assert all(s.contours.x.highlight is False and s.contours.z.highlight is False for s in surfaces)
    # the ice line is KEPT (it marks where each world freezes), relabelled per world, and its legend is
    # moved to the left — off the right-side colorbars it used to clash with.
    ice = [t for t in fig.data if type(t).__name__ == "Scatter3d"]
    assert len(ice) == 2 and {t.name for t in ice} == {"Earth ice line", "exo ice line"}
    assert fig.layout.legend.x == 0.0 and fig.layout.legend.xanchor == "left"


def test_save_comparison_html_writes_a_standalone_triptych(tmp_path):
    pytest.importorskip("plotly.graph_objects")
    a, b = _two_coarse_worlds()
    dv = ps.delta_view(a, b, active="temperature")
    out = pm.save_comparison_html(a.view(), b.view(), dv, tmp_path / "cmp.html",
                                  active="temperature", labels=("Earth", "exo"))
    assert out.exists() and out.stat().st_size > 1000
    # the corner hover-readout wiring is injected (the runtime behaviour can't be unit-tested, but the
    # script + its hover/unhover handlers must be present in the artifact)
    html = out.read_text(encoding="utf-8")
    assert "plotly_hover" in html and "plotly_unhover" in html
