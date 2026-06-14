"""Tests for the flow-field interchange seam (Planet R1, §9.3 / §11) — :mod:`planet.flow_serialize`.

R1 makes a :class:`~planet.flow_globe.FlowField` **serializable** through the planet-spec schema so a
future ocean engine (the spin-out) binds on the same interchange + renderer the emergent eddy band uses
today. The **load-bearing correctness** is the round-trip identity ``load(save(spec)) == spec`` on
**both** producers (the real eddy band, ``is_global=False``; a synthetic global field,
``is_global=True``) — that one serialization handles both **regardless of origin** *is* the
producer-agnosticism claim. The render is "reach, not correctness" (smoke only).

The honesty edge is pinned in the data, not a caption: the eddy band is embedded **NH-sector-only,
zeros everywhere else** — never mirrored to the south, never wrapped around the globe. And the win32
landmine (numpy scalars don't ``json.dumps``; surviving ones break the ``==``) is guarded directly.
"""
import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from planet import flow_globe as fg
from planet import flow_serialize as fs
from planet import planet_spec as ps
from planet import planetmap as pm

REPO_ROOT = Path(__file__).resolve().parents[2]

PROV_EDDY = "planet.eddy_flux (saturated frame)"
PROV_SYN = "synthetic analytic global flow"


def _band_flow_field() -> fg.FlowField:
    """A small NH-band :class:`FlowField` (``is_global=False``) with a scalar — fast, no sim.

    The same shape :func:`planet.flow_globe.flow_field_from_eddy` produces, built by hand so the
    serialization/embedding path runs in the fast lane (the eddy sim is slow-marked, below)."""
    lat = np.linspace(20.0, 60.0, 9)
    lon = np.linspace(-25.0, 25.0, 11)
    LON, LAT = np.meshgrid(np.radians(lon), np.radians(lat))
    u = 12.0 * np.cos(LAT)
    v = 4.0 * np.sin(2.0 * LON)
    theta = 5.0 * np.sin(LAT)
    cov = fg.Coverage(20.0, 60.0, -25.0, 25.0, is_global=False)
    return fg.FlowField(lat=lat, lon=lon, u=u, v=v, coverage=cov,
                        honesty="one NH band, ~50° sector — not global; the flux mostly reverses",
                        scalar=theta, scalar_label="θ (°C)", radius_m=6.371e6)


def _synthetic_eddy(ny=12, nx=20, nf=5):
    """A minimal EddyFlux+EddyFrames stand-in (fabricated frames — no sim), for the REAL producer path."""
    a = 6.371e6
    phi = np.linspace(30.0, 50.0, ny)
    y = a * np.radians(phi - phi.mean())
    x = np.linspace(0.0, 3.0e6, nx)
    rng = np.random.default_rng(0)
    frames = SimpleNamespace(times=np.linspace(0.0, 60.0, nf), u=rng.standard_normal((nf, ny, nx)),
                             v=rng.standard_normal((nf, ny, nx)), theta=rng.standard_normal((nf, ny, nx)),
                             phi=phi, x=x, y=y)
    return SimpleNamespace(frames=frames, saturation_period=30.0, irreversible_fraction=0.08)


def _vector_layer(spec: ps.PlanetSpec) -> pm.Layer:
    return spec.view().layer(fs.VECTOR_LAYER)


# --------------------------------------------------------------------------- #
# 0. Headless — the interchange seam runs on a bare core install (always-green)
# --------------------------------------------------------------------------- #
def test_importing_flow_serialize_stays_headless():
    code = ("import sys, planet.flow_serialize\n"
            "print(','.join(m for m in ('plotly', 'ipywidgets', 'matplotlib') if m in sys.modules))\n")
    out = subprocess.run([sys.executable, "-c", code], cwd=str(REPO_ROOT), capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip() == "", f"a heavy dep was pulled at import time: {out.stdout.strip()}"


# --------------------------------------------------------------------------- #
# 1. The round-trip identity on BOTH producers — the R1 proof (ADR 0004 #4)
# --------------------------------------------------------------------------- #
def test_round_trip_identity_of_the_eddy_band(tmp_path):
    spec = fs.vector_spec_from_flow_field(_band_flow_field(), provenance=PROV_EDDY)
    ps.save(spec, tmp_path / "eddy")
    assert ps.load(tmp_path / "eddy") == spec


def test_round_trip_identity_of_the_synthetic_global(tmp_path):
    spec = fs.vector_spec_from_flow_field(fs.synthetic_flow_field(), provenance=PROV_SYN)
    ps.save(spec, tmp_path / "syn")
    assert ps.load(tmp_path / "syn") == spec


def test_round_trip_identity_of_the_real_eddy_producer(tmp_path):
    # The real producer path: flow_field_from_eddy → embed → serialize → round-trip (fabricated frames,
    # so it stays fast; the live sim is the slow test below). NaN-gapped θ (outside the band) survives.
    field = fg.flow_field_from_eddy(_synthetic_eddy())
    spec = fs.vector_spec_from_flow_field(field, provenance=PROV_EDDY)
    ps.save(spec, tmp_path / "real")
    assert ps.load(tmp_path / "real") == spec


def test_round_trip_negative_control(tmp_path):
    # The pass above must not be a trivially-true __eq__: perturb one velocity cell after load → unequal.
    spec = fs.vector_spec_from_flow_field(fs.synthetic_flow_field(), provenance=PROV_SYN)
    ps.save(spec, tmp_path / "syn")
    loaded = ps.load(tmp_path / "syn")
    _vector_layer(loaded).data[0, 0, 0] += 1.0
    assert loaded != spec


# --------------------------------------------------------------------------- #
# 2. Coverage-extent + provenance survive — and is_global distinguishes the producers
# --------------------------------------------------------------------------- #
def test_coverage_and_provenance_survive_the_round_trip(tmp_path):
    spec = fs.vector_spec_from_flow_field(_band_flow_field(), provenance=PROV_EDDY)
    ps.save(spec, tmp_path / "eddy")
    style = _vector_layer(ps.load(tmp_path / "eddy")).style
    assert style["provenance"] == PROV_EDDY
    assert style["coverage"]["is_global"] is False
    assert style["coverage"]["lat_min"] == 20.0 and style["coverage"]["lat_max"] == 60.0
    assert style["radius_m"] == 6.371e6 and style["honesty"]


def test_is_global_distinguishes_the_two_producers():
    # The coverage machinery is only EXERCISED because one producer is non-global — the named R1 point.
    band = _vector_layer(fs.vector_spec_from_flow_field(_band_flow_field(), provenance=PROV_EDDY))
    syn = _vector_layer(fs.vector_spec_from_flow_field(fs.synthetic_flow_field(), provenance=PROV_SYN))
    assert band.style["coverage"]["is_global"] is False
    assert syn.style["coverage"]["is_global"] is True


# --------------------------------------------------------------------------- #
# 3. The honesty edge, carried IN THE DATA — NH sector only, zeros elsewhere
# --------------------------------------------------------------------------- #
def test_eddy_band_is_nh_sector_only_never_mirrored_or_wrapped():
    view = fs.vector_view_from_flow_field(_band_flow_field(), provenance=PROV_EDDY)
    grid = view.grid
    u, v = view.layer(fs.VECTOR_LAYER).data
    speed = np.hypot(u, v)

    south = grid.lat < 0.0                                       # the southern hemisphere
    assert np.all(speed[south, :] == 0.0)                       # NOT mirrored to the SH
    far_lon = np.abs(grid.lon) > 40.0                           # well outside the ~50° sector
    assert np.all(speed[:, far_lon] == 0.0)                     # NOT wrapped around the globe
    in_band = ((grid.lat[:, None] >= 20.0) & (grid.lat[:, None] <= 60.0)
               & (grid.lon[None, :] >= -25.0) & (grid.lon[None, :] <= 25.0))
    assert np.any(speed[in_band] > 0.0)                         # but real flow inside the band


def test_synthetic_global_fills_the_whole_grid():
    view = fs.vector_view_from_flow_field(fs.synthetic_flow_field(), provenance=PROV_SYN)
    u, v = view.layer(fs.VECTOR_LAYER).data
    lat = view.grid.lat
    speed = np.hypot(u, v)
    assert np.any(speed[lat < 0.0] > 0.0)                       # data in BOTH hemispheres (it is global)
    assert np.any(speed[lat > 0.0] > 0.0)


def test_speed_layer_is_the_velocity_magnitude():
    view = fs.vector_view_from_flow_field(fs.synthetic_flow_field(), provenance=PROV_SYN)
    u, v = view.layer(fs.VECTOR_LAYER).data
    assert np.allclose(view.layer(fs.SPEED_LAYER).data, np.hypot(u, v))


def test_full_globe_grid_spans_the_poles_so_the_renderer_does_not_smear():
    # The reason the band-embed is on a full-globe grid: render's _polecapped is a no-op (the repeated
    # rows are zeros), so no fabricated polar flow. Pin the grid spans ±90 / ±180 exactly.
    grid = fs.vector_view_from_flow_field(_band_flow_field(), provenance=PROV_EDDY).grid
    assert grid.lat[0] == -90.0 and grid.lat[-1] == 90.0
    assert grid.lon[0] == -180.0 and grid.lon[-1] == 180.0


# --------------------------------------------------------------------------- #
# 4. The win32 landmine — every style value is JSON-native (no numpy scalar)
# --------------------------------------------------------------------------- #
def test_style_values_are_json_native_python_types():
    # numpy ints don't json.dumps (and win32's default int is int32); a surviving np scalar breaks the
    # round-trip ==. Guard that the vector layer's style is plain str/float/bool/dict — no np.generic.
    style = fs.vector_view_from_flow_field(_band_flow_field(), provenance=PROV_EDDY).layer(fs.VECTOR_LAYER).style

    def assert_native(v):
        assert not isinstance(v, np.generic), f"numpy scalar leaked into style: {type(v)}"
        if isinstance(v, dict):
            for vv in v.values():
                assert_native(vv)
    assert_native(style)
    json.dumps(style)                                            # and it actually serializes


# --------------------------------------------------------------------------- #
# 5. Producer-agnosticism — both producers render through ONE generic renderer
# --------------------------------------------------------------------------- #
def test_both_producers_render_through_one_generic_renderer():
    pytest.importorskip("plotly")
    band = fs.vector_view_from_flow_field(_band_flow_field(), provenance=PROV_EDDY)
    syn = fs.vector_view_from_flow_field(fs.synthetic_flow_field(), provenance=PROV_SYN)
    # the same render(active="speed") consumes both — no per-producer special-casing (the claim)
    for view in (band, syn):
        fig = pm.render(view, active="speed", caption="test caption")
        assert fig is not None and len(fig.data) >= 2                  # a speed surface + ≥1 cone trace
    # structurally parallel: both carry the universal speed scalar + the circulation vector overlay
    for view in (band, syn):
        kinds = {ly.name: ly.kind for ly in view.layers}
        assert kinds[fs.SPEED_LAYER] is pm.LayerKind.SCALAR_FIELD
        assert kinds[fs.VECTOR_LAYER] is pm.LayerKind.VECTOR_OVERLAY


def test_bare_render_caption_of_a_flow_view_is_honest_by_default():
    # The honesty-by-DEFAULT guard: a flow view's VECTOR_OVERLAY carries `honesty`, so _field_caption uses
    # it WITHOUT the caller remembering caption= — a bare render() cannot silently publish the zonal-mean
    # band-vs-globe lie. (The single words below are token-safe against the caption's <br> word-wrapping.)
    view = fs.vector_view_from_flow_field(_band_flow_field(), provenance=PROV_EDDY)
    cap = pm._field_caption(view, fs.SPEED_LAYER).lower()
    assert "band" in cap and "reverses" in cap                    # the field's OWN honesty string
    assert "continents" not in cap                                # the zonal-mean clause is gone
    assert "emergent" not in cap                                  # the "emergent jet" clause is gone


def test_field_caption_without_honesty_keeps_the_zonal_mean_jet_clause():
    # The discriminator + coupler regression guard: a VECTOR_OVERLAY with NO `honesty` in its style (the
    # zonal-mean circulation_layer) keeps the emergent-jet + zonal-mean caption unchanged.
    grid = pm.Grid(np.linspace(-90.0, 90.0, 9), np.linspace(-180.0, 180.0, 5))
    view = pm.PlanetView(grid, (
        pm.Layer("speed", pm.LayerKind.SCALAR_FIELD, np.ones((9, 5)), "m/s", {}, z_order=0),
        pm.Layer("circulation", pm.LayerKind.VECTOR_OVERLAY,
                 np.stack([np.ones((9, 5)), np.zeros((9, 5))]), "m/s", {"colorscale": "RdBu_r"}, z_order=2),
    ))
    cap = pm._field_caption(view, "speed").lower()
    assert "emergent" in cap and "zonal-mean" in cap


def test_saved_html_carries_the_honest_caption_not_the_zonal_mean_one(tmp_path):
    pytest.importorskip("plotly")
    view = fs.vector_view_from_flow_field(fs.synthetic_flow_field(), provenance=PROV_SYN)
    honest = "Synthetic global flow — a producer-agnosticism probe, not a physical result."
    out = pm.save_html(view, tmp_path / "syn.html", active="speed", caption=honest)
    text = out.read_text(encoding="utf-8")
    assert "producer-agnosticism probe" in text                       # our honest caption is in
    assert "every point on a circle of latitude shares one climate" not in text   # the biome caption is out


# --------------------------------------------------------------------------- #
# 6. The real, live eddy producer end-to-end (slow — runs the short sim)
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_live_eddy_band_serializes_and_round_trips(tmp_path):
    from planet import demo_eddy_life
    field = fg.flow_field_from_eddy(demo_eddy_life.compute(nx=48, ny=48, n_frames=4).eddy)
    spec = fs.vector_spec_from_flow_field(field, provenance=PROV_EDDY)
    ps.save(spec, tmp_path / "live")
    assert ps.load(tmp_path / "live") == spec
    assert _vector_layer(spec).style["coverage"]["is_global"] is False
