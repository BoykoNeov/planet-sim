"""Tests for the O2 real-ocean producer (§9.6) — :mod:`planet.ocean_currents`.

The load-bearing claims, each pinned against **real data** (the committed 5° OSCAR subsample in
``fixtures/oscar_subsample.npz`` — see ``fixtures/README.md`` for provenance/regeneration; the raw
33 MB granule is never committed):

* the three O2-spike encodings — the 0–360→±180 longitude **rewrap** (with the column re-sort that
  keeps the axis monotone), the **NaN→mask→fill** rule (the O1 mask carries "no data", never a filled
  zero), and the ``(time, lon, lat)``→``(lat, lon)`` transpose (loader test, synthetic granule);
* the **R1 round-trip identity on real data** — the proof, now non-synthetic;
* the new **provenance honesty clause** ("real …; NOT computed by planet-sim's models"), machine-checked
  as a visible DOM element exactly like Rung C's reversibility clause — including on the banked artifact;
* the import discipline: the netCDF reader (``[ocean]`` extra) is lazy — the module and this suite's
  fixture path stay NumPy-only.
"""
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from planet import flow_globe as fg
from planet import ocean_currents as oc

REPO_ROOT = Path(__file__).resolve().parents[2]
FIXTURE = Path(__file__).parent / "fixtures" / "oscar_subsample.npz"
BANKED = REPO_ROOT / "docs" / "figures" / "planet-ocean-currents.html"
BANKED_SEASONAL = REPO_ROOT / "docs" / "figures" / "planet-ocean-currents-seasonal.html"


def _fixture_snapshot() -> oc.OceanSnapshot:
    """The committed real-data subsample, raw conventions intact (0–360 lon, NaN land)."""
    z = np.load(FIXTURE)
    return oc.OceanSnapshot(lat=z["lat"].astype(float), lon=z["lon"].astype(float),
                            u=z["u"].astype(float), v=z["v"].astype(float),
                            product=str(z["product"]), doi=str(z["doi"]), credit=str(z["credit"]),
                            date=str(z["date"]), depth_note=str(z["depth_note"]))


# --------------------------------------------------------------------------- #
# 0. Import discipline — NumPy-only at import; the reader dep is lazy
# --------------------------------------------------------------------------- #
def test_importing_ocean_currents_stays_headless_and_reader_free():
    # The [ocean] extra (h5netcdf/h5py, requests) must NOT be pulled at import time — the producer's
    # pure-array half runs on a bare core install (the fixture tests below prove it by running there).
    code = ("import sys, planet.ocean_currents\n"
            "print(','.join(m for m in ('plotly', 'ipywidgets', 'matplotlib', 'h5netcdf', 'h5py',"
            " 'requests') if m in sys.modules))\n")
    out = subprocess.run([sys.executable, "-c", code], cwd=str(REPO_ROOT), capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip() == "", f"an optional dep was pulled at import time: {out.stdout.strip()}"


# --------------------------------------------------------------------------- #
# 1. The fixture is genuinely OSCAR-shaped (raw conventions, real-data sanity)
# --------------------------------------------------------------------------- #
def test_fixture_keeps_the_raw_oscar_conventions():
    s = _fixture_snapshot()
    assert s.lon.min() >= 0.0 and s.lon.max() > 180.0          # 0–360 convention, NOT ±180
    assert s.lat.min() > -90.0 and s.lat.max() < 90.0          # cell-centred — never reaches the poles
    land_frac = float(np.isnan(s.u).mean())
    assert 0.3 < land_frac < 0.6                               # ~44% land/missing at native res
    assert float(np.nanmax(np.hypot(s.u, s.v))) < 3.5          # OSCAR's stated valid range is ±3 m/s
    assert "OSCAR" in s.product and s.date == "2020-06-01"


# --------------------------------------------------------------------------- #
# 2. The rewrap — 0–360 → ±180 with the column re-sort (the spike's monotone rule)
# --------------------------------------------------------------------------- #
def test_rewrap_makes_the_axis_monotone_and_moves_columns_with_it():
    lon = np.array([0.0, 90.0, 180.0, 270.0])
    u = np.array([[1.0, 2.0, 3.0, 4.0]])
    lon2, u2 = oc._rewrap_lon(lon, u)
    assert lon2.tolist() == [-180.0, -90.0, 0.0, 90.0]         # 180 → −180, 270 → −90
    assert u2[0].tolist() == [3.0, 4.0, 1.0, 2.0]              # the values travelled with their meridians
    # a grid already in ±180 passes through unchanged (the sort is the identity)
    lon3, u3 = oc._rewrap_lon(lon2, u2)
    assert np.array_equal(lon3, lon2) and np.array_equal(u3, u2)


def test_field_longitudes_are_monotone_pm180_on_real_data():
    field = oc.flow_field_from_ocean(_fixture_snapshot())
    assert bool(np.all(np.diff(field.lon) > 0))
    assert field.lon.min() >= -180.0 and field.lon.max() < 180.0


# --------------------------------------------------------------------------- #
# 3. NaN → mask → fill (the O1 rule: the mask carries "no data", never a filled zero)
# --------------------------------------------------------------------------- #
def test_mask_derives_from_nan_and_velocities_are_filled_finite():
    s = _fixture_snapshot()
    field = oc.flow_field_from_ocean(s)
    land = ~field.mask
    assert land.any() and field.mask.any()                     # both classes present (it IS an ocean planet)
    assert np.isfinite(field.u).all() and np.isfinite(field.v).all()   # no NaN reaches a bilinear sampler
    assert (field.u[land] == 0.0).all() and (field.v[land] == 0.0).all()
    # the mask is exactly the finiteness of the raw data (rewrapped alongside the velocities)
    _, raw_u, raw_v = oc._rewrap_lon(s.lon, s.u, s.v)
    assert np.array_equal(field.mask, np.isfinite(raw_u) & np.isfinite(raw_v))
    # valid cells keep their measured values, and the scalar is the speed of the (filled) field
    valid = field.mask
    assert np.array_equal(field.u[valid], raw_u[valid])
    assert np.array_equal(field.scalar, np.hypot(field.u, field.v))


def test_coverage_is_global_but_honestly_cell_centred():
    field = oc.flow_field_from_ocean(_fixture_snapshot())
    assert field.coverage.is_global is True                    # the contract's first true-global consumer
    assert field.coverage.lat_max < 90.0 and field.coverage.lat_min > -90.0   # OSCAR stops short of poles


# --------------------------------------------------------------------------- #
# 4. The R1 round-trip identity — the proof, now on REAL data (§9.6 O2)
# --------------------------------------------------------------------------- #
def test_round_trip_identity_on_real_ocean_data(tmp_path):
    from planet import planet_spec as ps
    from planet.flow_serialize import vector_spec_from_flow_field

    field = oc.flow_field_from_ocean(_fixture_snapshot())
    spec = vector_spec_from_flow_field(field, provenance="planet.ocean_currents (OSCAR, real data)")
    ps.save(spec, tmp_path / "ocean")
    assert ps.load(tmp_path / "ocean") == spec


def test_serialized_field_masks_the_poles_and_zeroes_the_land():
    # Integration of the O1 rules on real data: the interchange grid reaches ±90° exactly, OSCAR does
    # not — those rows must come out masked (no extrapolated "measured" polar ocean) — and every masked
    # cell of the embedded field must carry zero velocity and NaN scalar.
    from planet.flow_serialize import MASK_LAYER, SCALAR_LAYER, VECTOR_LAYER, vector_view_from_flow_field

    field = oc.flow_field_from_ocean(_fixture_snapshot())
    view = vector_view_from_flow_field(field, provenance="oscar")
    mask = view.layer(MASK_LAYER).data.astype(bool)
    assert not mask[0].any() and not mask[-1].any()            # ±90° rows: unobserved → masked
    uv = view.layer(VECTOR_LAYER).data
    assert (uv[:, ~mask] == 0.0).all()                         # no flow where there is no data
    assert np.isnan(view.layer(SCALAR_LAYER).data[~mask]).all()   # "no data", not "speed = 0"


# --------------------------------------------------------------------------- #
# 5. The provenance honesty clause — the new class, machine-checked in the DOM
# --------------------------------------------------------------------------- #
def test_honesty_leads_with_provenance_not_a_model_claim():
    h = oc.flow_field_from_ocean(_fixture_snapshot()).honesty
    assert "REAL data" in h and "OSCAR" in h and "PO.DAAC" in h and "DOI" in h
    assert "NOT computed by planet-sim" in h                   # the clause that flips the honesty class
    assert "illustrative" in h and "time-accelerated" in h     # streaming is still a stylization
    assert "no data and no particles" in h                     # land/pole coverage, stated on screen


def test_disclaimer_is_a_visible_dom_element_with_the_provenance_clause():
    html = fg.flow_globe_html(oc.flow_field_from_ocean(_fixture_snapshot()))
    m = re.search(r'<div class="disclaimer"[^>]*>(.*?)</div>', html, re.S)
    assert m, 'no visible <div class="disclaimer"> in the artifact'
    body = m.group(1)
    assert "OSCAR" in body and "NOT computed by planet-sim" in body
    assert "illustrative" in body


def test_banked_artifact_carries_the_provenance_clause():
    # The committed deliverable itself (docs/figures/planet-ocean-currents.html) is guarded: a future
    # regeneration cannot silently drop the provenance clause or the inlined three.js licence.
    text = BANKED.read_text(encoding="utf-8")
    assert 'class="disclaimer"' in text
    assert "OSCAR" in text and "NOT computed by planet-sim" in text
    assert "Three.js Authors" in text and "SPDX-License-Identifier: MIT" in text
    assert '"is_global":true' in text                          # the payload really is the global field


# --------------------------------------------------------------------------- #
# 6. The pacing knob — additive, default-off (the eddy artifact is byte-unchanged)
# --------------------------------------------------------------------------- #
def test_crossing_seconds_scales_the_pace_and_defaults_to_the_band_tuning():
    field = oc.flow_field_from_ocean(_fixture_snapshot())
    d_default = fg._build_data(field, 100, 0.03, 0.9, 0.5)
    d_band = fg._build_data(field, 100, 0.03, 0.9, 0.5, crossing_seconds=fg._BAND_CROSSING_SECONDS)
    assert d_default["accel"] == d_band["accel"]               # omitting the knob = the pre-O2 behaviour
    d_slow = fg._build_data(field, 100, 0.03, 0.9, 0.5, crossing_seconds=45.0)
    ratio = d_band["accel"] / d_slow["accel"]
    assert np.isclose(ratio, 45.0 / fg._BAND_CROSSING_SECONDS)   # pace scales exactly with the knob


# --------------------------------------------------------------------------- #
# 7. The loader — dim order, fill value, stride, the geostrophic knob (synthetic granule)
# --------------------------------------------------------------------------- #
def _write_synthetic_granule(path):
    """A tiny netCDF with OSCAR's exact layout: dims (time, longitude, latitude), 0–360 lon, _FillValue."""
    import h5netcdf

    lat = np.array([-60.0, 0.0, 60.0])
    lon = np.array([0.0, 90.0, 180.0, 270.0])
    u = np.arange(12, dtype=float).reshape(1, 4, 3)            # (time, lon, lat) — u[0, j, i]
    u[0, 1, 2] = -999.0                                        # one filled (no-data) cell
    with h5netcdf.File(str(path), "w") as f:
        f.dimensions = {"time": 1, "longitude": 4, "latitude": 3}
        f.create_variable("lat", ("latitude",), data=lat)
        f.create_variable("lon", ("longitude",), data=lon)
        for name, scale in (("u", 1.0), ("v", 2.0), ("ug", 10.0), ("vg", 20.0)):
            var = f.create_variable(name, ("time", "longitude", "latitude"), fillvalue=-999.0, data=u * scale)
            var.attrs["units"] = "m s-1"
        f.attrs["time_coverage_start"] = "2020-06-15T00:00:00"
        f.attrs["product_version"] = "v9.9"


# --------------------------------------------------------------------------- #
# 8. The seasonal series producer (§9.6 O4) — frames, the all-frames mask, honesty
# --------------------------------------------------------------------------- #
def _fixture_series(nt=12):
    """`nt` snapshots grounded in the real subsample: real OSCAR values, seasonally modulated. Raw
    conventions throughout (0–360 lon, NaN land) so the producer's rewrap/mask/stack runs on real shapes."""
    base = _fixture_snapshot()
    snaps = []
    for m in range(nt):
        s = float(np.cos(2 * np.pi * m / nt))
        snaps.append(oc.OceanSnapshot(lat=base.lat, lon=base.lon,
                                      u=base.u * (0.6 + 0.4 * s), v=base.v * (0.6 + 0.4 * s),
                                      product=base.product, doi=base.doi, credit=base.credit,
                                      date=f"2020-{m + 1:02d}-15", depth_note=base.depth_note))
    return snaps


def test_series_stacks_frames_and_defaults_to_month_labels():
    field = oc.flow_field_from_ocean_series(_fixture_series(12), period="monthly means for 2020")
    assert field.frames is not None
    nt, ny, nx = field.frames.u.shape
    assert nt == 12 and field.frames.v.shape == (12, ny, nx)
    assert field.frames.labels == ("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
    # the parent field is the representative (frame-0) snapshot (CPU fallback / land-ocean base still work)
    assert np.array_equal(field.u, field.frames.u[0]) and np.array_equal(field.v, field.frames.v[0])
    assert field.coverage.is_global is True


def test_series_mask_is_finite_in_every_frame_and_applied_across_the_stack():
    # The static mask is finite-in-ALL-frames (conservative — a cell measured only part of the year is
    # left bare, never blinking), and it is APPLIED to the whole stack (no filled-zero land leaks).
    ny, nx = 5, 8
    lat = np.linspace(-80.0, 80.0, ny)
    lon = np.linspace(0.0, 315.0, nx)                            # raw 0–360
    snaps = []
    for m in range(4):
        u = np.full((ny, nx), 0.2 + 0.1 * m)
        v = np.full((ny, nx), -0.1)
        u[0, 0] = np.nan                                         # land in EVERY frame
        if m == 2:
            u[2, 3] = np.nan                                     # a single-frame gap → must mask out in all
        snaps.append(oc.OceanSnapshot(lat=lat, lon=lon, u=u, v=v))
    field = oc.flow_field_from_ocean_series(snaps)
    assert int((~field.mask).sum()) == 2                         # the always-land cell AND the one-frame gap
    inv = ~field.mask
    assert np.all(field.frames.u[:, inv] == 0.0) and np.all(field.frames.v[:, inv] == 0.0)
    assert np.isfinite(field.frames.u).all() and np.isfinite(field.frames.v).all()   # no NaN reaches a sampler


def test_series_honesty_is_a_time_series_and_stays_acquisition_agnostic():
    # The honesty clause names it a TIME series and carries the caller's phrase VERBATIM — the producer
    # never invents "climatology"; the demo owns what the frames actually are (advisor's decouple call).
    field = oc.flow_field_from_ocean_series(_fixture_series(12), period="monthly means for 2020")
    h = field.honesty
    assert "REAL data" in h and "OSCAR" in h and "NOT computed by planet-sim" in h
    assert "TIME series" in h and "monthly means for 2020" in h
    assert "climatology" not in h                                # not claimed unless the caller's data is one
    assert "no data and no particles" in h


def test_framed_field_serializes_with_a_seasonal_frames_layer(tmp_path):
    from planet import planet_spec as ps
    from planet.flow_serialize import FRAMES_LAYER, vector_spec_from_flow_field

    field = oc.flow_field_from_ocean_series(_fixture_series(4), period="four snapshots, 2020")
    spec = vector_spec_from_flow_field(field, provenance="planet.ocean_currents (OSCAR series)")
    fl = spec.view().layer(FRAMES_LAYER)
    assert fl.data.ndim == 4 and fl.data.shape[0] == 4 and fl.data.shape[1] == 2
    assert fl.style["labels"] == list(field.frames.labels)
    ps.save(spec, tmp_path / "series")
    assert ps.load(tmp_path / "series") == spec                  # the R1 round-trip, now on a framed field


def test_seasonal_banked_artifact_carries_the_provenance_and_frames_when_present():
    # The O4 banked artifact is a token hand-off (12 granules); guard it IF it exists — a regeneration
    # can't silently drop the provenance clause, the frames payload, or the inlined three.js licence.
    if not BANKED_SEASONAL.exists():
        pytest.skip("seasonal artifact not banked yet (EARTHDATA_TOKEN hand-off — code + tests prove it)")
    text = BANKED_SEASONAL.read_text(encoding="utf-8")
    assert 'class="disclaimer"' in text
    assert "OSCAR" in text and "NOT computed by planet-sim" in text
    assert "Three.js Authors" in text and "SPDX-License-Identifier: MIT" in text
    assert '"frames":' in text and '"seconds_per_year":' in text   # the time axis really shipped
    assert 'id="timebadge"' in text                                # the month badge is in the DOM
    assert '"is_global":true' in text


def test_loader_transposes_fills_and_reads_attrs(tmp_path):
    h5netcdf = pytest.importorskip("h5netcdf")  # noqa: F841 — the [ocean] extra; skipped on a bare core install
    nc = tmp_path / "synthetic_oscar.nc"
    _write_synthetic_granule(nc)
    s = oc.load_oscar(nc)
    assert s.u.shape == (3, 4)                                 # (lat, lon) — the transpose happened
    assert s.u[0, 2] == 6.0                                    # u[0, lon=180 (j=2), lat=-60 (i=0)] = 2*3+0
    assert np.isnan(s.u[2, 1])                                 # the -999 fill cell (lat=60, lon=90) → NaN
    assert s.v[0, 2] == 12.0                                   # v = 2·u — the right variable was read
    assert s.date == "2020-06-15" and "v9.9" in s.product
    g = oc.load_oscar(nc, geostrophic_only=True)
    assert g.u[0, 2] == 60.0 and "geostrophic-only" in g.product
    strided = oc.load_oscar(nc, stride=2)
    assert strided.u.shape == (2, 2) and strided.lon.tolist() == [0.0, 180.0]
