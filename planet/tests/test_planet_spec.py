"""Tests for the planet-spec interchange schema (Planet §9.3 / ADR 0004 #3–4).

This is the one place in the deep-end map with a **real** correctness property. The map's renderer is
a reach layer (smoke-tested only), but **round-trip identity** ``load(save(spec)) == spec`` — array
identity within the schema — is a genuine invariant (ADR 0004 #4), so it gets a genuine, *always-green*
test (the schema imports neither Plotly nor ipywidgets — it runs on a bare core install).

The crown-jewel test (``test_round_trip_*``) is reinforced by a **negative control**
(``test_eq_detects_a_changed_array``) so the round-trip pass cannot be a trivially-true ``__eq__``,
and by a **mixed-dtype / NaN / inert / future-layer** case so the identity is exercised on exactly the
fields the schema must survive: int biome codes beside float fields, the ice-line NaN separators, the
inert geography seam, and an unknown future layer the versioned schema must carry forward.
"""
import dataclasses
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from planet import planet_spec as ps
from planet import planetmap as pm
from planet.albedo import EBMParams

REPO_ROOT = Path(__file__).resolve().parents[2]

COARSE = dict(n_cells=40, n_tau=0.25)


def _synthetic_spec() -> ps.PlanetSpec:
    """A small spec exercising every field the schema must survive — fast, no EBM solve.

    Float temperature + **int** biome codes + a NaN-gapped ice-line annotation + an **inert** geography
    layer + an **unknown future** layer (salinity) the versioned schema must carry forward.
    """
    lat = np.linspace(-90.0, 90.0, 9)
    lon = np.linspace(-180.0, 180.0, 5)
    grid = pm.Grid(lat, lon)
    layers = (
        pm.Layer("temperature", pm.LayerKind.SCALAR_FIELD, np.random.default_rng(0).random((9, 5)),
                 "°C", {"colorscale": "RdBu_r"}, z_order=0),
        pm.Layer("biome", pm.LayerKind.SCALAR_FIELD, np.arange(45).reshape(9, 5).astype(int),
                 "Whittaker biome code", {"categorical": True}, z_order=1),
        pm.Layer("ice_line", pm.LayerKind.ANNOTATION,
                 np.array([[60.0, 0.0], [np.nan, np.nan], [-60.0, 0.0]]), "degrees", {}, z_order=3),
        pm.Layer("elevation", pm.LayerKind.SCALAR_FIELD, np.zeros((9, 5)), "m", {}, z_order=-1, inert=True),
        pm.Layer("future_salinity", pm.LayerKind.SCALAR_FIELD, np.ones((9, 5)) * 35.0, "psu", {}, z_order=5),
    )
    knobs = dataclasses.asdict(EBMParams(n_cells=40))
    return ps.PlanetSpec(ps.SCHEMA_VERSION, grid, layers, knobs, dict(ps.KNOB_UNITS))


# --------------------------------------------------------------------------- #
# 0. The layering guard — the interchange schema is bare-core (no render deps)
# --------------------------------------------------------------------------- #
def test_importing_planet_spec_stays_headless():
    # The interchange schema must run on a bare core install (the round-trip test is always-green).
    # Checked in a clean subprocess (session-robust — see test_planetmap's guard for why).
    code = ("import sys, planet.planet_spec\n"
            "print(','.join(m for m in ('plotly', 'ipywidgets', 'matplotlib') if m in sys.modules))\n")
    out = subprocess.run([sys.executable, "-c", code], cwd=str(REPO_ROOT),
                         capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip() == "", f"a heavy dep was pulled at import time: {out.stdout.strip()}"


# --------------------------------------------------------------------------- #
# 1. The round-trip identity — the ADR 0004 #4 crown jewel
# --------------------------------------------------------------------------- #
def test_round_trip_identity_of_a_real_present_day_view(tmp_path):
    view = pm.climate_view(**COARSE)
    spec = ps.from_view(view, EBMParams(n_cells=COARSE["n_cells"]))
    ps.save(spec, tmp_path / "world")
    assert ps.load(tmp_path / "world") == spec


def test_round_trip_preserves_dtypes_nan_inert_and_future_layers(tmp_path):
    spec = _synthetic_spec()
    ps.save(spec, tmp_path / "w")
    loaded = ps.load(tmp_path / "w")
    assert loaded == spec
    by_name = {ly.name: ly for ly in loaded.layers}
    assert np.issubdtype(by_name["biome"].data.dtype, np.integer)     # int codes stay int
    assert np.issubdtype(by_name["temperature"].data.dtype, np.floating)
    assert by_name["elevation"].inert is True                         # the geography seam flag survives
    assert np.array_equal(by_name["ice_line"].data, spec.layers[2].data, equal_nan=True)  # NaN gap survives
    assert "future_salinity" in by_name                               # the versioned schema carries it forward


def test_save_writes_the_lean_json_plus_npz_encoding(tmp_path):
    # ADR 0004 #3: the v1 lean encoding is a JSON manifest + a .npz of the arrays, sharing a stem
    # (any passed suffix is dropped).
    json_path, npz_path = ps.save(_synthetic_spec(), tmp_path / "planet.ignored")
    assert json_path.name == "planet.json" and npz_path.name == "planet.npz"
    assert json_path.exists() and npz_path.exists()


def test_save_rejects_duplicate_layer_names(tmp_path):
    # Layer names key the .npz arrays (layer__<name>); duplicates would silently overwrite. Guard it.
    grid = pm.Grid(np.array([0.0]), np.array([0.0]))
    dup = (pm.Layer("x", pm.LayerKind.SCALAR_FIELD, np.zeros((1, 1)), "", {}, 0),
           pm.Layer("x", pm.LayerKind.SCALAR_FIELD, np.ones((1, 1)), "", {}, 1))
    spec = ps.PlanetSpec(ps.SCHEMA_VERSION, grid, dup, {}, {})
    with pytest.raises(ValueError, match="unique"):
        ps.save(spec, tmp_path / "dup")


def test_eq_detects_a_changed_array(tmp_path):
    # Negative control: the round-trip pass must not be a trivially-true __eq__. Perturb one array
    # after load and equality must fail (array-aware, not identity/metadata-only).
    spec = _synthetic_spec()
    ps.save(spec, tmp_path / "w")
    loaded = ps.load(tmp_path / "w")
    loaded.layers[0].data[0, 0] += 1.0
    assert loaded != spec
    # a changed knob is also caught
    other = ps.load(tmp_path / "w")
    other.knobs["S0"] += 1.0
    assert other != spec


# --------------------------------------------------------------------------- #
# 2. The schema is self-describing — explicit units, versioned (ADR 0004 #3)
# --------------------------------------------------------------------------- #
def test_manifest_is_self_describing(tmp_path):
    import json
    json_path, _ = ps.save(_synthetic_spec(), tmp_path / "w")
    manifest = json.loads(json_path.read_text(encoding="utf-8"))
    assert manifest["schema_version"] == ps.SCHEMA_VERSION
    assert manifest["grid"]["lat_units"] and manifest["grid"]["lon_units"]
    assert manifest["knob_units"]["S0"] == "W/m^2"
    # every layer carries its kind + units + the npz key that holds its array
    for m in manifest["layers"]:
        assert m["kind"] in {k.value for k in pm.LayerKind}
        assert "units" in m and m["npz_key"].startswith("layer__")


def test_knob_units_cover_every_ebmparams_field():
    # A guard: if a future EBMParams field is added, it must get a unit (nothing ambiguous, §7).
    assert set(ps.KNOB_UNITS) == set(dataclasses.asdict(EBMParams()))


# --------------------------------------------------------------------------- #
# 3. The world is re-runnable — knobs reconstruct EBMParams; the view re-renders
# --------------------------------------------------------------------------- #
def test_to_params_reconstructs_the_ebmparams():
    params = EBMParams(S0=1300.0, A=205.0, D=0.6, n_cells=40)
    view = pm.climate_view(S0=1300.0, A=205.0, D=0.6, **COARSE)
    spec = ps.from_view(view, params)
    assert spec.to_params() == params


def test_view_round_trips_to_a_renderable_planetview(tmp_path):
    spec = _synthetic_spec()
    ps.save(spec, tmp_path / "w")
    view = ps.load(tmp_path / "w").view()
    assert isinstance(view, pm.PlanetView)
    assert {ly.name for ly in view.layers} == {ly.name for ly in spec.layers}
