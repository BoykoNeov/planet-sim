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


# --------------------------------------------------------------------------- #
# 4. build_spec — the knobs → saveable-spec "design a world" entry point (§9.1 / ADR 0004 #4)
# --------------------------------------------------------------------------- #
def test_build_spec_defaults_recover_the_present_day_earth_and_round_trip(tmp_path):
    # build_spec() with default knobs is the Sun/Earth-size/Earth-tilt world: its params equal the live
    # map's climate_params, and it round-trips identically — the property the design bench rests on.
    spec = ps.build_spec(**COARSE)
    assert spec.to_params() == pm.climate_params(n_cells=COARSE["n_cells"])
    ps.save(spec, tmp_path / "earth")
    assert ps.load(tmp_path / "earth") == spec


def test_build_spec_captures_the_composed_knobs_not_the_base(tmp_path):
    # The trap guard (§9.1): a naive from_view(climate_view(...), EBMParams(S0, A, D)) would store the
    # UN-perturbed knobs for any non-default star/size/tilt (climate_view discards the real params), so
    # the exported world would re-run to the WRONG climate. build_spec must serialize the COMPOSED params
    # (ai←spectrum, D←size, s2←tilt) — and they must survive the round-trip.
    knobs = dict(T_star=3000.0, size=1.8, obliquity_deg=40.0)
    spec = ps.build_spec(**knobs, **COARSE)
    p = spec.to_params()
    assert p == pm.climate_params(**knobs, n_cells=COARSE["n_cells"])     # the composed params, not the base
    assert p != EBMParams(n_cells=COARSE["n_cells"])                      # genuinely perturbed off Earth
    ps.save(spec, tmp_path / "exo")
    assert ps.load(tmp_path / "exo") == spec


def test_build_spec_equals_the_explicit_from_view_path():
    # build_spec is exactly the one-call form of "relax the view + capture the params that made it":
    # the present-day spec equals the hand-composed from_view(climate_view(), climate_params()) — a guard
    # that build_spec cannot silently diverge from the live map it mirrors.
    spec = ps.build_spec(**COARSE)
    params = pm.climate_params(n_cells=COARSE["n_cells"])
    assert spec == ps.from_view(pm.climate_view(**COARSE), params)


# --------------------------------------------------------------------------- #
# 5. diff — the two-world "what changed from a to b" compare path (the design-bench follow-on)
# --------------------------------------------------------------------------- #
EXO = dict(T_star=3000.0, size=1.8, obliquity_deg=40.0)     # a world clearly off Earth (star/size/tilt)


def test_diff_self_is_empty_and_is_the_inverse_of_equality():
    # The reflexive property that mirrors the round-trip identity: a world diffed against itself reports
    # NO change on any axis, and bool(diff) is exactly the inverse of __eq__ (empty iff equal).
    spec = ps.build_spec(**COARSE)
    d = ps.diff(spec, spec)
    assert not d.knobs and not d.fields and not d.other_changed
    assert not d.only_in_a and not d.only_in_b and d.grids_compatible
    assert bool(d) is False and (spec == spec)                # empty diff  <=>  equal worlds


def test_diff_earth_vs_exoplanet_flags_the_composed_knobs_and_the_climate():
    # The headline case: design an exoplanet off Earth and ask "what changed?". The KNOB diff is exactly
    # the §9.1 composition (ai←spectrum, D←size, s2←tilt) — never the raw star/size/tilt levers; the
    # CLIMATE diff is the three scalar fields that moved (temperature/precip numeric, biome categorical),
    # the ice line shows in other_changed (an annotation, not cell-differenced), and bool == a != b.
    a = ps.build_spec(**COARSE)
    b = ps.build_spec(**EXO, **COARSE)
    d = ps.diff(a, b)

    assert set(d.knobs) == {"ai", "D", "s2"}                  # the composed knobs, not {T_star, size, obliquity}
    assert all(d.knobs[k][0] != d.knobs[k][1] for k in d.knobs)
    assert d.knob_units["D"] == "W/m^2/K"                     # the diff is self-describing (carries units)

    assert set(d.fields) == {"temperature", "precipitation", "biome"}   # elevation flat in both → unchanged
    temp = d.fields["temperature"]
    assert not temp.categorical and temp.max_abs_delta > 0 and temp.delta.shape == a.grid.lat.shape + (a.grid.lon.size,)
    biome = d.fields["biome"]
    assert biome.categorical and 0.0 < biome.changed_fraction < 1.0 and np.isnan(biome.max_abs_delta)

    assert d.other_changed == ("ice_line",)                  # the ice line moved (annotation, flagged not differenced)
    assert not d.only_in_a and not d.only_in_b and d.grids_compatible
    assert bool(d) is True and (a != b)


def test_diff_field_only_change_is_not_gated_on_a_knob_change():
    # The geography-seam case (and the guard that field comparison is NOT gated on knob differences):
    # import a heightmap into one world, leave every knob identical → an empty knob diff but a real field
    # diff on exactly the elevation layer. A diff that keyed fields off knob changes would miss this.
    a = ps.build_spec(**COARSE)
    heightmap = np.random.default_rng(1).random(a.view().layer("elevation").data.shape) * 1000.0
    b = ps.build_spec(elevation=heightmap, **COARSE)
    d = ps.diff(a, b)
    assert d.knobs == {}                                      # identical knobs — same climate
    assert set(d.fields) == {"elevation"} and d.fields["elevation"].changed_fraction > 0
    assert not d.other_changed and d.grids_compatible
    assert bool(d) is True and (a != b)


def test_diff_an_added_layer_shows_in_only_in_b():
    # The interchange edge: a world carrying an extra (future) layer the other lacks shows up in
    # only_in_b — and only_in_a under the reverse diff. Grids stay compatible (matching-shape layer).
    base = ps.build_spec(**COARSE)
    shape = base.view().layer("temperature").data.shape
    extra = pm.Layer("future_salinity", pm.LayerKind.SCALAR_FIELD, np.full(shape, 35.0), "psu", z_order=5)
    withlayer = dataclasses.replace(base, layers=base.layers + (extra,))

    d = ps.diff(base, withlayer)
    assert d.only_in_b == ("future_salinity",) and d.only_in_a == ()
    assert not d.fields and d.grids_compatible and bool(d) is True and (base != withlayer)
    assert ps.diff(withlayer, base).only_in_a == ("future_salinity",)    # reversed


def test_diff_is_symmetric_the_deltas_negate():
    # diff(a, b) and diff(b, a) describe the same difference from opposite ends: the same knobs changed
    # (values swapped), and the numeric field delta b−a is exactly the negation of a−b.
    a = ps.build_spec(**COARSE)
    b = ps.build_spec(**EXO, **COARSE)
    ab, ba = ps.diff(a, b), ps.diff(b, a)
    assert set(ab.knobs) == set(ba.knobs)
    assert all(ab.knobs[k] == (ba.knobs[k][1], ba.knobs[k][0]) for k in ab.knobs)
    assert np.array_equal(ab.fields["temperature"].delta, -ba.fields["temperature"].delta)
    assert np.array_equal(ab.fields["biome"].delta, ba.fields["biome"].delta)   # a≠b mask is symmetric


def test_diff_incompatible_grids_skips_the_fields_but_keeps_the_knobs():
    # Two worlds at different resolution (n_cells) have non-aligned grids: the per-cell field deltas are
    # meaningless and must be SKIPPED (grids_compatible False, empty fields) — but the knob diff (here
    # n_cells itself) still stands, so the comparison is not silently empty (the cell must surface this).
    a = ps.build_spec(**COARSE)
    b = ps.build_spec(n_cells=COARSE["n_cells"] + 20, n_tau=COARSE["n_tau"])
    d = ps.diff(a, b)
    assert d.grids_compatible is False and d.fields == {}
    assert "n_cells" in d.knobs and d.knobs["n_cells"] == (COARSE["n_cells"], COARSE["n_cells"] + 20)
    assert bool(d) is True
