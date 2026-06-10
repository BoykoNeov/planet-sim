"""The planet-spec interchange schema — export/import a planet's state (Planet §9.3, ADR 0004 #3–4).

The deep-end map (:mod:`planet.planetmap`) lets you author a climate; this module lets you
take that world **out** (share, inspect) and bring an externally-authored one **in** (a future
geography-editing app paints elevation/coastlines, the model imports them). ADR 0004 #3 fixes the
design call: **pin a *schema*, not a file format**. The portable artifact is a documented, versioned,
self-describing **planet-spec**:

* the **grid geometry** (lat/lon, in degrees);
* **explicit units** — each layer and each knob carries its unit string, so an external consumer
  cannot misread it (this program is unit-obsessive by §7 discipline);
* the **layer list** — *the* :class:`~planet.planetmap.PlanetView` registry **is** the export
  manifest (ADR 0004 #3: one structure serialized, not a second one invented — this module imports
  :class:`~planet.planetmap.Layer`/:class:`~planet.planetmap.Grid`, it does not
  reinvent them);
* the **knob values** (the :class:`~planet.albedo.EBMParams` that reproduce the climate); and
* a **``schema_version``** for forward/backward compatibility.

The **encoding** is chosen per consumer, *behind* the schema. v1 ships the **lean default — a JSON
manifest + a NumPy ``.npz``** of the arrays (the JSON names which ``.npz`` key holds each array). Two
further encodings are *documented future options behind the same schema*, deliberately not built in v1
(build the seam, not the machinery):

* **editable-geography heightmaps** (elevation / bathymetry / land-ocean masks) interchange as
  **16-bit grayscale PNG** — the native currency of paint tools and web canvases, the round-trip a
  future map-editor needs (8-bit is too coarse for elevation);
* **NetCDF** for climate-tool interop — domain-standard but browser-hostile, so *not* the v1 choice
  (the editor is the consumer; resolving the NetCDF-vs-browser tension by *consumer* is the whole
  point of pinning the schema rather than a format).

The one *real* correctness property here (ADR 0004 #4)
------------------------------------------------------
Unlike the map (whose only test is an execution smoke-test), **round-trip identity**
``load(save(spec)) == spec`` — array identity within the schema — is a genuine invariant and gets a
genuine test (``tests/test_planet_spec.py``). :meth:`PlanetSpec.__eq__` compares the metadata exactly
and every array with :func:`numpy.array_equal` (with ``equal_nan`` for the float fields — the ice-line
annotation carries NaN separators between its two latitude circles). This module imports neither Plotly
nor ipywidgets, so the round-trip test is **always-green** on a bare core install.

The geography seam is inert at v1 (the honesty flag)
----------------------------------------------------
A geography layer (elevation/bathymetry/mask) round-trips like any other — but **importing it and
running a model that *responds* to it is a separate, staircase-gated capability**. Until the consuming
physics exists (a lapse-rate diagnostic at the cheap tier; true 2-D orographic precip at rung 5), an
imported geometry is **inert**: carried, displayed, round-tripped, but not changing the climate. That
is the :attr:`~planet.planetmap.Layer.inert` flag — the round-trip guarantee is *array
identity*, not a changed climate. Named, not blurred.
"""
from __future__ import annotations

import dataclasses
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from planet import demo_biomes
from planet.albedo import EBMParams
from planet.ebm import A_OLR, D_TRANSPORT, S0_EARTH
from planet.exoplanet import T_SUN
from planet.obliquity import OBLIQUITY_EARTH
from planet.planetmap import (Grid, Layer, LayerKind, LIVE_N_TAU, N_LON, PlanetView,
                              build_view, climate_params)

SCHEMA_VERSION = 1

# Self-describing units for the knob bundle (EBMParams fields). The physical ones carry SI/climlab
# units; the dimensionless / structural ones ("" / a tag) are explicit too, so nothing is ambiguous.
KNOB_UNITS: dict[str, str] = {
    "S0": "W/m^2", "s2": "dimensionless",
    "A": "W/m^2", "B": "W/m^2/K", "D": "W/m^2/K",
    "T_freeze": "degC",
    "a0": "dimensionless", "a2": "dimensionless", "ai": "dimensionless",
    "water_depth": "m", "n_cells": "count", "face": "enum",
}


def _arrays_equal(a: np.ndarray, b: np.ndarray) -> bool:
    """Array equality for the round-trip invariant — ``equal_nan`` for floats (annotation NaN gaps).

    ``np.array_equal(equal_nan=True)`` is only valid for floating dtypes (``isnan`` is undefined on
    ints), so the NaN-tolerant path is taken only for float arrays; integer layers (the biome codes)
    compare exactly. Shapes and dtypes are preserved by the ``.npz`` round-trip, so an equality here
    is a genuine identity, not a lenient match.
    """
    a = np.asarray(a)
    b = np.asarray(b)
    if a.shape != b.shape:
        return False
    if np.issubdtype(a.dtype, np.floating) and np.issubdtype(b.dtype, np.floating):
        return np.array_equal(a, b, equal_nan=True)
    return np.array_equal(a, b)


def _layers_equal(x: Layer, y: Layer) -> bool:
    """Two layers are identical iff their metadata matches exactly and their data arrays are equal."""
    return (x.name == y.name and LayerKind(x.kind) == LayerKind(y.kind) and x.units == y.units
            and x.style == y.style and x.z_order == y.z_order and bool(x.inert) == bool(y.inert)
            and _arrays_equal(x.data, y.data))


@dataclass(eq=False)
class PlanetSpec:
    """A versioned, self-describing snapshot of a planet — the portable interchange artifact (ADR 0004 #3).

    Bundles the :class:`~planet.planetmap.PlanetView` registry (``grid`` + ``layers`` — the
    export manifest itself) with the ``knobs`` that reproduce the climate (the
    :class:`~planet.albedo.EBMParams` fields), their ``knob_units``, and a ``schema_version``.
    ``eq=False`` because equality must be array-aware: :meth:`__eq__` compares the metadata exactly and
    every array with :func:`numpy.array_equal` — the ADR 0004 #4 round-trip invariant.
    """

    schema_version: int
    grid: Grid
    layers: tuple[Layer, ...]
    knobs: dict
    knob_units: dict

    def __eq__(self, other) -> bool:
        if not isinstance(other, PlanetSpec):
            return NotImplemented
        if (self.schema_version != other.schema_version
                or self.knobs != other.knobs
                or self.knob_units != other.knob_units
                or self.grid.lat_units != other.grid.lat_units
                or self.grid.lon_units != other.grid.lon_units
                or not _arrays_equal(self.grid.lat, other.grid.lat)
                or not _arrays_equal(self.grid.lon, other.grid.lon)
                or len(self.layers) != len(other.layers)):
            return False
        return all(_layers_equal(x, y) for x, y in zip(self.layers, other.layers))

    def view(self) -> PlanetView:
        """The :class:`~planet.planetmap.PlanetView` (grid + layers) for the renderer."""
        return PlanetView(grid=self.grid, layers=tuple(self.layers))

    def to_params(self) -> EBMParams:
        """Reconstruct the :class:`~planet.albedo.EBMParams` from the knobs — the world made re-runnable.

        The ``knobs`` dict is exactly ``dataclasses.asdict(EBMParams)``, so this round-trips the
        forcing/machinery that produced the climate. (At v1 an imported *geography* layer is still
        inert — re-running reproduces the zonal-mean climate, not a geography-responsive one; §9.3.)
        """
        return EBMParams(**self.knobs)


def from_view(view: PlanetView, params: EBMParams, schema_version: int = SCHEMA_VERSION) -> PlanetSpec:
    """Build a :class:`PlanetSpec` from a rendered view + the knobs that produced it (the export side).

    The view's registry becomes the spec's manifest verbatim (ADR 0004 #3); the params become the
    self-describing ``knobs`` (with :data:`KNOB_UNITS`).
    """
    knobs = dataclasses.asdict(params)
    return PlanetSpec(schema_version=schema_version, grid=view.grid, layers=tuple(view.layers),
                      knobs=knobs, knob_units=dict(KNOB_UNITS))


def build_spec(S0: float = S0_EARTH, A: float = A_OLR, D: float = D_TRANSPORT, *,
               T_star: float = T_SUN, size: float = 1.0, obliquity_deg: float = OBLIQUITY_EARTH,
               n_cells: int = 180, n_tau: float = LIVE_N_TAU, n_lon: int = N_LON,
               elevation: np.ndarray | None = None) -> PlanetSpec:
    """Knobs → a fresh equilibrium climate → a saveable :class:`PlanetSpec` (the *design-a-world* entry point).

    The one-call path a "design your own planet" workflow needs: compose the knobs into the climate
    params (:func:`planet.planetmap.climate_params` — the single source of truth), relax to the
    biome-map :class:`~planet.planetmap.PlanetView` (:func:`planet.demo_biomes.compute` +
    :func:`planet.planetmap.build_view`), and bundle **that view with the params that produced it**
    into a spec ready for :func:`save`.

    Capturing the *composed* params is the whole point: a naive
    ``from_view(climate_view(...), EBMParams(S0, A, D))`` would store the **un-perturbed** knobs for any
    non-default star/size/tilt (``climate_view`` builds and discards the real params), so the exported
    world would be re-runnable to the *wrong* climate. Here :func:`climate_params` builds the params
    once and they are both solved *and* serialized, so :meth:`PlanetSpec.to_params` reconstructs the
    exact world (the §9.1 / ADR 0004 #4 re-runnability). The Sun/Earth-size/Earth-tilt defaults produce
    the present-day Earth spec, identical to ``from_view(climate_view(), EBMParams())``.
    """
    params = climate_params(S0, A, D, T_star=T_star, size=size, obliquity_deg=obliquity_deg, n_cells=n_cells)
    result = demo_biomes.compute(params, n_tau=n_tau)
    view = build_view(result, n_lon=n_lon, elevation=elevation)
    return from_view(view, params)


# --------------------------------------------------------------------------- #
# Two-world diff — "what changed between world a and world b" (the design-bench *compare* path)
# --------------------------------------------------------------------------- #
@dataclass(frozen=True, eq=False)
class FieldDelta:
    """How one scalar-field layer differs between two worlds — the render-ready per-cell delta.

    ``delta`` is the full-globe lat×lon difference, read according to ``categorical``: for a **numeric**
    field (temperature, precipitation, elevation) it is the signed ``b − a`` (a Δ-map paints straight
    from it); for a **categorical** field (the biome codes — code arithmetic is meaningless) it is instead
    the boolean ``a ≠ b`` *changed mask*. ``changed_fraction`` is the share of cells that differ — the
    honest metric for both, and *the* metric for the categorical field. ``max_abs_delta`` / ``mean_delta``
    summarise the numeric delta (both ``nan`` for a categorical field); ``mean_delta`` is an honest
    **area-weighted** mean because the grid is uniform in ``sin φ`` (equal-area latitude bands) and in
    longitude. ``eq=False`` for the same array reason as :class:`~planet.planetmap.Layer`.
    """

    name: str
    units: str
    categorical: bool
    delta: np.ndarray
    changed_fraction: float
    max_abs_delta: float
    mean_delta: float


@dataclass(frozen=True, eq=False)
class SpecDiff:
    """The structured difference between two planet-specs — world ``a`` → world ``b`` (the *compare* sibling of :meth:`PlanetSpec.__eq__`).

    Where ``__eq__`` answers *are these the same world*, :func:`diff` answers *how do they differ*, along
    the axes a "design two worlds and compare them" workflow cares about:

    * ``knobs`` — the design inputs that changed: ``name → (a_value, b_value)`` for each differing knob
      (an exoplanet off Earth differs in exactly ``ai``/``D``/``s2`` — the §9.1 composition);
      ``knob_units`` carries their unit strings so a diff renders without the original specs.
    * ``fields`` — the climate that changed: a :class:`FieldDelta` per **scalar-field** layer present in
      *both* worlds whose data differs (temperature, precipitation, biome, the inert elevation seam).
    * ``only_in_a`` / ``only_in_b`` — layer names present in one world but not the other (an imported
      future layer, a Phase-4 circulation overlay): the interchange edge.
    * ``other_changed`` — layers present in both whose data differs but that were **not** cell-differenced:
      non-scalar layers (the ice-line annotation, a vector overlay), *and* — when ``grids_compatible`` is
      ``False`` — any scalar layer too (a mismatched grid makes the per-cell ``b − a`` meaningless).
    * ``grids_compatible`` — ``False`` when the two grids do not match (e.g. different ``n_cells``); the
      per-cell ``fields`` deltas are then **skipped** (cell-to-cell subtraction is not meaningful) while
      the ``knobs`` diff still stands. A consumer must surface this rather than read an empty ``fields``
      as "no climate change".

    The diff is **data-scoped** — it compares knob values and layer *data*, not layer metadata
    (units/style/z_order/inert). So for specs that share the v1 layer metadata (everything
    :func:`build_spec` / :func:`from_view` produces), ``bool(diff(a, b))`` is exactly ``a != b``: the diff
    is empty iff the worlds are equal — the inverse of the round-trip invariant. ``eq=False`` for the same
    array reason as :class:`~planet.planetmap.Layer`.
    """

    knobs: dict
    knob_units: dict
    fields: dict
    only_in_a: tuple
    only_in_b: tuple
    other_changed: tuple
    grids_compatible: bool

    def __bool__(self) -> bool:
        """``True`` iff the two worlds differ on any tracked axis (knobs, fields, layer set, grid)."""
        return bool(self.knobs or self.fields or self.only_in_a or self.only_in_b
                    or self.other_changed or not self.grids_compatible)


def _grids_compatible(a: Grid, b: Grid) -> bool:
    """Two grids are field-comparable iff their coordinates *and* units match exactly (so a per-cell
    ``b − a`` aligns cell-for-cell)."""
    return (a.lat_units == b.lat_units and a.lon_units == b.lon_units
            and _arrays_equal(a.lat, b.lat) and _arrays_equal(a.lon, b.lon))


def _field_delta(a: Layer, b: Layer) -> FieldDelta:
    """The per-cell :class:`FieldDelta` for one scalar-field layer — numeric ``b − a``, or the categorical
    ``a ≠ b`` changed-mask for the biome codes (code subtraction is meaningless)."""
    categorical = bool(a.style.get("categorical") or b.style.get("categorical"))
    a_data, b_data = np.asarray(a.data), np.asarray(b.data)
    if categorical:
        changed = a_data != b_data
        return FieldDelta(a.name, a.units, True, changed, float(np.mean(changed)),
                          float("nan"), float("nan"))
    delta = b_data.astype(float) - a_data.astype(float)
    return FieldDelta(a.name, a.units, False, delta, float(np.mean(delta != 0)),
                      float(np.max(np.abs(delta))), float(delta.mean()))


def diff(a: PlanetSpec, b: PlanetSpec) -> SpecDiff:
    """Compare two worlds — *what changed from* ``a`` *to* ``b`` (the design-bench two-world diff).

    The :func:`diff` counterpart to :meth:`PlanetSpec.__eq__`: instead of a single same/different bit it
    returns a structured, render-ready :class:`SpecDiff` — the changed **knobs** (the design inputs), a
    per-cell :class:`FieldDelta` for each changed **scalar field** (the climate that moved), the layer-set
    edges (``only_in_a`` / ``only_in_b`` / ``other_changed``), and a ``grids_compatible`` flag. Pure and
    headless (no Plotly / ipywidgets), it reuses the schema's own :func:`_arrays_equal` (``equal_nan``) so
    "changed" here means exactly "not array-identical" there.

    The per-cell ``fields`` deltas need aligned grids; when the grids differ (e.g. two worlds at different
    ``n_cells``) ``grids_compatible`` is ``False`` and ``fields`` is left empty while the knob diff still
    stands — a consumer must surface that rather than read empty ``fields`` as "no climate change". The
    diff is **data-scoped** (knob values + layer *data*, not layer metadata), so for the v1-metadata specs
    :func:`build_spec` produces, ``bool(diff(a, b)) == (a != b)``.
    """
    knob_changes = {}
    for k in (*a.knobs, *(k for k in b.knobs if k not in a.knobs)):
        va, vb = a.knobs.get(k), b.knobs.get(k)
        if va != vb:
            knob_changes[k] = (va, vb)
    knob_units = {k: a.knob_units.get(k, b.knob_units.get(k)) for k in knob_changes}

    b_layers = {ly.name: ly for ly in b.layers}
    only_in_a = tuple(ly.name for ly in a.layers if ly.name not in b_layers)
    a_names = {ly.name for ly in a.layers}
    only_in_b = tuple(ly.name for ly in b.layers if ly.name not in a_names)

    grids_ok = _grids_compatible(a.grid, b.grid)
    fields: dict = {}
    other_changed = []
    for la in a.layers:                                   # layers in BOTH, in a's registration order
        lb = b_layers.get(la.name)
        if lb is None or _arrays_equal(la.data, lb.data):  # absent, or identical data → not a change
            continue
        if grids_ok and LayerKind(la.kind) is LayerKind.SCALAR_FIELD and LayerKind(lb.kind) is LayerKind.SCALAR_FIELD:
            fields[la.name] = _field_delta(la, lb)
        else:
            other_changed.append(la.name)
    return SpecDiff(knobs=knob_changes, knob_units=knob_units, fields=fields,
                    only_in_a=only_in_a, only_in_b=only_in_b,
                    other_changed=tuple(other_changed), grids_compatible=grids_ok)


def delta_view(spec_a: PlanetSpec, spec_b: PlanetSpec, active: str = "biome") -> PlanetView:
    """A single-layer :class:`~planet.planetmap.PlanetView` painting *what changed* in one field between
    two worlds — the **Δ-globe data** for :func:`planet.planetmap.render_comparison`.

    The headless data side of the deep-end two-world diff (it builds a :class:`~planet.planetmap.Layer`,
    pulls **no** Plotly — only the paint step needs the ``[webviz]`` extra). From :func:`diff` it takes the
    ``active`` field's :class:`FieldDelta` and wraps it as a renderable scalar layer on ``spec_a``'s grid:

    * a **continuous** field (temperature / precipitation / elevation) → the signed ``b − a`` with a
      *diverging* colorscale centred at zero (a ``cmid`` style hint the renderer honors); or
    * the **categorical** biome → the ``a ≠ b`` *changed-mask* (code arithmetic is meaningless), a 2-tone
      "unchanged / changed" layer answering *where did the biome flip*.

    Identical worlds give an all-zero Δ (a flat globe); the field must be a scalar layer present in both,
    and the two grids must match — a per-cell Δ across different ``n_cells`` is meaningless, so an
    incompatible grid raises (mirroring :attr:`SpecDiff.grids_compatible`). The style dicts are JSON-safe,
    so this module stays Plotly/ipywidgets-free.
    """
    d = diff(spec_a, spec_b)
    if not d.grids_compatible:
        raise ValueError("delta_view needs two worlds on the same grid — a per-cell Δ is meaningless across "
                         "different n_cells; build both at one resolution.")
    a_layer = spec_a.view().layer(active)              # validates `active`, and is the unchanged-case reference
    if LayerKind(a_layer.kind) is not LayerKind.SCALAR_FIELD:
        raise ValueError(f"delta_view differences a scalar field; {active!r} is "
                         f"{LayerKind(a_layer.kind).value}, not a SCALAR_FIELD")
    fd = d.fields.get(active)
    if a_layer.style.get("categorical"):
        mask = (np.asarray(fd.delta) if fd is not None
                else np.zeros_like(np.asarray(a_layer.data), dtype=bool)).astype(int)
        layer = Layer(f"{active} changed", LayerKind.SCALAR_FIELD, mask, "changed",
                      style={"categorical": True, "colorbar_title": "biome Δ",
                             "colors": {"0": "#e8e8e8", "1": "#6a3d9a"},     # purple ≠ the red ice line
                             "names": {"0": "same biome", "1": "biome changed"}}, z_order=0)
    else:
        delta = (np.asarray(fd.delta, dtype=float) if fd is not None
                 else np.zeros(np.asarray(a_layer.data).shape, dtype=float))
        layer = Layer(f"Δ {active}", LayerKind.SCALAR_FIELD, delta, a_layer.units,
                      style={"colorscale": "RdBu_r", "cmid": 0.0,
                             "colorbar_title": f"Δ {active} ({a_layer.units})"}, z_order=0)
    return PlanetView(grid=spec_a.grid, layers=(layer,))


def _stem_paths(path) -> tuple[Path, Path]:
    """The ``(manifest.json, arrays.npz)`` pair for a path stem (a suffix is stripped)."""
    stem = Path(path).with_suffix("")
    return stem.with_suffix(".json"), stem.with_suffix(".npz")


def save(spec: PlanetSpec, path) -> tuple[Path, Path]:
    """Write ``spec`` as the lean v1 encoding — a JSON manifest + a ``.npz`` of the arrays.

    ``path`` is a stem (any suffix is dropped): ``<stem>.json`` holds the self-describing manifest
    (schema version, grid units + array keys, knobs + units, and per-layer metadata naming each array's
    ``.npz`` key); ``<stem>.npz`` holds the grid and every layer's data array (dtypes preserved → the
    round-trip is a true identity). Returns the two written paths.
    """
    json_path, npz_path = _stem_paths(path)
    names = [ly.name for ly in spec.layers]
    if len(names) != len(set(names)):
        raise ValueError(f"layer names must be unique (they key the .npz arrays), got duplicates in {names}")
    arrays: dict[str, np.ndarray] = {
        "grid_lat": np.asarray(spec.grid.lat),
        "grid_lon": np.asarray(spec.grid.lon),
    }
    layer_meta = []
    for ly in spec.layers:
        key = f"layer__{ly.name}"
        arrays[key] = np.asarray(ly.data)
        layer_meta.append({
            "name": ly.name, "kind": LayerKind(ly.kind).value, "units": ly.units,
            "style": ly.style, "z_order": int(ly.z_order), "inert": bool(ly.inert),
            "npz_key": key, "dtype": str(arrays[key].dtype), "shape": list(arrays[key].shape),
        })
    manifest = {
        "schema_version": int(spec.schema_version),
        "grid": {"lat_units": spec.grid.lat_units, "lon_units": spec.grid.lon_units,
                 "lat_key": "grid_lat", "lon_key": "grid_lon"},
        "knobs": spec.knobs,
        "knob_units": spec.knob_units,
        "layers": layer_meta,
    }
    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    np.savez(npz_path, **arrays)
    return json_path, npz_path


def load(path) -> PlanetSpec:
    """Read a planet-spec back from its ``<stem>.json`` + ``<stem>.npz`` — the import side.

    Rebuilds the exact :class:`~planet.planetmap.Grid` / :class:`~planet.planetmap.Layer`
    objects (arrays from the ``.npz``, metadata from the JSON), so ``load(save(spec)) == spec``
    (the ADR 0004 #4 round-trip invariant). The arrays are copied out of the ``.npz`` handle before it
    closes.
    """
    json_path, npz_path = _stem_paths(path)
    manifest = json.loads(json_path.read_text(encoding="utf-8"))
    with np.load(npz_path) as npz:
        data = {k: np.asarray(npz[k]) for k in npz.files}

    g = manifest["grid"]
    grid = Grid(lat=data[g["lat_key"]], lon=data[g["lon_key"]],
                lat_units=g["lat_units"], lon_units=g["lon_units"])
    layers = tuple(
        Layer(name=m["name"], kind=LayerKind(m["kind"]), data=data[m["npz_key"]],
              units=m["units"], style=m["style"], z_order=int(m["z_order"]), inert=bool(m["inert"]))
        for m in manifest["layers"]
    )
    return PlanetSpec(schema_version=int(manifest["schema_version"]), grid=grid, layers=layers,
                      knobs=manifest["knobs"], knob_units=manifest["knob_units"])
