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

from planet.albedo import EBMParams
from planet.planetmap import Grid, Layer, LayerKind, PlanetView

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
