"""Serialize a vector flow field through the planet-spec schema — the viz/output seam (Planet R1, §9.3).

The §9.5 flow renderers consume a :class:`~planet.flow_globe.FlowField` — the renderer-agnostic
*vector-field-on-a-globe* contract (a lat×lon grid + per-cell ``(u, v)`` + a coverage extent + an
honesty label). R1 makes that field **serializable** so a *future* ocean engine (the spin-out,
plan §11) can hand a ``(u, v)`` + scalar field to the **same** interchange and renderer the emergent
eddy band uses today — proving the seam is **producer-agnostic**: the serialization and the renderer do
not care *what* generated the field.

The design call (ADR 0004 #3 — *one structure serialized, not a second one invented*)
----------------------------------------------------------------------------------------
A :class:`FlowField` is expressed as a ``VECTOR_OVERLAY`` :class:`~planet.planetmap.Layer` inside a
:class:`~planet.planetmap.PlanetView` / :class:`~planet.planet_spec.PlanetSpec` — **not** a new
file format. So :func:`planet.planet_spec.save` / :func:`~planet.planet_spec.load` round-trip it
**unchanged** (the ``(2, n_lat, n_lon)`` ``[u, v]`` rides the ``.npz`` like any layer's data), and
:func:`planet.planetmap.render` paints it with **no edit** (the kind-dispatching cone overlay the
``VECTOR_OVERLAY`` seam already built). The honest spatial **coverage-extent** and the **provenance**
(which producer made it) ride in the JSON-safe ``style`` dict beside the data — cast to native Python
scalars so the JSON bounce is loss-free and the round-trip ``==`` holds.

The embedding — the load-bearing honesty edge
----------------------------------------------
The field's ``(u, v)`` is laid onto a **full-globe** grid with **zeros outside its coverage box**. Two
reasons this is the right shape, not a compromise:

* **No smear, no disclaimer.** A full-globe grid trips :func:`planet.planetmap._polecapped`'s no-op
  guard, so the rows repeated to the poles are zeros — no fabricated polar flow. Zeros = *no flow*, so
  the result is honest-by-construction (the band-vs-globe edge is carried *in the data*, not a caption).
* **It is the spin-out's target shape.** An ECCO ocean field is a full-globe lat/lon ``(u, v)`` + SST
  with data everywhere; *"full-globe grid, ``is_global=False``, real data in a sub-box, coverage records
  the box"* is exactly what S1 binds on. ``is_global=False`` on a full-globe grid is not a contradiction
  — coverage describes *where the data is*, not the grid's extent.

**The band is placed honestly** (unlike the zonal-mean :func:`planet.planetmap.circulation_layer`,
which mirrors ``|lat|`` to both hemispheres and repeats across all longitude — valid only for a
zonally-symmetric jet): the emergent eddy band is the project's **only** longitudinally-structured
field, so it is laid down **NH-only, in its true ~55° longitude sector only, zeros everywhere else** —
never mirrored to the southern hemisphere, never wrapped around the globe. Mirroring or wrapping it
would fabricate flow the model never produced (the honesty edge the eddy memories police).

The mask increment (§9.6 O1) — validity inside the box
-------------------------------------------------------
:class:`Coverage` is a bounding *box*; a real ocean field has no-data cells **inside** it (land, and the
unobserved poles of a cell-centred product). O1 grows the contract once: ``FlowField.mask`` (``True`` =
valid; ``None`` = all-valid = the exact pre-mask path). Three resample rules come straight from the O2
OSCAR spike (plan §9.6): the mask resamples **nearest-neighbour**, never bilinearly (a thresholded
interpolated mask smears the coastline); destination latitudes **beyond the source's range come out
masked** (OSCAR stops at ±89.75° — extrapolated poles would be fabricated ocean); and the mask is
**applied** to the embedded arrays (zero ``(u, v)``, NaN scalar on invalid cells), so filled-zero land
can never read back as measured zero current. Serialized, the mask is one more (categorical 0/1) layer
in the same ``.npz`` — schema-additive, so a no-mask spec loads exactly as before.

Frames — a named, deferred increment
------------------------------------
The plan's R1 list names *frames* (a time axis of ``(u, v, scalar)``). R1 serializes a single
(saturated) **snapshot**: a time axis is **orthogonal** to producer-agnosticism (it is the §9.5
time-animation concern) and is trivially a stacked ``.npz`` array later — so it is a conscious deferral,
not a gap, revisited at S1 against a real ocean field's dimensions (the *Retarget-when-done* rule).

This module is **NumPy-only at import** (it imports the headless halves of :mod:`planet.planetmap`,
:mod:`planet.planet_spec`, and :mod:`planet.flow_globe` — no Plotly / matplotlib), so the round-trip
test is always-green on a bare core install.
"""
from __future__ import annotations

import numpy as np

from .flow_globe import Coverage, FlowField
from .planet_spec import SCHEMA_VERSION, PlanetSpec
from .planetmap import Grid, Layer, LayerKind, PlanetView

# The interchange target grid — a full-globe lat×lon mesh that **includes the poles and ±180°** so the
# renderer's pole-cap is a no-op (no smear). Coarse is fine: the field is a reach artifact, the cones
# sub-sample anyway, and the round-trip identity (the real proof) is resolution-independent.
GLOBE_N_LAT = 91                 # −90 … +90 inclusive (2° spacing)
GLOBE_N_LON = 181                # −180 … +180 inclusive (2° spacing)

VECTOR_LAYER = "circulation"     # the VECTOR_OVERLAY layer name (mirrors planetmap.circulation_layer)
SPEED_LAYER = "speed"            # the universal scalar: |(u, v)| — present for EVERY producer
SCALAR_LAYER = "scalar"          # the field's own optional scalar (θ for the eddy band; absent for synthetic)
MASK_LAYER = "mask"              # the optional per-cell validity mask (§9.6 O1) — absent when all-valid


def _globe_grid(n_lat: int = GLOBE_N_LAT, n_lon: int = GLOBE_N_LON) -> Grid:
    """A full-globe lat×lon :class:`~planet.planetmap.Grid` spanning the poles and ±180° exactly."""
    return Grid(lat=np.linspace(-90.0, 90.0, n_lat), lon=np.linspace(-180.0, 180.0, n_lon))


def _bilinear(src_lat: np.ndarray, src_lon: np.ndarray, f: np.ndarray,
              dst_lat: np.ndarray, dst_lon: np.ndarray) -> np.ndarray:
    """Bilinearly resample ``f`` (on the monotone ``src_lat`` × ``src_lon`` mesh) onto ``dst_lat`` × ``dst_lon``.

    Two 1-D :func:`numpy.interp` passes (over longitude per source row, then over latitude per
    destination column) — NumPy-only, no SciPy, so the schema stays bare-core. :func:`numpy.interp`
    clamps to the edge value outside the source range; the caller masks to zero outside the coverage box
    (a non-global field), so that edge-clamp never paints flow beyond the box.
    """
    f = np.asarray(f, dtype=float)
    over_lon = np.empty((src_lat.size, dst_lon.size))
    for j in range(src_lat.size):
        over_lon[j] = np.interp(dst_lon, src_lon, f[j])
    out = np.empty((dst_lat.size, dst_lon.size))
    for i in range(dst_lon.size):
        out[:, i] = np.interp(dst_lat, src_lat, over_lon[:, i])
    return out


def _nearest_mask(src_lat: np.ndarray, src_lon: np.ndarray, mask: np.ndarray,
                  dst_lat: np.ndarray, dst_lon: np.ndarray) -> np.ndarray:
    """Nearest-neighbour resample of the boolean validity ``mask`` onto ``dst_lat`` × ``dst_lon``.

    Deliberately **not** :func:`_bilinear` (the §9.6 O2-spike retarget): interpolating a 0/1 mask makes
    fractional coastline values, and thresholding those smears the land/ocean boundary by a source cell.
    Nearest-neighbour keeps the mask boolean and the coastline crisp.

    **Pole honesty** rides here too: destination *latitudes beyond the source's range* come out
    ``False`` — a cell-centred real product (OSCAR stops at ±89.75°) must not be edge-clamp-extrapolated
    into "measured" polar ocean. Longitude is periodic, so there the nearest edge cell is the honest
    neighbour (at worst half a source cell away) and stays open.
    """
    src_lat = np.asarray(src_lat, dtype=float)
    src_lon = np.asarray(src_lon, dtype=float)
    i = np.abs(dst_lat[:, None] - src_lat[None, :]).argmin(axis=1)
    j = np.abs(dst_lon[:, None] - src_lon[None, :]).argmin(axis=1)
    out = np.asarray(mask, dtype=bool)[np.ix_(i, j)]
    in_lat = (dst_lat >= src_lat.min()) & (dst_lat <= src_lat.max())
    return out & in_lat[:, None]


def _embed_on_globe(field: FlowField, grid: Grid):
    """Lay a :class:`FlowField` onto the full-globe ``grid`` — zeros outside its coverage box.

    Returns ``(u, v, scalar, mask)`` each ``(n_lat, n_lon)`` on the globe grid (``scalar``/``mask`` are
    ``None`` if the field carries none). A **non-global** field is masked to zero outside its coverage
    box — the band's own NH sector — so nothing is painted where the model resolves nothing (the honesty
    edge). A global field fills the whole grid.

    A field with a **validity mask** (§9.6 O1) gets it resampled nearest-neighbour (:func:`_nearest_mask`
    — boolean-preserving, coastline-crisp, poles-beyond-the-source masked) and then **applied**: an
    invalid cell carries zero ``(u, v)`` and NaN ``scalar``, so a filled-zero land cell in the source can
    never read back as "measured zero current" — the mask travels beside the data *and* the data agrees
    with it.
    """
    dst_lat, dst_lon = np.asarray(grid.lat, dtype=float), np.asarray(grid.lon, dtype=float)
    u = _bilinear(field.lat, field.lon, field.u, dst_lat, dst_lon)
    v = _bilinear(field.lat, field.lon, field.v, dst_lat, dst_lon)
    scalar = None if field.scalar is None else _bilinear(field.lat, field.lon, field.scalar, dst_lat, dst_lon)
    mask = None if field.mask is None else _nearest_mask(field.lat, field.lon, field.mask, dst_lat, dst_lon)

    if not field.coverage.is_global:
        cov = field.coverage
        inbox = ((dst_lat[:, None] >= cov.lat_min) & (dst_lat[:, None] <= cov.lat_max)
                 & (dst_lon[None, :] >= cov.lon_min) & (dst_lon[None, :] <= cov.lon_max))
        u = np.where(inbox, u, 0.0)
        v = np.where(inbox, v, 0.0)
        if scalar is not None:
            scalar = np.where(inbox, scalar, np.nan)      # NaN, not 0 — "no data here", not "θ = 0 °C"
        if mask is not None:
            mask = mask & inbox                           # outside the box is by definition not valid data

    if mask is not None:
        u = np.where(mask, u, 0.0)                        # no flow on invalid cells — zeros, like the box edge
        v = np.where(mask, v, 0.0)
        if scalar is not None:
            scalar = np.where(mask, scalar, np.nan)       # "no data", not "scalar = 0"
    return u, v, scalar, mask


def _coverage_style(cov: Coverage) -> dict:
    """The coverage extent as a **JSON-native** dict (plain ``float``/``bool`` — survives the save/load ``==``)."""
    return {"lat_min": float(cov.lat_min), "lat_max": float(cov.lat_max),
            "lon_min": float(cov.lon_min), "lon_max": float(cov.lon_max),
            "is_global": bool(cov.is_global)}


def vector_view_from_flow_field(field: FlowField, *, provenance: str,
                                n_lat: int = GLOBE_N_LAT, n_lon: int = GLOBE_N_LON) -> PlanetView:
    """A :class:`FlowField` → a renderable, serializable :class:`~planet.planetmap.PlanetView` (the export side).

    Embeds the field onto a full-globe grid (zeros outside its coverage box) and registers:

    * a **speed** ``SCALAR_FIELD`` — ``|(u, v)|``, the *universal* scalar present for **every** producer
      (a vector-only view cannot be rendered: :func:`planet.planetmap.render` needs a scalar surface, and
      speed is native to any vector field — so the eddy and the synthetic views are structurally parallel);
    * a **circulation** ``VECTOR_OVERLAY`` — the stacked ``(2, n_lat, n_lon)`` ``[u, v]`` (m/s), carrying
      the **coverage-extent**, the **provenance**, the planet ``radius_m``, and the field's **honesty**
      string in its JSON-safe ``style`` (all cast native so the round-trip ``==`` holds);
    * the field's own **scalar** ``SCALAR_FIELD`` (e.g. θ for the eddy band), *if* it carries one — NaN
      outside the coverage box ("no data", not "θ = 0");
    * the **mask** ``SCALAR_FIELD`` (§9.6 O1), *if* the field carries a validity mask — the 0/1
      per-cell "is there data here" map (categorical style, so ``render(active="mask")`` paints an
      honest coverage globe). It is one more array in the same ``.npz`` — the schema is **additive**:
      a no-mask spec has no such layer and loads exactly as before.

    Both producers (the real eddy band, a synthetic global field) flow through this one function, so the
    serialized artifact and the rendered globe are identical in structure regardless of origin.
    """
    grid = _globe_grid(n_lat, n_lon)
    u, v, scalar, mask = _embed_on_globe(field, grid)
    speed = np.hypot(u, v)

    cov = field.coverage
    extent = f"{cov.lat_max - cov.lat_min:.0f}°×{cov.lon_max - cov.lon_min:.0f}° sector" if not cov.is_global else "global"
    vec_style = {
        "colorscale": "RdBu_r",
        "arrow_color": "#1a1a1a",
        "coverage": _coverage_style(cov),
        "provenance": str(provenance),
        "radius_m": float(field.radius_m),
        "honesty": str(field.honesty),
        "label": f"flow — {provenance} ({extent})",
    }
    layers = [
        Layer(SPEED_LAYER, LayerKind.SCALAR_FIELD, speed, "m/s",
              style={"colorscale": "Viridis", "colorbar_title": "flow speed (m/s)"}, z_order=0),
        Layer(VECTOR_LAYER, LayerKind.VECTOR_OVERLAY, np.stack([u, v]), "m/s",
              style=vec_style, z_order=2),
    ]
    if scalar is not None:
        # `scalar_label` ("θ (°C)") is a display label, not a unit — it rides in the colorbar title; the
        # unit slot stays honest (FlowField carries no separate unit string, §7 is unit-obsessive).
        layers.insert(1, Layer(SCALAR_LAYER, LayerKind.SCALAR_FIELD, scalar, "",
                               style={"colorscale": "RdBu_r",
                                      "colorbar_title": str(field.scalar_label or "scalar")}, z_order=1))
    if mask is not None:
        # 0/1 int (not bool) — the delta_view precedent: Plotly's categorical surface path wants codes.
        # The unit slot is honest ("" — a validity flag has no unit); the names spell the semantics out.
        layers.append(Layer(MASK_LAYER, LayerKind.SCALAR_FIELD, mask.astype(int), "",
                            style={"categorical": True, "colorbar_title": "data coverage",
                                   "colors": {"0": "#4a4a52", "1": "#2e6fba"},
                                   "names": {"0": "no data (land / unobserved)", "1": "valid cell"}},
                            z_order=3))
    return PlanetView(grid=grid, layers=tuple(layers))


def vector_spec_from_flow_field(field: FlowField, *, provenance: str,
                                n_lat: int = GLOBE_N_LAT, n_lon: int = GLOBE_N_LON) -> PlanetSpec:
    """A :class:`FlowField` → a saveable :class:`~planet.planet_spec.PlanetSpec` (no climate knobs).

    A pure *field* interchange: the spec carries the vector view's grid + layers with **empty knobs**
    (a serialized flow field has no :class:`~planet.albedo.EBMParams` — it is not a re-runnable climate,
    it is a field handed across the seam). :func:`planet.planet_spec.save` / ``load`` round-trip it
    unchanged, and the round-trip identity ``load(save(spec)) == spec`` is the R1 proof.
    """
    view = vector_view_from_flow_field(field, provenance=provenance, n_lat=n_lat, n_lon=n_lon)
    return PlanetSpec(schema_version=SCHEMA_VERSION, grid=view.grid, layers=view.layers,
                      knobs={}, knob_units={})


def synthetic_flow_field(*, n_lat: int = GLOBE_N_LAT, n_lon: int = GLOBE_N_LON) -> FlowField:
    """A **synthetic global** ``(u, v)`` field — the second producer that proves producer-agnosticism.

    An analytic, clearly-not-a-model global flow on the full-globe grid: banded zonal winds
    ``u = U₀·cos(2φ)`` (easterly tropics, westerly mid-latitudes, easterly polar — the recognizable
    shape) plus a longitudinal wave ``v = V₀·cosφ·sin(3λ)`` giving genuine meridional flow (so the
    ``v`` component is exercised, not just ``u``). It carries ``is_global=True`` so it fills the globe —
    the deliberate contrast with the eddy band's ``is_global=False`` that *exercises the coverage
    machinery* (a named R1 deliverable: coverage only bites when one producer is non-global).

    Its honesty label states plainly that it is a synthetic probe, not a physical result — the
    interchange carries provenance/honesty for *every* producer, not only the model's own fields.
    """
    lat = np.linspace(-90.0, 90.0, n_lat)
    lon = np.linspace(-180.0, 180.0, n_lon)
    LON, LAT = np.meshgrid(np.radians(lon), np.radians(lat))     # (n_lat, n_lon)
    u0, v0 = 30.0, 8.0
    u = u0 * np.cos(2.0 * LAT)
    v = v0 * np.cos(LAT) * np.sin(3.0 * LON)
    honesty = (
        "Synthetic analytic global flow — NOT a physical model result. This is a producer-agnosticism "
        "probe: a fabricated global (u, v) field handed through the SAME serialization and renderer as "
        "the model's emergent eddy band, to show the seam does not care what produced the field. The "
        "pattern (banded zonal winds + a longitudinal wave) is illustrative only."
    )
    coverage = Coverage(lat_min=-90.0, lat_max=90.0, lon_min=-180.0, lon_max=180.0, is_global=True)
    return FlowField(lat=lat, lon=lon, u=u, v=v, coverage=coverage, honesty=honesty,
                     scalar=None, scalar_label="", radius_m=6.371e6)
