"""The real-ocean producer — OSCAR surface currents onto the flow-globe contract (§9.6 O2).

Every :class:`~planet.flow_globe.FlowField` consumer so far carried this project's *own* output (the
emergent eddy band) or a labelled synthetic probe. O2 feeds the seam its first **real-world** field: one
OSCAR global surface-current snapshot (NASA PO.DAAC, 0.25°, daily mean) → :func:`flow_field_from_ocean`
→ the **unchanged** Rung-C particle globe and the **unchanged** R1 serialization. Nothing here computes
ocean physics — this module *reads and reshapes measured-class data*; the spin-out boundary (plan §11)
is untouched.

The honesty class this producer introduces (the provenance clause)
------------------------------------------------------------------
The eddy showcase's disclaimer says "illustrative, mostly reversible" about *our* model. An OSCAR field
inverts the problem: the currents are *real* (reanalysis-class), but they are **not this project's
result** — a viewer who just watched the emergent eddy globe could reasonably assume planet-sim now
simulates oceans. So the honesty string leads with **provenance** ("real reanalysis-class ocean surface
currents — OSCAR …; NOT computed by this project's models") and it is machine-checked as a visible DOM
element exactly like Rung C's reversibility clause (``test_ocean_currents.py``).

The three O2-spike facts this module encodes (plan §9.6, spiked 2026-07-06)
---------------------------------------------------------------------------
* **Dim order**: the granule stores ``(time, longitude, latitude)`` — lon *before* lat, the transpose of
  the contract's ``(n_lat, n_lon)``. :func:`load_oscar` transposes; it is file-layout knowledge and lives
  only in the loader.
* **Longitude convention**: 0…359.75° (not ±180). The rewrap ``((lon+180) % 360) − 180`` un-sorts the
  axis, so :func:`flow_field_from_ocean` re-orders columns (`argsort`) before anything downstream sees
  the grid (``_bilinear``/the renderer's ``sample()`` both assume monotone axes).
* **NaN = land/missing (~44 % of cells)**: a NaN would poison the renderer's bilinear sampling and the
  interchange's :func:`~planet.flow_serialize._bilinear` at every coastline. So the producer (a) derives
  the O1 validity **mask** from finiteness, then (b) fills NaN→0 in ``(u, v)`` — and the mask is what
  keeps a filled zero from ever reading back as "measured zero current" (the O1 rule: the serializer
  re-applies the mask, the renderer seeds/keeps particles only on valid cells).

The grid is **cell-centred** (±89.75° — OSCAR never observes the true poles); the coverage box records
that honestly, and the O1 pole rule in :func:`~planet.flow_serialize._nearest_mask` masks the
interchange grid's ±90° rows rather than extrapolating them into "measured" polar ocean.

Import discipline: **NumPy-only at import** — the netCDF reader (`h5netcdf`, the ``[ocean]`` extra) is
imported lazily inside :func:`load_oscar`, so the producer's pure-array half (and its tests, which run
off a small committed fixture) needs no optional dependency. The raw 33 MB granule is **never
committed**; see :mod:`planet.demo_ocean_currents` for the token-authenticated download.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from .flow_globe import Coverage, FlowField, FlowFrames

# Calendar-month abbreviations — the default frame labels for a 12-snapshot monthly series (the O4 payload);
# any other length falls back to the snapshots' own dates, so the labels always say what the frames are.
_MONTH_ABBR = ("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

# The granule the spike settled on (one day, global, 0.25°, ~33 MB) — the demo's default target.
OSCAR_PRODUCT = "OSCAR L4 v2.0 total surface currents"
OSCAR_DOI = "10.5067/OSCAR-25F20"
OSCAR_CREDIT = "NASA PO.DAAC / Earth & Space Research"


@dataclass(frozen=True)
class OceanSnapshot:
    """One real ocean surface-current snapshot, loader-normalized but convention-raw.

    ``lat`` (n_lat,) / ``lon`` (n_lon,) cell-centre degrees, each monotone increasing but in the
    *product's own* longitude convention (OSCAR: 0…360); ``u``/``v`` (n_lat, n_lon) m/s with **NaN on
    land/unobserved cells** (fill values already converted). The provenance strings ride along so the
    honesty clause is composed from what the file actually says, not from hardcoded lore.
    """

    lat: np.ndarray
    lon: np.ndarray
    u: np.ndarray
    v: np.ndarray
    product: str = OSCAR_PRODUCT
    doi: str = OSCAR_DOI
    credit: str = OSCAR_CREDIT
    date: str = ""
    depth_note: str = "average over the top 30 m of the mixed layer"


def load_oscar(path, *, stride: int = 1, geostrophic_only: bool = False) -> OceanSnapshot:
    """Read one OSCAR v2.0 granule (netCDF4/HDF5) into an :class:`OceanSnapshot`.

    ``stride`` subsamples both axes (1 = native 0.25°; 2 = 0.5°, the banked-artifact default — the
    conscious size call, see :mod:`planet.demo_ocean_currents`). ``geostrophic_only`` reads ``ug``/``vg``
    instead of the total (geostrophic + Ekman) ``u``/``v`` — the named future knob from the spike.

    Lazy-imports ``h5netcdf`` (the ``[ocean]`` extra) so this module stays NumPy-only at import. Handles
    the granule's ``(time, longitude, latitude)`` dim order (transposed here, nowhere else) and converts
    the ``_FillValue`` to NaN so "no data" has exactly one spelling downstream.
    """
    import h5netcdf  # the [ocean] extra — deliberately not a module-level import

    uname, vname = ("ug", "vg") if geostrophic_only else ("u", "v")
    with h5netcdf.File(str(path), "r") as f:
        lat = np.asarray(f.variables["lat"][::stride], dtype=float)
        lon = np.asarray(f.variables["lon"][::stride], dtype=float)

        def read(name):
            var = f.variables[name]
            raw = np.asarray(var[0, ::stride, ::stride], dtype=float).T   # (lon, lat) → (lat, lon)
            fill = var.attrs.get("_FillValue")
            if fill is not None:
                raw = np.where(raw == float(fill), np.nan, raw)
            return raw

        u, v = read(uname), read(vname)
        attrs = {k: f.attrs[k] for k in f.attrs}

    def text(key, default=""):
        val = attrs.get(key, default)
        return val.decode() if isinstance(val, bytes) else str(val)

    date = text("time_coverage_start")[:10]
    version = text("product_version", "v2.0")
    kind = "geostrophic-only surface currents" if geostrophic_only else "total (geostrophic + Ekman) surface currents"
    return OceanSnapshot(lat=lat, lon=lon, u=u, v=v,
                         product=f"OSCAR L4 {version} {kind}",
                         doi=OSCAR_DOI, credit=OSCAR_CREDIT, date=date)


def _rewrap_lon(lon: np.ndarray, *fields: np.ndarray):
    """0…360 → ±180 longitude, with the columns re-ordered so the axis is monotone again.

    The rewrap alone leaves the array as two ascending runs (0…180, then −180…0 appended); ``argsort``
    (stable) is the roll that restores one monotone axis — required by every downstream bilinear sampler.
    A grid already in ±180 passes through unchanged (the sort is the identity).
    """
    lon180 = ((np.asarray(lon, dtype=float) + 180.0) % 360.0) - 180.0
    order = np.argsort(lon180, kind="stable")
    return (lon180[order], *[np.asarray(fld)[:, order] for fld in fields])


def flow_field_from_ocean(snap: OceanSnapshot) -> FlowField:
    """An :class:`OceanSnapshot` → the renderer/serializer contract — the O2 producer.

    Pure array work (no I/O, no optional deps): rewrap+re-sort the longitude axis to ±180, derive the O1
    validity ``mask`` from finiteness (NaN = land/unobserved), **then** fill NaN→0 in ``(u, v)`` so no
    NaN ever reaches a bilinear sampler — the mask, not the fill value, carries "no data here". The
    scalar is the current **speed** (over the filled arrays; masked cells are zero and are never sampled
    by the renderer, and the serializer re-masks them to NaN). ``is_global=True`` — the contract's first
    true-global consumer; the coverage box still records the honest *data* extent (cell-centred ±89.75°,
    short of the true poles).
    """
    lon, u, v = _rewrap_lon(snap.lon, snap.u, snap.v)
    lat = np.asarray(snap.lat, dtype=float)
    mask = np.isfinite(u) & np.isfinite(v)
    u = np.where(mask, u, 0.0)
    v = np.where(mask, v, 0.0)
    speed = np.hypot(u, v)

    ocean_pct = round(100.0 * float(mask.mean()))
    honesty = (
        f"REAL data, not this project's model: these are reanalysis-class ocean surface currents — "
        f"{snap.product}, daily mean for {snap.date or 'one day'}, {snap.credit}, DOI {snap.doi} — "
        f"measured-and-derived from satellite observations, NOT computed by planet-sim's models "
        f"(this project renders the field; it did not produce it). The velocity is the "
        f"{snap.depth_note}. The streaming particles are illustrative and heavily time-accelerated; "
        f"they trace the snapshot's velocity field, not actual water parcels. Land and unobserved "
        f"cells (~{100 - ocean_pct}% of the grid, including the poles beyond ±{abs(lat).max():.2f}°) "
        f"carry no data and no particles."
    )
    coverage = Coverage(lat_min=float(lat.min()), lat_max=float(lat.max()),
                        lon_min=float(lon.min()), lon_max=float(lon.max()), is_global=True)
    return FlowField(lat=lat, lon=lon, u=u, v=v, coverage=coverage, honesty=honesty,
                     scalar=speed, scalar_label="current speed (m/s)", mask=mask)


def flow_field_from_ocean_series(snaps, *, labels=None, period: str = "",
                                 pace_note: str = "one model year compressed to a few seconds") -> FlowField:
    """A time-ordered list of :class:`OceanSnapshot` → a **framed** :class:`FlowField` — the §9.6 O4 producer.

    The seasonal-currents deliverable: the same rewrap→mask→fill pipeline as :func:`flow_field_from_ocean`,
    run over every frame and stacked into a :class:`~planet.flow_globe.FlowFrames` time axis the renderer
    crossfades through the year (the Somali-Current monsoon reversal is the showpiece). All snapshots must
    share one grid (the OSCAR product's lat/lon); ``snaps`` is expected already time-ordered.

    **Acquisition-agnostic (advisor call).** This function does not know or claim *how* the frames were
    built — whether a true multi-year climatology or one year's 12 monthly means. The caller (the demo)
    owns acquisition and passes ``period`` — the honest phrase describing what the series actually is (e.g.
    ``"monthly means for 2020"``) — which rides verbatim into the provenance clause. Do **not** let the
    word "climatology" leak in unless the data is one.

    Two honesty rules carry from O1/O2, tightened for the time axis:

    * **Static mask = finite in *every* frame** (conservative): a cell is valid only where all frames carry
      data, so a cell measured for only part of the year is left bare rather than blinking in and out —
      never over-claiming. Sea-ice seasonality (the ice edge that genuinely moves) is a documented
      second-order edge folded into this "no data in any frame → bare" rule, not animated.
    * The mask is **applied to every frame** (zero ``(u, v)`` on invalid cells) so a filled-zero land cell
      can never read back as a measured zero current — the O1 rule, now across the stack.

    ``labels`` overrides the per-frame captions (the time badge); the default is calendar-month
    abbreviations for a 12-frame series, else each snapshot's own date. The parent field's
    ``u``/``v``/``scalar`` are the **representative (frame-0)** snapshot (the CPU fallback and land/ocean
    base still work); the per-frame *scalar* is not stacked (the renderer colours by in-shader speed).
    """
    snaps = list(snaps)
    if not snaps:
        raise ValueError("flow_field_from_ocean_series needs at least one snapshot")
    lat = np.asarray(snaps[0].lat, dtype=float)
    lon180 = ((np.asarray(snaps[0].lon, dtype=float) + 180.0) % 360.0) - 180.0
    order = np.argsort(lon180, kind="stable")      # the same monotone re-sort flow_field_from_ocean applies
    lon = lon180[order]
    ny, nx = np.asarray(snaps[0].u).shape

    # per-frame rewrap; the validity mask is finite-in-ALL-frames (the conservative static mask, above).
    us, vs, finite = [], [], np.ones((ny, nx), dtype=bool)
    for s in snaps:
        u = np.asarray(s.u, dtype=float)
        v = np.asarray(s.v, dtype=float)
        if u.shape != (ny, nx) or v.shape != (ny, nx):
            raise ValueError("all snapshots must share one grid (same OSCAR product / stride)")
        u, v = u[:, order], v[:, order]
        finite &= np.isfinite(u) & np.isfinite(v)
        us.append(u)
        vs.append(v)
    mask = finite
    uf = np.stack([np.where(mask, u, 0.0) for u in us])    # (nt, ny, nx) — mask applied across the whole stack
    vf = np.stack([np.where(mask, v, 0.0) for v in vs])

    nt = len(snaps)
    if labels is not None:
        frame_labels = tuple(str(x) for x in labels)
    elif nt == 12:
        frame_labels = _MONTH_ABBR
    else:
        frame_labels = tuple((s.date or f"frame {i + 1}") for i, s in enumerate(snaps))

    prov = snaps[0]
    ocean_pct = round(100.0 * float(mask.mean()))
    span = period or f"a {nt}-frame seasonal series"
    honesty = (
        f"REAL data, not this project's model: these are reanalysis-class ocean surface currents — "
        f"{prov.product}, {prov.credit}, DOI {prov.doi} — measured-and-derived from satellite "
        f"observations, NOT computed by planet-sim's models (this project renders the field; it did not "
        f"produce it). This is a TIME series — {span}, {nt} frames — and the globe crossfades between "
        f"frames, heavily time-accelerated ({pace_note}); the streaming particles are illustrative: they "
        f"trace each frame's snapshot velocity, not actual water parcels, and the smooth morph between "
        f"frames is interpolation, not resolved day-to-day flow. The velocity is the {prov.depth_note}. "
        f"Land and unobserved cells (~{100 - ocean_pct}% of the grid, including the poles beyond "
        f"±{abs(lat).max():.2f}° and anywhere lacking data in ANY frame) carry no data and no particles."
    )
    coverage = Coverage(lat_min=float(lat.min()), lat_max=float(lat.max()),
                        lon_min=float(lon.min()), lon_max=float(lon.max()), is_global=True)
    frames = FlowFrames(u=uf, v=vf, labels=frame_labels)
    return FlowField(lat=lat, lon=lon, u=uf[0], v=vf[0], coverage=coverage, honesty=honesty,
                     scalar=np.hypot(uf[0], vf[0]), scalar_label="current speed (m/s)",
                     mask=mask, frames=frames)
