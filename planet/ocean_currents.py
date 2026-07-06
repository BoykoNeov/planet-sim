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

from .flow_globe import Coverage, FlowField

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
