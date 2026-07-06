"""The O4 banked artifact: real OSCAR ocean currents through the year — the seasonal flow-globe (§9.6).

O2 streamed **one** OSCAR day; O4 gives the flow a **time axis**. Twelve OSCAR snapshots — one
representative day per month of 2020 — flow through :func:`planet.ocean_currents.flow_field_from_ocean_series`
into a :class:`~planet.flow_globe.FlowFrames`; the **unchanged** R1 serialization round-trips the framed
field, and the Rung-C three.js globe **crossfades** between months so the particles steer through the
seasonal cycle → ``docs/figures/planet-ocean-currents-seasonal.html``. The showpiece is the **Somali
Current monsoon reversal** — the one major current on Earth that flips direction with the monsoon
(SW-monsoon Jun–Sep vs NE-monsoon Nov–Feb) — legible as the time badge cycles the month.

**What these frames ARE (the honesty the label states, user-chosen 2026-07-06).** They are **twelve
monthly snapshots — one day per month of 2020**, NOT monthly means and NOT a multi-year climatology. The
Somali reversal is a large monsoon signal that reads clearly in a day-per-month series; the label says
exactly that (``"monthly snapshots (one day per month), 2020"``) so nothing over-claims. Averaging into
true means (one year) or a climatology (many years) is a heavier download left as a named alternative —
the producer is acquisition-agnostic, so only this demo's fetch + period phrase would change.

Two conscious calls (the O4 versions of O2's resolution/pace calls):

* **Frames coarser than the O2 still.** The animation advects at **1.5°** (``STRIDE = 6`` off native
  0.25°) — motion hides resolution, and twelve frames at 0.5° would be a ~30 MB HTML. The per-frame
  *scalar* is dropped (colour = in-shader current speed), so the payload carries velocity only; the file
  lands comparable to O2's single-frame 5 MB.
* **Seasonal pace.** ``seconds_per_year`` sets how fast the year cycles (a leisurely default so the
  reversal is watchable); ``crossing_seconds`` still sets the particle streaming pace (the O2 global value).

Data discipline (as O2): the raw granules are **never committed** — they live in gitignored ``outputs/``
and are fetched on demand with an Earthdata bearer token (``EARTHDATA_TOKEN``; get one free at
https://urs.earthdata.nasa.gov → *Generate Token*). Twelve granules (~400 MB) — a token hand-off.

Run headless (downloads the twelve granules if needed, saves the HTML, prints the summary)::

    python -m planet.demo_ocean_seasonal
"""
from __future__ import annotations

import os
from pathlib import Path

from .demo_ocean_currents import TOKEN_ENV, fetch_oscar
from .flow_globe import save_flow_globe_html
from .ocean_currents import flow_field_from_ocean_series, load_oscar

_REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_FIGURE = _REPO_ROOT / "docs" / "figures" / "planet-ocean-currents-seasonal.html"
OUTPUT_FIGURE = _REPO_ROOT / "outputs" / "planet-ocean-currents-seasonal.html"

YEAR = 2020
DAY = 15                    # a representative mid-month day (one snapshot per month)
GRANULE_BASE = ("https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-protected/"
                "OSCAR_L4_OC_FINAL_V2.0/")

STRIDE = 6                 # 1.5° animation grid — the conscious size call (module docstring); motion hides it
N_PARTICLES = 20000        # a full globe, same visual density as O2
CROSSING_SECONDS = 45.0    # particle streaming pace (fastest crosses 360° in ~45 s) — the O2 global value
SECONDS_PER_YEAR = 30.0    # one full trip through the twelve months in ~30 s (leisurely — the reversal reads)

# The honesty phrase for what the frames ACTUALLY are (user-chosen: day-per-month, NOT means/climatology).
PERIOD = f"twelve monthly snapshots (one day per month), {YEAR}"

TITLE = "planet-sim — real ocean surface currents through the year (OSCAR)"
SUBTITLE = ("twelve months of real reanalysis-class surface currents (OSCAR L4 v2.0, NASA PO.DAAC), one "
            "day per month of 2020, crossfaded on the same renderer as the eddy band (§9.6 O4) — watch the "
            "Somali Current reverse with the monsoon (time badge, top-right); not a planet-sim simulation")


def _granule(month: int) -> tuple[str, Path]:
    """The (URL, gitignored local path) for one month's representative-day granule."""
    name = f"oscar_currents_final_{YEAR}{month:02d}{DAY:02d}.nc"
    return GRANULE_BASE + name, _REPO_ROOT / "outputs" / name


def fetch_year(months=range(1, 13)) -> list[Path]:
    """Download the twelve monthly granules (bearer token), skipping any already present. Returns the paths."""
    paths = []
    for m in months:
        url, dest = _granule(m)
        paths.append(fetch_oscar(dest, url) if not dest.exists() else dest)
    return paths


def compute(nc_paths, stride: int = STRIDE):
    """Load each granule at the render stride → the framed contract field. Returns ``(snapshots, field)``."""
    snaps = [load_oscar(p, stride=stride) for p in nc_paths]
    field = flow_field_from_ocean_series(snaps, period=PERIOD, pace_note="one year compressed to ~30 s")
    return snaps, field


def verify_round_trip(field) -> None:
    """Assert the R1 round-trip identity on the *framed* field — the proof extends to the time axis (§9.6 O4)."""
    import tempfile

    from . import planet_spec as ps
    from .flow_serialize import FRAMES_LAYER, vector_spec_from_flow_field

    spec = vector_spec_from_flow_field(field, provenance="planet.ocean_currents (OSCAR series, real data)")
    assert any(ly.name == FRAMES_LAYER for ly in spec.layers), "the seasonal frames layer is missing"
    with tempfile.TemporaryDirectory() as td:
        ps.save(spec, Path(td) / "seasonal")
        assert ps.load(Path(td) / "seasonal") == spec, "round-trip identity failed on the framed ocean field"


def save_globe(field) -> Path:
    """Write the O4 showcase HTML (docs + outputs copies). Returns the docs path."""
    for target in (DOCS_FIGURE, OUTPUT_FIGURE):
        save_flow_globe_html(field, target, title=TITLE, subtitle=SUBTITLE,
                             n_particles=N_PARTICLES, crossing_seconds=CROSSING_SECONDS,
                             seconds_per_year=SECONDS_PER_YEAR,
                             colormap="speed",     # ocean speed is a 0→max field → the sequential ramp (§9.6 O3c)
                             trails=True)           # the Perpetual-Ocean motion-trail look (§9.6 O3b)
    return DOCS_FIGURE


def print_summary(snaps, field) -> None:
    import numpy as np

    nt, ny, nx = field.frames.u.shape
    ocean_pct = 100.0 * float(field.mask.mean())
    speed_max = float(np.hypot(field.frames.u, field.frames.v).max())
    print(f"OSCAR seasonal series: {nt} frames, {ny}x{nx} cells "
          f"({180.0 / max(1, ny - 1):.2g}-deg render grid), {ocean_pct:.0f}% valid-in-every-frame ocean")
    print(f"  frames: {', '.join(field.frames.labels)}  ({PERIOD})")
    print(f"  |current| max {speed_max:.2f} m/s — {snaps[0].product}, DOI {snaps[0].doi}")


def main() -> None:
    import sys
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")

    if not os.environ.get(TOKEN_ENV, "").strip() and not all(_granule(m)[1].exists() for m in range(1, 13)):
        print(f"Twelve OSCAR granules are needed and ${TOKEN_ENV} is not set.\n"
              f"  Get a free token at https://urs.earthdata.nasa.gov (Generate Token), then:\n"
              f"    set {TOKEN_ENV}=<token>  &&  python -m planet.demo_ocean_seasonal")
        raise SystemExit(1)

    paths = fetch_year()
    snaps, field = compute(paths)
    print_summary(snaps, field)
    verify_round_trip(field)
    print("  round-trip identity on the framed field: OK (the R1 proof, now with a time axis)")
    saved = save_globe(field)
    print(f"Seasonal ocean-currents flow-globe saved -> {saved.relative_to(_REPO_ROOT)}")
    print("  open it in a browser — the time badge (top-right) cycles the month; watch the Somali Current flip.")


if __name__ == "__main__":
    main()
