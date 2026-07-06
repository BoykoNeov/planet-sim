"""The O2 banked artifact: real OSCAR ocean surface currents on the particle flow-globe (§9.6).

One OSCAR global snapshot (0.25°, daily mean, NASA PO.DAAC) flows through the pipe built for exactly
this — :func:`planet.ocean_currents.flow_field_from_ocean` → the **unchanged** R1 serialization (the
round-trip identity is asserted here on the *real* field, no longer only on synthetic probes) → the
**unchanged** Rung-C three.js globe → ``docs/figures/planet-ocean-currents.html``.

Data discipline (plan §9.6): the raw 33 MB granule is **never committed** — it lives in the gitignored
``outputs/`` and is fetched on demand with an Earthdata bearer token (``EARTHDATA_TOKEN`` env var; get
one free at https://urs.earthdata.nasa.gov → *Generate Token*). The token is sent as one
``Authorization: Bearer`` header (the spike's finding — no ``.netrc``/URS redirect dance) and is never
logged or written anywhere.

The two conscious calls the plan flagged forward:

* **Resolution** — the *banked render* consumes the field at 0.5° (``STRIDE = 2`` off the native 0.25°),
  not the interchange's 2° proof grid ("beautiful dies at 2°"): 0.5° keeps the western-boundary currents
  (Gulf Stream, Kuroshio, Agulhas) sharp while the self-contained HTML stays committable; native 0.25°
  quadruples the payload for detail a 6500–20000-particle render cannot resolve anyway.
* **Pace** — the renderer's auto-accel was tuned so the fastest particle crosses the ~55° eddy band in
  ~6 s; on 360° of globe that sprint reads as a storm, so the global field passes
  ``crossing_seconds=45`` (the additive default-off knob; the eddy artifact is byte-unchanged).

Run headless (downloads if needed, saves the HTML, prints the summary)::

    python -m planet.demo_ocean_currents [path\\to\\oscar_granule.nc]
"""
from __future__ import annotations

import os
from pathlib import Path

from .flow_globe import save_flow_globe_html
from .ocean_currents import OceanSnapshot, flow_field_from_ocean, load_oscar

_REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_FIGURE = _REPO_ROOT / "docs" / "figures" / "planet-ocean-currents.html"
OUTPUT_FIGURE = _REPO_ROOT / "outputs" / "planet-ocean-currents.html"

GRANULE = "oscar_currents_final_20200601.nc"
GRANULE_URL = ("https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-protected/"
               "OSCAR_L4_OC_FINAL_V2.0/" + GRANULE)
DEFAULT_NC = _REPO_ROOT / "outputs" / GRANULE          # gitignored — the raw granule is never committed
TOKEN_ENV = "EARTHDATA_TOKEN"

STRIDE = 2                 # 0.5° render grid — the conscious size/beauty call (module docstring)
N_PARTICLES = 20000        # a full globe wants more than the band's 6500 (same visual density)
CROSSING_SECONDS = 45.0    # global pace: fastest particle crosses 360° in ~45 s (band default: 55° in 6 s)

OSCAR_PROVENANCE = "planet.ocean_currents (OSCAR L4 v2.0, PO.DAAC)"

TITLE = "planet-sim — real ocean surface currents (OSCAR)"
SUBTITLE = ("one day of real reanalysis-class surface currents (OSCAR L4 v2.0, NASA PO.DAAC) streamed "
            "through the same seam and renderer as the emergent eddy band (§9.6 O2) — not a planet-sim "
            "simulation")


def fetch_oscar(dest: Path = DEFAULT_NC, url: str = GRANULE_URL) -> Path:
    """Download the granule with an Earthdata bearer token (``EARTHDATA_TOKEN``). No-op if present.

    The token rides in one ``Authorization: Bearer`` header (the O2 spike's auth finding); it is read
    from the environment and never echoed, logged, or persisted. Raises a clear error naming the env var
    when it is absent — the granule cannot be fetched anonymously.
    """
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    token = os.environ.get(TOKEN_ENV, "").strip()
    if not token:
        raise RuntimeError(
            f"{dest.name} not found and ${TOKEN_ENV} is not set. Get a free token at "
            f"https://urs.earthdata.nasa.gov (Generate Token), then: set {TOKEN_ENV}=<token> "
            f"and re-run — or download the granule yourself and pass its path:\n  {url}"
        )
    import requests  # the [ocean] extra — lazy, like the reader

    dest.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(url, headers={"Authorization": f"Bearer {token}"}, stream=True, timeout=120) as r:
        r.raise_for_status()
        tmp = dest.with_suffix(dest.suffix + ".part")
        with open(tmp, "wb") as fh:
            for chunk in r.iter_content(chunk_size=1 << 20):
                fh.write(chunk)
        tmp.replace(dest)
    return dest


def compute(nc_path, stride: int = STRIDE):
    """Load the granule at the render stride and build the contract field. Returns ``(snapshot, field)``."""
    snap = load_oscar(nc_path, stride=stride)
    return snap, flow_field_from_ocean(snap)


def verify_round_trip(field) -> None:
    """Assert the R1 round-trip identity on the *real* field — the proof, now non-synthetic (§9.6 O2).

    Runs on the standard 2° interchange grid (the proof stays resolution-independent and cheap; the
    *render* is what consumes the fine grid). Raises ``AssertionError`` on any mismatch.
    """
    import tempfile

    from . import planet_spec as ps
    from .flow_serialize import vector_spec_from_flow_field

    spec = vector_spec_from_flow_field(field, provenance=f"{OSCAR_PROVENANCE} (real data)")
    with tempfile.TemporaryDirectory() as td:
        ps.save(spec, Path(td) / "ocean")
        assert ps.load(Path(td) / "ocean") == spec, "round-trip identity failed on the real ocean field"


def save_globe(field) -> Path:
    """Write the O2 showcase HTML (docs + outputs copies). Returns the docs path."""
    for target in (DOCS_FIGURE, OUTPUT_FIGURE):
        save_flow_globe_html(field, target, title=TITLE, subtitle=SUBTITLE,
                             n_particles=N_PARTICLES, crossing_seconds=CROSSING_SECONDS,
                             colormap="speed",     # ocean speed is a 0→max field → the sequential ramp (§9.6 O3c)
                             trails=True)           # the Perpetual-Ocean motion-trail look (§9.6 O3b)
    return DOCS_FIGURE


def print_summary(snap: OceanSnapshot, field) -> None:
    ny, nx = field.u.shape
    ocean_pct = 100.0 * float(field.mask.mean())
    print(f"OSCAR snapshot {snap.date}: {ny}x{nx} cells "
          f"({180.0 / max(1, ny - 1):.2g}-deg render grid), {ocean_pct:.0f}% valid ocean")
    print(f"  |current| max {float(field.scalar.max()):.2f} m/s — {snap.product}, DOI {snap.doi}")


def main() -> None:
    import sys
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")

    nc = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(os.environ.get("OSCAR_NC", DEFAULT_NC))
    nc = fetch_oscar(nc) if not nc.exists() else nc
    snap, field = compute(nc)
    print_summary(snap, field)
    verify_round_trip(field)
    print("  round-trip identity on the REAL field: OK (the R1 proof, now non-synthetic)")
    saved = save_globe(field)
    print(f"Ocean-currents flow-globe saved -> {saved.relative_to(_REPO_ROOT)}")
    print("  open it in a browser (works straight off disk — three.js is vendored inline).")


if __name__ == "__main__":
    main()
