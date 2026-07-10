"""Rung 5A: the mountain casts a rain shadow — the first step off the zonal mean.

Every rung up to here kept the climate **zonal-mean**: a single number per latitude circle, the same
at every longitude, and the elevation seam on the coarse globe was **inert** — carried in the schema but
driving nothing. This demo wakes it. A fine regional patch is placed on the sphere at a named range, a
north–south ridge is raised on it, the (purely zonal) mid-latitude westerly is read off the emergent jet
(:mod:`planet.coupler`), and the **Smith & Barstad (2004)** linear orographic model rains the air as it
is forced up the windward slope (:mod:`planet.orographic`). Then the wrung-out air descends the lee, and
an opt-in along-wind **moisture budget** (:mod:`planet.orographic_depletion`) draws the lee baseline
*below* the zonal mean — the real rain-shadow **desert**, not merely "no orographic bonus."

The payoff is the biome map itself: re-classify the patch under Phase 2's Whittaker rule with the
orographic precipitation, and a large fraction of cells change biome — the windward slope wettens toward
forest, the lee dries toward steppe/desert (:mod:`planet.orographic_scene`).

**At the honest altitude.** The flow is **prescribed**, not solved over the mountain: the linear theory
is *handed* the cross-mountain wind (the westerly at the patch latitude), it does not compute the flow
deflection. The absolute cm/yr amplitude rides ``OROGRAPHIC_HOURS_PER_YEAR`` — an effective annual
duration of active uplift, a loose-magnitude calibration knob in the spirit of :mod:`planet.precip`'s
band amplitudes. What is structurally solid is the **pattern** (windward-wet / lee-dry, the shadow in
the right place) and the two exact analytic anchors the engine tests (:mod:`planet.orographic`).

Run headless (saves the figure, prints the rain-shadow story):

    python -m planet.demo_orographic
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from . import orographic_scene as osc
from .orographic_scene import OrographicScene

RANGE_NAME = "cascades"        # the demo range (a N–S ridge under the westerlies — Cascades-like)

_REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_FIGURE = _REPO_ROOT / "docs" / "figures" / "planet-orographic.png"
OUTPUT_FIGURE = _REPO_ROOT / "outputs" / "planet-orographic.png"


def compute(range_name: str = RANGE_NAME, *, use_jet: bool = True) -> OrographicScene:
    """Build the rung-5A scene: place the ridge, source the wind, rain the windward slope, deplete the lee.

    ``use_jet`` reads the cross-mountain wind off the emergent coupled jet (the honest, but ~minutes,
    shallow-water spin-up); ``use_jet=False`` uses the Smith & Barstad reference westerly instead — a fast,
    deterministic fallback that gives essentially the same pattern (the jet and the reference westerly are
    the same ~15 m/s at this latitude). The lee moisture depletion (Rung 5A.3) is always on here — it is
    what turns the dry shadow into a real desert below the zonal-mean baseline.
    """
    return osc.demo_scene(range_name, use_jet=use_jet, deplete=True)


def print_summary(scene: OrographicScene) -> None:
    """Print the rain-shadow story — the payoff in text."""
    base = scene.baseline_precip_cm
    tot = scene.precip_cm
    lat, lon = scene.lat_deg, scene.lon_deg

    # The crest row (the wettest windward line) and its lee: read the windward peak and the lee trough.
    crest_row = int(np.argmax(scene.elevation_m.max(axis=1)))
    tot_row = tot[crest_row]
    windward_max = float(tot_row.max())
    lee_min = float(tot_row.min())
    base_mean = float(base.mean())

    print("\nOrographic rain shadow (rung 5A) — the mountain wakes the inert elevation seam\n")
    print(f"  range: '{RANGE_NAME}'  (patch {abs(lat.min()):.0f}–{abs(lat.max()):.0f}°N × "
          f"{abs(lon.max()):.0f}–{abs(lon.min()):.0f}°W, N–S ridge, crest {scene.elevation_m.max():.0f} m)")
    print(f"  cross-mountain wind: {scene.wind_speed:.1f} m/s from {scene.wind_direction_deg:.0f}° "
          f"(the westerly at {scene.lat_ref_deg:.0f}°N — prescribed, read off the jet)\n")

    print(f"  zonal-mean baseline precip:  {base_mean:5.1f} cm/yr  (the same at every longitude — rung 0)")
    print(f"  windward slope (forced ascent): peaks at {windward_max:5.1f} cm/yr  "
          f"(×{windward_max / base_mean:.1f} the baseline — the wet windward flank)")
    print(f"  lee (descent + upwind rainout): drops to {lee_min:5.1f} cm/yr  "
          f"(×{lee_min / base_mean:.2f} the baseline — BELOW it: the rain-shadow desert)\n")

    print(f"  the payoff — the mountain changes the biome map:")
    print(f"    biome reclassified on {100 * scene.biome_changed_fraction:4.0f}% of the patch "
          f"(windward → wetter, lee → drier than the zonal-mean map)")
    print(f"    lee-desert fraction {100 * scene.lee_desert_fraction:4.0f}% "
          f"(cells whose precip fell BELOW the zonal-mean baseline — the 5A.3 depletion at work,")
    print(f"      structurally impossible under enhancement-only)\n")


def save_figure(scene: OrographicScene) -> Path:
    """Render and save the orographic-scene artifact (needs the optional ``viz`` extra)."""
    import matplotlib
    matplotlib.use("Agg")                                    # headless
    from .plots import orographic_scene_figure

    fig = orographic_scene_figure(scene)
    for target in (DOCS_FIGURE, OUTPUT_FIGURE):
        target.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(target, dpi=130)
    return DOCS_FIGURE


def main() -> None:
    import sys
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")             # °, →, × on legacy codepages

    scene = compute()
    print_summary(scene)
    try:
        saved = save_figure(scene)
        print(f"Figure saved → {saved.relative_to(_REPO_ROOT)}")
    except ImportError:
        print("(matplotlib not installed — install the viz extra to render the figure: "
              "pip install -e .[viz])")


if __name__ == "__main__":
    main()
