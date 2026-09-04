"""Rung 5A.4: the mountain is *cold*, not only wet — elevation finally moves the thermometer.

Rung 5A gave the terrain its first real job: force the air up, rain the windward slope, dry the lee
(:mod:`planet.demo_orographic`). But the **temperature** underneath was still the zonal-mean climate
broadcast across the patch — the 2500 m crest and the valley beside it were handed the *same* number, so
the map could get wetter but never colder. This demo closes that: each cell is cooled by its own terrain
height (:mod:`planet.elevation_temperature`) before the Whittaker classifier runs, and the crest of the
Cascades-like ridge stops being forest and becomes **alpine tundra** — the vertical biome zonation of
Whittaker's own mountain transects, emerging on a planet whose climate was computed rather than drawn.

**The rung also carries a negative, and it is the more interesting half.** The repository's habit is to
retire a prescribed number by making it emergent (rung 4 retired the fixed ``B``). The obvious target here
is the 6.5 K/km lapse-rate constant, replaced by rung 4's own **moist adiabat** integrated up the column.
It does not survive its own benchmark. At the demo's mid-latitude range the emergent rate comes out
6.3 K/km — it *reproduces* the constant instead of retiring it. And in the tropics it is measurably worse:
the observed 0 °C isotherm sits near 4.5–5 km (Harris, Bowman & Shin 2000), the constant puts it at
~4.4 km — just below that band — and the moist adiabat puts it above 7 km, ~45 % high, because a
*saturated parcel* adiabat is not the *environmental* lapse rate of an unsaturated mean column. The
**ordering** is what decides it — one is close, the other is far — not either number landing inside the
band. So the constant stays the default, the emergent path ships opt-in as a diagnostic, and the failure
is printed here rather than hidden.

**At the honest altitude.** The correction is **diagnostic and one-way**: the cooled surface does not
re-solve the energy balance, does not grow snow (no albedo feedback), and does not re-run the orographic
model. And the Whittaker rule's cold bands ignore precipitation, so on the crest — wet *and* cold — the
rain shadow and the lapse rate act on the same cells but only the cooling can change the answer.

Run headless (saves the figure, prints the alpine story and the lapse-rate verdict):

    python -m planet.demo_alpine_biomes
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from . import biomes, elevation_temperature as elev, orographic_scene as osc
from .orographic_scene import OrographicScene

RANGE_NAME = "cascades"          # the demo range (a N–S ridge under the westerlies — Cascades-like)
RIDGE_M = 2500.0                 # m — the crest height of the demo ridge
SWEEP_LATITUDES = np.linspace(0.0, 89.0, 60)   # the latitude sweep the lapse-rate verdict panel plots

_REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_FIGURE = _REPO_ROOT / "docs" / "figures" / "planet-alpine-biomes.png"
OUTPUT_FIGURE = _REPO_ROOT / "outputs" / "planet-alpine-biomes.png"


def compute(range_name: str = RANGE_NAME, *, use_jet: bool = False
            ) -> tuple[OrographicScene, OrographicScene, elev.LapseDiagnostic]:
    """Build the alpine scene, its no-cooling control, and the lapse-rate diagnostic.

    Returns ``(scene, sea_level_scene, diagnostic)``: the same patch, the same wind and the same
    Rung-5A.3 rainfall, classified **with** and **without** the Rung-5A.4 terrain cooling — so the alpine
    cap in the figure is attributable to the lapse rate alone — plus the latitude sweep that sets the
    emergent moist adiabat against the pinned constant and the observed freezing level.

    ``use_jet`` reads the cross-mountain wind off the emergent coupled jet (honest, but a ~minutes
    shallow-water spin-up); the default uses the Smith & Barstad reference westerly, which is the same
    ~15 m/s at this latitude and keeps the demo fast. The rainfall is irrelevant to the cooling either way.
    """
    from . import demo_biomes

    scene = osc.demo_scene(range_name, use_jet=use_jet, deplete=True, lapse=True, amplitude_m=RIDGE_M)
    control = osc.demo_scene(range_name, use_jet=use_jet, deplete=True, lapse=False, amplitude_m=RIDGE_M)

    state = demo_biomes.compute().state
    T_sea = np.interp(SWEEP_LATITUDES, state.latitude_deg(), state.T)
    diagnostic = elev.lapse_diagnostic(SWEEP_LATITUDES, T_sea, elevation_m=RIDGE_M)
    return scene, control, diagnostic


def print_summary(scene: OrographicScene, control: OrographicScene,
                  diagnostic: elev.LapseDiagnostic) -> None:
    """Print the alpine story and the lapse-rate verdict — the payoff and the negative, in text."""
    mid = scene.temperature_C.shape[0] // 2
    crest = int(np.argmax(scene.elevation_m[mid, :]))
    T_sea = float(scene.sea_level_temperature_C[mid, crest])
    T_crest = float(scene.temperature_C[mid, crest])
    before = biomes.Biome(int(control.biome_codes[mid, crest]))
    after = biomes.Biome(int(scene.biome_codes[mid, crest]))

    print("\nThe mountain is COLD, not only wet (rung 5A.4) — elevation moves the thermometer\n")
    print(f"  range: '{RANGE_NAME}'  (N–S ridge, crest {scene.elevation_m.max():.0f} m, "
          f"patch centred {scene.lat_ref_deg:.0f}°)")
    print(f"  crest air:  sea-level zonal mean {T_sea:+6.2f} °C  →  terrain-cooled {T_crest:+6.2f} °C  "
          f"(a {T_sea - T_crest:.1f} K drop over {scene.elevation_m.max():.0f} m)\n")

    print("  the payoff — the crest changes biome because it got COLD, not because it got wet:")
    print(f"    crest biome:  {biomes.BIOME_NAMES[before.value]}  →  {biomes.BIOME_NAMES[after.value]}")
    print(f"    re-classified by the cooling alone: {100 * scene.alpine_fraction:4.0f}% of the patch")
    print(f"    re-classified vs the zonal-mean map overall: "
          f"{100 * control.biome_changed_fraction:3.0f}% (rain only) → "
          f"{100 * scene.biome_changed_fraction:3.0f}% (rain + cooling)")
    print("    (the two overlap: the Whittaker cold bands ignore precipitation, so on the cold crest the")
    print("     rain shadow can no longer change the answer — a named degeneracy, not a bug)\n")

    lat = diagnostic.latitude_deg
    i_mid = int(np.argmin(np.abs(lat - abs(scene.lat_ref_deg))))
    lo, hi = elev.OBSERVED_TROPICAL_FREEZING_LEVEL_M
    print("  the negative — making the lapse rate EMERGENT does not retire the pinned constant:")
    print(f"    at {lat[i_mid]:.0f}°: emergent moist adiabat {1e3 * diagnostic.gamma_moist[i_mid]:4.2f} K/km "
          f"vs the pinned {1e3 * diagnostic.gamma_constant[i_mid]:4.2f} K/km — it CONFIRMS it (~3%)")
    print(f"    in the deep tropics it is worse: freezing level "
          f"{diagnostic.freezing_constant_m[0] / 1e3:4.2f} km (constant) vs "
          f"{diagnostic.freezing_moist_m[0] / 1e3:4.2f} km (moist), observed {lo / 1e3:.1f}–{hi / 1e3:.1f} km")
    print("    a saturated PARCEL adiabat is not the ENVIRONMENTAL lapse rate of an unsaturated mean")
    print("    column — so the constant stays the default and the emergent path ships as a diagnostic.\n")


def save_figure(scene: OrographicScene, control: OrographicScene,
                diagnostic: elev.LapseDiagnostic) -> Path:
    """Render and save the alpine-biomes artifact (needs the optional ``viz`` extra)."""
    import matplotlib
    matplotlib.use("Agg")                                    # headless
    from .plots import alpine_biomes_figure

    fig = alpine_biomes_figure(scene, control, diagnostic)
    for target in (DOCS_FIGURE, OUTPUT_FIGURE):
        target.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(target, dpi=130)
    return DOCS_FIGURE


def main() -> None:
    import sys
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")             # °, →, × on legacy codepages

    scene, control, diagnostic = compute()
    print_summary(scene, control, diagnostic)
    try:
        saved = save_figure(scene, control, diagnostic)
        print(f"Figure saved → {saved.relative_to(_REPO_ROOT)}")
    except ImportError:
        print("(matplotlib not installed — install the viz extra to render the figure: "
              "pip install -e .[viz])")


if __name__ == "__main__":
    main()
