"""Rung 5B.1+: the seasonal ice-albedo feedback — a migrating ice edge, and continentality written in ice.

Rung 5B.1 woke the heat capacity and cast continentality with a **fixed** albedo. This demo lifts that
scope edge: it hands the marcher the rung-0 **step-function ice-albedo** (:func:`planet.seasonal.ice_coalbedo`,
the *same* nonlinearity Phase 1's Snowball rides), re-frozen on each tile's own temperature every
half-step. Now the model does something no equilibrium EBM below it can: it grows a **seasonal ice edge**
that forms in winter and melts in summer. And because the small-``C`` **land** tile plunges below freezing
over a wide winter band while the sluggish **ocean** tile barely does, the land grows a wide seasonal ice
zone the ocean lacks — **continentality, now written in ice** (:mod:`planet.seasonal`).

The feedback is *marcher-only* — a state-dependent albedo has no frequency-domain closed form, so the
spectral solver cannot carry it — and it restores Phase 1's **bistability**: a warm seed settles on a
finite-ice climate, a cold seed on a frozen snowball, at one and the same sun.

Run headless (saves the figure, prints the ice story):

    python -m planet.demo_seasonal_ice
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import seasonal as sea
from .seasonal import SeasonalEBM, SeasonalClimate

BAND_LAT_DEG = 45.0            # the midlatitude band the cycle panel reads (land freezes, ocean stays open)

_REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_FIGURE = _REPO_ROOT / "docs" / "figures" / "planet-seasonal-ice.png"
OUTPUT_FIGURE = _REPO_ROOT / "outputs" / "planet-seasonal-ice.png"


@dataclass(frozen=True)
class SeasonalIceResult:
    """The banked seasonal-ice artifact — the plain bundle the figure and summary consume."""

    model: SeasonalEBM
    warm: SeasonalClimate           # the finite-ice branch (warm seed) — the headline fields
    snowball: SeasonalClimate       # the frozen branch (cold seed) — the bistability partner
    band_lat_deg: float

    def band_index(self) -> int:
        return self.model.nearest_index(self.band_lat_deg)


def compute(land_fraction: float = sea.LAND_FRACTION, band_lat_deg: float = BAND_LAT_DEG,
            n_cells: int = 180, n_steps: int = 360) -> SeasonalIceResult:
    """Build the seasonal model, march both ice branches to their limit cycles, bundle (no plotting)."""
    m = SeasonalEBM(land_fraction=land_fraction, n_cells=n_cells, n_steps=n_steps)
    warm = m.march(coalbedo_fn=sea.ice_coalbedo, T_init=20.0, tol=1e-6, max_years=300)
    snowball = m.march(coalbedo_fn=sea.ice_coalbedo, T_init=-40.0, tol=1e-6, max_years=300)
    return SeasonalIceResult(m, warm, snowball, band_lat_deg)


def print_summary(r: SeasonalIceResult) -> None:
    """Print the ice story — the payoff in text."""
    m, c = r.model, r.warm
    i = r.band_index()
    lat = c.latitude_deg()
    print("\nSeasonal ice-albedo → a migrating ice edge, and continentality written in ice\n")

    edge_land = sea.ice_edge_latitude(m.x, c.T_land)
    edge_ocean = sea.ice_edge_latitude(m.x, c.T_ocean)
    print(f"  seasonal ice edge (equatorward-most winter ice):")
    print(f"    land  tile → {edge_land:4.1f}°   ocean tile → {edge_ocean:4.1f}°   "
          f"(land ice reaches {edge_ocean - edge_land:.0f}° further equatorward)\n")

    fL, fO = c.ice_fraction("land"), c.ice_fraction("ocean")
    print(f"  at {abs(lat[i]):.0f}°{'N' if lat[i] >= 0 else 'S'} — the same latitude, two surfaces:")
    print(f"    land  tile: frozen {100 * fL[i]:4.0f}% of the year")
    print(f"    ocean tile: frozen {100 * fO[i]:4.0f}% of the year")
    print(f"    → the land freezes a good part of the year where the ocean stays open — continentality,")
    print("      as an ice asymmetry.\n")

    land_perennial = sea.ice_edge_latitude(m.x, c.T_land, kind="perennial")
    ocean_perennial = sea.ice_edge_latitude(m.x, c.T_ocean, kind="perennial")
    print(f"  and heat capacity decides whether the ice MELTS:")
    print(f"    land  ice is purely seasonal — the tiny-C land climbs above freezing every summer "
          f"(perennial edge {land_perennial:.0f}°: ~no year-round land ice)")
    print(f"    ocean ice, once formed, ~persists — the huge-C ocean barely warms in summer "
          f"(perennial edge {ocean_perennial:.0f}° ≈ its winter edge {edge_ocean:.0f}°)\n")

    # Bistability — the ice feedback carries Phase 1's Snowball into the seasonal cycle.
    Tw = float(r.warm.annual_mean("mean").mean())
    Ts = float(r.snowball.annual_mean("mean").mean())
    edge_snow = sea.ice_edge_latitude(m.x, r.snowball.T_mean)
    print(f"  bistability (one sun, two climates — the seed picks the branch):")
    print(f"    warm start  → global mean {Tw:6.1f} °C  (finite ice, open midlatitudes)")
    print(f"    cold start  → global mean {Ts:6.1f} °C  (snowball — ice edge {edge_snow:.0f}°)\n")

    print(f"  (marcher converged: warm {r.warm.years} yr, snowball {r.snowball.years} yr; "
          f"the spectral solver cannot carry this — the albedo is state-dependent.)\n")


def save_figure(r: SeasonalIceResult) -> Path:
    """Render and save the seasonal-ice artifact (needs the optional ``viz`` extra)."""
    import matplotlib
    matplotlib.use("Agg")                                    # headless
    from .plots import seasonal_ice_figure

    fig = seasonal_ice_figure(r)
    for target in (DOCS_FIGURE, OUTPUT_FIGURE):
        target.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(target, dpi=130)
    return DOCS_FIGURE


def main() -> None:
    import sys
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")             # °C, →, ≈, × on legacy codepages

    r = compute()
    print_summary(r)
    try:
        saved = save_figure(r)
        print(f"Figure saved → {saved.relative_to(_REPO_ROOT)}")
    except ImportError:
        print("(matplotlib not installed — install the viz extra to render the figure: "
              "pip install -e .[viz])")


if __name__ == "__main__":
    main()
