"""Rung 5B.2: the seasons meet the continents — continentality becomes a *map*.

Rung 5B.1 woke the heat capacity and cast continentality, but as a zonal-mean caricature: one land tile
and one ocean tile per latitude, a single number for the whole latitude circle. This demo resolves the
**longitude** axis — a single temperature field ``T(φ, λ, t)`` over a real land mask, marched under the
axial-tilt seasons to a converged annual limit cycle (:mod:`planet.seasonal_map`, the North–Mengel–Short
1983 model). Now the heat capacity varies point-by-point, the transport diffuses heat between neighbouring
longitudes, and the payoff is a **map**: continental **interiors** reach Siberian-style seasonal extremes
while their **coasts** are moderated by the neighbouring ocean and the open ocean barely moves at all.

Two things the map makes visible that the zonal mean could not:
* **Continentality within a latitude.** At one midlatitude the seasonal range runs from a maritime few
  kelvin over the ocean up to tens of kelvin deep in a continent — the interior/coast/ocean gradient.
* **The annual mean is blind to the mask.** Average over the year and the heat capacity cancels: the
  annual-mean temperature map is *zonally flat*, the same aquaplanet climate at every longitude. The whole
  land/sea signal lives in the **seasonal amplitude**, nowhere in the mean (5B.1's ``⟨T_L⟩ = ⟨T_O⟩``, now a
  map). This is the NMS83 headline.

Run headless (saves the figure, prints the continentality-map story):

    python -m planet.demo_seasonal_map
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import seasonal_map as sm
from .seasonal_map import SeasonalMapEBM, SeasonalMapClimate

BAND_LAT_DEG = 45.0            # the midlatitude band the cross-section reads
INTERIOR_LON_DEG = 90.0       # deep inside the broad NH continent (Eurasia-like)
COAST_LON_DEG = 40.0          # just inside its west coast
OCEAN_LON_DEG = 250.0         # open ocean

_REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_FIGURE = _REPO_ROOT / "docs" / "figures" / "planet-seasonal-map.png"
OUTPUT_FIGURE = _REPO_ROOT / "outputs" / "planet-seasonal-map.png"


@dataclass(frozen=True)
class SeasonalMapResult:
    """The banked 2-D seasonal artifact — the plain bundle the figure and summary consume."""

    model: SeasonalMapEBM
    climate: SeasonalMapClimate
    band_lat_deg: float
    interior_lon_deg: float
    coast_lon_deg: float
    ocean_lon_deg: float

    def band_index(self) -> int:
        return self.model.nearest_index(self.band_lat_deg)

    def sample_columns(self) -> tuple[int, int, int]:
        """(interior, coast, ocean) longitude indices for the cross-section / cycle panels."""
        return (self.model.nearest_lon_index(self.interior_lon_deg),
                self.model.nearest_lon_index(self.coast_lon_deg),
                self.model.nearest_lon_index(self.ocean_lon_deg))


def compute(n_cells: int = 90, n_lon: int = 180, n_steps: int = 360,
            band_lat_deg: float = BAND_LAT_DEG) -> SeasonalMapResult:
    """Build the 2-D model on a coarse idealized-Earth mask, march to the limit cycle, bundle (no plotting)."""
    m = SeasonalMapEBM(n_cells=n_cells, n_lon=n_lon, n_steps=n_steps)
    mask = sm.earthlike_mask(m.x, m.lon)
    m = SeasonalMapEBM(land_mask=mask, n_cells=n_cells, n_lon=n_lon, n_steps=n_steps)
    c = m.march(tol=1e-6, max_years=60)
    return SeasonalMapResult(m, c, band_lat_deg, INTERIOR_LON_DEG, COAST_LON_DEG, OCEAN_LON_DEG)


def print_summary(r: SeasonalMapResult) -> None:
    """Print the continentality-map story — the payoff in text."""
    m, c = r.model, r.climate
    i = r.band_index()
    lat = c.latitude_deg()
    interior, coast, ocean = r.sample_columns()
    rng = c.seasonal_range()
    lon = c.longitude_deg()

    print("\n2-D seasonal EBM → continentality as a MAP (seasons × continents, NMS83)\n")
    print(f"  grid: {m.n_cells} lat × {m.n_lon} lon; land fraction {c.land_mask.mean():.2f}; "
          f"C_ocean/C_land = {m.C_ocean / m.C_land:.0f}×  (converged in {c.years} yr)\n")

    print(f"  at {abs(lat[i]):.0f}°{'N' if lat[i] >= 0 else 'S'}, three points at the SAME latitude:")
    for name, k in (("continental interior", interior), ("coast", coast), ("open ocean", ocean)):
        surf = "land " if c.land_mask[i, k] else "ocean"
        print(f"    {name:22s} (λ={lon[k]:5.0f}°, {surf}): seasonal range {rng[i, k]:5.1f} K")
    print(f"    → the interior swings {rng[i, interior] / rng[i, ocean]:.0f}× the ocean; the coast sits "
          f"between (the adjacent sea moderates it): continentality is a MAP, not one number.\n")

    # The NMS headline, verified: the annual-mean map is zonally flat — the mask is invisible in the mean.
    amean = c.annual_mean()
    zonal_spread = float(np.max(amean.max(axis=1) - amean.min(axis=1)))
    print(f"  the annual-mean map is zonally flat (max east–west spread {zonal_spread:.2f} K): average over")
    print("    the year and C cancels — the land/sea contrast lives ENTIRELY in the seasonal amplitude,")
    print("    nothing in the mean. The continentality map overlays an aquaplanet annual-mean climate.\n")


def save_figure(r: SeasonalMapResult) -> Path:
    """Render and save the continentality-map artifact (needs the optional ``viz`` extra)."""
    import matplotlib
    matplotlib.use("Agg")                                    # headless
    from .plots import seasonal_map_figure

    fig = seasonal_map_figure(r)
    for target in (DOCS_FIGURE, OUTPUT_FIGURE):
        target.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(target, dpi=130)
    return DOCS_FIGURE


def main() -> None:
    import sys
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")             # °C, →, × on legacy codepages

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
