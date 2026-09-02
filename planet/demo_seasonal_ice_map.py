"""Rung 5B.3: the seasons, the continents, and the snow — the seasonal ice map, and a mask the annual mean can see.

Rung 5B.2 marched the seasons over a land–sea map with a fixed albedo and proved a clean theorem: the
**annual mean is blind to the mask** — average over the year and the heat capacity cancels, so the
continents leave no trace in the mean climate. Rung 5B.1+ put the ice-albedo feedback on the zonal tiles.
This demo puts it **on the map** (:mod:`planet.seasonal_map`, ``march(coalbedo_fn=…)``): every grid point
freezes on its own temperature, so the small-``C`` continental interiors grow a **winter snow cover** that
the sluggish ocean at the same latitude never does, and the polar ocean grows sea ice that lingers. The
payoffs are a **map of seasonal ice** — where it forms, and how much of the year it lasts — and a
**broken theorem**: snow reflects the winter sun the ocean keeps absorbing, so the continents end
*colder in the annual mean* than the ocean at their latitude. The mask is now visible in the mean, by
**rectification** of the seasonal cycle through the albedo step — an effect no linear model can have.

Run headless (saves the figure, the monthly GIF and the month-slider globe, prints the ice story):

    python -m planet.demo_seasonal_ice_map
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import seasonal as sea
from . import seasonal_map as sm
from .seasonal_map import SeasonalMapEBM, SeasonalMapClimate

BAND_LAT_DEG = 55.0            # the midlatitude band read across (interior snow, open ocean)
INTERIOR_LON_DEG = 90.0        # deep inside the broad NH continent
OCEAN_LON_DEG = 200.0          # open ocean at the same latitude
N_CELLS, N_LON, N_STEPS = 45, 90, 360       # coarse enough to march in under a minute
FPS, SAVE_DPI = 6, 90

_REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_FIGURE = _REPO_ROOT / "docs" / "figures" / "planet-seasonal-ice-map.png"
OUTPUT_FIGURE = _REPO_ROOT / "outputs" / "planet-seasonal-ice-map.png"
DOCS_ANIMATION = _REPO_ROOT / "docs" / "figures" / "planet-seasonal-ice-map.gif"
OUTPUT_ANIMATION = _REPO_ROOT / "outputs" / "planet-seasonal-ice-map.gif"
DOCS_GLOBE = _REPO_ROOT / "docs" / "figures" / "planet-seasonal-ice-globe.html"


@dataclass(frozen=True)
class SeasonalIceMapResult:
    """The banked 2-D seasonal-ice artifact — the plain bundle the figure, animation and summary consume."""

    model: SeasonalMapEBM
    climate: SeasonalMapClimate          # the ice-albedo march (warm seed → the finite-ice branch)
    band_lat_deg: float
    interior_lon_deg: float
    ocean_lon_deg: float

    def band_index(self) -> int:
        return self.model.nearest_index(self.band_lat_deg)

    def sample_columns(self) -> tuple[int, int]:
        """(interior, ocean) longitude indices at the band."""
        return (self.model.nearest_lon_index(self.interior_lon_deg),
                self.model.nearest_lon_index(self.ocean_lon_deg))

    def month_steps(self) -> np.ndarray:
        """Time-sample index of the middle of each of the 12 months."""
        return ((np.arange(12) + 0.5) * self.climate.days.size / 12.0).astype(int)


def compute(n_cells: int = N_CELLS, n_lon: int = N_LON, n_steps: int = N_STEPS,
            band_lat_deg: float = BAND_LAT_DEG) -> SeasonalIceMapResult:
    """Build the 2-D model on the idealized-Earth mask, march the ice feedback to its limit cycle, bundle."""
    m = SeasonalMapEBM(n_cells=n_cells, n_lon=n_lon, n_steps=n_steps)
    mask = sm.earthlike_mask(m.x, m.lon)
    m = SeasonalMapEBM(land_mask=mask, n_cells=n_cells, n_lon=n_lon, n_steps=n_steps)
    c = m.march(coalbedo_fn=sea.ice_coalbedo, T_init=15.0, tol=1e-6, max_years=150)
    return SeasonalIceMapResult(m, c, band_lat_deg, INTERIOR_LON_DEG, OCEAN_LON_DEG)


def print_summary(r: SeasonalIceMapResult) -> None:
    """Print the seasonal-ice-map story — the payoff in text."""
    m, c = r.model, r.climate
    i = r.band_index()
    interior, ocean = r.sample_columns()
    lat, lon = c.latitude_deg(), c.longitude_deg()
    frac = c.ice_fraction()
    anom = c.zonal_anomaly()
    land = c.land_mask

    print("\n2-D seasonal EBM + ice-albedo → the seasonal ice MAP, and a mask the annual mean can see (rung 5B.3)\n")
    print(f"  grid {m.n_cells} lat × {m.n_lon} lon, land fraction {land.mean():.2f}; converged in {c.years} yr "
          f"(global annual mean {c.annual_mean().mean():.1f} °C)\n")

    print(f"  at {abs(lat[i]):.0f}°{'N' if lat[i] >= 0 else 'S'} — the same latitude, two surfaces:")
    print(f"    continental interior (λ={lon[interior]:3.0f}°): frozen {100 * frac[i, interior]:3.0f}% of the year"
          f"   — winter snow, melts every summer")
    print(f"    open ocean           (λ={lon[ocean]:3.0f}°): frozen {100 * frac[i, ocean]:3.0f}% of the year"
          f"   — stays open\n")

    # seasonal-ice extent by surface: the equatorward-most latitude that freezes at some point in the year
    nh = lat > 0.0
    def reach(surface_mask):
        cells = (frac > 0.0) & surface_mask & nh[:, None]
        rows = np.where(cells.any(axis=1))[0]
        return float(lat[rows.min()]) if rows.size else 90.0
    print(f"  seasonal-ice reach (NH): over land down to {reach(land):.0f}°, over ocean down to {reach(~land):.0f}°")
    peren = (c.T.max(axis=2) < sea.T_FREEZE)
    print(f"  perennial ice: {100 * peren[land].mean():.0f}% of land cells, {100 * peren[~land].mean():.0f}% of ocean cells"
          f" — land ice is seasonal, polar sea ice lingers\n")

    print("  the annual mean now SEES the mask (the 5B.2 theorem breaks — by design):")
    print(f"    interior ⟨T⟩ anomaly vs its zonal mean {anom[i, interior]:+.1f} K, open ocean {anom[i, ocean]:+.1f} K")
    print(f"    (max east–west spread of the annual-mean map {float(np.max(anom.max(axis=1) - anom.min(axis=1))):.1f} K; "
          f"exactly 0 with a fixed albedo)")
    print("    → winter snow reflects sun the ocean keeps absorbing: the seasonal cycle is RECTIFIED through the")
    print("      albedo step, and the continents end colder in the annual mean. A nonlinear effect — no linear")
    print("      (fixed-albedo) model can show it, whatever its heat capacities.\n")


def save_figure(r: SeasonalIceMapResult) -> Path:
    """Render and save the seasonal-ice-map still (needs the optional ``viz`` extra)."""
    import matplotlib
    matplotlib.use("Agg")                                    # headless
    from .plots import seasonal_ice_map_figure

    fig = seasonal_ice_map_figure(r)
    for target in (DOCS_FIGURE, OUTPUT_FIGURE):
        target.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(target, dpi=130)
    return DOCS_FIGURE


def save_animation(r: SeasonalIceMapResult, fps: int = FPS, dpi: int = SAVE_DPI) -> Path:
    """Render and save the month-by-month animation as a GIF (Pillow ships with matplotlib)."""
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.animation import PillowWriter
    from .plots import seasonal_ice_map_animation

    anim = seasonal_ice_map_animation(r)
    for target in (DOCS_ANIMATION, OUTPUT_ANIMATION):
        target.parent.mkdir(parents=True, exist_ok=True)
        anim.save(target, writer=PillowWriter(fps=fps), dpi=dpi)
    return DOCS_ANIMATION


def save_globe(r: SeasonalIceMapResult) -> Path:
    """Render and save the month-slider globe (needs the optional ``webviz`` extra — Plotly)."""
    from .seasonal_globe import save_seasonal_ice_globe
    return save_seasonal_ice_globe(r.climate, DOCS_GLOBE)


def main() -> None:
    import sys
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")             # °C, →, ⟨⟩ on legacy codepages

    r = compute()
    print_summary(r)
    try:
        saved = save_figure(r)
        print(f"Figure saved → {saved.relative_to(_REPO_ROOT)}")
        saved = save_animation(r)
        print(f"Animation saved → {saved.relative_to(_REPO_ROOT)}")
    except ImportError:
        print("(matplotlib not installed — install the viz extra to render the figure: "
              "pip install -e .[viz])")
    try:
        saved = save_globe(r)
        print(f"Interactive globe saved → {saved.relative_to(_REPO_ROOT)}")
    except ImportError:
        print("(plotly not installed — install the webviz extra to bank the interactive globe: "
              "pip install -e .[webviz])")


if __name__ == "__main__":
    main()
