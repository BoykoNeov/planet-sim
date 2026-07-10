"""Rung 5B.1: the seasonal cycle wakes heat capacity, and continentality appears.

Every EBM up the staircase solved for an **equilibrium**, where the heat capacity ``C`` cancels — so a
land column and an ocean column at the same latitude ended up **identical**, and continentality (the big
seasonal range of continental interiors) was *exactly zero*. This demo turns on the **seasons**: the
same diffusive transport + linear radiation, but marched forward under axial-tilt insolation ``S(x, t)``
to a converged **annual limit cycle**. Now ``C`` is load-bearing, and a small-``C`` **land** tile beside
a large-``C`` **ocean** tile at the same latitude tells the whole story — the land swings hard and
nearly in step with the sun; the ocean barely moves and lags ~2 months. *Continentality, from the heat-
capacity contrast alone* (:mod:`planet.seasonal`).

The result is checked two ways: a **time-marcher** (the engine-reuse method — Strang-split stepping to a
limit cycle) and an exact **frequency-domain** solve (the tight reference, whose ``n=0`` harmonic *is* the
old annual-mean EBM). They agree; the marcher's backward-Euler transport is not damping the swing.

Run headless (saves the figure, prints the continentality story):

    python -m planet.demo_seasonal
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import seasonal as sea
from .seasonal import SeasonalEBM, SeasonalClimate

BAND_LAT_DEG = 45.0            # the midlatitude band the cycle panel reads (a strong land/ocean contrast)

_REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_FIGURE = _REPO_ROOT / "docs" / "figures" / "planet-seasonal.png"
OUTPUT_FIGURE = _REPO_ROOT / "outputs" / "planet-seasonal.png"


@dataclass(frozen=True)
class SeasonalResult:
    """The banked seasonal artifact — the plain bundle the figure and summary consume."""

    model: SeasonalEBM
    climate: SeasonalClimate           # the exact spectral limit cycle (the headline fields)
    marched: SeasonalClimate           # the time-marcher's limit cycle (the cross-check)
    band_lat_deg: float

    def band_index(self) -> int:
        return self.model.nearest_index(self.band_lat_deg)


def compute(land_fraction: float = sea.LAND_FRACTION, band_lat_deg: float = BAND_LAT_DEG,
            n_cells: int = 180, n_steps: int = 360) -> SeasonalResult:
    """Build the seasonal model, solve the limit cycle both ways, bundle the result (no plotting)."""
    m = SeasonalEBM(land_fraction=land_fraction, n_cells=n_cells, n_steps=n_steps)
    spectral = m.spectral()                                  # the exact reference (the headline fields)
    marched = m.march(tol=1e-6, max_years=80)                # the engine-reuse cross-check
    return SeasonalResult(m, spectral, marched, band_lat_deg)


def print_summary(r: SeasonalResult) -> None:
    """Print the continentality story — the payoff in text."""
    m, c = r.model, r.climate
    i = r.band_index()
    lat = c.latitude_deg()
    print("\nSeasonal EBM → continentality (land vs ocean heat capacity)\n")
    print(f"  heat capacities: land C_L = {m.C_land:.3e}  ocean C_O = {m.C_ocean:.3e} J/m²/K "
          f"(ocean/land = {m.C_ocean / m.C_land:.0f}×)")
    print(f"  land fraction f_L = {m.f_land:.2f}\n")

    print(f"  at {abs(lat[i]):.0f}°{'N' if lat[i] >= 0 else 'S'} — the same latitude, two surfaces:")
    ampL, ampO = c.amplitude("land")[i], c.amplitude("ocean")[i]
    lagL, lagO = m.phase_lag_days(c.T_land)[i], m.phase_lag_days(c.T_ocean)[i]
    print(f"    land  tile: seasonal amplitude {ampL:5.1f} K (range {2*ampL:4.0f} K), lags the sun {lagL:4.1f} days")
    print(f"    ocean tile: seasonal amplitude {ampO:5.1f} K (range {2*ampO:4.0f} K), lags the sun {lagO:4.1f} days")
    print(f"    → the land tile swings {ampL/ampO:.0f}× harder, the ocean lags {lagO-lagL:.0f} days longer: "
          f"continentality, from C alone.\n")

    # The core insight, verified: the annual MEANS are identical — continentality is purely seasonal.
    mean_gap = float(np.max(np.abs(c.annual_mean("land") - c.annual_mean("ocean"))))
    print(f"  the annual means of land and ocean are identical (max gap {mean_gap:.1e} K): at the annual")
    print("    mean C cancels — continentality lives ENTIRELY in the seasonal amplitude, zero in the mean.\n")

    # The two solvers agree (the anti-damping cross-check).
    gap = max(np.max(np.abs(r.marched.T_land - c.T_land)), np.max(np.abs(r.marched.T_ocean - c.T_ocean)))
    print(f"  time-marcher vs exact spectral solve: max |ΔT| = {gap:.2e} K "
          f"(marcher converged in {r.marched.years} yr) — backward-Euler transport is not damping the swing.\n")


def save_figure(r: SeasonalResult) -> Path:
    """Render and save the seasonal-cycle / continentality artifact (needs the optional ``viz`` extra)."""
    import matplotlib
    matplotlib.use("Agg")                                    # headless
    from .plots import seasonal_figure

    fig = seasonal_figure(r)
    for target in (DOCS_FIGURE, OUTPUT_FIGURE):
        target.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(target, dpi=130)
    return DOCS_FIGURE


def main() -> None:
    import sys
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")             # °C, ₂, →, × on legacy codepages

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
