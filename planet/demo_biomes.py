"""The Phase-2 banked artifact: climate → habitability — the biome map (the dramatic end-to-end win).

The capstone's payoff, banked **early** (plan §3): the planet's climate state — the EBM temperature
``T(φ)`` (:mod:`projects.planet.ebm`) and the diagnostic precipitation ``P(φ)``
(:mod:`projects.planet.precip`) — is mapped through the **Whittaker classifier**
(:mod:`projects.planet.biomes`) into a **map of biomes**: tropical rain forest at the equator,
savanna and the great deserts in the subtropics, temperate forest in the midlatitudes, then boreal
forest and tundra toward the poles. *Planetary knobs in, bands of life out* — the planet analogue of
Steel's microstructure and Chip's device.

And the bands **migrate as the knobs turn**: this demo also warms the planet (a CO₂ increase = a lower
OLR offset ``A``) and shows the biomes shift **poleward** — the tropics expand, the ice retreats — the
"knob in, habitability out" story the interactive map (a later phase) makes live.

Run headless (saves the figure, prints the bands):

    python -m projects.planet.demo_biomes
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np

from . import biomes, precip
from .albedo import EBMParams, present_day_climate
from .biomes import Biome, BIOME_NAMES
from .ebm import A_OLR, ClimateState

PRESENT_N_TAU = 0.01           # finer relaxation step for the headline present-day equilibrium
CO2_WARMING_DELTA_A = 8.0      # W m⁻² — drop in the OLR offset A standing in for a CO₂ increase. A CO₂
                               #   doubling is ~4 W m⁻² of forcing, so ~8 ≈ a strong (multi-doubling /
                               #   high-emissions) scenario → a visible-but-credible ~+4–5 °C warming.

_REPO_ROOT = Path(__file__).resolve().parents[2]
DOCS_FIGURE = _REPO_ROOT / "docs" / "figures" / "planet-biomes.png"
OUTPUT_FIGURE = _REPO_ROOT / "outputs" / "planet-biomes.png"


@dataclass(frozen=True)
class BiomeResult:
    """The banked biome artifact — the plain bundle the figure and the summary consume."""

    params: EBMParams
    state: ClimateState
    precip_cm: np.ndarray
    codes: np.ndarray

    @property
    def area_fractions(self) -> dict[Biome, float]:
        return biomes.biome_area_fractions(self.codes)

    def fraction(self, biome: Biome) -> float:
        """The planetary area fraction of one biome (0 if absent)."""
        return self.area_fractions.get(biome, 0.0)

    def bands(self, latitudes=(0, 15, 30, 45, 60, 75, 90)):
        """Sample (latitude, T, P, biome) at representative latitudes — the equator→pole band story."""
        lat = self.state.latitude_deg()
        out = []
        for target in latitudes:
            i = int(np.argmin(np.abs(lat - target)))
            out.append((lat[i], float(self.state.T[i]), float(self.precip_cm[i]), Biome(int(self.codes[i]))))
        return out


def _classify(state: ClimateState) -> tuple[np.ndarray, np.ndarray]:
    """A climate state → (precip field cm/yr, biome codes) — the Phase-2 map composition."""
    precip_cm = precip.precip_field(state)
    return precip_cm, biomes.classify_field(state.T, precip_cm)


def compute(params: EBMParams | None = None) -> BiomeResult:
    """Present-day climate → precip → biomes → :class:`BiomeResult` (no plotting)."""
    if params is None:
        params = EBMParams()
    state = present_day_climate(params, n_tau=PRESENT_N_TAU)
    precip_cm, codes = _classify(state)
    return BiomeResult(params=params, state=state, precip_cm=precip_cm, codes=codes)


def warmed(params: EBMParams | None = None, delta_A: float = CO2_WARMING_DELTA_A) -> BiomeResult:
    """A warmer planet (CO₂ up ≈ lower OLR offset ``A``) — to show the bands migrate poleward."""
    if params is None:
        params = EBMParams()
    warm_params = replace(params, A=params.A - delta_A)
    state = present_day_climate(warm_params, n_tau=PRESENT_N_TAU)
    precip_cm, codes = _classify(state)
    return BiomeResult(params=warm_params, state=state, precip_cm=precip_cm, codes=codes)


def print_summary(r: BiomeResult, warm: BiomeResult) -> None:
    """Print the biome-band story and the poleward migration under warming — the demo's payoff in text."""
    print("\nClimate → habitability: the biome map (Whittaker classifier on the EBM climate)\n")
    print(f"  present-day (S₀ = {r.params.S0:.1f} W/m², global mean T̄ = {r.state.global_mean_T:.2f} °C):")
    print("    latitude   T (°C)   P (cm/yr)   biome")
    for lat, T, P, b in r.bands():
        print(f"      {lat:5.1f}°   {T:6.2f}    {P:6.1f}     {BIOME_NAMES[b]}")
    rainforest = r.fraction(Biome.TROPICAL_RAIN_FOREST)
    tundra = r.fraction(Biome.TUNDRA)
    print(f"\n  area fractions sum to {sum(r.area_fractions.values()):.6f} (the map tiles the planet — "
          f"the consistency check)")
    print(f"  tropical rain forest covers {100*rainforest:.1f}% of the planet; tundra {100*tundra:.1f}%")
    print()
    dA = r.params.A - warm.params.A
    print(f"  Turn a knob — warm the planet (CO₂ up ≈ OLR offset A: {r.params.A:.0f} → {warm.params.A:.0f} W/m², "
          f"−{dA:.0f}):")
    print(f"    global mean T̄: {r.state.global_mean_T:.2f} → {warm.state.global_mean_T:.2f} °C   "
          f"(ice line {r.state.ice_line_lat:.1f}° → {warm.state.ice_line_lat:.1f}°)")
    print(f"    tropical rain forest: {100*rainforest:.1f}% → {100*warm.fraction(Biome.TROPICAL_RAIN_FOREST):.1f}%   "
          f"tundra: {100*tundra:.1f}% → {100*warm.fraction(Biome.TUNDRA):.1f}%")
    print("    → the bands migrate poleward: the tropics expand, the ice (and tundra) retreat.\n")


def save_figure(r: BiomeResult) -> Path:
    """Render and save the present-day biome-map artifact (needs the optional ``viz`` extra)."""
    import matplotlib
    matplotlib.use("Agg")                            # headless
    from .plots import biomes_figure

    fig = biomes_figure(r.state, r.precip_cm, r.codes)
    for target in (DOCS_FIGURE, OUTPUT_FIGURE):
        target.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(target, dpi=130)
    return DOCS_FIGURE


def main() -> None:
    import sys
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")     # °C, ₂, →, ° on legacy codepages

    r = compute()
    warm = warmed()
    print_summary(r, warm)
    try:
        saved = save_figure(r)
        print(f"Figure saved → {saved.relative_to(_REPO_ROOT)}")
    except ImportError:
        print("(matplotlib not installed — install the viz extra to render the figure: "
              "pip install -e .[viz])")


if __name__ == "__main__":
    main()
