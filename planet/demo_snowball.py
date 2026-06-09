"""The Phase-1 banked artifact: Snowball-Earth hysteresis — one knob, two climates, a jump.

The planet's first complete, demonstrable artifact and the counterpart of Steel's four-curves
and Chip's process→device flow. It drives the latitudinal energy-balance model
(:mod:`projects.planet.ebm`) with the ice-albedo feedback (:mod:`projects.planet.albedo`) and
traces the **Snowball-Earth hysteresis** by a parameter-continuation sweep in the solar constant:

  1. **Present-day climate** — the temperate planet with a polar ice cap (ice line ~70°). This is the
     **finite-cap branch** (Earth's), and at today's sun it is *distinct from and slightly colder
     than* the near-ice-free branch the dimming sweep traces — today's S₀ already admits multiple
     stable climates (the bistability seen before the sun is even touched).
  2. **Dim the sun** (sweep S₀ down, warm-started) — the ice line creeps equatorward, then past a
     threshold the ice-albedo feedback runs away and the planet **freezes over in a catastrophic
     jump** (Snowball Earth).
  3. **Brighten it back** (sweep S₀ up) — the white planet reflects so much that it **stays frozen
     far past where it froze**: a wide hysteresis loop, two stable climates for one sun.

*Knob in, frozen-or-temperate planet out* — the dramatic, counter-intuitive bifurcation that is
the cheapest, most teachable payoff in the capstone, and simultaneously the integration test of
every Phase-1 module (transport + the Strang-split radiation + the feedback).

Run headless (saves the figure, prints the loop):

    python -m projects.planet.demo_snowball
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from . import ebm
from .albedo import EBMParams, HysteresisLoop, present_day_climate, snowball_hysteresis
from .ebm import ClimateState

# --- The sweep recipe ------------------------------------------------------- #
SWEEP_S0_MIN = 1000.0          # W m⁻² — dim enough to be deep in the Snowball trap
SWEEP_S0_MAX = 1900.0          # W m⁻² — bright enough to re-melt the frozen planet
SWEEP_STEPS = 60              # continuation points per branch
SWEEP_N_TAU = 0.05           # relaxation step (× τ_rad) — gentle, tracks the ice line smoothly
PRESENT_N_TAU = 0.01         # finer step for the headline present-day equilibrium

_REPO_ROOT = Path(__file__).resolve().parents[2]
DOCS_FIGURE = _REPO_ROOT / "docs" / "figures" / "planet-snowball.png"
OUTPUT_FIGURE = _REPO_ROOT / "outputs" / "planet-snowball.png"


@dataclass(frozen=True)
class SnowballResult:
    """The banked Snowball artifact — the plain bundle the figure and the summary consume."""

    params: EBMParams
    present: ClimateState
    loop: HysteresisLoop

    @property
    def present_ice_line(self) -> float:
        return self.present.ice_line_lat

    @property
    def freeze_S0(self) -> float:
        return self.loop.freeze_S0

    @property
    def melt_S0(self) -> float:
        return self.loop.melt_S0

    @property
    def hysteresis_width(self) -> float:
        return self.loop.hysteresis_width

    @property
    def freeze_dimming_pct(self) -> float:
        """The dimming (% below present S₀) at which the planet snowball-jumps."""
        return 100.0 * (self.params.S0 - self.freeze_S0) / self.params.S0

    @property
    def snowball_Tbar(self) -> float:
        """The global-mean temperature of the frozen (Snowball) state — the cold branch's floor."""
        return float(self.loop.Tbar_down.min())


def compute(params: EBMParams | None = None) -> SnowballResult:
    """Run the present-day equilibrium + the full hysteresis sweep → :class:`SnowballResult` (no plotting)."""
    if params is None:
        params = EBMParams()
    present = present_day_climate(params, n_tau=PRESENT_N_TAU)
    loop = snowball_hysteresis(
        params, S0_min=SWEEP_S0_MIN, S0_max=SWEEP_S0_MAX, n_steps=SWEEP_STEPS, n_tau=SWEEP_N_TAU,
    )
    return SnowballResult(params=params, present=present, loop=loop)


def print_summary(r: SnowballResult) -> None:
    """Print the Snowball story — the demo's payoff in text."""
    print("\nSnowball Earth: one knob (the solar constant), two climates, a catastrophic jump\n")
    print(f"  present-day (S₀ = {r.params.S0:.1f} W/m²): "
          f"global mean T̄ = {r.present.global_mean_T:.2f} °C, "
          f"ice line at {r.present_ice_line:.1f}° latitude (the finite-cap branch — Earth's)")
    print("    (a near-ice-free climate is ALSO stable at today's sun — present-day Earth is the "
          "colder, capped one)")
    print(f"  net top-of-atmosphere imbalance at equilibrium: {r.present.net_toa:+.2e} W/m²  "
          f"(absorbed solar = OLR — the conservation check)")
    print()
    print("  Hysteresis loop (continuation sweep, dim → freeze, brighten → melt):")
    print(f"    freezes over (Snowball jump) when dimmed to  S₀ ≈ {r.freeze_S0:.0f} W/m²  "
          f"(−{r.freeze_dimming_pct:.1f} %)")
    print(f"    re-melts only when brightened back up to     S₀ ≈ {r.melt_S0:.0f} W/m²")
    print(f"    hysteresis width = {r.hysteresis_width:.0f} W/m²  "
          f"(the Snowball is a deep trap — stays frozen far past where it froze)")
    print(f"    Snowball global-mean temperature: {r.snowball_Tbar:.1f} °C (a frozen white planet)\n")


def save_figure(r: SnowballResult) -> Path:
    """Render and save the Snowball hysteresis artifact (needs the optional ``viz`` extra)."""
    import matplotlib
    matplotlib.use("Agg")                            # headless
    from .plots import snowball_figure

    fig = snowball_figure(r.loop, r.present)
    for target in (DOCS_FIGURE, OUTPUT_FIGURE):
        target.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(target, dpi=130)
    return DOCS_FIGURE


def main() -> None:
    import sys
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")     # °C, ₂, →, ° on legacy codepages

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
