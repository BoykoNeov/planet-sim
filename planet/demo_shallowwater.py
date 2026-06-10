"""The Phase-3 banked artifact: the rotating shallow-water engine, exercised with planetary numbers.

The capstone's new shared engine (:mod:`engines.fluid`) made visible. It drives the planetary
β-plane (:mod:`planet.circulation`) through the two classic rotating-fluid demonstrations:

  1. **Geostrophic adjustment** — an unbalanced height anomaly (a pressure bump at rest) radiates
     inertia-gravity waves and **settles into a geostrophically-balanced vortex of scale ``L_R``**
     (the deformation radius). Most of the bump's height is shed as waves; the balanced remnant *is*
     the analytic Helmholtz-adjusted state. The conservation diagnostics (mass machine-exact;
     energy / potential enstrophy bounded) hold flat throughout.
  2. **A westward-propagating Rossby wave** — a balanced large-scale undulation drifts *west*, the
     planetary-vorticity-gradient (β) signature, at the analytic phase speed.

*The rotating fluid engine, validated and banked* — the structural counterpart of Steel freezing
the diffusion spine, now for the program's second engine. Phase 4 will couple the EBM to it so a
midlatitude jet emerges.

Run headless (saves the figure, prints the summary):

    python -m planet.demo_shallowwater
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .circulation import AdjustmentRun, RossbyRun, geostrophic_adjustment, rossby_wave

_REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_FIGURE = _REPO_ROOT / "docs" / "figures" / "planet-shallowwater.png"
OUTPUT_FIGURE = _REPO_ROOT / "outputs" / "planet-shallowwater.png"


@dataclass(frozen=True)
class ShallowWaterResult:
    """The banked Phase-3 artifact — the plain bundle the figure and the summary consume."""

    adjustment: AdjustmentRun
    rossby: RossbyRun

    @property
    def helmholtz_rel_error(self) -> float:
        a = self.adjustment
        return float(np.max(np.abs(a.eta_balanced - a.eta_helmholtz)) / np.max(np.abs(a.eta_helmholtz)))

    @property
    def drawdown_fraction(self) -> float:
        """Fraction of the initial bump height shed as radiated gravity waves."""
        a = self.adjustment
        return float(1.0 - a.eta_balanced.max() / a.eta_init.max())


def compute(nx: int = 96, ny: int = 96) -> ShallowWaterResult:
    """Run the geostrophic-adjustment + Rossby-wave demos → :class:`ShallowWaterResult` (no plotting)."""
    return ShallowWaterResult(
        adjustment=geostrophic_adjustment(nx=nx, ny=ny, n_periods=25.0),
        rossby=rossby_wave(nx=nx, ny=ny, mk=1, ml=1, frac_period=0.5),
    )


def print_summary(r: ShallowWaterResult) -> None:
    """Print the Phase-3 story — the engine's behaviour in text."""
    a, ros = r.adjustment, r.rossby
    print("\nThe rotating shallow-water engine (engines/fluid) — geostrophic adjustment & Rossby waves\n")
    print(f"  deformation radius  L_R = {a.L_R/1e3:.0f} km   (√(gH)/f₀, the midlatitude scale)")
    print("  Geostrophic adjustment (a pressure bump relaxes to balance):")
    print(f"    {100*r.drawdown_fraction:.0f} % of the bump height radiated away as gravity waves;")
    print(f"    the balanced remnant matches the analytic Helmholtz state to "
          f"{100*r.helmholtz_rel_error:.1f} %  (scale ≈ L_R)")
    print(f"    mass drift  {np.abs(a.mass).max():.1e} (machine-exact),  "
          f"energy {np.abs(a.energy).max():.1e},  enstrophy {np.abs(a.enstrophy).max():.1e}")
    print("  Rossby wave (β-plane):")
    print(f"    phase speed  c = {ros.c_measured:.1f} m/s  (WESTWARD; analytic {ros.c_analytic:.1f} m/s)\n")


def save_figure(r: ShallowWaterResult) -> Path:
    """Render and save the Phase-3 artifact (needs the optional ``viz`` extra)."""
    import matplotlib
    matplotlib.use("Agg")                            # headless
    from .plots import shallowwater_figure

    fig = shallowwater_figure(r.adjustment, r.rossby)
    for target in (DOCS_FIGURE, OUTPUT_FIGURE):
        target.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(target, dpi=130)
    return DOCS_FIGURE


def main() -> None:
    import sys
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")     # ², ₀, →, ° on legacy codepages

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
