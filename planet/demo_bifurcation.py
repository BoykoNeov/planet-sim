"""Rung 0+: the complete equilibrium diagram — every climate the sun allows, and the second cliff.

Phase 1's Snowball demo *swept* the sun and watched the planet jump. This demo asks the inverse question
— *for an ice line here, what sun holds it?* — and answers it with one linear solve per ice line
(:mod:`planet.bifurcation`), which traces **every** equilibrium at once: the stable branches the sweep
rode *and* the unstable ones it could never sit on. The picture is the classic S-curve, and it has
**two** folds, not one:

* the **Snowball fold** (the freeze catastrophe Phase 1 found, near 33° / −8 % sun), and
* the **small-ice-cap fold** near 80°: a polar cap smaller than a critical radius ``θ_c ≈ 10°`` cannot
  be held by *any* sun. A brightening planet's cap does not shrink gracefully to nothing — it reaches
  ``θ_c`` and **vanishes in a jump** (North 1984's small-ice-cap instability, the Snowball's mirror).

Between the folds lies the whole **finite-cap window** — the only band of suns in which a planet with a
polar ice cap can exist. Today's sun sits inside it, a couple of W/m² below its upper edge. The stability
of every point is read off the curve's *slope* (Cahalan & North 1979) and *checked* by marching; the
Phase-1 sweep is overlaid to show it riding exactly the stable branches and jumping at the folds; and the
marcher's O(Δt) fixed-point bias — the reason Phase 1 needed a tiny step — is quantified against the
exact curve and retired.

Run headless (saves the figure, prints the diagram):

    python -m planet.demo_bifurcation
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import bifurcation as bf
from .albedo import EBMParams, HysteresisLoop, snowball_hysteresis
from .ebm import S0_EARTH

D_SWEEP = (0.2, 0.3, 0.4, 0.555, 0.7, 0.85, 1.0, 1.15, 1.3, 1.5)   # W m⁻² K⁻¹ — the transport sweep
N_TAU_SWEEP = (0.5, 0.2, 0.1, 0.05, 0.02)                            # the relaxation steps (× τ_rad)
SWEEP_N_TAU = 0.05                                                    # the Phase-1 demo's sweep step
CURVE_N_CELLS = 720            # the exact curve on a fine grid (the polar fold needs it); the sweep stays on Phase 1's

_REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_FIGURE = _REPO_ROOT / "docs" / "figures" / "planet-bifurcation.png"
OUTPUT_FIGURE = _REPO_ROOT / "outputs" / "planet-bifurcation.png"


@dataclass(frozen=True)
class BifurcationResult:
    """The banked equilibrium-diagram artifact — the plain bundle the figure and summary consume."""

    curve: bf.EquilibriumCurve
    loop: HysteresisLoop                 # the Phase-1 continuation sweep, overlaid
    D_values: np.ndarray
    theta_c: np.ndarray                  # critical cap radius per D (NaN = no finite-cap branch)
    window_lo: np.ndarray
    window_hi: np.ndarray
    n_taus: np.ndarray
    relaxed_ice_line: np.ndarray         # Phase-1 relaxation ice line per n_tau (present sun)

    @property
    def present(self) -> bf.Equilibrium:
        return self.curve.stable_finite_cap_at(S0_EARTH)


def compute(params: EBMParams | None = None) -> BifurcationResult:
    """Trace the curve, the Phase-1 sweep, the D-sweep and the step-bias sweep; bundle (no plotting)."""
    from dataclasses import replace
    p = params if params is not None else EBMParams()
    curve = bf.equilibrium_curve(replace(p, n_cells=CURVE_N_CELLS))
    loop = snowball_hysteresis(p, n_steps=60, n_tau=SWEEP_N_TAU)
    D_values = np.array(D_SWEEP)
    theta_c, lo, hi = bf.critical_cap_sweep(D_values, p)
    n_taus = np.array(N_TAU_SWEEP)
    relaxed = bf.relaxation_bias_sweep(n_taus, p, tol=1e-11, max_iter=400000)
    return BifurcationResult(curve, loop, D_values, theta_c, lo, hi, n_taus, relaxed)


def print_summary(r: BifurcationResult) -> None:
    """Print the diagram — the payoff in text."""
    c = r.curve
    lo, hi = c.snowball_fold, c.small_ice_cap_fold
    print("\nThe complete equilibrium diagram of the ice-albedo EBM (every branch, both folds)\n")
    print(f"  Snowball fold        : ice line {lo.latitude_deg:5.1f}°, S₀ = {lo.S0:7.1f} W/m² "
          f"({100 * (lo.S0 / S0_EARTH - 1):+.1f} % of today) — a dimming planet freezes over here")
    print(f"  small-ice-cap fold   : ice line {hi.latitude_deg:5.1f}°, S₀ = {hi.S0:7.1f} W/m² "
          f"({100 * (hi.S0 / S0_EARTH - 1):+.2f} %) — the cap cannot shrink below θ_c = {hi.cap_radius_deg:.1f}°")
    print(f"  finite-cap window    : {lo.S0:.0f} … {hi.S0:.0f} W/m² — the only suns that hold a polar cap\n")
    print(f"  ice-free branch from : S₀ ≥ {c.ice_free_threshold_S0:.1f} (below it the pole must freeze)")
    print(f"  Snowball branch to   : S₀ ≤ {c.snowball_threshold_S0:.1f} (above it the equator must thaw)\n")

    eqs = c.equilibria_at(S0_EARTH)
    print(f"  at today's sun ({S0_EARTH} W/m²) the model holds {len(eqs)} equilibria:")
    for e in sorted(eqs, key=lambda e: -e.x_ice):
        tag = "stable  " if e.stable else "UNSTABLE"
        print(f"    {tag}  {e.kind:10s}  ice line {e.latitude_deg:5.1f}°")
    pres = r.present
    print(f"  → Earth's branch is the stable finite cap at {pres.latitude_deg:.1f}° (benchmark ~70°), only "
          f"{hi.S0 - S0_EARTH:.1f} W/m² ({100 * (hi.S0 - S0_EARTH) / S0_EARTH:.2f} %) below the small-ice-cap cliff.\n")

    print("  the Phase-1 continuation sweep (an independent method) against the exact folds:")
    print(f"    sweep freezes at {r.loop.freeze_S0:.0f}  vs fold {lo.S0:.0f};  re-melts at {r.loop.melt_S0:.0f}  "
          f"vs branch end {c.snowball_threshold_S0:.0f}  (sweep step {abs(r.loop.S0_up[1] - r.loop.S0_up[0]):.0f})\n")

    print("  critical cap radius θ_c vs transport D (~10° at weak transport, growing once D ≳ 0.4):")
    for D, th, w0, w1 in zip(r.D_values, r.theta_c, r.window_lo, r.window_hi):
        if np.isnan(th):
            print(f"    D = {D:5.3f}:  no stable finite cap at all — ice-free or Snowball only")
        else:
            print(f"    D = {D:5.3f}:  θ_c = {th:5.1f}°   window {w0:6.0f} … {w1:6.0f} W/m²  (width {w1 - w0:5.0f})")

    print("\n  Phase 1's relaxation vs the exact curve at today's sun (the O(Δt) fixed-point bias, quantified):")
    for nt, lat in zip(r.n_taus, r.relaxed_ice_line):
        print(f"    n_tau = {nt:4.2f}: relaxed ice line {lat:5.1f}°   (exact {pres.latitude_deg:.1f}°, "
              f"gap {lat - pres.latitude_deg:+.1f}°)")
    print("    → the bias shrinks with the step; the last ~1° is the marcher's cell-quantized ice edge on its\n"
          "      180-cell grid (it halves per grid doubling). The exact diagram needs no step at all.\n")


def save_figure(r: BifurcationResult) -> Path:
    """Render and save the equilibrium-diagram artifact (needs the optional ``viz`` extra)."""
    import matplotlib
    matplotlib.use("Agg")                                    # headless
    from .plots import bifurcation_figure

    fig = bifurcation_figure(r)
    for target in (DOCS_FIGURE, OUTPUT_FIGURE):
        target.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(target, dpi=130)
    return DOCS_FIGURE


def main() -> None:
    import sys
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")             # °, ₀, →, θ on legacy codepages

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
