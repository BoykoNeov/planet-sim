"""Rung-4 completion: wiring the emergent gray OLR per latitude → tropical amplification.

Rung 4 (:mod:`planet.radiation`) derived where climlab's ``B = 2`` comes from, but in a *single
global-mean column*. This demo banks the figure for the named completion (:mod:`planet.radiative_ebm`):
letting that emergent ``OLR(Ts)`` *drive* the climate **per latitude**. The headline:

    The OLR slope is NOT one number. It is small at the warm equator (water vapour) and large at the
    cold pole — so under a uniform forcing the warming concentrates in the TROPICS.

That is the **mirror image** of rung-2.5's moisture-*transport* polar amplification: the water-vapour
*radiative* feedback alone favours the tropics, while it is transport (and the lapse-rate and ice
feedbacks held out of scope) that make Earth's poles amplify. Two rungs, two mechanisms, opposite signs.

**At the honest altitude:** the **sign was measured, not assumed** (the spike's discriminator — whichever
latitude has the smallest local slope warms most). The model runs at the *climlab-matched* water-vapour
loading (global-mean ``B = 2``, sub-runaway), solved by a coupled Newton iteration (the nonlinear
generalisation of rung-0's direct solve — the Strang relaxation carries a splitting error and goes unstable
near the warm-equator runaway edge). Direction banked; magnitude loose (it rides the water-vapour loading).

Three panels:
  1. The discriminator: the **local** OLR slope ``B_loc(φ)`` across the planet vs rung-0's flat ``B = 2`` —
     smallest at the warm equator (the whole story in one curve).
  2. Present-day climate: gray vs rung-0 — a near-uniform **Jensen warm shift** (concave OLR), the
     equator-to-pole contrast essentially unchanged.
  3. The warming response to a uniform forcing: gray (**tropical-amplified**) vs a uniform-slope null
     (warms uniformly) — the headline, the mirror of rung-2.5's polar amplification.

Run headless:  python -m planet.demo_radiative_ebm
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import radiative_ebm as rebm
from .albedo import EBMParams
from .ebm import B_OLR, EnergyBalanceModel
from .moist_ebm import constant_albedo_absorbed

_REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_FIGURE = _REPO_ROOT / "docs" / "figures" / "planet-radiative-ebm.png"
OUTPUT_FIGURE = _REPO_ROOT / "outputs" / "planet-radiative-ebm.png"


@dataclass(frozen=True)
class RadiativeEBMResult:
    """The banked rung-4-completion evidence — the discriminator, the present climates, and the warming."""

    ta: rebm.TropicalAmplification        # the headline experiment (gray + uniform-B null, warmed)
    rung0_present_T: np.ndarray           # rung-0's dt-free present climate (steady_linear) for comparison
    wv_fraction: float                    # the climlab-matched water-vapour loading


def compute() -> RadiativeEBMResult:
    """Run the headline experiment at the climlab-matched loading and the rung-0 reference beside it."""
    params = EBMParams()
    column = rebm.climlab_matched_column()
    ta = rebm.tropical_amplification(column, params, forcing=10.0)
    rung0 = EnergyBalanceModel(D=params.D, n_cells=params.n_cells,
                               face="harmonic").steady_linear(constant_albedo_absorbed(params))
    return RadiativeEBMResult(ta=ta, rung0_present_T=rung0.T, wv_fraction=column.wv_fraction)


def print_summary(r: RadiativeEBMResult) -> None:
    ta = r.ta
    print("\nGray OLR per latitude → tropical amplification (rung-4 completion)\n")
    print(f"  water-vapour loading     : {r.wv_fraction:.3f}  (climlab-matched: global-mean B = 2)")
    print(f"  local slope B_loc(equator/mid/pole) : "
          f"{ta.B_loc_present[0]:.2f} / {ta.B_loc_present[len(ta.B_loc_present)//2]:.2f} / "
          f"{ta.B_loc_present[-1]:.2f}  W/m²/K   (rung-0 assumes a flat {B_OLR:.0f})")
    print(f"  present mean  : gray {ta.gray_present.global_mean_T:.2f} °C  vs rung-0 "
          f"{r.rung0_present_T.mean():.2f} °C   (Jensen warm shift)")
    print(f"  present contrast (eq−pole): gray {ta.gray_present.T[0]-ta.gray_present.T[-1]:.1f} °C  vs rung-0 "
          f"{r.rung0_present_T[0]-r.rung0_present_T[-1]:.1f} °C   (≈ unchanged)")
    print(f"  warming ΔA = {ta.forcing:.0f} W/m²:  ⟨δT⟩ = {ta.mean_delta_T_gray:.2f} °C  "
          f"(naive ΔA/B_tan = {ta.dA_over_B_tan:.2f} — NOT pinned; Jensen + WV feedback amplify the mean)")
    print(f"  amplification δT(pole)/δT(equator): gray {ta.amp_gray:.2f}  (TROPICAL, < 1)  "
          f"vs uniform-B null {ta.amp_null:.2f}")
    print(f"  band ratio mean(≥60°)/mean(≤30°)  : gray {ta.amp_gray_band:.2f}")
    print(f"  => the mirror of rung-2.5's polar amplification (~1.4): WV radiative feedback favours the "
          f"tropics.\n")


def save_figure(r: RadiativeEBMResult) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ta = r.ta
    phi = ta.phi
    fig, ax = plt.subplots(1, 3, figsize=(16.5, 5.0))

    # Panel 1 — the discriminator: B_loc(φ) vs rung-0's flat B = 2.
    ax[0].plot(phi, ta.B_loc_present, "-", color="C3", lw=2,
               label="gray local slope B_loc(φ)  (computed)")
    ax[0].axhline(B_OLR, color="0.45", lw=1.6, ls="--", label=f"rung-0 prescribed flat B = {B_OLR:.0f}")
    imin = int(np.argmin(ta.B_loc_present))
    ax[0].plot([phi[imin]], [ta.B_loc_present[imin]], "o", color="C3", ms=8, zorder=5)
    ax[0].annotate("smallest damping\n→ warms most", (phi[imin], ta.B_loc_present[imin]),
                   textcoords="offset points", xytext=(28, 6), fontsize=8, color="C3")
    ax[0].set_xlabel("latitude φ (°)"); ax[0].set_ylabel("local OLR slope B_loc (W/m²/K)")
    ax[0].set_title("The OLR slope is NOT one number\nsmall at the warm equator (water vapour), large at the pole")
    ax[0].legend(fontsize=8, loc="upper right"); ax[0].grid(alpha=0.3)

    # Panel 2 — present-day climate: gray vs rung-0 (the Jensen warm shift, ~same contrast).
    ax[1].plot(phi, ta.gray_present.T, "-", color="C0", lw=2,
               label=f"gray  (mean {ta.gray_present.global_mean_T:.1f} °C)")
    ax[1].plot(phi, r.rung0_present_T, "--", color="0.45",
               label=f"rung-0  (mean {r.rung0_present_T.mean():.1f} °C)")
    ax[1].set_xlabel("latitude φ (°)"); ax[1].set_ylabel("temperature (°C)")
    ax[1].set_title("Present climate: a near-uniform Jensen warm shift\n(concave OLR → warmer mean, contrast ≈ unchanged)")
    ax[1].legend(fontsize=8, loc="upper right"); ax[1].grid(alpha=0.3)

    # Panel 3 — the headline: warming δT(φ), gray (tropical) vs the uniform-B null (uniform).
    ax[2].plot(phi, ta.delta_T_gray, "-", color="C3", lw=2,
               label=f"gray δT  (amp {ta.amp_gray:.2f} — TROPICAL)")
    ax[2].plot(phi, ta.delta_T_null, "--", color="0.45",
               label=f"uniform-B null  (amp {ta.amp_null:.2f})")
    ax[2].axhline(ta.dA_over_B_tan, color="C2", lw=1.0, ls=":",
                  label=f"naive ΔA/B_tan = {ta.dA_over_B_tan:.1f} (not the mean)")
    ax[2].set_xlabel("latitude φ (°)"); ax[2].set_ylabel("warming δT under uniform ΔA (°C)")
    ax[2].set_title("Warming is TROPICALLY amplified\n(mirror of rung-2.5's transport-driven polar amplification)")
    ax[2].legend(fontsize=8, loc="upper right"); ax[2].grid(alpha=0.3)

    fig.suptitle("Gray OLR per latitude (rung-4 completion): the latitude-varying OLR slope makes warming "
                 "tropical — the mirror of rung-2.5's polar amplification", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    for target in (DOCS_FIGURE, OUTPUT_FIGURE):
        target.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(target, dpi=110)
    plt.close("all")
    return DOCS_FIGURE


def main() -> None:
    import sys
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    r = compute()
    print_summary(r)
    try:
        saved = save_figure(r)
        print(f"Figure saved → {saved.relative_to(_REPO_ROOT)}")
    except ImportError:
        print("(matplotlib not installed — install the viz extra: pip install -e .[viz])")


if __name__ == "__main__":
    main()
