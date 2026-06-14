"""Rung-4 evidence: the emergent lapse-rate feedback (moist adiabat) — and how it OVERTURNED the scope.

Rung 4's gray column (:mod:`planet.radiation`) explained climlab's ``B = 2`` as ``Planck − water-vapour
+ the lapse-rate feedback the gray column omits``, with the lapse-rate term ``λ_LR ≈ 0.84`` *imported*
from Soden & Held (2006). This slice makes that term **emergent**: swap the fixed convective ``Γ`` for a
**moist adiabat** that flattens as it warms, so surface warming amplifies in the upper troposphere and the
outgoing longwave steepens with ``Ts``.

**The §12 scoping guess was OVERTURNED.** The expectation was "supplies ``λ_LR ≈ 0.84``, closing the gap to
climlab's 2." It does not — the emergent value is **``≈ +1.5``** and the moist-adiabat column **overshoots**:
its with-water-vapour ``B ≈ 3.1`` sits *above* climlab's 2, not at it. Banked: the **sign and kind**
(upper-troposphere amplification ⇒ a positive contribution to ``B``; the kernel closes; resolution-
converged). Loose: the **magnitude** — a single *global* moist-adiabat column captures only the **tropical**
branch (the deep tropics are moist-adiabatic, the extratropics are not), so it overshoots the *global-mean*
Soden & Held ``0.84``, and it rides the prescribed ``τ`` shape + water-vapour loading (the wall).

Three panels:
  1. The warming amplification ``ΔT(p)/ΔTs`` with height: fixed ``Γ`` warms uniformly (≈ 1 in the
     troposphere, falling to 0 at the pinned stratosphere); the moist adiabat warms the **upper
     troposphere more** (> 1) — the mechanism, measured.
  2. The Soden & Held **kernel waterfall** on the moist-adiabat column: ``Planck − water-vapour
     + lapse-rate = B_total``, the lapse-rate term now **emergent** (≈ 1.5), pushing ``B`` *above* 2.
  3. The **overturn**: fixed ``Γ`` (no feedback) leaves ``B`` below 2; the emergent (tropical) lapse-rate
     feedback overshoots to ``B`` above 2 — with the emergent ``λ_LR`` against the global-mean ``0.84``.

Run headless:  python -m planet.demo_lapse_rate
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import radiation as rad
from .ebm import B_OLR

_REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_FIGURE = _REPO_ROOT / "docs" / "figures" / "planet-lapse-rate.png"
OUTPUT_FIGURE = _REPO_ROOT / "outputs" / "planet-lapse-rate.png"


@dataclass(frozen=True)
class LapseRateResult:
    """The banked lapse-rate-feedback evidence — the warming amplification and the emergent kernel split."""

    pressure_hpa: np.ndarray     # pressure levels (hPa), TOA→surface
    amp_fixed: np.ndarray        # ΔT(p)/ΔTs for the fixed-Γ column (≈ 1 in the troposphere)
    amp_moist: np.ndarray        # ΔT(p)/ΔTs for the moist-adiabat column (> 1 aloft)
    planck: float                # kernel Planck term on the moist-adiabat column (W m⁻²/K)
    water_vapour: float          # kernel water-vapour term (< 0)
    lapse_rate_feedback: float   # kernel lapse-rate term — EMERGENT, ≈ +1.5 (the headline)
    b_total_moist: float         # with-WV slope of the moist-adiabat column (overshoots 2)
    b_total_fixed: float         # with-WV slope of the fixed-Γ column (below 2)
    kernel_residual: float       # closure residual (≈ 0)
    sh_lapse_rate: float         # Soden & Held global-mean λ_LR (0.84) — the overshot touchstone


def _amplification(col, dT=4.0):
    """ΔT(p)/ΔTs at each level — the vertical warming structure (1 = uniform, > 1 = amplified aloft)."""
    _, hi = col._profile(rad.PRESENT_SURFACE_T + dT)
    _, lo = col._profile(rad.PRESENT_SURFACE_T - dT)
    return (hi - lo) / (2.0 * dT)


def compute() -> LapseRateResult:
    """Calibrate both columns and read off the warming amplification and the emergent kernel split."""
    fixed = rad.calibrate_column()
    moist = rad.calibrate_column(moist_adiabat=True)

    p, _ = moist._profile(rad.PRESENT_SURFACE_T)
    k = moist.feedback_kernel()

    return LapseRateResult(
        pressure_hpa=p / 100.0,
        amp_fixed=_amplification(fixed),
        amp_moist=_amplification(moist),
        planck=k.planck,
        water_vapour=k.water_vapour,
        lapse_rate_feedback=k.lapse_rate,
        b_total_moist=k.total,
        b_total_fixed=fixed.feedback_kernel().total,
        kernel_residual=k.closure_residual,
        sh_lapse_rate=rad.SH_LAPSE_RATE,
    )


def print_summary(r: LapseRateResult) -> None:
    print("\nEmergent lapse-rate feedback: the moist adiabat (rung 4 within-rung upgrade)\n")
    print(f"  kernel split (moist-adiabat column):")
    print(f"    Planck      = {r.planck:+.2f} W/m²/K")
    print(f"    water vapour= {r.water_vapour:+.2f} W/m²/K")
    print(f"    lapse rate  = {r.lapse_rate_feedback:+.2f} W/m²/K   <-- EMERGENT (target was 0.84)")
    print(f"    -----------------------------------")
    print(f"    B_total     = {r.b_total_moist:+.2f} W/m²/K   (closure residual {r.kernel_residual:+.1e})")
    print(f"\n  OVERTURN: the scoped λ_LR ≈ 0.84 was expected to close the gap to climlab's 2.")
    print(f"    emergent λ_LR = {r.lapse_rate_feedback:.2f}  >>  global-mean Soden & Held 0.84  → OVERSHOOT")
    print(f"    fixed-Γ B = {r.b_total_fixed:.2f} (below 2, omits LR);  moist-adiabat B = {r.b_total_moist:.2f} "
          f"(ABOVE climlab's {B_OLR:.0f})")
    print(f"  WHY loose: a single GLOBAL moist-adiabat column = the TROPICAL branch only (extratropics are")
    print(f"    not moist-adiabatic); + it rides the τ shape & water-vapour loading (the wall).\n")


def save_figure(r: LapseRateResult) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(16.5, 5.2))

    # Panel 1 — the warming amplification with height (the mechanism).
    ax[0].plot(r.amp_fixed, r.pressure_hpa, "-", color="0.5", lw=2, label="fixed Γ (uniform shift)")
    ax[0].plot(r.amp_moist, r.pressure_hpa, "-", color="C3", lw=2, label="moist adiabat")
    ax[0].axvline(1.0, color="0.3", lw=1.0, ls=":")
    ax[0].annotate("upper troposphere\nwarms MORE\n(ΔT/ΔTs > 1)", (1.25, 320), fontsize=8, color="C3")
    ax[0].invert_yaxis()
    ax[0].set_xlabel("warming amplification  ΔT(p) / ΔTs"); ax[0].set_ylabel("pressure (hPa)")
    ax[0].set_title("The mechanism: a moist adiabat flattens as it warms\n→ the upper troposphere warms more than the surface")
    ax[0].legend(fontsize=8, loc="lower right"); ax[0].grid(alpha=0.3)

    # Panel 2 — the emergent Soden & Held kernel waterfall on the moist-adiabat column.
    steps = [
        ("Planck\n(uniform warming)", r.planck, 0.0, "C1"),
        ("− water vapour", r.water_vapour, r.planck, "C0"),
        (f"+ lapse rate\n(EMERGENT {r.lapse_rate_feedback:.2f})", r.lapse_rate_feedback,
         r.planck + r.water_vapour, "C2"),
    ]
    for i, (lab, delta, base, col) in enumerate(steps):
        ax[1].bar(i, delta, bottom=base, color=col, width=0.62, edgecolor="0.3")
        top = base + delta
        ax[1].text(i, top + (0.06 if delta >= 0 else -0.18), f"{top:.2f}", ha="center", fontsize=9)
        if i > 0:
            ax[1].plot([i - 1 + 0.31, i - 0.31], [base, base], color="0.4", lw=0.9, ls=":")
    ax[1].axhline(B_OLR, color="0.35", lw=1.6, ls="--", label=f"climlab prescribed B = {B_OLR:.0f}")
    ax[1].set_xticks(range(3)); ax[1].set_xticklabels([s[0] for s in steps], fontsize=7.5)
    ax[1].set_ylim(0, max(r.planck, r.b_total_moist) * 1.15)
    ax[1].set_ylabel("OLR slope contribution B (W/m²/K)")
    ax[1].set_title(f"Emergent kernel split: B = {r.b_total_moist:.2f}\n"
                    "the lapse-rate term is now COMPUTED, not imported")
    ax[1].legend(fontsize=8, loc="upper right"); ax[1].grid(alpha=0.3, axis="y")

    # Panel 3 — the overturn: fixed Γ below 2, emergent (tropical) LR overshoots above 2.
    labels = ["fixed Γ\n(no LR feedback)", "moist adiabat\n(emergent LR)"]
    vals = [r.b_total_fixed, r.b_total_moist]
    colors = ["0.6", "C3"]
    bars = ax[2].bar(labels, vals, color=colors, width=0.55, edgecolor="0.3")
    ax[2].axhline(B_OLR, color="0.35", lw=1.6, ls="--", label=f"climlab B = {B_OLR:.0f}")
    for b, v in zip(bars, vals):
        ax[2].text(b.get_x() + b.get_width() / 2, v + 0.05, f"{v:.2f}", ha="center", fontsize=10)
    ax[2].annotate(f"emergent λ_LR = {r.lapse_rate_feedback:.2f}\nvs global-mean S&H = {r.sh_lapse_rate:.2f}\n→ OVERSHOOT",
                   (0.5, max(vals) * 0.62), ha="center", fontsize=8.5, color="C2",
                   bbox=dict(boxstyle="round", fc="0.95", ec="C2"))
    ax[2].set_ylabel("with-water-vapour OLR slope B (W/m²/K)")
    ax[2].set_title("The OVERTURN: emergent LR OVERSHOOTS\n(tropical branch only → above 2, not at it)")
    ax[2].legend(fontsize=8, loc="upper left"); ax[2].grid(alpha=0.3, axis="y")

    fig.suptitle("Emergent lapse-rate feedback (rung 4): a moist adiabat supplies a POSITIVE, upper-trop-"
                 "amplified feedback (≈ 1.5) — the right sign & kind, but it OVERSHOOTS the global-mean 0.84",
                 fontsize=10)
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
