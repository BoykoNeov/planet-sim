"""Rung-4 evidence: where the prescribed OLR ``A + B·T`` comes from (gray radiative transfer).

The rung-0 EBM *prescribes* the outgoing longwave ``OLR = A + B·T`` from the climlab constants
(``A = 210``, ``B = 2``) — it cannot say where ``B = 2`` comes from, cannot produce a water-vapour
feedback, and cannot compute a CO₂ forcing. This demo banks the figure for the **emergent** OLR computed
by a gray radiative–convective column (:mod:`planet.radiation`), and the rung's headline:

    climlab's B = 2  ≈  the gray emission-level Planck slope (≈ 3.4 W m⁻² K⁻¹)  −  the water-vapour feedback.

**At the honest altitude:** the present operating point (33 K greenhouse, OLR = 239) is matched **by
construction** (the optical depth is calibrated to it), so the emergent OLR and the linear ``A + B·T``
agree in *value* at present — the **slope** is the finding. The no-WV slope sits near the ``4σTₑ³`` Planck
touchstone (tight); the water-vapour feedback (``τ`` rising with ``Ts`` through Clausius–Clapeyron) drops
it through climlab's 2 — **direction banked, magnitude loose** (it rides on the water-vapour optical-depth
loading, the calibrated wall). And the gray CO₂ forcing is **saturating, not logarithmic** — the named
limitation that motivates the within-rung spectral-band upgrade.

Three panels:
  1. The emergent ``OLR(Ts)`` (with water vapour), the present operating point, and rung-0's ``A + B·T``
     tangent — the linear OLR is the *local tangent* of the computed curve.
  2. The decomposition waterfall: climlab's ``B = 2 ≈ Planck slope − water-vapour feedback + the
     lapse-rate feedback the gray column omits`` — every term at a Soden & Held (2006) order, not tuned.
  3. The CO₂ forcing per doubling: it **decreases** at high CO₂ (saturating gray band), not the constant
     per-doubling logarithmic law.

Run headless:  python -m planet.demo_radiation
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import radiation as rad
from .ebm import A_OLR, B_OLR

_REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_FIGURE = _REPO_ROOT / "docs" / "figures" / "planet-radiation.png"
OUTPUT_FIGURE = _REPO_ROOT / "outputs" / "planet-radiation.png"


@dataclass(frozen=True)
class RadiationResult:
    """The banked rung-4 evidence — the emergent OLR, the slope decomposition, and the CO₂ forcing."""

    Te: float                    # emission temperature (K)
    ts_curve: np.ndarray         # surface temperatures for the OLR curve (K)
    olr_curve: np.ndarray        # emergent OLR with water vapour (W m⁻²)
    olr_climlab: np.ndarray      # rung-0's prescribed A + B·T over the same Ts (W m⁻²)
    planck_slope: float          # 4σTₑ³ — the emission-level Planck touchstone
    b_dry: float                 # emergent slope, no water-vapour feedback (≈ Planck)
    b_moist: float               # emergent slope with water vapour (drops through 2)
    wv_feedback: float           # b_dry − b_moist > 0
    wv_fractions: np.ndarray     # the water-vapour optical-depth loading sweep
    b_moist_sweep: np.ndarray    # b_moist at each loading
    wv_fraction_at_climlab: float  # the loading at which b_moist = climlab's B (the recovered wall)
    co2_mult: np.ndarray         # CO₂ optical-depth multiples
    co2_per_doubling: np.ndarray  # forcing added by each doubling (W m⁻²) — decreasing = saturating


def compute() -> RadiationResult:
    """Calibrate the present-day column and read off the emergent OLR, slope decomposition, and forcing."""
    col = rad.calibrate_column()
    Te = rad.emission_temperature()

    ts_curve = np.linspace(273.0, 303.0, 31)
    olr_curve = np.array([col.outgoing_longwave(t, water_vapour=True) for t in ts_curve])
    olr_climlab = A_OLR + B_OLR * (ts_curve - 273.15)

    b_dry, b_moist, wv_feedback = col.feedback_decomposition()
    planck_slope = 4.0 * rad.STEFAN_BOLTZMANN * Te ** 3

    # The water-vapour loading sweep — how much the feedback subtracts, and where climlab's 2 is recovered.
    wv_fractions = np.linspace(0.2, 0.8, 13)
    b_moist_sweep = np.array([rad.calibrate_column(wv_fraction=f).feedback_slope(water_vapour=True)
                              for f in wv_fractions])
    # b_moist DECREASES with loading; find the crossing b_moist = B_OLR by interpolation on the reversed arrays
    order = np.argsort(b_moist_sweep)
    wv_at_climlab = float(np.interp(B_OLR, b_moist_sweep[order], wv_fractions[order]))

    co2_mult = np.array([2.0, 4.0, 8.0, 16.0, 32.0])
    olr_at = {1.0: col.outgoing_longwave(rad.PRESENT_SURFACE_T, co2_factor=1.0, water_vapour=False)}
    for m in co2_mult:
        olr_at[m] = col.outgoing_longwave(rad.PRESENT_SURFACE_T, co2_factor=m, water_vapour=False)
    prev = 1.0
    per_doubling = []
    for m in co2_mult:
        per_doubling.append(olr_at[prev] - olr_at[m])
        prev = m

    return RadiationResult(
        Te=Te, ts_curve=ts_curve, olr_curve=olr_curve, olr_climlab=olr_climlab,
        planck_slope=planck_slope, b_dry=b_dry, b_moist=b_moist, wv_feedback=wv_feedback,
        wv_fractions=wv_fractions, b_moist_sweep=b_moist_sweep,
        wv_fraction_at_climlab=wv_at_climlab,
        co2_mult=co2_mult, co2_per_doubling=np.array(per_doubling),
    )


def print_summary(r: RadiationResult) -> None:
    print("\nGray radiative transfer: where A + B·T comes from (rung 4)\n")
    print(f"  emission temperature Tₑ = {r.Te:.1f} K   (greenhouse = {rad.PRESENT_SURFACE_T - r.Te:.1f} K)")
    print(f"  present operating point  : OLR({rad.PRESENT_SURFACE_T:.0f} K) = {rad.PRESENT_OLR:.0f} W/m² "
          f"(matched by construction)")
    print(f"  Planck touchstone 4σTₑ³  : {r.planck_slope:.2f} W/m²/K")
    print(f"  emergent slope, no WV    : B = {r.b_dry:.2f} W/m²/K  (≈ Planck λ₀ {rad.SH_PLANCK}, ABOVE 2)")
    print(f"  emergent slope, with WV  : B = {r.b_moist:.2f} W/m²/K  (WV subtracts {r.wv_feedback:.2f} "
          f"≈ λ_wv {rad.SH_WATER_VAPOUR})")
    print(f"  => climlab's B = {B_OLR:.1f}  ≈  Planck {r.b_dry:.2f} − WV {r.wv_feedback:.2f} + lapse-rate "
          f"{rad.SH_LAPSE_RATE:.2f}  = {r.b_dry - r.wv_feedback + rad.SH_LAPSE_RATE:.2f}")
    print(f"     (every term at a Soden & Held 2006 order — ORDER-VALIDATED, not tuned; the lapse-rate")
    print(f"      feedback is the one a FIXED gray lapse rate omits. Planck+WV alone = 2 at loading "
          f"{r.wv_fraction_at_climlab:.2f}.)")
    print(f"  CO₂ forcing per doubling : {np.array2string(r.co2_per_doubling, precision=1)} W/m² "
          f"(DECREASING → saturating, NOT logarithmic)\n")


def save_figure(r: RadiationResult) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(16.5, 5.0))

    # Panel 1 — the emergent OLR(Ts) and rung-0's A + B·T tangent.
    ax[0].plot(r.ts_curve, r.olr_curve, "-", color="C0", lw=2,
               label="emergent OLR (computed, with water vapour)")
    ax[0].plot(r.ts_curve, r.olr_climlab, "--", color="0.45",
               label=f"rung-0 prescribed A + B·T  (B={B_OLR:.0f})")
    ax[0].plot([rad.PRESENT_SURFACE_T], [rad.PRESENT_OLR], "o", color="C3", ms=8, zorder=5)
    ax[0].annotate("present\n(288 K, 239 W/m²)\nmatched by construction",
                   (rad.PRESENT_SURFACE_T, rad.PRESENT_OLR), textcoords="offset points",
                   xytext=(-105, -8), fontsize=8, color="C3")
    ax[0].set_xlabel("surface temperature Ts (K)"); ax[0].set_ylabel("outgoing longwave OLR (W/m²)")
    ax[0].set_title("The OLR is COMPUTED, not prescribed\nrung-0's A + B·T shares its present operating point")
    ax[0].legend(fontsize=8, loc="upper left"); ax[0].grid(alpha=0.3)

    # Panel 2 — the decomposition waterfall: climlab's B = Planck − water-vapour + lapse-rate (gray omits).
    net = r.b_dry - r.wv_feedback + rad.SH_LAPSE_RATE
    steps = [
        ("Planck slope\n(gray, no WV)", r.b_dry, 0.0, "C1"),
        (f"− water vapour\n(≈ λ_wv {rad.SH_WATER_VAPOUR})", -r.wv_feedback, r.b_dry, "C0"),
        (f"+ lapse rate\n(λ_LR {rad.SH_LAPSE_RATE}, gray OMITS)", rad.SH_LAPSE_RATE, r.b_dry - r.wv_feedback, "C2"),
    ]
    tops = []
    for i, (lab, delta, base, col) in enumerate(steps):
        ax[1].bar(i, delta, bottom=base, color=col, width=0.62, edgecolor="0.3")
        top = base + delta
        tops.append(top)
        ax[1].text(i, top + (0.06 if delta >= 0 else -0.16), f"{top:.2f}", ha="center", fontsize=9)
        if i > 0:                                            # waterfall connector from the previous top
            ax[1].plot([i - 1 + 0.31, i - 0.31], [base, base], color="0.4", lw=0.9, ls=":")
    ax[1].axhline(B_OLR, color="0.35", lw=1.6, ls="--", label=f"climlab prescribed B = {B_OLR:.0f}")
    ax[1].set_xticks(range(3)); ax[1].set_xticklabels([s[0] for s in steps], fontsize=7.5)
    ax[1].set_ylim(0, r.b_dry * 1.12)
    ax[1].set_ylabel("OLR slope contribution B (W/m²/K)")
    ax[1].set_title(f"climlab's B = Planck − WV + lapse rate = {net:.2f} ≈ {B_OLR:.0f}\n"
                    "(every term at a Soden–Held 2006 order — not tuned)")
    ax[1].legend(fontsize=8, loc="upper right"); ax[1].grid(alpha=0.3, axis="y")

    # Panel 3 — the saturating (non-logarithmic) CO₂ forcing.
    doublings = np.arange(1, len(r.co2_mult) + 1)
    ax[2].bar(doublings, r.co2_per_doubling, color="C0", width=0.6)
    ax[2].axhline(r.co2_per_doubling[0], color="0.45", lw=1.2, ls="--",
                  label="a LOGARITHMIC law would be flat\n(constant per doubling)")
    for d, v in zip(doublings, r.co2_per_doubling):
        ax[2].text(d, v + 0.8, f"{v:.0f}", ha="center", fontsize=8)
    ax[2].set_xticks(doublings)
    ax[2].set_xticklabels([f"{int(2**d)}×" for d in doublings])
    ax[2].set_xlabel("CO₂ optical depth (× present)")
    ax[2].set_ylabel("forcing added by this doubling ΔF (W/m²)")
    ax[2].set_title("Gray CO₂ forcing SATURATES\n(per-doubling falls → not the log law)")
    ax[2].legend(fontsize=8, loc="upper right"); ax[2].grid(alpha=0.3, axis="y")

    fig.suptitle("Gray radiative transfer (rung 4): the emergent OLR explains climlab's B = 2 as Planck − "
                 "water-vapour + the lapse-rate feedback gray omits (Soden–Held orders, not tuned)", fontsize=10)
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
