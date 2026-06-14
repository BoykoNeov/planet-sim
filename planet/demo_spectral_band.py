"""Rung-4 within-rung slice: the spectral-band log law — why CO₂ forcing is LOGARITHMIC.

The gray column (:mod:`planet.demo_radiation`, panel 3) gets the CO₂ forcing's *sign* right but the
*shape* wrong: a band-independent (gray) absorption **saturates** — each doubling adds *less* (per
doubling ``ΔF`` ≈ 48 → 53 → 41 → 25 → 20 W m⁻², decreasing) — because adding CO₂ pushes the *whole*
emission spectrum to the cold upper atmosphere. The observed law is **logarithmic** (Myhre+ 1998:
``ΔF = 5.35·ln(C/C₀)`` ≈ 3.7 W m⁻² *per doubling, constant*).

This demo banks the figure for the fix (:class:`planet.radiation.SpectralCO2Band`): resolve the CO₂
15-µm band into spectral bins whose absorption falls off **exponentially in the wings**. The band
centre is saturated, so the forcing comes only from the wings — and an exponential wing's ``τ = 1``
emission level spreads by a *constant* spectral width per doubling, so ``ΔF`` is **constant per
doubling** = the log law. The exponential wing is the whole ingredient: flatten it (uniform ``k``) and
the forcing saturates again.

**At the honest altitude:** the *functional form* (logarithmic) is the win; the *magnitude* (~4 W m⁻²
per doubling, the Myhre band) rides the band parameters (wing scale, band-centre τ, half-width) —
calibrated to **order**, not line-by-line spectroscopy (the wall, the same cited-closure status as the
gray ``τ ↔ greenhouse`` mapping). And the log law is **range-limited**: below where the band centre
saturates it is linear/√, above where the wings exhaust it saturates again — the realistic 0.5×–8×
window sits in the flat middle.

Two panels:
  1. ``ΔF`` per doubling — gray (decreasing, saturating) vs the spectral band (flat, logarithmic) vs
     the Myhre constant. The mirror of the gray saturation panel.
  2. Cumulative forcing ``F(C)`` on a log CO₂ axis — the spectral band is a **straight line** (= the
     logarithm) where gray **bends over** (saturates); Myhre's ``5.35·ln(C/C₀)`` for reference.

Run headless:  python -m planet.demo_spectral_band
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import radiation as rad

_REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_FIGURE = _REPO_ROOT / "docs" / "figures" / "planet-spectral-band.png"
OUTPUT_FIGURE = _REPO_ROOT / "outputs" / "planet-spectral-band.png"


@dataclass(frozen=True)
class SpectralResult:
    """The banked spectral-band evidence — the per-doubling forcing and the cumulative log law."""

    doubling_factors: np.ndarray   # the successive CO₂ multiples (0.5,1,2,4,8,16,32)
    gray_per_doubling: np.ndarray  # gray ΔF per doubling — DECREASING (saturating)
    band_per_doubling: np.ndarray  # spectral ΔF per doubling — CONSTANT (logarithmic)
    band_per_doubling_mid: np.ndarray  # the flat-middle subset (0.5×–8×) that is pinned constant
    co2_curve: np.ndarray          # CO₂ multiples for the cumulative-forcing curves (log-spaced)
    gray_forcing: np.ndarray       # cumulative gray F(C) (W m⁻²) — bends over
    band_forcing: np.ndarray       # cumulative spectral F(C) (W m⁻²) — straight in ln C
    myhre_forcing: np.ndarray      # Myhre 5.35·ln(C/C₀) reference (W m⁻²)
    band_per_doubling_mean: float  # the constant (W m⁻²/doubling), in the Myhre band
    myhre_per_doubling: float      # 5.35·ln2 ≈ 3.71 W m⁻²
    analytic_coefficient: float    # the τ=1 wing-formula dF/dlnC (the sharp-limit derivation)


def compute() -> SpectralResult:
    """Read off the gray (saturating) and spectral (logarithmic) CO₂ forcings for the figure."""
    col = rad.calibrate_column()
    band = rad.SpectralCO2Band()
    Ts = rad.PRESENT_SURFACE_T

    factors = np.array([0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0])
    n_doublings = np.log2(factors[1:] / factors[:-1])

    def gray_olr(m):
        return col.outgoing_longwave(Ts, co2_factor=m, water_vapour=False)
    gray_olr_at = np.array([gray_olr(m) for m in factors])
    gray_per = -np.diff(gray_olr_at) / n_doublings           # OLR falls → positive forcing

    band_per = band.forcing_per_doubling(factors)
    mid = band.forcing_per_doubling((0.5, 1, 2, 4, 8))       # the flat-middle subset that is pinned

    # Cumulative forcing F(C) = OLR(1×) − OLR(C×); on a log CO₂ axis the log law is a straight line.
    co2_curve = np.geomspace(0.5, 32.0, 25)
    gray_base = gray_olr(1.0)
    gray_F = np.array([gray_base - gray_olr(m) for m in co2_curve])
    band_F = np.array([band.co2_forcing(m) for m in co2_curve])
    myhre_F = rad.MYHRE_COEFFICIENT * np.log(co2_curve)

    return SpectralResult(
        doubling_factors=factors, gray_per_doubling=gray_per, band_per_doubling=band_per,
        band_per_doubling_mid=mid, co2_curve=co2_curve, gray_forcing=gray_F, band_forcing=band_F,
        myhre_forcing=myhre_F, band_per_doubling_mean=float(mid.mean()),
        myhre_per_doubling=rad.MYHRE_PER_DOUBLING, analytic_coefficient=band.log_law_coefficient(),
    )


def print_summary(r: SpectralResult) -> None:
    print("\nSpectral-band log law: why CO₂ forcing is logarithmic, not saturating (rung 4 slice)\n")
    labels = [f"{a:g}→{b:g}×" for a, b in zip(r.doubling_factors[:-1], r.doubling_factors[1:])]
    print("  per-doubling ΔF (W/m²):")
    print(f"    gray (band-indep.) : " + "  ".join(f"{lab}:{v:5.1f}" for lab, v in zip(labels, r.gray_per_doubling))
          + "   → DECREASING (saturating)")
    print(f"    spectral (exp wings): " + "  ".join(f"{lab}:{v:5.2f}" for lab, v in zip(labels, r.band_per_doubling))
          + "   → CONSTANT (logarithmic)")
    spread = r.band_per_doubling_mid.max() / r.band_per_doubling_mid.min()
    print(f"\n  spectral flat-middle 0.5×–8×: {r.band_per_doubling_mean:.2f} W/m²/doubling "
          f"(spread {spread:.3f}×)")
    print(f"  Myhre 5.35·ln2             : {r.myhre_per_doubling:.2f} W/m²/doubling (order-validated, not tuned)")
    print(f"  analytic τ=1 wing formula  : {r.analytic_coefficient * np.log(2):.2f} W/m²/doubling "
          f"(sharp limit; column sits ~20–30% above)\n")


def save_figure(r: SpectralResult) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(13.0, 5.2))

    # Panel 1 — ΔF per doubling: gray saturates, the spectral band is flat (the log law).
    labels = [f"{int(b)}×" if b >= 1 else f"{b:g}×" for b in r.doubling_factors[1:]]
    x = np.arange(len(labels))
    ax[0].bar(x - 0.2, r.gray_per_doubling, width=0.38, color="0.6", edgecolor="0.3",
              label="gray (band-independent) — saturates")
    ax[0].bar(x + 0.2, r.band_per_doubling, width=0.38, color="C0", edgecolor="0.3",
              label="spectral band (exp. wings) — logarithmic")
    ax[0].axhline(r.myhre_per_doubling, color="C3", lw=1.6, ls="--",
                  label=f"Myhre 5.35·ln2 ≈ {r.myhre_per_doubling:.1f} (constant)")
    for xi, v in zip(x - 0.2, r.gray_per_doubling):
        ax[0].text(xi, v + 0.8, f"{v:.0f}", ha="center", fontsize=7.5, color="0.3")
    for xi, v in zip(x + 0.2, r.band_per_doubling):
        ax[0].text(xi, v + 0.8, f"{v:.1f}", ha="center", fontsize=7.5, color="C0")
    ax[0].set_xticks(x); ax[0].set_xticklabels(labels)
    ax[0].set_xlabel("CO₂ optical depth (× present)")
    ax[0].set_ylabel("forcing added by this doubling ΔF (W/m²)")
    ax[0].set_title("Gray SATURATES (ΔF falls); exponential wings give the\nCONSTANT-per-doubling Myhre log law")
    ax[0].legend(fontsize=8, loc="upper right"); ax[0].grid(alpha=0.3, axis="y")

    # Panel 2 — cumulative F(C) on a log CO₂ axis: the log law is a straight line.
    ax[1].semilogx(r.co2_curve, r.band_forcing, "-", color="C0", lw=2.2,
                   label="spectral band — straight (= logarithmic)")
    ax[1].semilogx(r.co2_curve, r.gray_forcing, "-", color="0.55", lw=2.0,
                   label="gray — bends over (saturates)")
    ax[1].semilogx(r.co2_curve, r.myhre_forcing, "--", color="C3", lw=1.6,
                   label="Myhre 5.35·ln(C/C₀)")
    ax[1].axvline(1.0, color="0.7", lw=0.9, ls=":")
    ax[1].set_xlabel("CO₂ optical depth C/C₀  (log axis)")
    ax[1].set_ylabel("cumulative forcing F(C) = OLR(1×) − OLR(C) (W/m²)")
    ax[1].set_title("On a log CO₂ axis the log law is a STRAIGHT LINE\n(spectral straight; gray curves over)")
    ax[1].legend(fontsize=8, loc="upper left"); ax[1].grid(alpha=0.3, which="both")

    fig.suptitle("Spectral-band log law (rung-4 slice): exponential band wings turn the gray SATURATING "
                 "CO₂ forcing into the logarithmic Myhre law — form banked, magnitude calibrated to order",
                 fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
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
