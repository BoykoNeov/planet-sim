"""Rung-1 evidence: wet-get-wetter, dry-get-drier — the thermodynamic contrast sharpening.

:mod:`planet.precip` scales the *whole* precip pattern by **one** Clausius–Clapeyron amplitude
``CC(T̄)`` (~7 %/K), so under warming **everywhere wettens — including the deserts** (the opposite of
what is observed). The Held & Soden 2006 "rich-get-richer" pattern instead **sharpens the contrast**:
the global *mean* precipitation is energy-constrained to a slower ~2.5 %/K, while the wet−dry *anomaly*
intensifies at the faster ~7 %/K moisture rate. This demo banks the before/after figure for the opt-in
contrast split (:func:`planet.moist.wet_get_wetter_precip_field`).

**At the honest altitude** (the ``precip.py`` class — a better *prescribed* parameterization, not a
derived field): the split *direction* (the contrast grows faster than the mean) is the **prescription**;
the two rates are **cited closures** (the energy-constrained slope is the named sub-grid wall, the C–C
7 %/K the moisture-capacity rate). What is structurally exact is the **mean-zero anomaly split** (so the
area-mean precip scales at the energy rate) and the **reduction** to both rung-0 and the
energy-constrained field when the rates coincide. The *dynamic* pattern shift (changing circulation) is
the moisture-convergence / rung-3+ route — named, not here.

**The named overreach (pinned, not papered): it dries the POLES, the wrong sign.** Scaling *every*
below-mean anomaly at the C–C rate dries all latitudes below the mean — including the high latitudes
(`P ≈ 20 ≪ ⟨P⟩`), which drop *harder* than the deserts (the poleward orange band in panel 2). But
observed high-latitude precipitation **increases** under warming (poleward moisture transport — a
*dynamic* effect). "Wet-get-wetter, dry-get-drier" is a **tropical/subtropical** idealization: right in
the subtropics (deserts dry), **wrong at the poles** (rung 3+) — annotated on the figure, not hidden.

Two panels:
  1. ``P(φ)``: the present field vs a warmed planet under the rung-0 *uniform* C–C vs the wet-get-wetter
     split — the ITCZ/storm-track bands intensify under both, but only the split **dries the deserts**.
  2. The warming response ``ΔP(φ) = warmed − present``: rung-0 uniform is positive *everywhere* (deserts
     wetten — the flaw), while wet-get-wetter is **negative in the subtropics** (dry-get-drier) — and also
     at the **poles** (the named overreach: high-lat P should *increase*, annotated as wrong).

The committed figure + the ``slow`` test assertion (the desert sign, the conserving mean, the reduction)
let a fresh clone reproduce the headline, not just read the conclusion.

Run headless:  python -m planet.demo_wet_get_wetter
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from . import moist, precip
from .albedo import EBMParams, present_day_climate

_REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_FIGURE = _REPO_ROOT / "docs" / "figures" / "planet-wet-get-wetter.png"
OUTPUT_FIGURE = _REPO_ROOT / "outputs" / "planet-wet-get-wetter.png"

WARMING_K = 6.0          # the warming used to expose the contrast split (unclamped — the floor is inactive)


def compute() -> moist.WetGetWetter:
    """The wet-get-wetter trade for the present-day climate warmed by :data:`WARMING_K`."""
    return moist.wet_get_wetter(present_day_climate(EBMParams()), delta_T=WARMING_K)


def print_summary(w: moist.WetGetWetter) -> None:
    pbar = float(np.mean(w.precip_present))
    print("\nWet-get-wetter, dry-get-drier (rung 1, the thermodynamic contrast split)\n")
    print(f"  warming:           +{w.delta_T:.0f} K above present")
    print(f"  ITCZ peak:         {w.itcz_change_wgw:+7.1f} cm/yr  (wet-get-wetter — the wet band intensifies)")
    print(f"  desert @ {w.subtropics_lat:.0f}°:     wet-get-wetter {w.subtropics_change_wgw:+6.1f} cm/yr  (DRY-get-drier)")
    print(f"                     rung-0 uniform {w.subtropics_change_uniform:+6.1f} cm/yr  (WETTENS — the flaw this fixes)")
    print(f"  global mean ⟨P⟩:    {pbar:.1f} → {float(np.mean(w.precip_wgw)):.1f} cm/yr "
          f"(scales at the energy rate ~{moist.energy_constrained_rate()*100:.1f} %/K, not C–C 7 %/K)")
    print("    => the SPLIT direction is prescribed; what is exact is the mean-zero anomaly + the reductions.\n")


def save_figure(w: moist.WetGetWetter) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    lat = w.phi
    fig, ax = plt.subplots(1, 2, figsize=(13, 5.0))

    # Panel 1 — the precip fields (present vs warmed-uniform vs warmed-split).
    ax[0].plot(lat, w.precip_present, color="0.55", lw=1.4, label="present")
    ax[0].plot(lat, w.precip_uniform, "--", color="C3",
               label=f"+{w.delta_T:.0f} K, rung-0 uniform C–C (deserts wetten)")
    ax[0].plot(lat, w.precip_wgw, "-", color="C0",
               label=f"+{w.delta_T:.0f} K, wet-get-wetter (deserts dry)")
    ax[0].axvspan(20, 35, color="0.85", alpha=0.5, zorder=0)
    ax[0].annotate("subtropical\ndeserts", (27.5, 0.82 * float(w.precip_present.max())),
                   ha="center", fontsize=8, color="0.35")
    ax[0].set_xlim(0, 90)
    ax[0].set_xlabel("latitude (°)"); ax[0].set_ylabel("P (cm/yr)")
    ax[0].set_title("Warming sharpens the pattern: the ITCZ/storm-track\nbands intensify, the deserts dry")
    ax[0].legend(fontsize=8, loc="upper right"); ax[0].grid(alpha=0.3)

    # Panel 2 — the warming response ΔP(φ): uniform positive everywhere vs split negative in the subtropics.
    d_uniform = w.precip_uniform - w.precip_present
    d_wgw = w.precip_wgw - w.precip_present
    ax[1].axhline(0, color="0.6", lw=0.8)
    ax[1].plot(lat, d_uniform, "--", color="C3", label="rung-0 uniform C–C")
    ax[1].plot(lat, d_wgw, "-", color="C0", label="wet-get-wetter")
    ax[1].fill_between(lat, d_wgw, 0, where=(d_wgw < 0), color="C1", alpha=0.35, label="dry-get-drier")
    ax[1].axvspan(20, 35, color="0.85", alpha=0.5, zorder=0)
    ax[1].plot([w.subtropics_lat], [w.subtropics_change_wgw], "o", color="C0", ms=6)
    ax[1].annotate(f"desert: {w.subtropics_change_wgw:+.0f} vs uniform {w.subtropics_change_uniform:+.0f}",
                   (w.subtropics_lat, w.subtropics_change_wgw), textcoords="offset points",
                   xytext=(10, -2), fontsize=8, color="C0")
    # The named overreach: this idealization also dries the poles, where observed high-lat P INCREASES.
    pole = int(np.argmax(lat))
    ax[1].annotate("OVERREACH: dries the poles\n(obs. high-lat P increases —\ntropical/subtropical idealization)",
                   (lat[pole], d_wgw[pole]), textcoords="offset points", xytext=(-148, 62),
                   fontsize=7.5, color="C3", ha="left",
                   arrowprops=dict(arrowstyle="->", color="C3", lw=0.8))
    ax[1].set_xlim(0, 90)
    ax[1].set_xlabel("latitude (°)"); ax[1].set_ylabel("ΔP under warming (cm/yr)")
    ax[1].set_title("The rung-0 uniform amplitude wettens EVERYWHERE;\nwet-get-wetter dries the subtropics")
    ax[1].legend(fontsize=8, loc="upper right"); ax[1].grid(alpha=0.3)

    fig.suptitle("Wet-get-wetter, dry-get-drier (rung 1): the global mean scales at the slow energy rate, "
                 "the wet−dry contrast at the faster C–C rate (Held & Soden 2006)", fontsize=10)
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
    w = compute()
    print_summary(w)
    try:
        saved = save_figure(w)
        print(f"Figure saved → {saved.relative_to(_REPO_ROOT)}")
    except ImportError:
        print("(matplotlib not installed — install the viz extra: pip install -e .[viz])")


if __name__ == "__main__":
    main()
