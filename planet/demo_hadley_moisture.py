"""Rung-2 evidence: the deferred Hadley moisture-convergence fix (the deep-tropical ITCZ).

:mod:`planet.moist`'s default ``P − E`` is pure down-gradient **eddy** diffusion of moisture, which is
**backwards at the moist equator** — it *exports* moisture from the ITCZ (there is no diffusive way to
converge moisture at a maximum). The real ITCZ convergence is the **mean Hadley circulation**: its
low-level, moist branch flows equatorward, carrying water toward the ascent. This demo banks the
before/after figure for the opt-in fix (:func:`planet.moist.hadley_moisture_convergence`).

**At the honest altitude** (the advisor's framing): the convergence-at-the-ITCZ / divergence-in-the-
subtropics *structure* is **by construction** (the cell strength ``HADLEY_STRENGTH`` is the prescribed
sub-grid wall, calibrated to observed *order*, not derived). What is genuinely **emergent** — and bankable
— is the *amplitude*: ``q(T)`` is carried from the EBM, so the ITCZ convergence **intensifies at the ~C–C
moisture rate** under warming (the "rich-get-richer" P−E scaling, Held & Soden 2006), *faster* than the
energy-constrained global mean. And it is a conserving **budget** (``∫(P−E) = 0`` — the ITCZ convergence is
paid for by subtropical divergence), not a painted band. The *fully emergent* mean circulation (a resolved
ascent, not an imposed Ψ) is the gross-moist-stability / overturning route at **rung 3+** — named, not built.

Two panels:
  1. ``P − E(φ)``: the eddy-only default (equator EXPORTS — backwards) vs eddy + Hadley (equator CONVERGES,
     the subtropical desert emerges as the descending branch) — the fix, made visible.
  2. The warming response of the ITCZ convergence at FIXED cell strength: it intensifies at the ~C–C rate
     (emergent, from ``q(T)``), against the slower energy-constrained global-mean rate.

The committed figure + the ``slow`` test assertion (the sign flip, conservation, the C–C-order emergent
rate) let a fresh clone reproduce the headline, not just read the conclusion.

Run headless:  python -m planet.demo_hadley_moisture
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np

from . import moist
from .albedo import EBMParams, present_day_climate

_REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_FIGURE = _REPO_ROOT / "docs" / "figures" / "planet-hadley-moisture.png"
OUTPUT_FIGURE = _REPO_ROOT / "outputs" / "planet-hadley-moisture.png"

WARMING_K = 4.0          # the warming used to expose the emergent C–C-rate intensification


@dataclass(frozen=True)
class HadleyResult:
    """The banked rung-2 Hadley-fix evidence — the before/after budget and the emergent warming rate."""

    lat: np.ndarray              # latitude (deg)
    pme_eddy: np.ndarray         # eddy-only P − E (cm/yr) — the backwards default
    pme_full: np.ndarray         # eddy + Hadley P − E (cm/yr) — the fix
    eq_eddy: float               # equatorial P − E, eddy-only (< 0, export)
    eq_full: float               # equatorial P − E, with Hadley (> 0, convergence)
    subtropics_eddy: float       # canonical 25–35° mean, eddy-only (> 0, mislocated)
    subtropics_full: float       # canonical 25–35° mean, with Hadley (still > 0 — desert NOT relocated)
    dry_belt_min: float          # the Hadley dry-belt minimum P − E equatorward of 35° (< 0)
    dry_belt_lat: float          # latitude of that minimum (deg — equatorward of the canonical subtropics)
    net_full: float              # ∫(P − E) dx with Hadley (≈ 0 — conserved)
    itcz_rate: float             # ITCZ-convergence intensification per K (emergent, ~C–C)
    cc_rate: float               # the C–C moisture rate at the equatorial T (per K)
    energy_rate: float           # the energy-constrained global-mean rate (per K, slower)


def compute() -> HadleyResult:
    """Run the present-day climate's eddy-only and eddy+Hadley moisture budgets and the warming response."""
    st = present_day_climate(EBMParams())
    lat = st.latitude_deg()

    eddy = moist.moisture_budget(st)
    full = moist.moisture_budget(st, hadley=True)

    # The Hadley dry belt (the descending branch) — equatorward of the canonical 25–35° subtropics.
    inner = lat < 35.0
    i_min = int(np.argmin(np.where(inner, full.p_minus_e, np.inf)))

    # The emergent warming response: warm the climate, hold the cell strength FIXED — only q(T) moves.
    warm = replace(st, T=st.T + WARMING_K, global_mean_T=st.global_mean_T + WARMING_K)
    eq_now = moist.hadley_moisture_convergence(st)[0]
    eq_warm = moist.hadley_moisture_convergence(warm)[0]
    itcz_rate = float((eq_warm / eq_now - 1.0) / WARMING_K)
    cc_rate = float(moist.L_VAPOR / (moist.R_VAPOR * (st.T[0] + moist.T0_KELVIN) ** 2))

    return HadleyResult(
        lat=lat, pme_eddy=eddy.p_minus_e, pme_full=full.p_minus_e,
        eq_eddy=eddy.equatorial_export, eq_full=full.equatorial_export,
        subtropics_eddy=eddy.subtropical_balance, subtropics_full=full.subtropical_balance,
        dry_belt_min=float(full.p_minus_e[i_min]), dry_belt_lat=float(lat[i_min]),
        net_full=full.net_p_minus_e, itcz_rate=itcz_rate, cc_rate=cc_rate,
        energy_rate=float(moist.energy_constrained_rate()),
    )


def print_summary(r: HadleyResult) -> None:
    print("\nThe Hadley moisture-convergence fix (rung 2, deep-tropical ITCZ)\n")
    print(f"  equator P-E:     eddy-only {r.eq_eddy:+8.1f} cm/yr (EXPORT, backwards)")
    print(f"                   + Hadley  {r.eq_full:+8.1f} cm/yr (CONVERGENCE — the fix)")
    print(f"  Hadley dry belt:  min {r.dry_belt_min:+7.1f} cm/yr at {r.dry_belt_lat:.0f}° "
          f"(E>P — but EQUATORWARD of the canonical 25-35°)")
    print(f"  canonical 25-35:  eddy-only {r.subtropics_eddy:+7.1f}, + Hadley {r.subtropics_full:+7.1f} cm/yr "
          f"(P>E on BOTH — desert NOT relocated; hyper-peaked q, rung 3+)")
    print(f"  water conserved:  ∫(P-E) = {r.net_full:+.2e} cm/yr")
    print(f"  EMERGENT amplitude: ITCZ convergence intensifies {r.itcz_rate*100:.1f} %/K "
          f"(~C-C {r.cc_rate*100:.1f} %/K)")
    print(f"                      vs the energy-constrained global mean {r.energy_rate*100:.1f} %/K (slower)")
    print(f"    => the cell is PRESCRIBED (strength = the wall); the amplitude is EMERGENT from q(T).\n")


def save_figure(r: HadleyResult) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(13, 5.0))

    ax[0].axhline(0, color="0.6", lw=0.8)
    ax[0].plot(r.lat, r.pme_eddy, "--", color="C3", label="eddy only (default): equator EXPORTS")
    ax[0].plot(r.lat, r.pme_full, "-", color="C0", label="eddy + Hadley (opt-in): equator CONVERGES")
    ax[0].plot([0.0], [r.eq_full], "o", color="C0", ms=6)
    ax[0].plot([0.0], [r.eq_eddy], "o", color="C3", ms=6)
    ax[0].annotate("ITCZ\nconvergence", (0.0, r.eq_full), textcoords="offset points", xytext=(8, -2),
                   fontsize=8, color="C0")
    ax[0].annotate(f"dry belt (E>P)\nat ~{r.dry_belt_lat:.0f}°", (r.dry_belt_lat, r.dry_belt_min),
                   textcoords="offset points", xytext=(6, 6), fontsize=8, color="C0")
    ax[0].axvspan(25, 35, color="0.85", alpha=0.5, zorder=0)
    ax[0].annotate("canonical subtropics\n(25–35°): still P>E\n→ desert not relocated", (30.0, 150.0),
                   ha="center", fontsize=7, color="0.35")
    ax[0].set_xlim(0, 90)
    ax[0].set_xlabel("latitude (°)"); ax[0].set_ylabel("P − E (cm/yr)")
    ax[0].set_title("The deep-tropical fix: the mean Hadley cell\nconverges moisture the eddy diffusion cannot")
    ax[0].legend(fontsize=8, loc="lower right"); ax[0].grid(alpha=0.3)

    rates = np.array([r.energy_rate, r.itcz_rate, r.cc_rate]) * 100.0
    labels = ["energy-constrained\nglobal mean ⟨P⟩", "ITCZ convergence\n(emergent, this work)", "C–C moisture\ncapacity"]
    colors = ["0.6", "C0", "C2"]
    ax[1].bar(labels, rates, color=colors)
    for i, v in enumerate(rates):
        ax[1].text(i, v + 0.1, f"{v:.1f} %/K", ha="center", fontsize=9)
    ax[1].set_ylabel("intensification rate (%/K)")
    ax[1].set_title("The ITCZ convergence intensifies at the ~C–C rate\n(emergent from q(T)) — faster than ⟨P⟩")
    ax[1].grid(alpha=0.3, axis="y")

    fig.suptitle("The Hadley moisture-convergence fix (rung 2): a PRESCRIBED cell with EMERGENT moisture "
                 "content flips the backwards ITCZ — a conserving budget, not a painted band", fontsize=10)
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
