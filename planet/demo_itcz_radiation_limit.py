"""Rung-2.x refinement: why NO transport tightens the ITCZ sensitivity — it is a RADIATION limit.

The backlog asked to *tighten* the full-sphere ITCZ-migration sensitivity (:mod:`planet.sphere_ebm`) — the
dry EBM gives ``≈ −6.3 deg/PW``, a factor ~1.5–2 above the observed ``~3`` (Donohoe et al. 2013) — by
**re-deriving the transport ``D``**. This demo banks the answer, which is a **negative with a mechanism**:

At the symmetric steady state the equatorial energy balance pins ``D·T̄ₓₓ(0) = −NEI(0)``, so the sensitivity
is *identically*

    δ / AHT_eq = −1 / (2π a² · NEI(0))            (Bischoff & Schneider 2014)

a function of the **net radiative input at the equator ``NEI(0)``** — a *radiation* quantity that the
transport ``D`` cancels out of. Turning ``D`` only slides the equatorial temperature (and hence ``NEI(0)``)
**along** this one curve: from the radiative-equilibrium limit (``D → 0``, ``NEI → 0``, infinite
sensitivity) up to the **isothermal ceiling** (``D → ∞``, equator at the global mean, ``NEI ≈ 57 W/m²``,
sensitivity floored at ``≈ −3.9 deg/PW``). The moist-MSE upgrade (:mod:`planet.sphere_moist_ebm`) is *also*
a transport intervention — it only nudges ``NEI(0)`` ~10 % via the ~1.7 K cooler moist equator, landing at
``≈ −5.7``. **Observed ``−3 deg/PW`` needs ``NEI(0) ≈ 75 W/m²`` — beyond even the isothermal ceiling — so no
transport, dry or moist, reaches it.** The lever is a stronger equatorial radiative surplus (better
radiation → rung 4), or the gross-moist-stability dynamics the diffusive closure omits (rung 3+).

Two panels:
  1. The master curve ``deg/PW = −1/(2π a² NEI(0))``: the transport-reachable band ``NEI ∈ (0, 57]`` shaded,
     the dry / moist / isothermal-ceiling points on it, and the observed target off the reachable end.
  2. The cancellation that hides the identity: dry → moist, ``D_eff(0)`` rises but ``|T̄ₓₓ(0)|`` flattens in
     lockstep, so their product ``|D_eff·T̄ₓₓ| = NEI(0)`` barely moves — the sensitivity with it.

Run headless:  python -m planet.demo_itcz_radiation_limit
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import ebm
from .sphere_ebm import AREA_FACTOR, PW, SphereEBM, itcz_sensitivity_from_nei
from . import sphere_moist_ebm as sm

_REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_FIGURE = _REPO_ROOT / "docs" / "figures" / "planet-itcz-radiation-limit.png"
OUTPUT_FIGURE = _REPO_ROOT / "outputs" / "planet-itcz-radiation-limit.png"

OBSERVED_DEG_PER_PW = -3.0          # Donohoe et al. 2013


def _absorbed(x, T):
    return ebm.insolation(x) * (1.0 - ebm.ALBEDO_A0)


def _nei_for(deg_per_pw: float) -> float:
    """Invert the identity: the NEI(0) that a given deg/PW sensitivity requires (W m⁻²)."""
    return -(180.0 / math.pi) * PW / (AREA_FACTOR * deg_per_pw)


@dataclass(frozen=True)
class RadiationLimit:
    """The banked evidence: the sensitivity is a function of NEI(0), transport-floored above observed."""

    nei_dry: float
    nei_moist: float
    nei_ceiling: float          # isothermal (D→∞) limit
    nei_observed: float         # what observed −3 deg/PW would require
    slope_dry: float
    slope_moist: float
    slope_ceiling: float
    d_sweep_nei: np.ndarray     # NEI(0) as D is turned up (transport slides you along the curve)
    d_sweep_slope: np.ndarray
    Deff0_dry: float
    Deff0_moist: float
    Txx_dry: float
    Txx_moist: float


def compute(n_cells: int = 360) -> RadiationLimit:
    r = sm.moist_vs_dry_sensitivity(n_cells=n_cells)

    # isothermal (D→∞) ceiling: equator relaxes to the global mean ⟹ NEI(0) at its maximum.
    m = SphereEBM(n_cells=n_cells)
    absorbed = _absorbed(m.x, np.zeros_like(m.x))
    Tbar = (m.global_mean(absorbed) - m.A) / m.B
    absorbed0 = float(np.interp(0.0, m.x, absorbed))
    nei_ceiling = absorbed0 - m.A - m.B * Tbar
    slope_ceiling = itcz_sensitivity_from_nei(nei_ceiling)

    # the dry D-sweep: turning D up slides NEI(0) up the curve toward the ceiling, never to observed.
    d_facs = np.array([0.25, 0.5, 1.0, 2.0, 5.0, 20.0, 100.0])
    nei_sweep, slope_sweep = [], []
    for f in d_facs:
        md = SphereEBM(n_cells=n_cells, D=ebm.D_TRANSPORT * f)
        c = md.steady_linear(_absorbed)
        nei = md.net_radiative_input_equator(_absorbed, c.T)
        nei_sweep.append(nei); slope_sweep.append(itcz_sensitivity_from_nei(nei))

    Deff0_dry = ebm.D_TRANSPORT
    Deff0_moist = float(sm.effective_diffusivity(r.Teq_moist, r.D_s, r.RH))
    cdry = SphereEBM(n_cells=n_cells).steady_linear(_absorbed)
    Txx_dry = SphereEBM(n_cells=n_cells).equatorial_curvature(cdry.T)

    return RadiationLimit(
        nei_dry=r.nei_dry, nei_moist=r.nei_moist, nei_ceiling=float(nei_ceiling),
        nei_observed=_nei_for(OBSERVED_DEG_PER_PW),
        slope_dry=r.slope_dry, slope_moist=r.slope_moist, slope_ceiling=float(slope_ceiling),
        d_sweep_nei=np.array(nei_sweep), d_sweep_slope=np.array(slope_sweep),
        Deff0_dry=Deff0_dry, Deff0_moist=Deff0_moist, Txx_dry=float(Txx_dry), Txx_moist=r.Txx_moist,
    )


def print_summary(r: RadiationLimit) -> None:
    print("\nThe ITCZ sensitivity is a RADIATION limit — no transport tightens it (rung 2.x)\n")
    print("  sensitivity  δ/AHT = −1/(2π a² NEI(0))   [Bischoff & Schneider 2014]\n")
    print(f"  dry EBM:            NEI(0) = {r.nei_dry:5.1f} W/m²  →  {r.slope_dry:5.2f} deg/PW")
    print(f"  moist MSE upgrade:  NEI(0) = {r.nei_moist:5.1f} W/m²  →  {r.slope_moist:5.2f} deg/PW  "
          f"(transport moved it only ~{abs((r.slope_moist-r.slope_dry)/r.slope_dry)*100:.0f}%)")
    print(f"  isothermal ceiling: NEI(0) = {r.nei_ceiling:5.1f} W/m²  →  {r.slope_ceiling:5.2f} deg/PW  "
          f"(the D→∞ FLOOR)")
    print(f"  observed ~3 deg/PW  NEEDS   NEI(0) = {r.nei_observed:5.1f} W/m²  — ABOVE the ceiling ⟹ "
          f"unreachable by any transport")
    print(f"\n  the cancellation:  D_eff(0) {r.Deff0_dry:.2f} → {r.Deff0_moist:.2f} (×{r.Deff0_moist/r.Deff0_dry:.2f}),  "
          f"|T̄ₓₓ(0)| {abs(r.Txx_dry):.0f} → {abs(r.Txx_moist):.0f} (×{abs(r.Txx_moist/r.Txx_dry):.2f}),  "
          f"product ≈ const = NEI\n")


def save_figure(r: RadiationLimit) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(13.5, 5.2))

    # -- panel 1: the master curve deg/PW vs NEI(0) -------------------------- #
    nei = np.linspace(6.0, 90.0, 400)
    curve = -(180.0 / math.pi) * PW / (AREA_FACTOR * nei)      # the identity, vectorized for the curve
    ax[0].axvspan(0.0, r.nei_ceiling, color="#9ecae1", alpha=0.30,
                  label="reachable by transport  (NEI ≤ isothermal ceiling)")
    ax[0].plot(nei, curve, color="0.25", lw=2, label="δ/AHT = −1/(2π a² NEI)")
    ax[0].plot([r.nei_dry], [r.slope_dry], "o", ms=9, color="#d95f02", label=f"dry EBM ({r.slope_dry:.1f}°/PW)")
    ax[0].plot([r.nei_moist], [r.slope_moist], "s", ms=9, color="#1b9e77",
               label=f"moist MSE ({r.slope_moist:.1f}°/PW)")
    ax[0].plot([r.nei_ceiling], [r.slope_ceiling], "D", ms=8, color="#7570b3",
               label=f"isothermal floor ({r.slope_ceiling:.1f}°/PW)")
    ax[0].plot(r.d_sweep_nei, r.d_sweep_slope, ".", color="0.45", ms=7,
               label="turning D up (slides along ↑)")
    ax[0].axhline(OBSERVED_DEG_PER_PW, color="#e7298a", ls="--", lw=1.4)
    ax[0].plot([r.nei_observed], [OBSERVED_DEG_PER_PW], "*", ms=15, color="#e7298a",
               label=f"observed ~3°/PW needs NEI={r.nei_observed:.0f} (off-range)")
    ax[0].annotate("no transport reaches here\n→ needs stronger equatorial\nradiative surplus (rung 4)",
                   xy=(r.nei_observed, OBSERVED_DEG_PER_PW), xytext=(r.nei_observed - 2, -2.1),
                   ha="right", va="center", fontsize=8.5, color="#e7298a")
    ax[0].set_xlim(0, 90); ax[0].set_ylim(-9, -1.5)
    ax[0].set_xlabel("equatorial net radiative input  NEI(0)  (W m⁻²)")
    ax[0].set_ylabel("ITCZ-migration sensitivity  (deg / PW)")
    ax[0].set_title("The sensitivity is a RADIATION limit, not a transport one:\n"
                    "transport only slides you along the curve — never to observed", fontsize=10)
    ax[0].legend(fontsize=7.6, loc="lower right"); ax[0].grid(alpha=0.3)

    # -- panel 2: the cancellation that hides the identity ------------------- #
    groups = ["D_eff(0)\n(W/m²/K)", "|T̄ₓₓ(0)|\n(°C)", "|D_eff·T̄ₓₓ| = NEI(0)\n(W/m²)", "|sensitivity|\n(deg/PW)"]
    dry_vals = [r.Deff0_dry, abs(r.Txx_dry), abs(r.nei_dry), abs(r.slope_dry)]
    moist_vals = [r.Deff0_moist, abs(r.Txx_moist), abs(r.nei_moist), abs(r.slope_moist)]
    # normalize each pair to the dry value so the (non-)changes are visible on one axis.
    dry_norm = [1.0 for _ in dry_vals]
    moist_norm = [m / d for m, d in zip(moist_vals, dry_vals)]
    xi = np.arange(len(groups)); w = 0.36
    ax[1].axhline(1.0, color="0.7", lw=0.8)
    ax[1].bar(xi - w / 2, dry_norm, w, color="#d95f02", label="dry")
    ax[1].bar(xi + w / 2, moist_norm, w, color="#1b9e77", label="moist")
    for i, (mv, dv) in enumerate(zip(moist_vals, dry_vals)):
        ax[1].text(xi[i] + w / 2, mv / dv + 0.02, f"×{mv/dv:.2f}", ha="center", fontsize=8.5)
    ax[1].set_xticks(xi); ax[1].set_xticklabels(groups, fontsize=8.2)
    ax[1].set_ylabel("moist ÷ dry")
    ax[1].set_title("Why moisture barely moves it: D_eff(0) rises but |T̄ₓₓ(0)|\n"
                    "flattens in lockstep → their product (= NEI) is pinned", fontsize=10)
    ax[1].legend(fontsize=8.5); ax[1].grid(alpha=0.3, axis="y")

    fig.suptitle("Tightening the ITCZ sensitivity? The identity δ/AHT = −1/(2π a² NEI(0)) says no transport "
                 "can — it is floored above observed by the isothermal radiative surplus", fontsize=10.5)
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
