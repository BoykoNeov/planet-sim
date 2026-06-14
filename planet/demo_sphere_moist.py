"""Rung-2.x evidence: emergent ITCZ rain from a conserving moisture budget, co-located with the EFE.

Rung 2.x's precip wire (:func:`planet.sphere_ebm.itcz_informed_precip`) **relocates a prescribed Gaussian
band** to the energy-flux equator — a dry model painting a rain belt. :mod:`planet.sphere_moist` instead
lets the ITCZ rain **emerge from a conserving ``P − E`` budget** (eddy diffusion + a two-cell mean Hadley
circulation whose ascent is anchored on the EFE), so the rain is *rained*, not *painted*. This demo banks
the headline figure.

**At the honest altitude** (the advisor's framing): the genuinely-new content is largely rung-2.x's EFE
displacement recombined with the rung-2 Hadley convergence into a **conserving budget** (the architectural
"what"). Two findings survive that are not free: (1) the **net** ``P − E`` peaks **on** the EFE only because
the prescribed cell beats the down-gradient eddy *export* there (a falsifiable check, ~2–3× margin); and
(2) the displaced-ITCZ peak intensification is **geometric** (the pinned-edge near cell narrows), **not**
emergent ``q`` (a clean negative result). The wet-NH/dry-SH dipole that accompanies the shift is
displacement-driven and its direction (toward the warm hemisphere) is by-construction; ``∫(P − E) = 0`` to
machine precision pays for the rain belt latitude by latitude.

Three panels:
  1. The emergent budget (eddy + Hadley = net ``P − E``) for a Q-displaced climate, beside the rung-2.x
     **painted** band relocated to the same EFE — *rained* vs *painted*.
  2. Zoom on the EFE: the eddy term EXPORTS, the Hadley term CONVERGES and wins (~2–3×), so the net rain
     belt lands on the EFE — co-location as a check, not an assumption.
  3. The conserving wet/dry dipole vs the imposed displacement (warm hemisphere wets, cold dries, ∫ = 0),
     with the geometric-not-``q`` negative result annotated.

The committed figure + the ``slow`` test let a fresh clone reproduce the headline, not just read it.

Run headless:  python -m planet.demo_sphere_moist
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np

from . import ebm, precip
from . import sphere_moist as sm
from .sphere_ebm import SphereEBM

_REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_FIGURE = _REPO_ROOT / "docs" / "figures" / "planet-sphere-moist.png"
OUTPUT_FIGURE = _REPO_ROOT / "outputs" / "planet-sphere-moist.png"

Q_DRIVE = 6.0       # imposed cross-equatorial energy flux (W m⁻²) that displaces the EFE/ITCZ


def _absorbed_const(x, T):
    """Constant-albedo absorbed shortwave — the splitting-free linear reference forcing."""
    return ebm.insolation(x) * (1.0 - ebm.ALBEDO_A0)


@dataclass(frozen=True)
class SphereMoistResult:
    """The banked evidence — the emergent budget, the co-location margin, the conserving dipole."""

    lat: np.ndarray              # latitude (deg), full sphere
    pme_net: np.ndarray          # net emergent P − E (cm/yr), Q-displaced climate
    pme_eddy: np.ndarray         # eddy component (backwards at the ITCZ)
    pme_hadley: np.ndarray       # mean-Hadley component (anchored on the EFE)
    painted_band: np.ndarray     # the rung-2.x prescribed Gaussian band relocated to the EFE (cm/yr)
    phi_efe: float               # the energy-flux equator (deg)
    rain_max_lat: float          # latitude of the net P − E maximum (co-located with phi_efe)
    eddy_at_efe: float           # eddy P − E at the EFE (< 0, export)
    hadley_at_efe: float         # Hadley P − E at the EFE (> 0, convergence — beats the export)
    net_p_minus_e: float         # ∫(P − E) dx (≈ 0 — conserved)
    drive_grid: np.ndarray       # imposed Q-flux samples (W m⁻²)
    efe_grid: np.ndarray         # the resulting EFE latitudes (deg)
    nh_grid: np.ndarray          # NH ∫P − E per displacement (cm/yr, > 0 — warm hemisphere wets)
    sh_grid: np.ndarray          # SH ∫P − E per displacement (cm/yr, < 0 — cold hemisphere dries)
    peak_real: float             # ITCZ peak with real q (cm/yr)
    peak_symq: float             # ITCZ peak with hemispherically-symmetric q (≈ peak_real → geometric)


def compute(nx: int = 360, q_drive: float = Q_DRIVE) -> SphereMoistResult:
    """Run the symmetric + Q-displaced full-sphere moisture budgets and the displacement sweep."""
    m = SphereEBM(n_cells=nx)
    lat = np.degrees(np.arcsin(np.clip(m.x, -1, 1)))

    c = m.steady_linear(_absorbed_const, Q=q_drive * m.x)
    b = sm.sphere_moisture_budget(c)
    i_efe = int(np.argmin(np.abs(lat - b.phi_efe)))

    # The rung-2.x "painted band" for the same EFE — a prescribed Gaussian relocated, not a budget.
    painted = precip.precip_pattern(lat, itcz_center_deg=b.phi_efe)

    # The geometric-not-q negative result: same cell position, symmetric q → essentially the same peak.
    b_symq = sm.sphere_moisture_budget(replace(c, T=0.5 * (c.T + c.T[::-1])))

    # The conserving dipole vs the imposed displacement.
    drives = np.array([0.0, 2.0, 4.0, 6.0, 8.0])
    efes, nhs, shs = [], [], []
    for qd in drives:
        bj = sm.sphere_moisture_budget(m.steady_linear(_absorbed_const, Q=qd * m.x))
        efes.append(bj.phi_efe); nhs.append(bj.nh_p_minus_e); shs.append(bj.sh_p_minus_e)

    return SphereMoistResult(
        lat=lat, pme_net=b.p_minus_e, pme_eddy=b.p_minus_e_eddy, pme_hadley=b.p_minus_e_hadley,
        painted_band=painted, phi_efe=b.phi_efe, rain_max_lat=b.rain_max_lat,
        eddy_at_efe=float(b.p_minus_e_eddy[i_efe]), hadley_at_efe=float(b.p_minus_e_hadley[i_efe]),
        net_p_minus_e=b.net_p_minus_e, drive_grid=drives, efe_grid=np.array(efes),
        nh_grid=np.array(nhs), sh_grid=np.array(shs),
        peak_real=float(b.p_minus_e.max()), peak_symq=float(b_symq.p_minus_e.max()),
    )


def print_summary(r: SphereMoistResult) -> None:
    print("\nEmergent ITCZ rain — the full-sphere moisture budget co-located with the EFE (rung 2.x)\n")
    print(f"  imposed Q-flux displaces the EFE to {r.phi_efe:+.2f} deg")
    print(f"  net P-E rain belt at {r.rain_max_lat:+.2f} deg  (EFE {r.phi_efe:+.2f}) — co-located, a check")
    print(f"    at the EFE: eddy {r.eddy_at_efe:+.1f} (EXPORT) vs Hadley {r.hadley_at_efe:+.1f} (CONVERGES, "
          f"x{r.hadley_at_efe / abs(r.eddy_at_efe):.1f})")
    print(f"  conserved budget: ∫(P-E) dx = {r.net_p_minus_e:+.2e} cm/yr")
    print(f"  conserving dipole at this displacement: warm(NH) {r.nh_grid[-2]:+.2f} / cold(SH) "
          f"{r.sh_grid[-2]:+.2f} cm/yr")
    print(f"  ITCZ peak: real-q {r.peak_real:.1f} vs symmetric-q {r.peak_symq:.1f} cm/yr "
          f"(=> intensity is GEOMETRIC, not emergent q)\n")


def save_figure(r: SphereMoistResult) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(16.5, 4.7))

    # Panel 1 — rained vs painted.
    ax[0].axhline(0, color="0.6", lw=0.8)
    ax[0].plot(r.lat, r.pme_eddy, ":", color="C3", lw=1.3, label="eddy (exports at ITCZ)")
    ax[0].plot(r.lat, r.pme_hadley, "--", color="C2", lw=1.3, label="mean Hadley (ascent on EFE)")
    ax[0].plot(r.lat, r.pme_net, "-", color="C0", lw=2.0, label="net P − E (emergent rain)")
    ax[0].plot(r.lat, r.painted_band - r.painted_band.min(), "-", color="0.6", lw=1.2,
               label="rung-2.x painted band (precip, not P−E)")
    ax[0].axvline(r.phi_efe, color="r", lw=0.8, alpha=0.6)
    ax[0].set_xlim(-50, 50)
    ax[0].set_xlabel("latitude (°)"); ax[0].set_ylabel("P − E  (cm/yr)")
    ax[0].set_title("Rained, not painted: rain emerges from a\nconserving budget (∫(P−E)=0), not a band")
    ax[0].legend(fontsize=7.5, loc="upper right"); ax[0].grid(alpha=0.3)

    # Panel 2 — co-location as a check (the eddy-vs-Hadley margin at the EFE).
    ax[1].axhline(0, color="0.6", lw=0.8)
    ax[1].plot(r.lat, r.pme_eddy, ":", color="C3", lw=1.5, label="eddy: EXPORT")
    ax[1].plot(r.lat, r.pme_hadley, "--", color="C2", lw=1.5, label="Hadley: CONVERGES")
    ax[1].plot(r.lat, r.pme_net, "-", color="C0", lw=2.2, label="net")
    ax[1].axvline(r.phi_efe, color="r", lw=1.0)
    ax[1].plot([r.rain_max_lat], [r.pme_net.max()], "v", color="C0", ms=9)
    ax[1].annotate(f"EFE = {r.phi_efe:+.2f}°\nrain max = {r.rain_max_lat:+.2f}°",
                   (r.phi_efe, r.pme_net.max()), textcoords="offset points", xytext=(10, -6), fontsize=8)
    ax[1].annotate(f"Hadley beats export\n×{r.hadley_at_efe / abs(r.eddy_at_efe):.1f}",
                   (r.phi_efe, r.hadley_at_efe), textcoords="offset points", xytext=(8, 4),
                   fontsize=8, color="C2")
    ax[1].set_xlim(-12, 12)
    ax[1].set_xlabel("latitude (°)"); ax[1].set_ylabel("P − E  (cm/yr)")
    ax[1].set_title("Co-location is a CHECK: net P − E peaks on the\nEFE because the cell beats the eddy export")
    ax[1].legend(fontsize=8, loc="lower right"); ax[1].grid(alpha=0.3)

    # Panel 3 — the conserving dipole vs displacement + the geometric-not-q negative.
    ax[2].axhline(0, color="0.6", lw=0.8)
    ax[2].plot(r.efe_grid, r.nh_grid, "o-", color="C0", label="warm (NH) ∫P − E  (> 0, wets)")
    ax[2].plot(r.efe_grid, r.sh_grid, "o-", color="C3", label="cold (SH) ∫P − E  (< 0, dries)")
    ax[2].set_xlabel("EFE displacement (°)"); ax[2].set_ylabel("hemispheric ∫(P − E)  (cm/yr)")
    ax[2].set_title("Conserving wet/dry dipole (∫ = 0): direction is\nby-construction; peak intensity is GEOMETRIC, not q")
    ax[2].legend(fontsize=8, loc="upper left"); ax[2].grid(alpha=0.3)
    ax[2].text(0.5, 0.04, f"ITCZ peak  real-q {r.peak_real:.0f}  ≈  sym-q {r.peak_symq:.0f} cm/yr",
               transform=ax[2].transAxes, ha="center", fontsize=8, color="0.3")

    fig.suptitle("Emergent ITCZ rain (rung 2.x): a conserving moisture budget puts the rain belt ON the "
                 "energy-flux equator — rained, not painted; co-location a check, intensity geometric",
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
