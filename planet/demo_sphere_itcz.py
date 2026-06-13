"""Rung-2.x evidence: the full-sphere energy-flux equator and the energetic ITCZ migration.

:mod:`planet.sphere_ebm` lifts rung 0 to the full sphere and locates the ITCZ at the **energy-flux
equator** (the zero of the atmospheric energy transport ``H(x)``). This demo banks the headline figure —
**at the honest, lower altitude** the rung is banked: the ITCZ-migration sensitivity is a **closed-form
consequence of the already-calibrated transport ``D`` and the mean-state curvature**, of the observed
*order* (~3 deg/PW) but a **factor ~1.5–2 high** — a corroboration that ``D`` is realistic, **not** a
prediction of the ITCZ. The shift *direction* (toward the warm hemisphere) is by-construction; the
forcing-independence is a linear-operator identity.

Three panels:
  1. ``H(x)`` for a symmetric climate (EFE at the equator) and a climate with an imposed cross-equatorial
     Q-flux (EFE displaced into the warmed hemisphere) — the migration, made visible.
  2. ``φ_EFE`` vs the imposed cross-equatorial energy transport: the engine points fall on the **closed-form
     line** ``δ = AHT/(2π a² D T̄ₓₓ(0))``, beside the observed ``~3 deg/PW`` reference slope (the factor ~2).
  3. The precipitation pattern with the ITCZ band relocated to ``φ_EFE`` (the opt-in `itcz_center_deg`
     wire) — a *dry* model **relocating a prescribed band**, not emergent rainfall.

The committed figure + the ``slow`` test assertion (the migration slope matches the closed form, of the
observed order and the right sign) let a fresh clone reproduce the headline, not just read the conclusion.

Run headless:  python -m planet.demo_sphere_itcz
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import ebm, precip
from .sphere_ebm import SphereEBM, itcz_sensitivity_closed

_REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_FIGURE = _REPO_ROOT / "docs" / "figures" / "planet-sphere-itcz.png"
OUTPUT_FIGURE = _REPO_ROOT / "outputs" / "planet-sphere-itcz.png"

OBSERVED_DEG_PER_PW = -3.0          # Donohoe et al. 2013 (the cross-equatorial-AHT ITCZ sensitivity)


def _absorbed_const(x, T):
    """Constant-albedo absorbed shortwave — the splitting-free linear reference forcing."""
    return ebm.insolation(x) * (1.0 - ebm.ALBEDO_A0)


@dataclass(frozen=True)
class ITCZResult:
    """The banked rung-2.x evidence — the EFE migration, its closed-form sensitivity, and the precip wire."""

    lat: np.ndarray              # latitude (deg), full sphere
    H_sym: np.ndarray            # atmospheric transport H(x) for the symmetric climate (PW)
    H_asym: np.ndarray           # H(x) for the Q-forced climate (PW)
    phi_efe_asym: float          # the displaced EFE/ITCZ latitude (deg)
    aht_grid: np.ndarray         # imposed cross-equatorial AHT samples (PW)
    phi_grid: np.ndarray         # the resulting EFE latitudes (deg)
    slope: float                 # measured deg/PW
    slope_closed: float          # closed-form deg/PW (1/(2π a² D Txx0))
    precip_sym: np.ndarray       # rung-0 (symmetric) precip pattern (cm/yr)
    precip_shift: np.ndarray     # precip with the ITCZ relocated to phi_efe_asym (cm/yr)


def compute(nx: int = 360, q_drive: float = 4.0) -> ITCZResult:
    """Run the symmetric + Q-forced climates and the migration sweep (all splitting-free `steady_linear`)."""
    m = SphereEBM(n_cells=nx)
    lat = np.degrees(np.arcsin(np.clip(m.x, -1, 1)))

    c_sym = m.steady_linear(_absorbed_const)
    c_asym = m.steady_linear(_absorbed_const, Q=q_drive * m.x)

    q_amps = np.array([0.0, 1.0, 2.0, 4.0, 6.0])
    ahts, phis = [], []
    for q in q_amps:
        c = m.steady_linear(_absorbed_const, Q=q * m.x)
        ahts.append(c.aht_eq); phis.append(c.phi_efe)
    ahts, phis = np.array(ahts), np.array(phis)
    slope = float(np.polyfit(ahts, phis, 1)[0])
    slope_closed = itcz_sensitivity_closed(m.equatorial_curvature(c_sym.T), m.D)

    precip_sym = precip.precipitation(lat, c_sym.global_mean_T, itcz_center_deg=0.0)
    precip_shift = precip.precipitation(lat, c_asym.global_mean_T, itcz_center_deg=c_asym.phi_efe)

    return ITCZResult(
        lat=lat, H_sym=m.atmospheric_transport(c_sym.T), H_asym=m.atmospheric_transport(c_asym.T),
        phi_efe_asym=c_asym.phi_efe, aht_grid=ahts, phi_grid=phis, slope=slope, slope_closed=slope_closed,
        precip_sym=precip_sym, precip_shift=precip_shift,
    )


def print_summary(r: ITCZResult) -> None:
    print("\nFull-sphere EBM + the energetic ITCZ (rung 2.x)\n")
    print(f"  symmetric climate:  EFE at the equator (H crosses 0 at 0°)")
    print(f"  Q-forced climate:   EFE displaced to {r.phi_efe_asym:+.2f}°  (toward the warmed hemisphere)")
    print(f"  ITCZ sensitivity:   measured {r.slope:.2f} deg/PW,  closed form {r.slope_closed:.2f} deg/PW")
    print(f"    => a CLOSED-FORM consequence of the calibrated D (not a prediction); same ORDER as the")
    print(f"       observed ~{abs(OBSERVED_DEG_PER_PW):.0f} deg/PW but a factor ~{abs(r.slope)/abs(OBSERVED_DEG_PER_PW):.1f} high.")
    print(f"  precip wire: ITCZ rain belt relocated to {r.phi_efe_asym:+.2f}° (a prescribed band moved, "
          f"NOT emergent rain)\n")


def save_figure(r: ITCZResult) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))

    ax[0].axhline(0, color="0.6", lw=0.8); ax[0].axvline(0, color="0.85", lw=0.8)
    ax[0].plot(r.lat, r.H_sym, label="symmetric (EFE = 0°)")
    ax[0].plot(r.lat, r.H_asym, label=f"Q-forced (EFE = {r.phi_efe_asym:+.1f}°)")
    ax[0].plot([r.phi_efe_asym], [0], "ro", ms=6)
    ax[0].set_xlabel("latitude (°)"); ax[0].set_ylabel("northward atmospheric transport H (PW)")
    ax[0].set_title("Energy transport H(x): the ITCZ sits where H = 0"); ax[0].legend(); ax[0].grid(alpha=0.3)

    ax[1].plot(r.aht_grid, r.phi_grid, "o", ms=5, label="engine (EFE)")
    aline = np.linspace(r.aht_grid.min(), r.aht_grid.max(), 50)
    ax[1].plot(aline, r.slope_closed * aline, "-", label=f"closed form  {r.slope_closed:.1f}°/PW")
    ax[1].plot(aline, OBSERVED_DEG_PER_PW * aline, "--", color="0.5",
               label=f"observed ~{OBSERVED_DEG_PER_PW:.0f}°/PW (Donohoe 2013)")
    ax[1].set_xlabel("cross-equatorial energy transport AHT_eq (PW)")
    ax[1].set_ylabel("ITCZ latitude φ_EFE (°)")
    ax[1].set_title(f"Migration sensitivity = {r.slope:.1f}°/PW\n(closed-form consequence of D — a factor "
                    f"~{abs(r.slope)/abs(OBSERVED_DEG_PER_PW):.1f} above observed)")
    ax[1].legend(); ax[1].grid(alpha=0.3)

    ax[2].plot(r.lat, r.precip_sym, label="symmetric ITCZ (rung 0)")
    ax[2].plot(r.lat, r.precip_shift, label=f"ITCZ → {r.phi_efe_asym:+.1f}° (wired)")
    ax[2].set_xlim(-40, 40)
    ax[2].set_xlabel("latitude (°)"); ax[2].set_ylabel("precipitation (cm/yr)")
    ax[2].set_title("Precip wire: a DRY model relocates a\nprescribed rain band (not emergent rainfall)")
    ax[2].legend(); ax[2].grid(alpha=0.3)

    fig.suptitle("Full-sphere EBM + the energetic ITCZ (rung 2.x): the ITCZ migrates with the energy-flux "
                 "equator — a closed-form consequence of the calibrated transport D, of the observed order",
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
