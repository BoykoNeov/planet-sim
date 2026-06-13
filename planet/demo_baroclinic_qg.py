"""Rung-3 Phase-B evidence: the saturated two-layer QG state is genuine TURBULENCE (the bet's headline).

:mod:`planet.baroclinic_qg` produces the saturated baroclinic eddy thickness flux that makes the
rung-1 reduction-to-EBM non-vacuous. But — the load-bearing subtlety — a *down-gradient,
irreversible* flux is guaranteed for **any** sustained baroclinic state (it is the APE→eddy energy
conversion), so the flux sign alone cannot distinguish developed turbulence from a quasi-steady
finite-amplitude wave. This demo banks the **decisive evidence** that it is turbulence:

  1. **The PV fields** — coherent vortices, rolled-up filaments, and vortex stripping across scales
     (a turbulent field, not a clean wave train).
  2. **The isotropic KE spectrum** — the energy peak has migrated **below the injection wavenumber
     ``k*``** (an **inverse energy cascade** to the box/Rhines scale, ``v'_rms ≫ U_s``), broadband and
     continuous. A steady wave would peak *at* ``k*`` with discrete harmonics; the upscale peak is the
     signature only turbulence produces.

This is the figure the slow-test assertion (``k_peak < k*``) guards in CI — committed so a fresh
clone can reproduce the headline, not just read the conclusion. The dimensional ``κ`` is *not* the
banked quantity (it is box/drag/resolution-dependent — the magnitude shifts ~45 % nx96→128); the
banked claim is **dimensionless + qualitative** (down-gradient, irreversible, turbulent, ``κ/(v'L_d)
~ O(1)`` vs rung-1's ~1e-3).

Run headless (saves the figure, prints the summary):

    python -m planet.demo_baroclinic_qg
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .baroclinic_qg import TwoLayerQG

_REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_FIGURE = _REPO_ROOT / "docs" / "figures" / "planet-baroclinic-qg-turbulence.png"
OUTPUT_FIGURE = _REPO_ROOT / "outputs" / "planet-baroclinic-qg-turbulence.png"

# Idealized rung-3 parameters (resolvable L_d, modest speeds — the saturation spike's).
PARAMS = dict(f0=1.0e-4, gp=2.0, H1=400.0, H2=400.0, beta=1.6e-11)
US = 4.0


@dataclass(frozen=True)
class QGTurbulenceResult:
    """The banked rung-3 Phase-B evidence — one saturated turbulent run, its diagnostics + figure data."""

    q: np.ndarray            # final PV snapshot (2, ny, nx)
    K_over_kstar: np.ndarray  # spectrum wavenumbers / k*
    E_norm: np.ndarray       # time-mean KE spectrum, normalised to its peak
    L: float                 # box size (m)
    kstar: float             # most-unstable (injection) wavenumber (1/m)
    Ld: float                # deformation radius (m)
    kappa: float             # emergent eddy thickness diffusivity (m^2/s) — NOT banked (dimensional)
    vrms: float              # eddy v rms (m/s)
    irr: float               # irreversible fraction
    ratio: float             # kappa / (v'_rms * L_d) — the dimensionless banked quantity
    condensate_frac: float   # fraction of KE below 0.7 k* (the inverse-cascade condensate)
    k_peak_over_kstar: float  # spectral peak / k* (< 1 ⟺ inverse cascade ⟺ turbulence)
    r_fac: float             # bottom-drag fraction of sigma


def compute(nx: int = 96, n_lam: int = 3, r_fac: float = 0.5,
            n_efold_total: float = 20.0, n_efold_avg: float = 8.0,
            seed: int = 0) -> QGTurbulenceResult:
    """Run the unstable shear to saturated turbulence (with dissipation), banking the time-mean KE
    spectrum, the final PV snapshot, and the flux diagnostics. Defaults to the weak-drag (``r=0.5σ``)
    condensate case — the clearest turbulence; a smaller ``nx``/fewer e-folds is the smoke-test path.
    """
    disp = TwoLayerQG.symmetric(8, 8, 1.0, 1.0, Us=US, **PARAMS)
    kstar, sig = disp.most_unstable()
    L = n_lam * 2 * np.pi / kstar
    dx = L / nx
    m = TwoLayerQG.symmetric(nx, nx, L, L, Us=US, nu4=0.1 * US * dx ** 3,
                             r_drag=r_fac * sig, **PARAMS)
    s = m.random_state(amplitude=1e-3, seed=seed)
    t_end, t_avg = n_efold_total / sig, (n_efold_total - n_efold_avg) / sig
    t, n = 0.0, 0
    fluxes, vrms, spec, Kc = [], [], None, None
    while t < t_end:
        dt = m.max_dt(s, 0.3)
        s = m.step(s, dt)
        t += dt
        n += 1
        if t >= t_avg and n % 20 == 0:
            f1, f2 = m.bulk_eddy_flux(s)
            fluxes.append(0.5 * (f1 + f2))
            vrms.append(m.v_rms(s))
            Kc, E = m.ke_spectrum(s)
            spec = E if spec is None else spec + E
    fluxes = np.array(fluxes)
    kappa = float(fluxes.mean() / m.Us)
    vr = float(np.mean(vrms))
    irr = float(abs(fluxes.mean()) / np.abs(fluxes).mean())
    Kn = Kc / kstar
    Esum = spec.sum()
    return QGTurbulenceResult(
        q=s.q, K_over_kstar=Kn, E_norm=spec / spec.max(), L=L, kstar=kstar, Ld=m.Ld,
        kappa=kappa, vrms=vr, irr=irr, ratio=kappa / (vr * m.Ld),
        condensate_frac=float(spec[Kn < 0.7].sum() / Esum),
        k_peak_over_kstar=float(Kn[np.argmax(spec)]), r_fac=r_fac,
    )


def print_summary(r: QGTurbulenceResult) -> None:
    """Print the rung-3 Phase-B story — saturated turbulence with a non-vacuous down-gradient flux."""
    print("\nTwo-layer QG turbulence — the saturated baroclinic eddy flux (rung 3 Phase B)\n")
    print(f"  box {r.L/1e3:.0f} km, L_d {r.Ld/1e3:.0f} km, bottom drag r = {r.r_fac:g}σ")
    print(f"  v'_rms = {r.vrms:.1f} m/s  (≫ U_s = {US} → an inverse-cascade condensate)")
    print(f"  KE spectrum peak at K/k* = {r.k_peak_over_kstar:.2f}  "
          f"({'< 1 ⟹ inverse cascade ⟹ TURBULENCE' if r.k_peak_over_kstar < 1 else 'AT k* ⟹ wave-like'})")
    print(f"  condensate fraction (KE below 0.7 k*) = {r.condensate_frac:.2f}")
    print(f"  emergent flux: down-gradient (κ = {r.kappa:.2e} m²/s, NOT banked — dimensional),")
    print(f"    irreversible fraction = {r.irr:.2f}  (vs rung-1's ~0.1 reversible barotropic flux)")
    print(f"  κ/(v'_rms·L_d) = {r.ratio:.2f}  (the banked DIMENSIONLESS efficiency; rung-1 ~1e-3)\n")


def save_figure(r: QGTurbulenceResult) -> Path:
    """Render and save the PV-fields + KE-spectrum figure (needs the optional ``viz`` extra)."""
    import matplotlib
    matplotlib.use("Agg")                                # headless
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.7))
    Lkm = r.L / 1e3
    for k, lab in [(0, "upper-layer PV  q₁"), (1, "lower-layer PV  q₂")]:
        im = ax[k].imshow(r.q[k], origin="lower", cmap="RdBu_r", extent=[0, Lkm, 0, Lkm])
        ax[k].set_title(f"{lab}\n(coherent vortices + filaments = turbulence)")
        ax[k].set_xlabel("x (km)")
        ax[k].set_ylabel("y (km)")
        fig.colorbar(im, ax=ax[k], fraction=0.046, pad=0.04)
    ax[2].loglog(r.K_over_kstar[1:], r.E_norm[1:], "-o", ms=3)
    ax[2].axvline(1.0, color="k", ls="--", lw=0.9, label="injection k*")
    ax[2].set_title(f"isotropic KE spectrum\n(peak at K/k*={r.k_peak_over_kstar:.2f} < 1 = inverse cascade)")
    ax[2].set_xlabel("K / k*  (wavenumber relative to the most-unstable mode)")
    ax[2].set_ylabel("E(K) / E_max")
    ax[2].legend()
    ax[2].grid(alpha=0.3, which="both")
    fig.suptitle("Two-layer QG turbulence (rung 3 Phase B): the energy peak migrates BELOW the "
                 "injection scale — an inverse cascade only turbulence makes, not a steady wave",
                 fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    for target in (DOCS_FIGURE, OUTPUT_FIGURE):
        target.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(target, dpi=110)
    plt.close("all")
    return DOCS_FIGURE


def main() -> None:
    import sys
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")         # °, ₂, ≫, ⟹, ² on legacy codepages

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
