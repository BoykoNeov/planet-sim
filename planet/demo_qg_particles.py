"""§9.6 O5 — the QG condensate as a **particle flow-globe**: the second *emergent* producer.

The ocean rungs O2–O4 stream *real* OSCAR data through the flow-globe; O5 is the independent bonus —
the project's own **rung-3 Phase-B** two-layer QG turbulence (:mod:`planet.baroclinic_qg`) rendered as
particles on the globe. The saturated baroclinic field is an inverse-cascade **condensate** — coherent
vortices and rolled-up potential-vorticity filaments — and here it streams on a real, rotatable 3-D
planet (the same three.js renderer the eddy band and the ocean currents use).

It is the *direct analog of the eddy-particles showcase*: an emergent, single-snapshot, **signed**
scalar field (the upper-layer PV anomaly, coloured on the diverging RdBu_r ramp), a bounded box patch
placed at an illustrative display latitude, honest-by-disclosure. The whole rung-3 win is that unlike
the ~90 %-reversible eddy slosh, this saturated turbulent transport is *genuinely persistent* — so the
disclaimer (fresh, not the eddy's) says so while still naming the idealized-box / not-real-data edges.

Architecturally O5 is the **third producer** through the R1 interchange (:mod:`planet.flow_serialize`),
re-tripping the §9.4 rule-of-three for the globe-geometry helpers — and *re-affirming the hold*: the QG
box carries no latitude metric to share with the eddy band's ``_band_geometry``, so its sector embedding
is inlined in :func:`planet.flow_globe.flow_field_from_qg` (see that function's rule-of-three note).

Run headless (runs a short turbulence sim, saves the HTML, prints the summary)::

    python -m planet.demo_qg_particles
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .baroclinic_qg import QGState, TwoLayerQG
from .demo_baroclinic_qg import PARAMS, US
from .flow_globe import flow_field_from_qg, save_flow_globe_html

_REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_FIGURE = _REPO_ROOT / "docs" / "figures" / "planet-qg-particles.html"
OUTPUT_FIGURE = _REPO_ROOT / "outputs" / "planet-qg-particles.html"

NX = 96                   # matches the banked rung-3 turbulence figure (demo_baroclinic_qg)
N_LAM = 3                 # box = 3 most-unstable wavelengths
R_FAC = 0.5               # weak bottom drag — the clearest inverse-cascade condensate
N_EFOLD = 20.0            # e-folds of linear growth to run (well past saturation)


@dataclass(frozen=True)
class QGParticleResult:
    """One saturated turbulent run: the model + final state the producer consumes, plus a scale summary."""

    model: TwoLayerQG
    state: QGState
    L: float                 # box size (m)
    Ld: float                # deformation radius (m)
    vrms: float              # eddy v rms (m/s) — ≫ U_s ⟹ an inverse-cascade condensate


def compute(nx: int = NX, n_lam: int = N_LAM, r_fac: float = R_FAC,
            n_efold: float = N_EFOLD, seed: int = 0) -> QGParticleResult:
    """Run the unstable shear to saturated turbulence and return the model + final PV state.

    Mirrors :func:`planet.demo_baroclinic_qg.compute` (same idealized params, box sizing, dissipation)
    but keeps the **model** (the renderer needs the grid metric to embed the box) and the final
    :class:`~planet.baroclinic_qg.QGState`, rather than only the banked diagnostics.
    """
    disp = TwoLayerQG.symmetric(8, 8, 1.0, 1.0, Us=US, **PARAMS)
    kstar, sig = disp.most_unstable()
    L = n_lam * 2 * np.pi / kstar
    dx = L / nx
    m = TwoLayerQG.symmetric(nx, nx, L, L, Us=US, nu4=0.1 * US * dx ** 3,
                             r_drag=r_fac * sig, **PARAMS)
    s = m.random_state(amplitude=1e-3, seed=seed)
    t_end = n_efold / sig
    t = 0.0
    while t < t_end:
        dt = m.max_dt(s, 0.3)
        s = m.step(s, dt)
        t += dt
    return QGParticleResult(model=m, state=s, L=L, Ld=m.Ld, vrms=m.v_rms(s))


def print_summary(r: QGParticleResult) -> None:
    """Print the box scales and the condensate signature (the showcase's one-line story)."""
    print("\nTwo-layer QG condensate — the saturated baroclinic eddy field as a particle globe (§9.6 O5)\n")
    print(f"  box {r.L / 1e3:.0f} km, L_d {r.Ld / 1e3:.0f} km")
    print(f"  v'_rms = {r.vrms:.1f} m/s  (≫ U_s = {US} → an inverse-cascade condensate)")
    print("  colour = upper-layer PV anomaly (the vortex-filament field); placed at an illustrative 45° N\n")


def save_globe(r: QGParticleResult) -> Path:
    """Build the generic flow field from the QG state and write the O5 HTML showcase (PV colour, RdBu_r)."""
    field = flow_field_from_qg(r.model, r.state)
    for target in (DOCS_FIGURE, OUTPUT_FIGURE):
        save_flow_globe_html(
            field, target,
            title="planet-sim — QG turbulence flow-globe (showcase)",
            subtitle="the emergent rung-3 two-layer QG condensate — coherent vortices and PV filaments — "
                     "as a particle flow on a rotatable planet (§9.6 O5, an idealized box, not the sea)",
        )
    return DOCS_FIGURE


def main() -> None:
    import sys
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")     # °, ₂, ≫, →, κ, β on legacy codepages

    r = compute()
    print_summary(r)
    saved = save_globe(r)
    print(f"Particle flow-globe saved → {saved.relative_to(_REPO_ROOT)}")
    print("  open it in a browser (works straight off disk — three.js is vendored inline).")


if __name__ == "__main__":
    main()
