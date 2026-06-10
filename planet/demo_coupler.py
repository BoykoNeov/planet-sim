"""The Phase-4 banked artifact: the EBM gradient forces the shallow-water → an emergent jet.

The capstone's closing payoff — and the moment the program's **two shared engines are coupled**.
The frozen diffusion spine produces an equilibrium climate (:mod:`planet.ebm`); its
**meridional temperature gradient** is handed to the frozen rotating shallow-water engine
(:mod:`planet.circulation`); and a **geostrophically-balanced midlatitude westerly jet
emerges** (:mod:`planet.coupler`). *Climate in, circulation out.*

  1. **The forcing chain** — a warm equator / cold pole sets a high/low target height field; the
     coupler relaxes the flow toward it (thermal relaxation + weak drag, split around the bare engine).
  2. **The emergent jet** — rotation turns the unbalanced height field into a **westerly jet** at the
     midlatitude baroclinic zone, *flanked by a compensating easterly return* (the doubly-periodic
     channel's near-zero net zonal momentum). The east–west–east sign banding resembles the general
     circulation, but the single-layer periodic channel does **not** reproduce the observed
     westerly-dominant magnitudes (a named scope edge). The jet sits at the EBM gradient maximum, not
     the channel centre: it is **emergent**.
  3. **Geostrophic balance** — the steady jet satisfies ``f·u ≈ −g·∂h/∂y`` in its core (the anchor).
  4. **Conservation, reframed** — mass is machine-exact under forcing; switch the forcing **off** and
     the bare engine conserves mass / energy / enstrophy *and the jet persists* (the release test).

One-way only (climate → circulation); two-way is rung 1 of the GCM climb (plan §3–4).

Run headless (saves the figure, prints the summary):

    python -m planet.demo_coupler
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .albedo import present_day_climate
from .coupler import CoupledJet, couple_jet
from .ebm import ClimateState

_REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_FIGURE = _REPO_ROOT / "docs" / "figures" / "planet-coupler.png"
OUTPUT_FIGURE = _REPO_ROOT / "outputs" / "planet-coupler.png"
DOCS_GLOBE = _REPO_ROOT / "docs" / "figures" / "planet-coupler-map.html"
OUTPUT_GLOBE = _REPO_ROOT / "outputs" / "planet-coupler-map.html"


@dataclass(frozen=True)
class CouplerResult:
    """The banked Phase-4 artifact — the climate state and the emergent jet it forces."""

    state: ClimateState
    jet: CoupledJet


def compute(nx: int = 96, ny: int = 96) -> CouplerResult:
    """Relax the present-day EBM, force the shallow-water engine with it → :class:`CouplerResult`."""
    state = present_day_climate(n_tau=0.05)
    jet = couple_jet(state=state, nx=nx, ny=ny)
    return CouplerResult(state=state, jet=jet)


def print_summary(r: CouplerResult) -> None:
    """Print the Phase-4 story — the two engines coupled, the jet emergent and balanced."""
    j = r.jet
    print("\nOne-way EBM → shallow-water coupler: the climate gradient forces an emergent jet\n")
    print(f"  EBM climate: global mean T̄ = {r.state.global_mean_T:.2f} °C, ice line {r.state.ice_line_lat:.1f}°")
    print(f"  channel: midlatitude β-plane, deformation radius L_R = {j.L_R/1e3:.0f} km, "
          f"latitudes [{j.phi[0]:.0f}°, {j.phi[-1]:.0f}°]")
    print(f"  spin-up: {'converged' if j.converged else 'capped'} after {j.iterations} steps")
    print("\n  The emergent jet:")
    print(f"    westerly jet  {j.jet_speed:.1f} m/s  @ {j.jet_lat:.1f}°   "
          f"(EBM ∂T/∂φ maximum at {j.gradient_peak_lat:.1f}° — the jet sits there, NOT the channel centre)")
    print(f"    flanking easterly return  {j.u_profile.min():.1f} m/s  "
          f"(the doubly-periodic channel's near-zero net zonal momentum)")
    print(f"    geostrophic balance: core residual {100*j.core_balance_residual:.1f}% "
          f"(f·u ≈ −g ∂h/∂y — the analytic anchor)")
    print("\n  Conservation (reframed — forced–dissipative):")
    print(f"    forced:  mass drift {np.abs(j.mass).max():.1e} (machine-exact);  "
          f"energy {np.abs(j.energy).max():.2f}, enstrophy {np.abs(j.enstrophy).max():.2f} (NOT conserved — "
          f"forcing–drag balance selects the jet)")
    print(f"    RELEASE: forcing off → mass {np.abs(j.mass_release).max():.1e}, "
          f"energy {np.abs(j.energy_release).max():.1e}, enstrophy {np.abs(j.enstrophy_release).max():.1e} "
          f"(engine invariants re-confirmed)")
    print(f"             jet persists: {j.jet_speed:.1f} → {j.u_profile_release.max():.1f} m/s "
          f"(a genuine balanced state, not forcing-propped)\n")


def save_figure(r: CouplerResult) -> Path:
    """Render and save the Phase-4 artifact (needs the optional ``viz`` extra)."""
    import matplotlib
    matplotlib.use("Agg")                            # headless
    from .plots import coupler_figure

    fig = coupler_figure(r.jet, r.state)
    for target in (DOCS_FIGURE, OUTPUT_FIGURE):
        target.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(target, dpi=130)
    return DOCS_FIGURE


def save_circulation_globe(r: CouplerResult) -> Path:
    """Bank the interactive-map artifact: the emergent jet **drawn on the globe** (needs ``[webviz]``).

    The plan's Phase-4 map artifact (§9): the deep-end interactive map *registers* the coupled jet as a
    ``circulation`` ``vector_overlay`` over the temperature field. Recomposes the biome climate
    (:func:`planet.demo_biomes.compute`) into a :class:`~planet.planetmap.PlanetView`
    with the jet overlay and writes a standalone Plotly HTML globe — the seam the renderer paints with no
    restructuring (ADR 0004 #1). Returns the written path.
    """
    from . import demo_biomes, planetmap
    from .albedo import EBMParams

    biome = demo_biomes.compute(EBMParams(), n_tau=0.05)        # the climate the jet was forced from
    view = planetmap.build_view(biome, jet=r.jet)
    for target in (DOCS_GLOBE, OUTPUT_GLOBE):
        planetmap.save_html(view, target, active="temperature")  # circulation cones over T(φ)
    return DOCS_GLOBE


def main() -> None:
    import sys
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")     # °C, ₂, →, ∂, ° on legacy codepages

    r = compute()
    print_summary(r)
    try:
        saved = save_figure(r)
        print(f"Figure saved → {saved.relative_to(_REPO_ROOT)}")
    except ImportError:
        print("(matplotlib not installed — install the viz extra to render the figure: "
              "pip install -e .[viz])")
    try:
        globe = save_circulation_globe(r)
        print(f"Interactive globe (jet over temperature) saved → {globe.relative_to(_REPO_ROOT)}")
    except ImportError:
        print("(plotly not installed — install the webviz extra for the circulation globe: "
              "pip install -e .[webviz])")


if __name__ == "__main__":
    main()
