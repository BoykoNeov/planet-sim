"""The rung-A banked artifact: the emergent eddy life cycle, **animated** (the mechanism, §9.5).

Rung 1's Phase-B module (:mod:`planet.eddy_flux`) released a barotropically-unstable jet and let a
passive temperature tracer stir — diagnosing an *emergent* eddy diffusivity ``κ_eff``. That life
cycle is the only genuinely time-varying, longitudinally-structured 2-D flow the project produces.
This demo **animates** it (the user's forward requirement, referencing NASA's *Perpetual Ocean* /
Ventusky as *broad* visual references) — the program's **first time-animation primitive**.

It is deliberately the **honest** rung (plan §9.5: build A first). Two panels:

  1. **The stirring** — the tracer ``θ`` advected by the released eddies on the midlatitude β-plane
     channel, with the eddy velocity overlaid. A *band*, not a globe (the named scope edge).
  2. **The transport budget** — the cumulative meridional flux: the **throughput** ``Σ∫|F̄|dt`` raging
     upward while the **net** ``Σ|∫F̄dt|`` stays small. This makes :mod:`planet.eddy_flux`'s headline
     finding — the instantaneous flux is **~90 % reversible** — *visible*: without it, a stirring movie
     would silently overclaim "ocean currents carrying heat", the two things this single-layer periodic
     channel lacks (genuine net transport, and a globe).

The banked frames are the **diagnostic-pure** ``n_frames`` side-channel: ``n_frames=0`` leaves the κ
result bit-for-bit unchanged (the inert-seam discipline, §9.3, applied to motion).

Run headless (saves the GIF, prints the summary):

    python -m planet.demo_eddy_life
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from . import eddy_flux as ef
from .albedo import EBMParams, present_day_climate

_REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_FIGURE = _REPO_ROOT / "docs" / "figures" / "planet-eddy-life.gif"
OUTPUT_FIGURE = _REPO_ROOT / "outputs" / "planet-eddy-life.gif"

N_FRAMES = 48              # frames banked over the full release (modest → a small committed GIF)
FPS = 12                  # GIF playback rate
SAVE_DPI = 72             # modest dpi so the committed GIF stays lean (~1–2 MB, on par with the HTML globes)


@dataclass(frozen=True)
class EddyLifeResult:
    """The banked rung-A artifact — one released eddy life cycle with its viz frames."""

    eddy: ef.EddyFlux


def compute(nx: int = 80, ny: int = 80, n_frames: int = N_FRAMES) -> EddyLifeResult:
    """Run the present-day (steep-gradient) eddy life cycle with the frame side-channel banked."""
    eddy = ef.eddy_life_cycle(present_day_climate(EBMParams(s2=-0.48)), nx=nx, ny=ny, n_frames=n_frames)
    return EddyLifeResult(eddy=eddy)


def print_summary(r: EddyLifeResult) -> None:
    """Print the rung-A story — the emergent, mostly-reversible eddy life cycle on the band."""
    e = r.eddy
    print("\nThe emergent eddy life cycle — released on the barotropically-unstable jet\n")
    print(f"  jet: {e.jet_speed:.1f} m/s @ {e.jet_lat:.1f}°   "
          f"(Rayleigh–Kuo unstable: {e.rayleigh_kuo})")
    print(f"  eddy-KE saturates at {e.saturation_period:.0f} inertial periods (growth → saturation)")
    print(f"  emergent diffusivity: κ_bulk = {e.kappa_bulk:.3e} m²/s  →  D_eff = {e.D_eff:.5f} W m⁻² K⁻¹")
    print(f"  irreversible fraction = {e.irreversible_fraction:.2f}  "
          f"(~{100*(1-e.irreversible_fraction):.0f}% of the instantaneous flux just sloshes)")
    if e.frames is not None:
        fr = e.frames
        print(f"  banked frames: {fr.times.size} over [0, {fr.times[-1]:.0f}] periods "
              f"({fr.theta.shape[1]}×{fr.theta.shape[2]} grid)\n")


def save_animation(r: EddyLifeResult, fps: int = FPS, dpi: int = SAVE_DPI) -> Path:
    """Render and save the rung-A GIF (needs the optional ``viz`` extra; Pillow ships with matplotlib)."""
    import matplotlib
    matplotlib.use("Agg")                            # headless
    from matplotlib.animation import PillowWriter
    import matplotlib.pyplot as plt

    from .plots import eddy_life_animation

    anim = eddy_life_animation(r.eddy)
    for target in (DOCS_FIGURE, OUTPUT_FIGURE):
        target.parent.mkdir(parents=True, exist_ok=True)
        anim.save(target, writer=PillowWriter(fps=fps), dpi=dpi)
    plt.close("all")
    return DOCS_FIGURE


def main() -> None:
    import sys
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")     # °C, ₂, →, ∂, ², ⁻ on legacy codepages

    r = compute()
    print_summary(r)
    try:
        saved = save_animation(r)
        print(f"Animation saved → {saved.relative_to(_REPO_ROOT)}")
    except ImportError:
        print("(matplotlib not installed — install the viz extra to render the animation: "
              "pip install -e .[viz])")


if __name__ == "__main__":
    main()
