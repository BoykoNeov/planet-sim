"""The Rung-C banked artifact: the emergent eddy life cycle as a **particle flow-globe** (§9.5).

Rung A banked the flat two-panel GIF; Rung B banked the scalar-field Plotly globe. This demo banks the
**showcase**: the *same* released eddy life cycle (:func:`planet.eddy_flux.eddy_life_cycle`'s
diagnostic-pure ``n_frames`` side-channel) rendered as streaming particles on a real, rotatable three.js
planet (:mod:`planet.flow_globe`) — a self-contained HTML page (three.js vendored inline) that opens
straight off the filesystem.

It is the one renderer governed by the **honest-by-disclosure** carve-out (ADR 0002, 2026-06-12): the
particles imply persistent currents, but the flux is ~90 % reversible, so a visible on-screen disclaimer
documents the departure (see :mod:`planet.flow_globe`). The flow is confined to the band's *true*
coverage — we do not fabricate a global field.

Run headless (saves the HTML, prints the summary)::

    python -m planet.demo_eddy_particles
"""
from __future__ import annotations

from pathlib import Path

from . import demo_eddy_life
from .flow_globe import flow_field_from_eddy, save_flow_globe_html

_REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_FIGURE = _REPO_ROOT / "docs" / "figures" / "planet-eddy-particles.html"
OUTPUT_FIGURE = _REPO_ROOT / "outputs" / "planet-eddy-particles.html"

NX = NY = 80              # the eddy-channel resolution (matches the Rung-A/B banked artifacts)
N_FRAMES = 40             # frames over the full release (Rung C streams the saturated frame; B animates all)


def compute(nx: int = NX, ny: int = NY, n_frames: int = N_FRAMES) -> demo_eddy_life.EddyLifeResult:
    """Run the present-day (steep-gradient) eddy life cycle with the frame side-channel banked.

    Reuses :func:`planet.demo_eddy_life.compute` — Rungs A, B and C are three *views* of one identical
    life cycle (only the renderer differs), so the physics call is shared, never duplicated.
    """
    return demo_eddy_life.compute(nx=nx, ny=ny, n_frames=n_frames)


def save_globe(r: demo_eddy_life.EddyLifeResult) -> Path:
    """Build the generic flow field from the eddy and write the Rung-C HTML showcase."""
    field = flow_field_from_eddy(r.eddy)
    for target in (DOCS_FIGURE, OUTPUT_FIGURE):
        save_flow_globe_html(field, target)
    return DOCS_FIGURE


def main() -> None:
    import sys
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")     # °C, ₂, →, κ, β on legacy codepages

    r = compute()
    demo_eddy_life.print_summary(r)
    saved = save_globe(r)
    print(f"Particle flow-globe saved → {saved.relative_to(_REPO_ROOT)}")
    print("  open it in a browser (works straight off disk — three.js is vendored inline).")


if __name__ == "__main__":
    main()
