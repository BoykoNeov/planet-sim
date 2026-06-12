"""The Rung-B banked artifact: the emergent eddy life cycle, animated **on the globe** (§9.5).

Rung A (:mod:`planet.demo_eddy_life`) banked the flat two-panel GIF. This demo banks the **globe**
form: the *same* released eddy life cycle (:func:`planet.eddy_flux.eddy_life_cycle`'s diagnostic-pure
``n_frames`` side-channel) lifted onto the existing Plotly planet (:mod:`planet.eddy_globe`), with the
flux-budget panel beside it — a self-contained, play/slider HTML globe. *No new stack* (it reuses the
``[webviz]`` Plotly the biome map depends on, plan §9.5).

It carries the two honesty edges geometrically: the band is rendered at its **true ~55° width** in a
**single hemisphere** (not a 360° wrap, not mirrored), and the throughput-vs-net panel keeps the
**~90 %-reversible** finding on screen. See :mod:`planet.eddy_globe` for the full rationale.

Run headless (saves the HTML, prints the summary):

    python -m planet.demo_eddy_globe
"""
from __future__ import annotations

from pathlib import Path

from . import demo_eddy_life
from .eddy_globe import save_eddy_globe_html

_REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_FIGURE = _REPO_ROOT / "docs" / "figures" / "planet-eddy-globe.html"
OUTPUT_FIGURE = _REPO_ROOT / "outputs" / "planet-eddy-globe.html"

NX = NY = 80              # the eddy-channel resolution (matches the Rung-A banked GIF)
N_FRAMES = 40             # frames over the full release (surfacecolor-only updates → a lean HTML)
FRAME_MS = 120            # play-back rate (ms per frame)


def compute(nx: int = NX, ny: int = NY, n_frames: int = N_FRAMES) -> demo_eddy_life.EddyLifeResult:
    """Run the present-day (steep-gradient) eddy life cycle with the frame side-channel banked.

    Reuses :func:`planet.demo_eddy_life.compute` — Rung A and Rung B are two *views* of one identical
    life cycle (the only difference is the renderer), so the physics call is shared, not duplicated.
    """
    return demo_eddy_life.compute(nx=nx, ny=ny, n_frames=n_frames)


def save_globe(r: demo_eddy_life.EddyLifeResult, frame_ms: int = FRAME_MS) -> Path:
    """Render and save the Rung-B HTML globe (needs the optional ``[webviz]`` extra — Plotly)."""
    for target in (DOCS_FIGURE, OUTPUT_FIGURE):
        save_eddy_globe_html(r.eddy, target, frame_ms=frame_ms)
    return DOCS_FIGURE


def main() -> None:
    import sys
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")     # °C, ₂, →, ∂, ², ⁻ on legacy codepages

    r = compute()
    demo_eddy_life.print_summary(r)
    try:
        saved = save_globe(r)
        print(f"Globe animation saved → {saved.relative_to(_REPO_ROOT)}")
    except ImportError:
        print("(plotly not installed — install the webviz extra to render the globe: "
              "pip install -e .[webviz])")


if __name__ == "__main__":
    main()
