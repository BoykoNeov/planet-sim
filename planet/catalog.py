"""The one catalogue of planet-sim's demos — the single source of truth shared by every surface.

The interactive launcher (:mod:`planet.__main__`), the generated landing page
(:mod:`planet.site`), and their drift-guard tests all read :data:`DEMOS` from here. Add a demo
to this tuple and it appears in the menu, the CLI, *and* the webpage at once — there is no second
list to keep in sync. Kept in its own module (not in ``__main__``) so importing it never risks
re-executing the launcher: under ``python -m planet`` the launcher *is* ``__main__``, and importing
``planet.__main__`` a second time would mint a duplicate ``Demo`` class and ``DEMOS`` tuple.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

# Repo root (…/planet/catalog.py → parents[1]); the artifact paths below are relative to it.
_REPO_ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class Demo:
    """One catalogue entry — a runnable ``planet.<module>`` exposing ``main()``."""

    key: str                 # the short name the user types
    module: str              # the importable module whose main() we call
    title: str               # one-line headline for the menu / page card
    blurb: str               # a sentence of what it shows
    extras: tuple[str, ...]  # the pip extras its figure needs ("viz" / "webviz")
    artifact: str            # repo-relative primary artifact (a figure or an interactive globe)
    section: str             # grouping header (menu sections / page gallery groups)
    sim: bool = False        # True ⇒ runs a multi-second fluid simulation (warn first)
    interactive: str = ""    # repo-relative interactive globe HTML, when the primary artifact is a still


# The single source of truth. The menu, the CLI dispatch, the landing page, and the catalogue
# tests all read this one list — add a demo here and it shows up in every surface at once.
DEMOS: tuple[Demo, ...] = (
    Demo("snowball", "planet.demo_snowball", "Snowball-Earth hysteresis",
         "one knob (the solar constant), two stable climates, a catastrophic freeze",
         ("viz",), "docs/figures/planet-snowball.png", "Climate — energy balance"),
    Demo("biomes", "planet.demo_biomes", "Climate → biome map",
         "the Whittaker (temperature, rainfall) classifier; warming migrates the bands poleward",
         ("viz",), "docs/figures/planet-biomes.png", "Climate — energy balance"),
    Demo("exoplanet", "planet.demo_exoplanet", "Exoplanet knobs",
         "a redder star and a bigger planet reshape the climate and the ice line",
         ("viz",), "docs/figures/planet-exoplanet.png", "Climate — energy balance"),
    Demo("obliquity", "planet.demo_obliquity", "Axial tilt (obliquity)",
         "how the planet's tilt reshapes the pole-to-equator sunlight",
         ("viz",), "docs/figures/planet-obliquity.png", "Climate — energy balance"),
    Demo("shallowwater", "planet.demo_shallowwater", "Rotating shallow-water atmosphere",
         "geostrophic adjustment on the sphere — the circulation engine on its own",
         ("viz",), "docs/figures/planet-shallowwater.png",
         "Circulation — shallow-water (runs a short sim)", sim=True),
    Demo("coupler", "planet.demo_coupler", "EBM → circulation coupler",
         "an emergent jet grows from the pole-to-equator temperature gradient (one-way coupling)",
         ("viz", "webviz"), "docs/figures/planet-coupler.png",
         "Circulation — shallow-water (runs a short sim)", sim=True,
         interactive="docs/figures/planet-coupler-map.html"),
    Demo("eddy_life", "planet.demo_eddy_life", "Eddy life cycle — GIF",
         "the emergent eddy stirring the temperature, animated as a two-panel GIF",
         ("viz",), "docs/figures/planet-eddy-life.gif",
         "Circulation — shallow-water (runs a short sim)", sim=True),
    Demo("eddy_globe", "planet.demo_eddy_globe", "Eddy life cycle — globe",
         "the same eddy life cycle, animated on the interactive globe",
         ("webviz",), "docs/figures/planet-eddy-globe.html",
         "Circulation — shallow-water (runs a short sim)", sim=True),
    Demo("eddy_particles", "planet.demo_eddy_particles", "Eddy life cycle — particle globe",
         "the showcase: the eddy flow streamed as particles on a real, rotatable 3-D planet",
         (), "docs/figures/planet-eddy-particles.html",
         "Circulation — shallow-water (runs a short sim)", sim=True),
    Demo("map", "planet.planetmap", "Interactive biome-map globe",
         "the present-day globe — rotate / zoom / hover (the live sliders run in the notebook)",
         ("webviz",), "docs/figures/planet-map.html", "Interactive globes"),
)


def globe_href(demo: Demo) -> str | None:
    """The interactive globe HTML for this demo, if it has one (primary artifact or the secondary)."""
    if demo.artifact.endswith(".html"):
        return demo.artifact
    return demo.interactive or None
