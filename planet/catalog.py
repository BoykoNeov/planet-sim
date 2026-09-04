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
    Demo("bifurcation", "planet.demo_bifurcation", "The complete equilibrium diagram",
         "every climate the sun allows, stable and unstable — and the second cliff: a polar cap smaller "
         "than a critical size cannot exist (the small-ice-cap instability)",
         ("viz",), "docs/figures/planet-bifurcation.png", "Climate — energy balance"),
    Demo("biomes", "planet.demo_biomes", "Climate → biome map",
         "the Whittaker (temperature, rainfall) classifier; warming migrates the bands poleward",
         ("viz",), "docs/figures/planet-biomes.png", "Climate — energy balance"),
    Demo("exoplanet", "planet.demo_exoplanet", "Exoplanet knobs",
         "a redder star and a bigger planet reshape the climate and the ice line",
         ("viz",), "docs/figures/planet-exoplanet.png", "Climate — energy balance"),
    Demo("obliquity", "planet.demo_obliquity", "Axial tilt (obliquity)",
         "how the planet's tilt reshapes the pole-to-equator sunlight",
         ("viz",), "docs/figures/planet-obliquity.png", "Climate — energy balance"),
    Demo("seasonal", "planet.demo_seasonal", "Seasonal cycle & continentality",
         "turn on the seasons and heat capacity wakes: a land tile swings far more than the ocean tile at "
         "the same latitude — continentality from the C contrast alone",
         ("viz",), "docs/figures/planet-seasonal.png", "Climate — energy balance"),
    Demo("seasonal_sici", "planet.demo_seasonal_sici", "Do the seasons dissolve the second cliff?",
         "the annual-mean model says a polar cap smaller than ~10° cannot exist; switch the seasons on and "
         "it grows one grid cell at a time with no jump — and a deep, sluggish ocean brings the cliff back",
         ("viz",), "docs/figures/planet-seasonal-sici.png", "Climate — energy balance", sim=True),
    Demo("seasonal_ice_map", "planet.demo_seasonal_ice_map", "Seasonal ice map — animated",
         "the seasons, the continents and the snow: winter snow spreads over the land, sea ice lingers at the "
         "poles, and the continents end colder in the annual mean (a month-by-month GIF)",
         ("viz", "webviz"), "docs/figures/planet-seasonal-ice-map.gif", "Climate — energy balance",
         interactive="docs/figures/planet-seasonal-ice-globe.html"),
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
    Demo("flow_serialize", "planet.demo_flow_serialize", "Vector-field interchange — globe",
         "any producer's (u, v) through one schema + one renderer: a synthetic global flow, the seam the spin-out binds on",
         ("webviz",), "docs/figures/planet-flow-serialize.html",
         "Circulation — shallow-water (runs a short sim)", sim=True),
    Demo("qg_particles", "planet.demo_qg_particles", "QG turbulence — particle globe",
         "the emergent rung-3 two-layer QG condensate (coherent vortices + PV filaments) streamed as "
         "particles through the same seam and renderer as the eddy band — an idealized box, not the sea",
         (), "docs/figures/planet-qg-particles.html",
         "Circulation — shallow-water (runs a short sim)", sim=True),
    Demo("map", "planet.planetmap", "Interactive biome-map globe",
         "the present-day globe — rotate / zoom / hover (the live sliders run in the notebook)",
         ("webviz",), "docs/figures/planet-map.html", "Interactive globes"),
    Demo("ocean_currents", "planet.demo_ocean_currents", "Real ocean currents — particle globe",
         "one day of real OSCAR surface currents (NASA PO.DAAC) streamed through the same seam and "
         "renderer as the eddy band — not a planet-sim simulation",
         ("ocean",), "docs/figures/planet-ocean-currents.html", "Interactive globes"),
    Demo("ocean_seasonal", "planet.demo_ocean_seasonal", "Seasonal ocean currents — animated globe",
         "twelve monthly snapshots of real OSCAR currents (2020) crossfaded on a time axis — watch the "
         "month badge cycle and the Somali Current reverse with the monsoon",
         ("ocean",), "docs/figures/planet-ocean-currents-seasonal.html", "Interactive globes"),
)


def globe_href(demo: Demo) -> str | None:
    """The interactive globe HTML for this demo, if it has one (primary artifact or the secondary)."""
    if demo.artifact.endswith(".html"):
        return demo.artifact
    return demo.interactive or None
