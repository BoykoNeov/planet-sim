"""The §9.1 banked artifact: the exoplanet knobs — other-world climate from two pinned levers.

The capstone's *exoplanet sandbox* made into a demonstration. Two knobs (:mod:`planet.exoplanet`)
turn the present-day Earth model into an other-world climate, each a **parameter derivation**, not new
physics:

* **Stellar spectrum → ice albedo.** A redder, cooler host star (an M-dwarf) emits more near-IR light,
  where snow and ice are dark, so its broadband ice albedo is lower and the **ice-albedo feedback
  weakens**. The headline: an M-dwarf planet's Snowball hysteresis loop is **much narrower and shifted
  to a lower freeze threshold** — *a redder star is harder to snowball* (Joshi & Haberle 2012).
* **Planet size → meridional transport.** A bigger planet transports heat less effectively per unit
  area (``D ∝ 1/size²``), so the **equator-to-pole gradient sharpens** and the ice cap reaches
  equatorward — while the 0-D global mean is (nearly) unmoved (size enters only the mean-preserving
  transport).

Run headless (saves the figure, prints the story):

    python -m planet.demo_exoplanet
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np

from . import exoplanet as ex
from .albedo import EBMParams, HysteresisLoop, present_day_climate, snowball_hysteresis
from .ebm import ClimateState, equilibrium_temperature_0d

MDWARF_TEFF = 3050.0           # K — a cool M-dwarf (≈ M5V; Pecaut & Mamajek 2013) — the redder host star
SIZES = (0.5, 1.0, 2.0)        # Earth radii — a small world, Earth, a super-Earth (transport-only effect)
SIZE_N_TAU = 0.02              # relaxation step for the size profiles (fine enough for a steady T(φ))

_REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_FIGURE = _REPO_ROOT / "docs" / "figures" / "planet-exoplanet.png"
OUTPUT_FIGURE = _REPO_ROOT / "outputs" / "planet-exoplanet.png"


@dataclass(frozen=True)
class ExoplanetResult:
    """The banked exoplanet artifact — the plain bundle the figure and the summary consume.

    ``stellar_ai`` maps each labelled stellar type to its effective ice albedo; ``sun_loop`` /
    ``mdwarf_loop`` are the Snowball hysteresis loops of a Sun-like and an M-dwarf planet (the
    stellar headline); ``sizes`` / ``size_states`` are the relaxed present-day climates at a range of
    planet sizes (the size headline). Plain dataclasses/arrays — the loose-coupling currency.
    """

    stellar_ai: dict[str, float]
    mdwarf_label: str
    sun_loop: HysteresisLoop
    mdwarf_loop: HysteresisLoop
    sizes: tuple[float, ...]
    size_states: tuple[ClimateState, ...]


def compute(mdwarf_teff: float = MDWARF_TEFF, sizes: tuple[float, ...] = SIZES,
            n_cells: int = 90, n_steps: int = 22, sweep_n_tau: float = 0.1,
            size_n_tau: float = SIZE_N_TAU) -> ExoplanetResult:
    """Run both knobs → :class:`ExoplanetResult` (no plotting).

    The stellar headline traces the Snowball hysteresis (:func:`~planet.albedo.snowball_hysteresis`)
    for a Sun-like and an ``mdwarf_teff`` host star (only the ice albedo differs); the size headline
    relaxes the present-day climate (:func:`~planet.albedo.present_day_climate`) at each of
    ``sizes`` (only the transport ``D`` differs). ``n_cells`` / ``n_steps`` / ``sweep_n_tau`` are the
    sweep's resolution knobs (kept modest — this is a banked demo, not the gate); the structural physics
    is sealed fast by :mod:`tests.test_exoplanet`.
    """
    base = EBMParams(n_cells=n_cells)
    # Stellar headline: Sun vs M-dwarf Snowball loops (ai is the only difference).
    sun_loop = snowball_hysteresis(base, S0_min=900.0, S0_max=1900.0, n_steps=n_steps, n_tau=sweep_n_tau)
    mdwarf_params = replace(base, ai=ex.stellar_ice_albedo(mdwarf_teff))
    mdwarf_loop = snowball_hysteresis(mdwarf_params, S0_min=900.0, S0_max=1900.0,
                                      n_steps=n_steps, n_tau=sweep_n_tau)
    stellar_ai = {name: ex.stellar_ice_albedo(T) for name, T in ex.STELLAR_TYPES.items()}

    # Size headline: relaxed present-day climate at each size (D ∝ 1/size² is the only difference).
    size_states = tuple(
        present_day_climate(replace(base, D=ex.transport_for_size(size)), n_tau=size_n_tau)
        for size in sizes
    )
    mdwarf_label = _stellar_label(mdwarf_teff)
    return ExoplanetResult(stellar_ai=stellar_ai, mdwarf_label=mdwarf_label,
                           sun_loop=sun_loop, mdwarf_loop=mdwarf_loop,
                           sizes=tuple(sizes), size_states=size_states)


def _stellar_label(teff: float) -> str:
    """The nearest labelled stellar type for ``teff`` (for the figure/summary), e.g. ``'M5V (3050 K)'``."""
    name = min(ex.STELLAR_TYPES, key=lambda n: abs(ex.STELLAR_TYPES[n] - teff))
    return f"{name.split(' ')[0]} ({teff:.0f} K)"


def print_summary(r: ExoplanetResult) -> None:
    """Print the two exoplanet-knob stories — the redder-star snowball shift and the size gradient."""
    print("\nExoplanet knobs: other-world climate from two pinned levers\n")
    print("  Knob 1 — stellar spectrum → ice albedo (a redder star is harder to snowball):")
    print("    stellar type        effective ice albedo a_ice")
    for name, ai in r.stellar_ai.items():
        bar = "█" * int(round(ai * 40))
        print(f"      {name:14s}   {ai:.3f}  {bar}")
    sw, mw = r.sun_loop.hysteresis_width, r.mdwarf_loop.hysteresis_width
    print(f"\n    Snowball hysteresis loop width:  Sun (G2V) {sw:.0f} W/m²  →  {r.mdwarf_label} {mw:.0f} W/m²"
          f"  ({100*(1-mw/sw):.0f}% narrower)")
    print(f"    freeze threshold S₀:             Sun {r.sun_loop.freeze_S0:.0f}  →  "
          f"{r.mdwarf_label} {r.mdwarf_loop.freeze_S0:.0f} W/m²  (must dim further to snowball)")
    print()
    T0 = equilibrium_temperature_0d()      # the analytic constant-albedo 0-D mean — size-INVARIANT
    print("  Knob 2 — planet size → meridional transport (a bigger planet sharpens the gradient):")
    print("    size (R⊕)   D (W m⁻² K⁻¹)   analytic T̄₀   relaxed T̄   ice line")
    for size, st in zip(r.sizes, r.size_states):
        print(f"      {size:4.1f}      {ex.transport_for_size(size):7.3f}      {T0:6.2f} °C   {st.global_mean_T:6.2f} °C   "
              f"{st.ice_line_lat:5.1f}°")
    print(f"    → the analytic 0-D mean T̄₀ ≈ {T0:.1f} °C is size-INVARIANT (size enters only the transport),")
    print("      so the relaxed mean barely moves from ice-free to a small cap (0.5 → 1.0 R⊕); it then drops")
    print("      sharply at 2 R⊕ where the enlarged ice cap's albedo feedback cools the planet. Either way")
    print("      the gradient steepens and the ice cap creeps equatorward.\n")


def save_figure(r: ExoplanetResult) -> Path:
    """Render and save the §9.1 exoplanet artifact (needs the optional ``viz`` extra)."""
    import matplotlib
    matplotlib.use("Agg")                            # headless
    from .plots import exoplanet_figure

    fig = exoplanet_figure(r.sun_loop, r.mdwarf_loop, r.mdwarf_label, r.stellar_ai,
                           r.size_states, r.sizes)
    for target in (DOCS_FIGURE, OUTPUT_FIGURE):
        target.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(target, dpi=130)
    return DOCS_FIGURE


def main() -> None:
    import sys
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")     # °C, ₂, →, ° on legacy codepages

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
