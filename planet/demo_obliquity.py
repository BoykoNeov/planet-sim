"""The §9.1 banked artifact: the obliquity knob — axial tilt sets the insolation gradient.

The capstone's *obliquity* lever made into a demonstration. Like the two exoplanet knobs
(:mod:`planet.demo_exoplanet`), it is a **parameter derivation**, not new physics: a planet's
**axial tilt** ``ε`` fixes the second-Legendre coefficient ``s₂`` of its annual-mean insolation
(:mod:`planet.obliquity`), and that single number sets how steeply the year's sunlight is
graded from equator to pole.

* **Less tilt → a steeper gradient.** At zero tilt the sun is pinned over the equator
  (``s₂ = −5/8`` exactly) — a hot equator, a frozen pole.
* **More tilt → a flatter planet.** As ``ε`` grows the summer sun reaches higher latitudes, ``s₂``
  rises toward zero, the pole warms and the **ice cap retreats** (Earth's 23.44° lands the ice line at
  ~70°, the climlab benchmark; by ~40° the cap is gone).
* **Past ≈55° the gradient reverses** — the poles receive *more* annual sunlight than the equator
  (the high-obliquity world). Real and surfaced, but the EBM's single-P₂-mode insolation only
  *approximates* the strongly-flattened profile there (the named scope edge).

Run headless (saves the figure, prints the story):

    python -m planet.demo_obliquity
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import obliquity as ob
from .albedo import EBMParams, present_day_climate
from .ebm import ClimateState, S2_INSOLATION

CLIMATE_TILTS = (0.0, 10.0, ob.OBLIQUITY_EARTH, 40.0)   # ° — the relaxed-climate panel (no tilt → Earth → flat)
CLIMATE_N_TAU = 0.02           # relaxation step for the T(φ) profiles (fine enough for a steady climate)
CURVE_N = 46                   # samples of the s₂(ε) curve over 0–90° (every 2°)

_REPO_ROOT = Path(__file__).resolve().parents[2]
DOCS_FIGURE = _REPO_ROOT / "docs" / "figures" / "planet-obliquity.png"
OUTPUT_FIGURE = _REPO_ROOT / "outputs" / "planet-obliquity.png"


@dataclass(frozen=True)
class ObliquityResult:
    """The banked obliquity artifact — the plain bundle the figure and the summary consume.

    ``eps_grid`` / ``s2_grid`` trace the geometric ``s₂(ε)`` curve (the mechanism); ``s2_earth`` is the
    value at Earth's tilt (the climlab cross-check); ``tilts`` / ``climate_states`` are the relaxed
    present-day climates at a range of obliquities (the climate response). Plain arrays/dataclasses.
    """

    eps_grid: np.ndarray
    s2_grid: np.ndarray
    s2_earth: float
    tilts: tuple[float, ...]
    climate_states: tuple[ClimateState, ...]


def compute(tilts: tuple[float, ...] = CLIMATE_TILTS, n_cells: int = 90,
            curve_n: int = CURVE_N, climate_n_tau: float = CLIMATE_N_TAU) -> ObliquityResult:
    """Run the obliquity knob → :class:`ObliquityResult` (no plotting).

    Traces the geometric ``s₂(ε)`` curve over 0–90° (:func:`~planet.obliquity.insolation_p2_coefficient`)
    and relaxes the present-day climate (:func:`~planet.albedo.present_day_climate`) at each of
    ``tilts`` (only ``s₂`` differs — the knob feeds it through
    :func:`~planet.obliquity.insolation_s2`). ``curve_n`` / ``n_cells`` / ``climate_n_tau`` are
    resolution knobs (kept modest — this is a banked demo, not the gate; the physics is sealed fast by
    :mod:`tests.test_obliquity`).
    """
    eps_grid = np.linspace(ob.OBLIQUITY_MIN, ob.OBLIQUITY_MAX, curve_n)
    s2_grid = np.array([ob.insolation_p2_coefficient(e) for e in eps_grid])
    climate_states = tuple(
        present_day_climate(EBMParams(s2=ob.insolation_s2(e), n_cells=n_cells), n_tau=climate_n_tau)
        for e in tilts
    )
    return ObliquityResult(eps_grid=eps_grid, s2_grid=s2_grid,
                           s2_earth=ob.insolation_p2_coefficient(ob.OBLIQUITY_EARTH),
                           tilts=tuple(tilts), climate_states=climate_states)


def print_summary(r: ObliquityResult) -> None:
    """Print the two obliquity stories — the s₂(ε) mechanism and the relaxed-climate response."""
    print("\nObliquity knob: axial tilt → the annual-mean-insolation gradient\n")
    print("  Mechanism — the insolation P₂ coefficient s₂(ε) (from the daily-insolation geometry):")
    print("    tilt ε     s₂(ε)     note")
    notes = {0.0: "−5/8 exactly (sun pinned at the equator)",
             ob.OBLIQUITY_EARTH: f"≈ climlab −0.48 (cross-check: {S2_INSOLATION:+.2f})",
             54.0: "near the sign flip — poles warmer than equator above here"}
    for e in (0.0, 10.0, ob.OBLIQUITY_EARTH, 45.0, 54.0, 65.0, 90.0):
        s2 = ob.insolation_p2_coefficient(e)
        note = next((v for k, v in notes.items() if abs(k - e) < 0.6), "")
        flip = "  ←gradient reversed" if s2 > 0 and not note else ""
        print(f"      {e:5.1f}°   {s2:+.4f}   {note}{flip}")
    print()
    print("  Climate response — the relaxed present-day climate at each tilt (only s₂ differs):")
    print("    tilt ε    s₂ (knob)   global mean T̄    ice line")
    for e, st in zip(r.tilts, r.climate_states):
        print(f"      {e:5.1f}°   {ob.insolation_s2(e):+.4f}    {st.global_mean_T:6.2f} °C     {st.ice_line_lat:5.1f}°")
    print("    → more tilt spreads the year's sunlight poleward: the pole warms and the ice cap retreats")
    print("      (Earth's 23.44° lands the ice line at ~70°, the climlab benchmark; by ~40° the cap is gone).\n")


def save_figure(r: ObliquityResult) -> Path:
    """Render and save the §9.1 obliquity artifact (needs the optional ``viz`` extra)."""
    import matplotlib
    matplotlib.use("Agg")                            # headless
    from .plots import obliquity_figure

    fig = obliquity_figure(r.eps_grid, r.s2_grid, r.s2_earth, r.climate_states, r.tilts)
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
