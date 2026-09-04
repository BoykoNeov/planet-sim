"""Rung 5B.4: the seasons dissolve the second cliff — a small ice cap the annual mean says cannot exist.

Rung 0+ (:mod:`planet.demo_bifurcation`) found the ice-albedo EBM's *second* fold: a polar cap smaller
than a critical radius ``θ_c ≈ 10°`` is held by **no** sun, so a brightening planet's cap does not shrink
gracefully — it reaches ``θ_c`` and **vanishes in a jump**, leaving a hysteresis loop of several W/m² in
the solar constant. That is North's small-ice-cap instability, and it is a statement about a planet at
**equilibrium**, where the heat capacity has cancelled out.

Rung 5B put the same ice-albedo feedback on the *seasonal* marcher, where a polar cap melts every summer
and re-freezes every winter. This demo asks whether the fold survives that — and the answer is no. At
Earth's tilt with a 50 m ocean mixed layer the perennial cap grows **one grid cell at a time**, straight
through the radius the annual-mean model forbids, and the dimming and brightening sweeps retrace each
other exactly: **no jump, no hysteresis**. Deepen the mixed layer — damping the seasonal swing without
touching the annual-mean sunlight — and the instability comes **back**: a planted cap survives at a sun
where a warm start grows none, which is the bistability the fold is made of.

The knob is depth, not tilt, and that is deliberate: obliquity would move the annual-mean reference *and*
the seasonal swing together (at zero tilt the annual-mean model has no fold at all), so the comparison
would attribute the difference to neither. Depth moves only the seasons.

Guarding against the grid, which matters here more than usual: the polar cell is 4.3° wide at 720 cells
and shrinks only as ``√Δx``, so *"the cap shrank smoothly to nothing"* and *"the critical cap fell below
one cell"* would look identical if read off the interpolated cap **radius**. So the verdict is read off
the perennial ice **cell count** (which the albedo feedback actually sees), off the hysteresis **loop
width** (which converges under refinement where quantization does not), and repeated at three
resolutions.

Run headless (saves the figure, prints the diagram) — this marches many limit cycles, ~5-10 minutes:

    python -m planet.demo_seasonal_sici
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import seasonal_sici as sici
from .bifurcation import EquilibriumCurve
from .seasonal import OCEAN_MIXED_DEPTH

# The main sweep: fine enough in S₀ that the perennial cap can be caught growing one cell at a time.
SWEEP_N_CELLS = 720
SWEEP_N_STEPS = 180
SWEEP_S0 = np.arange(1366.0, 1382.01, 1.0)          # W m⁻² — brackets where the cap appears at Earth tilt

# The depth axis — the seasonal swing damped toward the annual-mean limit (~geometric in h_ml).
DEPTHS = (OCEAN_MIXED_DEPTH, 200.0, 800.0)          # m
# 720, NOT a coarser grid, even though these marches are long (τ = C/B grows with h). At 360 cells a
# θ_c-sized cap spans a single cell, so the deep-ocean positive control — "a planted cap survives where a
# warm start grows none" — would come down to one cell of ice versus none, on the very grid the resolution
# panel below prints as "too coarse to tell". The control has to stand on a grid the verdict trusts.
DEPTH_N_CELLS = 720

# The resolution check — the same verdict on three grids.
RES_N_CELLS = (360, 720, 1440)
RES_S0 = np.arange(1372.0, 1380.01, 1.0)

_REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_FIGURE = _REPO_ROOT / "docs" / "figures" / "planet-seasonal-sici.png"
OUTPUT_FIGURE = _REPO_ROOT / "outputs" / "planet-seasonal-sici.png"


@dataclass(frozen=True)
class SeasonalSICIResult:
    """The banked artifact: the annual-mean parent, the seasonal sweep, and both robustness checks."""

    cfg: sici.SICIConfig
    parent: EquilibriumCurve                    # the annual-mean equilibrium diagram of this same model
    loop: sici.HysteresisLoop                   # the seasonal continuation, down and up
    depths: tuple[float, ...]
    seed_points: tuple[sici.SeedDependence, ...]   # the planted-cap test at each depth
    res_n_cells: tuple[int, ...]
    res_loop_width: tuple[float, ...]
    res_max_jump: tuple[int, ...]
    res_polar_cell: tuple[float, ...]
    res_fold_cells: tuple[int, ...]      # cells a θ_c cap covers — the jump a FOLD would make

    @property
    def theta_c(self) -> float:
        """The annual-mean parent's critical cap radius (degrees) — the size the seasons are asked about."""
        fold = self.parent.small_ice_cap_fold
        return float("nan") if fold is None else fold.cap_radius_deg

    @property
    def parent_loop_width(self) -> float:
        """The parent's hysteresis loop (W m⁻²) — the fold's width, what the seasonal sweep is compared to."""
        fold = self.parent.small_ice_cap_fold
        return float("nan") if fold is None else fold.S0 - self.parent.ice_free_threshold_S0

    @property
    def fold_cells(self) -> int:
        """Cells a ``θ_c`` cap covers on the main sweep's grid — the step a genuine fold would make."""
        return self.cfg.cells_in_cap(self.theta_c)

    @property
    def forbidden_caps(self) -> np.ndarray:
        """Seasonal cap radii that fall strictly **inside** the annual-mean model's forbidden band."""
        caps = self.loop.down.perennial_cap_deg
        return caps[(caps > 0.0) & (caps < self.theta_c)]


def compute() -> SeasonalSICIResult:
    """Run the parent solve, the seasonal continuation, the depth axis, and the resolution check."""
    cfg = sici.SICIConfig(n_cells=SWEEP_N_CELLS, n_steps=SWEEP_N_STEPS)
    parent = sici.annual_mean_curve(cfg)
    loop = sici.hysteresis_loop(cfg, SWEEP_S0)

    depth_cfg = sici.SICIConfig(n_cells=DEPTH_N_CELLS, n_steps=SWEEP_N_STEPS)
    depth_parent = sici.annual_mean_curve(depth_cfg)
    fold = depth_parent.small_ice_cap_fold
    S0_seed = fold.S0 if fold is not None else depth_parent.ice_free_threshold_S0
    seeds = tuple(sici.seed_dependence(depth_cfg, S0_seed, curve=depth_parent,
                                       ocean_mixed_depth=d, max_years=8000) for d in DEPTHS)

    widths, jumps, cells, fold_cells = [], [], [], []
    for n in RES_N_CELLS:
        c = sici.SICIConfig(n_cells=n, n_steps=SWEEP_N_STEPS)
        res_parent = sici.annual_mean_curve(c)
        res_fold = res_parent.small_ice_cap_fold
        lp = sici.hysteresis_loop(c, RES_S0)
        widths.append(lp.width)
        jumps.append(max(lp.down.max_cell_jump, lp.up.max_cell_jump))
        cells.append(c.polar_cell_deg)
        fold_cells.append(0 if res_fold is None else c.cells_in_cap(res_fold.cap_radius_deg))

    return SeasonalSICIResult(cfg, parent, loop, DEPTHS, seeds,
                              RES_N_CELLS, tuple(widths), tuple(jumps), tuple(cells), tuple(fold_cells))


def print_summary(r: SeasonalSICIResult) -> None:
    """Print the two diagrams side by side and the verdict."""
    cfg = r.cfg
    print("Rung 5B.4 — does the small-ice-cap instability survive the seasons?")
    print("=" * 78)
    print(f"Grid: {cfg.n_cells} cells on x = sin φ (polar cell {cfg.polar_cell_deg:.2f}° wide), "
          f"{cfg.n_steps} steps/yr, tilt {cfg.obliquity_deg:.2f}°, all ocean.")
    print()
    print("The ANNUAL-MEAN parent of this very model (dt-free, exact):")
    print(f"  critical cap radius θ_c   = {r.theta_c:6.2f}°   (no smaller cap is held by any sun)")
    print(f"  hysteresis loop           = {r.parent_loop_width:6.2f} W/m²  "
          f"(the jump's width in the solar constant)")
    print()
    print("The SEASONAL model, same operator, same annual-mean sunlight:")
    print(f"  loop width                = {r.loop.width:6.2f} W/m²  "
          f"(sampled every {r.loop.dS0:.1f} W/m²; detection threshold "
          f"{r.loop.threshold_deg:.2f}° = ½ polar cell)")
    print(f"  largest cap gap down↔up   = {np.max(r.loop.cap_gap_deg):6.2f}°")
    print(f"  perennial cells per step  = {r.loop.down.max_cell_jump:6d}     "
          f"(a fold at θ_c would switch on {r.fold_cells} cells in ONE step)")
    print()
    print("  S₀ (W/m²)   perennial cells   cap radius (°)   [down-sweep]")
    for p in r.loop.down.points:
        mark = "  ← below θ_c, forbidden in the annual mean" if 0.0 < p.perennial_cap_deg < r.theta_c else ""
        print(f"   {p.S0:7.1f}        {p.n_perennial_cells:4d}            {p.perennial_cap_deg:6.2f}{mark}")
    print()
    print("Planting a cap of θ_c at the parent's fold, versus a warm start (the bistability test):")
    print("   depth (m)   seasonal swing @80°   warm cap   planted cap survives as   two climates?")
    for d, sd in zip(r.depths, r.seed_points):
        print(f"   {d:8.0f}        {sd.warm_polar_amplitude_K:6.2f} K          "
              f"{sd.warm_cap_deg:6.2f}°           {sd.survived_cap_deg:6.2f}°            "
              f"{'YES' if sd.bistable else 'no'}")
    print()
    print("Resolution check — the verdict must not move with the grid:")
    print("   cells    polar cell (°)   loop width (W/m²)   max cells/step   a fold would flip")
    for n, pc, w, j, fc in zip(r.res_n_cells, r.res_polar_cell, r.res_loop_width,
                               r.res_max_jump, r.res_fold_cells):
        note = "   ← too coarse to tell" if fc <= 2 else ""
        print(f"   {n:5d}       {pc:6.2f}           {w:6.2f}              {j:3d}"
              f"              {fc:3d} cells{note}")
    print("   (a step of 1-2 cells where a fold would flip 5-11 is an order of magnitude below one;")
    print("    on the coarsest grid a θ_c cap is barely a cell, so that row cannot settle it either way.)")
    print()
    n_forbidden = r.forbidden_caps.size
    print(f"VERDICT: the seasonal planet holds {n_forbidden} stable cap sizes below θ_c = {r.theta_c:.2f}° "
          f"({', '.join(f'{c:.1f}°' for c in r.forbidden_caps)}),")
    print("         grows them one cell at a time, and shows no hysteresis — the fold is gone at 50 m.")
    print(f"         It returns as the mixed layer deepens: bistability at "
          f"{', '.join(f'{d:.0f} m' for d, s in zip(r.depths, r.seed_points) if s.bistable) or 'no depth tested'}.")
    print("         Banked as a mechanism (seasonal amplitude vs. an equilibrium fold), not as a general")
    print("         claim: one tilt, a step-function albedo, no sea-ice thermodynamics.")


def save_figure(r: SeasonalSICIResult) -> Path:
    """Render and save the figure (needs the optional ``viz`` extra)."""
    import matplotlib
    matplotlib.use("Agg")                                    # headless
    from .plots import seasonal_sici_figure

    fig = seasonal_sici_figure(r)
    for target in (DOCS_FIGURE, OUTPUT_FIGURE):
        target.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(target, dpi=130)
    return DOCS_FIGURE


def main() -> None:
    import sys
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")             # °, ₀, →, θ on legacy codepages

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
