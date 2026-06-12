"""``engines.fluid`` — the program's second shared engine: a rotating shallow-water solver.

A hyperbolic, **explicit** (CFL-limited, C-grid, SSP-RK3) shallow-water solver on a
doubly-periodic β-plane — deliberately sharing *no* machinery with the parabolic-implicit
:mod:`engines.diffusion`. Built here (Planet Phase 3) so the circulation
coupler (Phase 4) and the documented GCM climb reuse it behind ``engines/fluid/CONTRACT.md``.

Load the **one-page contract** (`engines/fluid/CONTRACT.md`), not this package's internals.
"""
from .shallowwater import (
    Grid2D,
    ShallowWater,
    SWState,
    uniform_grid,
)
from .layered import (
    LayeredShallowWater,
    LayeredState,
    ThermalWindBackground,
)
from .stability import TwoLayerStability

__all__ = [
    "Grid2D",
    "ShallowWater",
    "SWState",
    "uniform_grid",
    # N-layer baroclinic extension (rung 3)
    "LayeredShallowWater",
    "LayeredState",
    "ThermalWindBackground",
    "TwoLayerStability",
]
