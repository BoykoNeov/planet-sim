"""Does the **seasonal cycle** survive the small-ice-cap instability — or dissolve it? (rung 5B.4)

Rung 0+ (:mod:`planet.bifurcation`) traced the *complete* equilibrium diagram of the **annual-mean**
ice-albedo EBM and found its second fold: the **small-ice-cap instability** (SICI, North 1984). Read off
the curve ``S₀(x_s)``, a polar cap smaller than a critical angular radius ``θ_c`` is held by **no** sun —
so a brightening planet's cap cannot shrink smoothly to nothing. It **jumps**: the last stable cap has a
finite size, and past the fold the planet flips ice-free. Dimming back, the ice-free branch survives to a
*lower* sun before a cap re-appears — a **hysteresis loop** in the solar constant, ``S₀(fold) − S₀(threshold)``,
narrow but real.

That whole diagram lives at ``∂T/∂t = 0``, where the heat capacity ``C`` cancels. Rung 5B.1 woke ``C`` by
marching the **seasons** (:mod:`planet.seasonal`), and 5B.1+ hung the same rung-0 step-function ice-albedo
on the marcher — which grew a *migrating* seasonal ice edge and put Phase 1's Snowball bistability inside
the annual cycle, but explicitly **named-and-deferred** the bifurcation study. This module is that study,
and it asks the question the seasonal literature asks:

    A polar cap that **melts every summer and re-freezes every winter** is not the object the annual-mean
    fold is a statement about. Does the instability survive being seasonally cycled?

The experiment — one axis, one pinned reference
-----------------------------------------------
The knob is the **ocean mixed-layer depth** ``h_ml`` at **fixed Earth obliquity**, and that choice is
load-bearing. Obliquity looks like the natural axis (``ε = 0`` is "no seasons") but it is *confounded*:
tilt sets the **annual-mean** insolation gradient as well as the seasonal swing, so the annual-mean parent
moves with it — at ``ε = 0`` the parent has **no SICI at all**, at ``ε = 15°`` its ``θ_c`` is still shrinking
at 2880 cells (resolution-limited, not converged), and only at Earth's ``ε = 23.44°`` does it converge
(``θ_c → 9.87°``, loop ``→ 6.76 W m⁻²``). Sweeping tilt would move the reference and the effect together and
attribute the difference to neither. Depth moves **only** the seasonal amplitude: the annual-mean forcing
is untouched, so :func:`annual_mean_curve` is one fixed, converged reference for the whole sweep, and

* ``h_ml → ∞``  ⟹ the seasonal swing → 0 and the marcher must **recover the annual-mean parent**, folds and
  all — the reduction (approached as ``1/(ωC)``, so *monotone approach*, never equality: see below);
* ``h_ml = 50 m`` (Earth) ⟹ the payoff.

Reading the verdict without being fooled by the grid
----------------------------------------------------
On the uniform-``x`` grid the polar cell's *latitude* width is ``≈ √(2Δx)`` — **4.3° at 720 cells, 3.0° at
1440**, and it falls only as ``√Δx``, so quadrupling the cells merely halves it. That matters more than it
looks: *"the cap shrank smoothly to nothing"* and *"θ_c fell below the polar cell"* are the **same data** on
a coarse grid. Three guards, in increasing order of how much weight they carry:

1. **The hysteresis loop width in ``S₀``** (:attr:`HysteresisLoop.width`) is the resolution-robust signature —
   a genuine fold gives a finite interval where the up- and down-sweeps disagree, and that interval
   converges under refinement; quantization gives one that shrinks with ``Δx``. It is defined against an
   explicit detection threshold (half a polar cell, :attr:`SICIConfig.cap_resolution_deg`) so that
   *"no loop"* and *"a loop below detection"* stay distinguishable.
2. **Seed dependence at one sun** (:func:`seed_dependence`) is the *second*, independent observable, and it
   can dissociate from the first: a warm-started continuation can walk continuously **through** a region
   where two states coexist, because the warm start keeps it pinned to one branch. So a cap of the parent's
   critical size is **planted** into a converged ice-free limit cycle (:func:`plant_cap`) and the march
   decides whether it survives — the seasonal translation of :func:`planet.bifurcation.relax_from_curve`.
   Two outcomes at one ``S₀``, reached without reference to any sweep direction, *is* the bistability.
3. **The perennial ice-cell count** (:attr:`SweepPoint.n_perennial_cells`) is the load-bearing one. The cap
   *radius* is interpolated (:func:`planet.ebm.ice_line_latitude` reads the temperature crossing, so it
   returns sub-cell values), but the albedo feedback only ever sees **whole cells flip**. A radius curve can
   therefore look smooth over a state that is stepping. Counting the cells that stay frozen all year settles
   it — but the count must be read against a scale, because a step of "one cell" means different things on
   different grids and at different ``S₀`` samplings. That scale is :meth:`SICIConfig.cells_in_cap` at
   ``θ_c``: **a fold turns a whole ``θ_c``-sized cap on in a single step**, so what disqualifies the fold is
   a step *an order of magnitude smaller than that*, not the literal value 1.

What is banked, and what is not
--------------------------------
**The finding (loose, calibrated).** At Earth's tilt with a 50 m mixed layer the perennial cap grows
**one cell at a time** — 0 → 1 → 2 → … , against the ~5 cells at 720 (~10 at 1440) that a fold at ``θ_c``
would switch on at once — straight through the annual-mean parent's critical radius, with
**no hysteresis** detectable at either 720 or 1440 cells, where the parent's own loop is a clear ~6.8 W m⁻².
Deepening the mixed layer brings the instability **back**. So in this model the SICI is not a robust property
of the ice-albedo feedback but of the *annual-mean idealisation*: it is the seasonal cycle's amplitude that
dissolves it, exactly the mechanism Huang & Bowman (1992) isolate ("the temperature gradient and the
**amplitude of the seasonal cycle**") and the direction Wagner & Eisenman (2015) report ("the stability of
the sea ice cover vastly increases with the inclusion of ... a seasonal cycle in solar forcing").

**What this is NOT.** It is not a replication of either paper, and the agreement must not be read as one.
Wagner & Eisenman carry sea-ice **thermodynamics** (thickness, latent heat, a separate ice surface); this
model has a step-function albedo on temperature alone, so "perennial ice" here means only *"this cell never
rose above ``T_f`` this year"*. Huang & Bowman ran realistic 2-D geography and found the instability present
in one hemisphere and absent in the other; this model has a uniform two-tile land/ocean mix and is run at
``f_L = 0`` (all ocean) for the headline, because the SICI is a statement about a **perennial** cap and rung
5B.1+ already established there is **no year-round land ice** at Earth insolation (the small-``C`` land tile
climbs above freezing every summer) — so an all-land sweep is a different question, not the mirror image.
One tilt, one albedo law, no ice thermodynamics: what is banked is the **mechanism**, not a general claim
that "the seasonal cycle removes the SICI".

A SIBLING module — the marcher is untouched
--------------------------------------------
Every rung here is a sibling (ADR 0005). This module *consumes* :class:`planet.seasonal.SeasonalEBM`
unchanged and reuses :class:`planet.bifurcation.EquilibriumCurve` for its fold algebra (the slope-stability
theorem, the parabola-refined turning points, the branch-end thresholds) rather than re-deriving it. The
single edit it needed upstream is :meth:`~planet.seasonal.SeasonalEBM.march`'s ``T_init`` growing a
per-tile **continuation seed** — a new branch beside the untouched scalar/``None`` paths, guarded by 5B.1+'s
two bit-identical reductions.

Units — SI, climlab-conventional: ``S₀`` in W m⁻², ``T`` in °C, ``x = sin φ`` on [−1, 1], cap radius
``θ = 90° − φ_ice`` in degrees, mixed-layer depth in m. Sources ([[seasonal-ebm-source]],
[[ebm-radiation-source]]): North 1984 (*J. Atmos. Sci.* 41, 3390 — the SICI in the annual-mean model);
Huang & Bowman 1992 (*Climate Dynamics* 7, 205 — the SICI in **seasonal** EBMs, the question this asks);
Wagner & Eisenman 2015 (*J. Climate* 28, 3998 — seasonal cycle and meridional transport stabilizing ice);
Cahalan & North 1979 (the slope-stability theorem the parent curve reads).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, replace
from typing import Optional, Sequence

import numpy as np
from scipy.linalg import solve_banded

from .albedo import EBMParams
from .bifurcation import EquilibriumCurve
from .ebm import A_OLR, B_OLR, D_TRANSPORT, S0_EARTH, T_FREEZE, ice_line_latitude
from .obliquity import OBLIQUITY_EARTH, insolation_s2
from .seasonal import (
    LAND_SOIL_DEPTH, OCEAN_MIXED_DEPTH, SeasonalClimate, SeasonalEBM,
    ice_coalbedo, ice_edge_latitude,
)

__all__ = [
    "SICIConfig", "SweepPoint", "Continuation", "HysteresisLoop", "SeedDependence", "DepthPoint",
    "annual_mean_curve", "perennial_ice_cells", "continuation_sweep", "hysteresis_loop", "plant_cap",
    "seed_dependence", "depth_sweep",
]


# --------------------------------------------------------------------------- #
# The experiment's configuration — everything the sweep holds fixed.
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class SICIConfig:
    """The fixed half of the experiment: grid, year, tilt, tiles, and the OLR/transport constants.

    Only ``S₀`` (swept by the continuation) and ``ocean_mixed_depth`` (the seasonal-amplitude axis) vary
    around this. ``land_fraction`` defaults to **0** — an all-ocean planet — because the instability is a
    statement about the **perennial** cap and there is no year-round land ice at these insolations (module
    docstring); the Earth-like ``0.30`` is the reference mix, not the headline.

    ``n_cells`` is the full-sphere cell count on ``x = sin φ ∈ [−1, 1]``, so the equivalent hemisphere
    resolution is ``n_cells / 2`` — :func:`planet.bifurcation.critical_cap_sweep` runs 720 *hemisphere*
    cells for this fold, i.e. ``n_cells = 1440`` here. The default 720 is the sweep-cost compromise; every
    headline claim is checked at 1440.
    """

    n_cells: int = 720
    n_steps: int = 180
    obliquity_deg: float = OBLIQUITY_EARTH
    land_fraction: float = 0.0
    ocean_mixed_depth: float = OCEAN_MIXED_DEPTH
    land_soil_depth: float = LAND_SOIL_DEPTH
    A: float = A_OLR
    B: float = B_OLR
    D: float = D_TRANSPORT
    tile: str = "ocean"

    def __post_init__(self):
        if self.tile not in ("land", "ocean", "mean"):
            raise ValueError(f"tile must be 'land', 'ocean' or 'mean', got {self.tile!r}")
        if self.n_cells % 2 != 0:
            raise ValueError(f"n_cells must be even (the equator must fall on a face), got {self.n_cells}")

    def model(self, S0: float = S0_EARTH, ocean_mixed_depth: Optional[float] = None) -> SeasonalEBM:
        """The :class:`~planet.seasonal.SeasonalEBM` at this configuration and solar constant."""
        return SeasonalEBM(A=self.A, B=self.B, D=self.D, land_fraction=self.land_fraction,
                           land_soil_depth=self.land_soil_depth,
                           ocean_mixed_depth=(self.ocean_mixed_depth if ocean_mixed_depth is None
                                              else float(ocean_mixed_depth)),
                           obliquity_deg=self.obliquity_deg, S0=float(S0),
                           n_cells=self.n_cells, n_steps=self.n_steps)

    @property
    def polar_cell_deg(self) -> float:
        """Latitude width of the **polar** cell (degrees) — the grid's cap-size quantum, ``≈ √(2Δx)``.

        The uniform-``x`` grid is coarsest in latitude exactly at the pole, where this fold lives: with
        ``Δx = 2/n_cells`` the outermost cell spans ``90° − asin(1 − Δx)``. It shrinks only as ``√Δx``, so
        it is quoted beside every cap radius — a cap smaller than this is a *sub-cell* read of a state the
        albedo feedback cannot actually resolve.
        """
        dx = 2.0 / self.n_cells
        return math.degrees(math.acos(max(-1.0, min(1.0, 1.0 - dx))))

    @property
    def cap_resolution_deg(self) -> float:
        """Half a polar cell — the detection threshold two cap radii must differ by to count as different.

        Below it, a difference between the up- and down-sweeps is grid/tolerance noise, not hysteresis;
        :attr:`HysteresisLoop.width` is defined against it so that *"no loop"* and *"a loop below
        detection"* never blur together.
        """
        return 0.5 * self.polar_cell_deg

    def cells_in_cap(self, cap_deg: float) -> int:
        """How many cells (one hemisphere) a cap of angular radius ``cap_deg`` covers on this grid.

        The scale a cell-count *step* has to be read against. A genuine fold at ``θ_c`` flips a whole
        ``θ_c``-sized cap on in **one** step of the solar constant, so ``cells_in_cap(θ_c)`` is the jump a
        fold would produce — and it grows as the grid refines. A step of one or two cells on a fine grid is
        therefore not "nearly a fold": it is an order of magnitude below one, and this is the number that
        says so. (Without it, "the count stepped by 2" would look worse at 1440 cells than "stepped by 1"
        at 720, when it is the *same physical cap growth* sampled on cells half as wide.)
        """
        x_edge = math.sin(math.radians(90.0 - float(cap_deg)))
        grid_x = self.model().x
        return int(((grid_x >= x_edge) & (grid_x >= 0.0)).sum())

    @property
    def params(self) -> EBMParams:
        """A rung-0 :class:`~planet.albedo.EBMParams` mirroring these constants — carried on the parent curve.

        ``s2`` is filled from the pinned :func:`planet.obliquity.insolation_s2` for this tilt so the record
        is self-describing; the curve itself never uses it (it integrates the model's **own** annual-mean
        insolation, not a P₂ truncation — the 5B.1 lesson).
        """
        return replace(EBMParams(), A=self.A, B=self.B, D=self.D, n_cells=self.n_cells // 2,
                       s2=insolation_s2(self.obliquity_deg))


# --------------------------------------------------------------------------- #
# The annual-mean parent — the dt-free reference the seasons are read against.
# --------------------------------------------------------------------------- #
def annual_mean_curve(cfg: SICIConfig) -> EquilibriumCurve:
    """The **annual-mean** equilibrium diagram of this seasonal model — the parent, solved exactly.

    Rung 0+'s inverse solve (:func:`planet.bifurcation.equilibrium_curve`), transplanted onto the
    *seasonal* model's own full-sphere operator and its own **annual-mean insolation**: prescribe the ice
    line at grid face ``k`` (a cap in **both** hemispheres, since the annual mean is symmetric), which makes
    the step-function albedo a known field and the EBM linear; one tridiagonal solve gives the profile per
    unit solar constant ``u(x)``; the ice-line condition ``T(x_s) = T_f`` closes it as
    ``S₀(x_s) = (T_f + A/B)/u(x_s)``. Sweeping the northern faces from the equator to the pole traces every
    equilibrium with an ice line, stable and unstable alike.

    Two deliberate choices make this the *exact parent of this marcher* rather than an approximation of it:

    * the transport operator is the model's **own** ``L_T`` (:meth:`~planet.seasonal.SeasonalEBM._transport_tridiag`,
      itself the engine's assembly), so the operator cannot drift between parent and child; and
    * the forcing is the **time-mean of the marcher's own seasonal insolation series**, not the P₂-truncated
      ``insolation()``. That is the 5B.1 reduction lesson: the true annual-mean insolation carries moments
      beyond ``P₂``, and matching against the truncation puts a ~1e-2 floor under a reduction that should be
      exact. (It agrees with the pinned :func:`planet.obliquity.annual_mean_insolation` kernel to machine
      precision at ``n_steps = 720`` — the same uniform-in-orbit sampling.)

    Returned as a :class:`~planet.bifurcation.EquilibriumCurve` so the whole fold apparatus — the
    slope-stability verdict ``dS₀/dx_s > 0``, the parabola-refined turning points, ``small_ice_cap_fold``,
    ``finite_cap_window``, ``stable_finite_cap_at`` — comes along unchanged.
    """
    m = cfg.model()
    s_unit = m.insolation_series().mean(axis=1) / m.S0          # ⟨S⟩/S₀ — the annual-mean profile
    ice_free_albedo = 1.0 - m.coalbedo()                        # a₀ + a₂P₂, the model's own
    a_ice = float(EBMParams().ai)
    sub, diag, sup = m._LT
    n = m.n_cells
    ab = np.zeros((3, n))
    ab[0, 1:] = sup[:-1]
    ab[1, :] = diag - m.B
    ab[2, :-1] = sub[1:]
    faces = np.asarray(m.grid.edges, dtype=float)
    Tf_shift = T_FREEZE + m.A / m.B

    ks = np.arange(n // 2, n + 1)                               # northern faces: equator (x=0) → pole (x=1)
    S0 = np.empty(ks.size)
    xs = np.empty(ks.size)
    T = np.empty((ks.size, n))
    Tbar = np.empty(ks.size)
    for j, k in enumerate(ks):
        alb = ice_free_albedo.copy()
        alb[k:] = a_ice                                          # the northern cap
        alb[:n - k] = a_ice                                      # its southern mirror
        u = solve_banded((1, 1), ab, -s_unit * (1.0 - alb))
        u_face = u[-1] if k == n else 0.5 * (u[k - 1] + u[k])    # 2nd order at the face; the pole by symmetry
        S0[j] = Tf_shift / u_face
        xs[j] = faces[k]
        T[j] = S0[j] * u - m.A / m.B
        Tbar[j] = float(np.mean(T[j]))                           # equal-area cells ⟹ mean = area mean
    stable = np.gradient(S0, xs) > 0.0                           # the slope-stability theorem
    return EquilibriumCurve(xs, S0, stable, T, Tbar, m.x, cfg.params)


# --------------------------------------------------------------------------- #
# Diagnostics on a converged limit cycle.
# --------------------------------------------------------------------------- #
def _tile_field(climate: SeasonalClimate, tile: str) -> np.ndarray:
    return {"land": climate.T_land, "ocean": climate.T_ocean, "mean": climate.T_mean}[tile]


def perennial_ice_cells(model: SeasonalEBM, climate: SeasonalClimate, tile: str = "ocean",
                        T_freeze: float = T_FREEZE) -> int:
    """Number of **northern** cells that stay below freezing for the *whole* year — the un-interpolated cap.

    The cap *radius* (:func:`planet.seasonal.ice_edge_latitude`) interpolates the temperature crossing and
    so returns sub-cell values; the albedo feedback, by contrast, only ever flips **whole cells**. This
    integer is therefore the honest read of the cap the model actually carries, and its step size between
    adjacent solar constants is what separates a genuinely continuous cap (``+1`` at a time) from a fold the
    interpolated radius was smoothing over (several cells at once). Counted on the northern hemisphere only
    (the seasonal cycle is hemispherically antisymmetric, so the two caps are mirror images half a year
    apart).
    """
    T = _tile_field(climate, tile)
    return int(((T.max(axis=1) < float(T_freeze)) & (model.x >= 0.0)).sum())


@dataclass(frozen=True)
class SweepPoint:
    """One converged annual limit cycle in a continuation sweep, reduced to its ice diagnostics.

    ``perennial_cap_deg`` is the angular radius ``90° − φ`` of the ice that survives the warmest month
    (0 = none); ``seasonal_cap_deg`` the radius of ice present at *some* point in the year (always the
    larger). ``n_perennial_cells`` is the same perennial cap counted in whole cells (the un-interpolated
    read). ``converged`` records whether the march reached its limit cycle — a fold must **never** be read
    from a point where it is ``False``, since a still-drifting cap is exactly what an unresolved jump looks
    like.
    """

    S0: float
    perennial_cap_deg: float
    seasonal_cap_deg: float
    n_perennial_cells: int
    global_mean_T: float
    polar_amplitude_K: float
    converged: bool
    years: int


@dataclass(frozen=True)
class Continuation:
    """A warm-started sweep of ``S₀`` in one direction, plus the state it ended on.

    ``points`` are in sweep order; ``final_state`` is the ``(T_land, T_ocean)`` day-0 pair the last march
    settled on, ready to seed the return leg. ``direction`` is ``"down"`` (dimming) or ``"up"`` (brightening).
    """

    cfg: SICIConfig
    points: tuple[SweepPoint, ...]
    direction: str
    final_state: tuple[np.ndarray, np.ndarray]

    @property
    def S0(self) -> np.ndarray:
        return np.array([p.S0 for p in self.points])

    @property
    def perennial_cap_deg(self) -> np.ndarray:
        return np.array([p.perennial_cap_deg for p in self.points])

    @property
    def seasonal_cap_deg(self) -> np.ndarray:
        return np.array([p.seasonal_cap_deg for p in self.points])

    @property
    def n_perennial_cells(self) -> np.ndarray:
        return np.array([p.n_perennial_cells for p in self.points])

    @property
    def all_converged(self) -> bool:
        return all(p.converged for p in self.points)

    @property
    def max_cell_jump(self) -> int:
        """Largest change in the perennial **cell count** between adjacent solar constants.

        ``1`` means the cap grew one cell at a time over the whole sweep — continuous at the resolution the
        feedback actually sees. Anything larger is a candidate fold (or too coarse an ``S₀`` step: check the
        sampling before reading it as physics).
        """
        n = self.n_perennial_cells
        return 0 if n.size < 2 else int(np.max(np.abs(np.diff(n))))


def continuation_sweep(cfg: SICIConfig, S0_values: Sequence[float],
                       seed: float | np.ndarray | tuple[np.ndarray, np.ndarray] = 25.0,
                       ocean_mixed_depth: Optional[float] = None,
                       tol: float = 1e-7, max_years: int = 800) -> Continuation:
    """March the ice-albedo seasonal EBM along ``S0_values``, each solve **warm-started** from the last.

    The continuation is the whole method: each solar constant resumes from its neighbour's converged day-0
    state rather than re-spinning up, which is both far cheaper *and* the thing that keeps the sweep **on a
    branch** — under the ice feedback the seed selects the climate, so a sweep that re-seeded uniformly at
    every point would sample whichever branch the seed happened to fall into instead of tracing the one the
    planet is actually on. Falling *off* a branch (the jump) is then the signal, not an artifact.

    ``seed`` is the first point's :meth:`~planet.seasonal.SeasonalEBM.march` seed — a scalar (uniform), a
    profile, or a ``(T_land, T_ocean)`` pair. ``tol``/``max_years`` are deliberately tighter and longer than
    the 5B.1+ defaults: critical slowing near a fold is exactly where a march quits early, and a
    still-drifting cap reads as a smooth shrink — i.e. the one systematic error that biases this experiment
    *toward* finding no instability. Every point carries its own ``converged`` flag for that reason.
    """
    S0_values = [float(s) for s in S0_values]
    state = seed
    points = []
    for S0 in S0_values:
        m = cfg.model(S0, ocean_mixed_depth=ocean_mixed_depth)
        c = m.march(coalbedo_fn=ice_coalbedo, T_init=state, tol=tol, max_years=max_years)
        state = (c.T_land[:, 0].copy(), c.T_ocean[:, 0].copy())
        T = _tile_field(c, cfg.tile)
        i_polar = m.nearest_index(80.0)
        points.append(SweepPoint(
            S0=S0,
            perennial_cap_deg=90.0 - ice_edge_latitude(m.x, T, kind="perennial"),
            seasonal_cap_deg=90.0 - ice_edge_latitude(m.x, T, kind="seasonal"),
            n_perennial_cells=perennial_ice_cells(m, c, cfg.tile),
            global_mean_T=float(T.mean()),
            polar_amplitude_K=float(c.amplitude(cfg.tile)[i_polar]),
            converged=bool(c.converged),
            years=int(c.years),
        ))
    direction = "down" if len(S0_values) > 1 and S0_values[1] < S0_values[0] else "up"
    return Continuation(cfg, tuple(points), direction, state)


# --------------------------------------------------------------------------- #
# Observable 1 — the hysteresis loop (the resolution-robust signature).
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class HysteresisLoop:
    """A down-sweep and the up-sweep that retraces it — and the ``S₀`` interval where they disagree.

    A fold makes the two legs differ over a **finite** band of solar constants: dimming, the ice-free planet
    holds past the sun at which a cap is already stable; brightening, the cap holds past the sun at which
    ice-free is. That band is the loop, and unlike "did the cap reach zero smoothly" it does not shrink as
    the grid refines — which is why it, not the cap radius, carries the verdict.
    """

    down: Continuation
    up: Continuation
    threshold_deg: float

    @property
    def S0(self) -> np.ndarray:
        """The common ``S₀`` samples, ascending."""
        return np.sort(self.down.S0)

    @property
    def dS0(self) -> float:
        """The (uniform) sampling interval — the finest loop this sweep could resolve."""
        s = self.S0
        return float(np.min(np.diff(s))) if s.size > 1 else 0.0

    @property
    def cap_gap_deg(self) -> np.ndarray:
        """``|cap_down − cap_up|`` at each ``S₀`` (ascending) — the raw disagreement between the legs."""
        d = dict(zip(np.round(self.down.S0, 9), self.down.perennial_cap_deg))
        u = dict(zip(np.round(self.up.S0, 9), self.up.perennial_cap_deg))
        return np.array([abs(d[k] - u[k]) for k in np.round(self.S0, 9)])

    @property
    def detected(self) -> bool:
        """Whether any sample's disagreement clears the detection threshold (half a polar cell)."""
        return bool(np.any(self.cap_gap_deg > self.threshold_deg))

    @property
    def interval(self) -> tuple[float, float] | None:
        """``(S₀_lo, S₀_hi)`` spanned by the disagreeing samples, or ``None`` if none clear the threshold."""
        mask = self.cap_gap_deg > self.threshold_deg
        if not mask.any():
            return None
        s = self.S0[mask]
        return (float(s.min()), float(s.max()))

    @property
    def width(self) -> float:
        """Loop width in W m⁻²: the span of disagreeing samples **plus one sampling interval**.

        The ``+ dS0`` is not padding — a single disagreeing sample already means the two branches part
        somewhere inside the neighbouring gaps, so the loop is at least one sample wide. Returns ``0.0``
        when nothing clears the threshold, which reads as *"no loop wider than* ``dS0`` *and no cap gap
        above half a polar cell"* — **not** as a proven zero. Quote it against :attr:`dS0` and the parent
        curve's own loop.
        """
        iv = self.interval
        return 0.0 if iv is None else float(iv[1] - iv[0] + self.dS0)

    @property
    def all_converged(self) -> bool:
        return self.down.all_converged and self.up.all_converged


def hysteresis_loop(cfg: SICIConfig, S0_values: Sequence[float],
                    warm_seed: float = 25.0, ocean_mixed_depth: Optional[float] = None,
                    tol: float = 1e-7, max_years: int = 800) -> HysteresisLoop:
    """Sweep ``S0_values`` **down** from an ice-free start, then **up** again from where that ended.

    The down leg starts bright and warm (no ice) and dims; the up leg resumes from the coldest, most-capped
    state the down leg reached and brightens back. If a fold exists the two legs cannot retrace each other
    across it. ``S0_values`` may be given in either order and must be uniformly spaced (the loop width is
    quoted against that spacing).
    """
    S0_asc = np.sort(np.asarray([float(s) for s in S0_values]))
    if S0_asc.size > 2:
        steps = np.diff(S0_asc)
        if not np.allclose(steps, steps[0], rtol=1e-6, atol=1e-9):
            raise ValueError("S0_values must be uniformly spaced (the loop width is read against the step)")
    down = continuation_sweep(cfg, S0_asc[::-1], seed=warm_seed, ocean_mixed_depth=ocean_mixed_depth,
                              tol=tol, max_years=max_years)
    up = continuation_sweep(cfg, S0_asc, seed=down.final_state, ocean_mixed_depth=ocean_mixed_depth,
                            tol=tol, max_years=max_years)
    return HysteresisLoop(down, up, cfg.cap_resolution_deg)


# --------------------------------------------------------------------------- #
# Observable 2 — seed dependence at one sun (bistability, independent of the sweep).
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class SeedDependence:
    """Plant a polar cap of a given size at one sun, and see whether the seasonal planet **keeps** it.

    The direct bistability test, and — unlike :class:`HysteresisLoop` — it never consults a sweep direction.
    ``warm_cap_deg`` is what an ice-free start settles on at this ``S₀``; ``survived_cap_deg`` is what a start
    with a cap of ``planted_cap_deg`` **already there** settles on. If the parent's fold is real, a planted
    cap at (or just above) the critical radius is a second stable climate and survives where the warm start
    has none. If the seasons have dissolved the fold, the planted cap simply melts back to whatever the warm
    start found, and the two agree.

    This matters because it can **dissociate** from the loop width: a warm-started continuation can walk
    continuously through a band where two climates coexist, pinned to one branch the whole way. Two
    observables, two ways to be fooled, and they are not fooled by the same thing.
    """

    S0: float
    planted_cap_deg: float
    warm_cap_deg: float
    survived_cap_deg: float
    threshold_deg: float
    converged: bool
    warm_polar_amplitude_K: float = float("nan")   # the seasonal swing at 80° the warm run carries

    @property
    def gap_deg(self) -> float:
        return abs(self.warm_cap_deg - self.survived_cap_deg)

    @property
    def survived(self) -> bool:
        """Whether the planted cap is still there (not melted away, not run away to a snowball)."""
        return self.converged and 0.0 < self.survived_cap_deg < 89.0

    @property
    def bistable(self) -> bool:
        """Two distinct climates at one sun — the planted cap survived *and* differs from the warm start."""
        return self.survived and self.gap_deg > self.threshold_deg


def plant_cap(model: SeasonalEBM, state: tuple[np.ndarray, np.ndarray], cap_deg: float,
              margin_K: float = 5.0, T_freeze: float = T_FREEZE) -> tuple[np.ndarray, np.ndarray]:
    """Cool both poles below freezing out to angular radius ``cap_deg`` — a seed with a cap already in it.

    Takes a converged day-0 ``(T_land, T_ocean)`` pair and drops every cell poleward of ``90° − cap_deg``
    (in **both** hemispheres) to ``T_f − margin_K``, leaving the rest of the state untouched. The point of
    perturbing an *already converged seasonal* state rather than seeding a fresh profile is that everything
    except the cap stays on the attractor — so what the march then does is a statement about the cap, not
    about a transient. (Seeding instead with the annual-mean parent's own profile does **not** work and is
    not a subtle failure: a seasonal planet handed an annual-mean field is effectively cold-started, because
    winter at high latitude runs tens of K below the annual mean, and the ice-albedo feedback takes it
    straight to a snowball.)
    """
    x_edge = math.sin(math.radians(90.0 - float(cap_deg)))
    cold = np.abs(model.x) >= x_edge
    seeded = []
    for T in state:
        T = np.array(T, dtype=float, copy=True)
        T[cold] = float(T_freeze) - float(margin_K)
        seeded.append(T)
    return seeded[0], seeded[1]


def seed_dependence(cfg: SICIConfig, S0: float, cap_deg: Optional[float] = None,
                    curve: Optional[EquilibriumCurve] = None, warm_seed: float = 25.0,
                    ocean_mixed_depth: Optional[float] = None, margin_K: float = 5.0,
                    tol: float = 1e-7, max_years: int = 2000) -> SeedDependence:
    """Does a cap of ``cap_deg`` survive at this sun, when a warm start would not have grown one?

    The seasonal counterpart of :func:`planet.bifurcation.relax_from_curve` — rung 0+ checked the
    slope-stability theorem by seeding the nonlinear relaxation *at* a curve equilibrium and watching
    whether it stayed. Here the cap is planted into an already-converged seasonal limit cycle
    (:func:`plant_cap`) and the march decides. ``cap_deg`` defaults to the parent curve's critical radius
    ``θ_c`` — the smallest cap the **annual-mean** model can hold — which is precisely the size the seasonal
    model is being asked about. ``curve`` defaults to this configuration's own :func:`annual_mean_curve`
    (pass a precomputed one to avoid re-solving it).
    """
    c = annual_mean_curve(cfg) if curve is None else curve
    if cap_deg is None:
        fold = c.small_ice_cap_fold
        cap_deg = fold.cap_radius_deg if fold is not None else 2.0 * cfg.polar_cell_deg
    m = cfg.model(S0, ocean_mixed_depth=ocean_mixed_depth)
    warm = m.march(coalbedo_fn=ice_coalbedo, T_init=float(warm_seed), tol=tol, max_years=max_years)
    planted = plant_cap(m, (warm.T_land[:, 0], warm.T_ocean[:, 0]), float(cap_deg), margin_K=margin_K)
    capped = m.march(coalbedo_fn=ice_coalbedo, T_init=planted, tol=tol, max_years=max_years)
    return SeedDependence(
        S0=float(S0),
        planted_cap_deg=float(cap_deg),
        warm_cap_deg=90.0 - ice_edge_latitude(m.x, _tile_field(warm, cfg.tile), kind="perennial"),
        survived_cap_deg=90.0 - ice_edge_latitude(m.x, _tile_field(capped, cfg.tile), kind="perennial"),
        threshold_deg=cfg.cap_resolution_deg,
        converged=bool(warm.converged and capped.converged),
        warm_polar_amplitude_K=float(warm.amplitude(cfg.tile)[m.nearest_index(80.0)]),
    )


# --------------------------------------------------------------------------- #
# The headline — both observables against the seasonal amplitude (mixed-layer depth).
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class DepthPoint:
    """Both observables at one ocean mixed-layer depth, plus the seasonal amplitude that depth buys."""

    depth_m: float
    loop_width: float
    loop_detected: bool
    seed_gap_deg: float
    bistable: bool
    polar_amplitude_K: float
    max_cell_jump: int
    converged: bool


def depth_sweep(depths: Sequence[float], cfg: SICIConfig, S0_values: Sequence[float],
                seed_S0: Optional[float] = None, curve: Optional[EquilibriumCurve] = None,
                tol: float = 1e-7, max_years: int = 4000) -> tuple[DepthPoint, ...]:
    """``(loop width, seed dependence)`` versus mixed-layer depth — the experiment's headline.

    Deepening the mixed layer damps the seasonal swing (``amplitude ∝ 1/√(B² + ω²C²)``) without touching the
    annual-mean forcing, so this sweep walks from the seasonal planet toward its own annual-mean parent
    along a single axis. **The approach is asymptotic, not attained**: the swing falls only as ``1/(ωC)``
    while the spin-up cost grows *linearly* in ``C``, so the deepest tractable point still carries a residual
    seasonal cycle and a residual gap from the parent's loop. Read the trend and quote the gap; do not
    expect the parent's numbers to be reproduced exactly at any finite depth.

    ``seed_S0`` (default: the parent's SICI fold, the sun where bistability is widest) is where
    :func:`seed_dependence` is evaluated. ``max_years`` is generous because a deep mixed layer has a
    relaxation time ``τ = C/B`` of decades — a deep point that quits early would look like a *stable* one,
    the same bias the loop's tolerance guards against.
    """
    c = annual_mean_curve(cfg) if curve is None else curve
    fold = c.small_ice_cap_fold
    S0_seed = float(seed_S0) if seed_S0 is not None else (
        fold.S0 if fold is not None else float(c.ice_free_threshold_S0))
    out = []
    for depth in [float(d) for d in depths]:
        loop = hysteresis_loop(cfg, S0_values, ocean_mixed_depth=depth, tol=tol, max_years=max_years)
        sd = seed_dependence(cfg, S0_seed, curve=c, ocean_mixed_depth=depth, tol=tol, max_years=max_years)
        amp = float(np.median([p.polar_amplitude_K for p in loop.down.points]))
        out.append(DepthPoint(
            depth_m=depth,
            loop_width=loop.width,
            loop_detected=loop.detected,
            seed_gap_deg=sd.gap_deg,
            bistable=sd.bistable,
            polar_amplitude_K=amp,
            max_cell_jump=max(loop.down.max_cell_jump, loop.up.max_cell_jump),
            converged=loop.all_converged and sd.converged,
        ))
    return tuple(out)
