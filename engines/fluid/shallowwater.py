"""2-D rotating shallow-water solver — the program's second shared engine.

Solves the **rotating-frame shallow-water equations** on a doubly-periodic
β-plane, in **vector-invariant** form,

    ∂u/∂t = (f + ζ) v − ∂_x B            (Coriolis + vorticity flux)
    ∂v/∂t = −(f + ζ) u − ∂_y B
    ∂h/∂t = −∇·(h u)                      (mass continuity, flux form)

with  ζ = ∂_x v − ∂_y u  the relative vorticity,
      B = g(h + h_b) + ½(u² + v²)         the Bernoulli function,
      f = f₀ + β·(y − y_ref)              the β-plane Coriolis parameter,

where ``h`` is the fluid-layer thickness, ``(u, v)`` the depth-independent
horizontal velocity, ``g`` gravity (or reduced gravity for a baroclinic mode),
and ``h_b`` an optional bottom topography. The absolute vorticity ``f + ζ`` and
the **potential vorticity** ``q = (f + ζ)/h`` are the rotation-aware quantities
the validation suite watches.

This is a **hyperbolic, explicit** solver — it shares *no* machinery with the
parabolic-implicit :mod:`engines.diffusion`. That is the point: the portfolio's
second engine is a genuinely different solver class (a CFL-limited wave solver,
not an unconditionally-stable tridiagonal one), built and frozen here so Planet's
circulation (Phase 4) and the documented GCM climb can reuse it behind a contract.

Discretization
--------------
* **Arakawa C-grid (staggered).** Thickness ``h`` (and ``h_b``) live at cell
  *centers*; ``u`` on east–west *faces*; ``v`` on north–south *faces*; relative
  vorticity ``ζ``, the Coriolis parameter ``f``, and the potential vorticity ``q``
  at cell *corners*. The staggering is what lets a centered scheme represent
  geostrophic balance and the Rossby radius without a computational checkerboard
  mode — and it is exactly where a *wrong Coriolis averaging* hides, which is why
  the PV / potential-enstrophy seal (a finite-amplitude test) is the discriminating
  validation leg, not the gravity-wave speed.
* **Vector-invariant form with a symmetric corner-PV Coriolis flux.** The rotation
  term is written as the PV ``q`` (at corners) times the mass flux, averaged
  symmetrically to the velocity points (Sadourny-1975 lineage). This engine does
  **not** claim Sadourny's *exact* semi-discrete enstrophy conservation (that needs
  a specific averaging; Arakawa–Lamb conserves both energy and enstrophy at a large
  step up in complexity). What it *guarantees and the suite measures* is: **mass
  conserved to machine precision** (flux-form continuity telescopes on the periodic
  domain), and **energy and potential enstrophy conserved to a bounded, dt→0
  convergent drift** — the honest claim for a centered explicit scheme.
* **SSP-RK3 (Shu–Osher) explicit time stepping.** Three-stage, third-order, with a
  stability region that includes part of the imaginary axis — so the *non-dissipative*
  centered advection + gravity-wave + Coriolis operator is stable up to a CFL based on
  the fastest signal ``√(gH) + |u|`` (gravity-wave CFL binds; ``f·dt ≪ 1`` is automatic,
  so the explicit Coriolis needs no semi-implicit treatment).

The frozen data boundary (ADR 0001)
-----------------------------------
The ``state`` is a :class:`SWState` — three plain 2-D ``ndarray`` fields
``(h, u, v)`` of identical shape ``(ny, nx)``, *stacked fields* and nothing more.
That bundle is the frozen data contract: :meth:`ShallowWater.step` / :meth:`solve`
consume and return exactly it, the diagnostics consume it, and no live object
crosses the per-step boundary. The grid, ``g``, ``H``, ``f₀``, ``β``, and topography
are **construction-time configuration** (they reduce to numbers / arrays during RHS
assembly), exactly as the diffusion engine treats its ``Grid``/``D``/BCs. The
stacked-field shape is the GCM-climb seam: **N vertical layers** is a leading-axis
extension, and an **advected tracer** is one more field — both contract *extensions*
that do not change the per-step array boundary. The tracer slot is declared on
:class:`SWState` (``tracer``) but **not advected in v1** (:meth:`step` raises if it
is set), naming rung 1 of the §5 staircase — build the seam, not the machinery.

Units & sign convention
-----------------------
**SI throughout** (the engine is unit-agnostic in the diffusion sense — the consumer
supplies the constants — but the natural units are SI): ``h, h_b, H`` in m, ``u, v``
in m/s, ``g`` in m/s², ``f₀`` in 1/s, ``β`` in 1/(m·s), lengths in m, time in s. The
domain is doubly periodic; ``f = f₀ + β·(y − y_ref)`` (β-plane). ``y_ref`` defaults to
the domain center so ``f₀`` is the central-latitude value. Rigid channel walls in y
are a *named, unbuilt* BC extension (the periodic β-plane is the v1 geometry).
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Optional

import numpy as np


# --------------------------------------------------------------------------- #
# Grid (plain numeric data — not a stateful object), mirroring engines.diffusion
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Grid2D:
    """Uniform, doubly-periodic 2-D finite-difference grid for the C-grid solver.

    ``nx``/``ny`` cells of size ``dx``/``dy`` span a domain ``Lx = nx·dx`` by
    ``Ly = ny·dy``. Cell *centers* (where ``h`` lives) are at
    ``x = (i+½)dx, y = (j+½)dy``; cell *corners* (where ``ζ``/``f``/``q`` live) at
    ``x = i·dx, y = j·dy`` (the SW corner of cell ``[j, i]``). Velocity faces are at
    the obvious half-offsets (see the module docstring). Arrays are indexed
    ``[j, i]`` (row = y, col = x); periodicity is realized with ``np.roll``.
    """

    nx: int
    ny: int
    dx: float
    dy: float

    @property
    def Lx(self) -> float:
        return self.nx * self.dx

    @property
    def Ly(self) -> float:
        return self.ny * self.dy

    @property
    def cell_area(self) -> float:
        return self.dx * self.dy

    def x_centers(self) -> np.ndarray:
        return (np.arange(self.nx) + 0.5) * self.dx

    def y_centers(self) -> np.ndarray:
        return (np.arange(self.ny) + 0.5) * self.dy

    def x_corners(self) -> np.ndarray:
        return np.arange(self.nx) * self.dx

    def y_corners(self) -> np.ndarray:
        return np.arange(self.ny) * self.dy

    def center_mesh(self):
        """``(X, Y)`` meshgrids of cell-center coordinates, shape ``(ny, nx)``."""
        return np.meshgrid(self.x_centers(), self.y_centers())


def uniform_grid(Lx: float, Ly: float, nx: int, ny: int) -> Grid2D:
    """A uniform doubly-periodic grid spanning ``[0, Lx] × [0, Ly]`` with ``nx×ny`` cells."""
    if nx < 2 or ny < 2:
        raise ValueError("need at least 2 cells in each direction")
    return Grid2D(nx=int(nx), ny=int(ny), dx=float(Lx) / nx, dy=float(Ly) / ny)


# --------------------------------------------------------------------------- #
# State — the frozen data boundary: stacked plain 2-D fields
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class SWState:
    """Shallow-water prognostic state: stacked plain 2-D fields on the C-grid.

    ``h`` (layer thickness, cell centers), ``u`` (east–west faces), ``v`` (north–south
    faces) — all shape ``(ny, nx)``, the frozen array boundary (ADR 0001). ``tracer``
    is the **declared-but-unbuilt** advected-scalar slot (rung 1 of the GCM climb,
    §5): carrying it is a contract extension; :meth:`ShallowWater.step` *raises* if it
    is set, naming the deferral rather than silently failing to advect it (the same
    "seam, not machinery" idiom as the map's unpainted ``vector_overlay``).
    """

    h: np.ndarray
    u: np.ndarray
    v: np.ndarray
    tracer: Optional[np.ndarray] = None

    def copy(self) -> "SWState":
        return SWState(
            h=np.array(self.h, dtype=float),
            u=np.array(self.u, dtype=float),
            v=np.array(self.v, dtype=float),
            tracer=None if self.tracer is None else np.array(self.tracer, dtype=float),
        )


# --------------------------------------------------------------------------- #
# Finite-difference shifts (periodic neighbours via np.roll)
# --------------------------------------------------------------------------- #
# Convention: arrays are [j, i] = [y, x]. A field value at cell/face/corner index
# (j, i) and the named neighbour:
def _xp(a: np.ndarray) -> np.ndarray:   # value at i+1 (east)
    return np.roll(a, -1, axis=1)


def _xm(a: np.ndarray) -> np.ndarray:   # value at i-1 (west)
    return np.roll(a, 1, axis=1)


def _yp(a: np.ndarray) -> np.ndarray:   # value at j+1 (north)
    return np.roll(a, -1, axis=0)


def _ym(a: np.ndarray) -> np.ndarray:   # value at j-1 (south)
    return np.roll(a, 1, axis=0)


# --------------------------------------------------------------------------- #
# Solver
# --------------------------------------------------------------------------- #
class ShallowWater:
    """Rotating shallow-water solver on a doubly-periodic β-plane; see the module docstring.

    Parameters
    ----------
    grid : Grid2D
        The uniform doubly-periodic C-grid (use :func:`uniform_grid`).
    g : float
        Gravitational acceleration (m/s²); a *reduced* gravity ``g'`` for a baroclinic
        layer mode. With ``mean_depth`` it fixes the gravity-wave speed ``c = √(gH)``.
    mean_depth : float
        Resting layer depth ``H`` (m). Sets the gravity-wave speed and (with ``f₀``)
        the Rossby radius ``L_R = √(gH)/f₀``.
    f0 : float
        Coriolis parameter at ``y_ref`` (1/s). ``0.0`` gives an f-plane with no rotation.
    beta : float
        Meridional Coriolis gradient ``β = df/dy`` (1/(m·s)); ``0.0`` is an f-plane.
    y_ref : float, optional
        Reference latitude for ``f = f₀ + β(y − y_ref)``; defaults to the domain center.
    bottom : ndarray, optional
        Bottom topography ``h_b`` at cell centers (m); ``None`` is a flat bottom.
    """

    def __init__(
        self,
        grid: Grid2D,
        g: float,
        mean_depth: float,
        f0: float = 0.0,
        beta: float = 0.0,
        y_ref: Optional[float] = None,
        bottom: Optional[np.ndarray] = None,
    ) -> None:
        if g <= 0.0:
            raise ValueError(f"g must be positive, got {g}")
        if mean_depth <= 0.0:
            raise ValueError(f"mean_depth H must be positive, got {mean_depth}")
        self.grid = grid
        self.g = float(g)
        self.H = float(mean_depth)
        self.f0 = float(f0)
        self.beta = float(beta)
        self.y_ref = float(0.5 * grid.Ly if y_ref is None else y_ref)
        if bottom is None:
            self.h_b = np.zeros((grid.ny, grid.nx))
        else:
            self.h_b = np.asarray(bottom, dtype=float)
            if self.h_b.shape != (grid.ny, grid.nx):
                raise ValueError(f"bottom must have shape {(grid.ny, grid.nx)}, got {self.h_b.shape}")
        # Coriolis parameter at CORNERS (where ζ and q live): f = f0 + β(y_corner − y_ref).
        y_corner = grid.y_corners()                       # length ny, at j·dy
        self.f_corner = (self.f0 + self.beta * (y_corner - self.y_ref))[:, None] * np.ones((1, grid.nx))

    # -- derived physical scales -------------------------------------------- #
    @property
    def gravity_wave_speed(self) -> float:
        """The (long) gravity-wave speed ``c = √(gH)`` (m/s)."""
        return float(np.sqrt(self.g * self.H))

    @property
    def rossby_radius(self) -> float:
        """The deformation radius ``L_R = √(gH)/|f₀|`` (m); ``inf`` on an f-plane with f₀=0."""
        if self.f0 == 0.0:
            return float("inf")
        return self.gravity_wave_speed / abs(self.f0)

    def f_at_corners(self) -> np.ndarray:
        """The β-plane Coriolis parameter at cell corners, shape ``(ny, nx)`` (1/s)."""
        return self.f_corner

    # -- diagnostics on the corner / center grids --------------------------- #
    def relative_vorticity(self, state: SWState) -> np.ndarray:
        """Relative vorticity ``ζ = ∂_x v − ∂_y u`` at cell corners (1/s).

        Corner ``[j, i]`` (SW corner of cell ``[j, i]``): ``(v[j,i]−v[j,i−1])/dx −
        (u[j,i]−u[j−1,i])/dy``. ``∫ζ dA = 0`` to machine precision on the periodic
        domain (a discrete curl), the structural fact behind circulation conservation.
        """
        u, v = state.u, state.v
        return (v - _xm(v)) / self.grid.dx - (u - _ym(u)) / self.grid.dy

    def _h_corner(self, h: np.ndarray) -> np.ndarray:
        """Thickness averaged to cell corners (4-point), shape ``(ny, nx)``."""
        return 0.25 * (h + _xm(h) + _ym(h) + _xm(_ym(h)))

    def potential_vorticity(self, state: SWState) -> np.ndarray:
        """Potential vorticity ``q = (f + ζ)/h`` at cell corners (1/(m·s)).

        The materially-conserved invariant of inviscid shallow water (``Dq/Dt = 0``):
        the discriminating, rotation-aware quantity. A wrong ``f``/``β`` or a wrong
        Coriolis averaging shows up as spurious PV extrema or gross potential-enstrophy
        drift in the finite-amplitude seal.
        """
        return (self.f_corner + self.relative_vorticity(state)) / self._h_corner(state.h)

    def mass(self, state: SWState) -> float:
        """Total mass ``∫h dA = Σ hᵢⱼ ΔxΔy`` — conserved to machine precision (flux form)."""
        return float(np.sum(state.h) * self.grid.cell_area)

    def _kinetic_energy_density(self, state: SWState) -> np.ndarray:
        """Kinetic energy per unit area ``½ h(u²+v²)`` at cell centers (J/m² with ρ=1)."""
        u, v = state.u, state.v
        u2c = 0.5 * (u ** 2 + _xp(u) ** 2)        # u² averaged to centers (its two x-faces)
        v2c = 0.5 * (v ** 2 + _yp(v) ** 2)        # v² averaged to centers (its two y-faces)
        return 0.5 * state.h * (u2c + v2c)

    def energy(self, state: SWState) -> float:
        """Total energy ``∫[½h(u²+v²) + ½g h²] dA`` (KE + PE; ρ=1) — bounded, dt→0 convergent drift.

        Not a machine-precision invariant of the centered scheme (see the module
        docstring): the suite asserts the drift is small and *converges* as dt→0, the
        honest conservation claim for an explicit centered shallow-water solver.

        Note this is the **total** energy, whose PE carries the large resting ``½gH²`` background;
        a *relative* drift ``|ΔE|/E`` is therefore conservative w.r.t. the (much smaller) dynamical
        energy. The dt³-convergence ratio is unaffected (the constant background cancels).
        """
        pe = 0.5 * self.g * state.h ** 2
        return float(np.sum(self._kinetic_energy_density(state) + pe) * self.grid.cell_area)

    def potential_enstrophy(self, state: SWState) -> float:
        """Total potential enstrophy ``∫ ½ q² h dA`` (corner grid) — bounded, dt→0 convergent drift.

        The rotation-aware integral invariant; its finite-amplitude conservation is the
        Coriolis seal. ``h`` and ``q`` are evaluated at corners (``h_corner`` is the
        4-point average), consistent with where ``q`` is defined.
        """
        q = self.potential_vorticity(state)
        return float(np.sum(0.5 * q ** 2 * self._h_corner(state.h)) * self.grid.cell_area)

    # -- the right-hand side (method of lines) ------------------------------ #
    def _rhs(self, h: np.ndarray, u: np.ndarray, v: np.ndarray):
        """Tendencies ``(dh/dt, du/dt, dv/dt)`` of the vector-invariant shallow-water system.

        All quantities live on their C-grid positions; ``np.roll`` realizes the periodic
        neighbour stencils described in the module docstring.
        """
        dx, dy = self.grid.dx, self.grid.dy

        # Mass fluxes at velocity faces: U = h_u·u (u-points), V = h_v·v (v-points).
        h_u = 0.5 * (h + _xm(h))                  # h at u-point (i−½,j): cells i, i−1
        h_v = 0.5 * (h + _ym(h))                  # h at v-point (i,j−½): cells j, j−1
        U = h_u * u
        V = h_v * v

        # Continuity (flux form): ∂h/∂t = −[(U_east−U_west)/dx + (V_north−V_south)/dy].
        # East face of cell i is the u-point of cell i+1, hence _xp(U); likewise _yp(V).
        dh = -((_xp(U) - U) / dx + (_yp(V) - V) / dy)

        # Potential vorticity at corners.
        zeta = (v - _xm(v)) / dx - (u - _ym(u)) / dy
        h_corner = 0.25 * (h + _xm(h) + _ym(h) + _xm(_ym(h)))
        q = (self.f_corner + zeta) / h_corner

        # Bernoulli function at centers: B = g(h + h_b) + ½(u²+v²).
        u2c = 0.5 * (u ** 2 + _xp(u) ** 2)
        v2c = 0.5 * (v ** 2 + _yp(v) ** 2)
        B = self.g * (h + self.h_b) + 0.5 * (u2c + v2c)

        # Coriolis + vorticity flux, symmetric corner-PV averaging to the velocity points.
        # u-point (i−½,j): + ¼[ q_S·(V + V_W) + q_N·(V_N + V_NW) ], q_S=q[j,i], q_N=q[j+1,i].
        cor_u = 0.25 * (
            q * (V + _xm(V)) + _yp(q) * (_yp(V) + _xm(_yp(V)))
        )
        # v-point (i,j−½): − ¼[ q_W·(U + U_S) + q_E·(U_E + U_ES) ], q_W=q[j,i], q_E=q[j,i+1].
        cor_v = -0.25 * (
            q * (U + _ym(U)) + _xp(q) * (_xp(U) + _ym(_xp(U)))
        )

        du = cor_u - (B - _xm(B)) / dx
        dv = cor_v - (B - _ym(B)) / dy
        return dh, du, dv

    # -- time stepping (SSP-RK3) -------------------------------------------- #
    def step(self, state: SWState, dt: float) -> SWState:
        """Advance ``state`` by one SSP-RK3 step ``dt``; returns the new state (no mutation).

        Raises if ``dt`` is non-positive, if a tracer is set (the unbuilt rung-1 slot),
        or if ``dt`` exceeds the CFL limit (a wrong-by-design large step would blow up —
        the explicit-solver analogue of the diffusion engine's stability guarantee, but
        here *conditional*, so it is enforced rather than promised).
        """
        if dt <= 0.0:
            raise ValueError("dt must be positive")
        if state.tracer is not None:
            raise NotImplementedError(
                "tracer advection is the rung-1 extension (GCM climb §5); v1 carries the "
                "slot on SWState but does not advect it"
            )
        # Guard the *true* stability threshold (Courant ≈ 1 on the fastest signal), not the
        # conservative recommended step (safety≈0.3) — so passing max_dt(0.3) never false-trips.
        cfl_limit = self.max_dt(state, safety=1.0)
        if dt > cfl_limit:
            raise ValueError(
                f"dt={dt:g} exceeds the CFL stability limit {cfl_limit:g} "
                "(= min(dx,dy)/(√(gH)+|u|)); reduce dt (the recommended step is max_dt())"
            )
        h, u, v = state.h, state.u, state.v

        # SSP-RK3 (Shu–Osher): three convex-combination stages.
        k1h, k1u, k1v = self._rhs(h, u, v)
        h1, u1, v1 = h + dt * k1h, u + dt * k1u, v + dt * k1v

        k2h, k2u, k2v = self._rhs(h1, u1, v1)
        h2 = 0.75 * h + 0.25 * (h1 + dt * k2h)
        u2 = 0.75 * u + 0.25 * (u1 + dt * k2u)
        v2 = 0.75 * v + 0.25 * (v1 + dt * k2v)

        k3h, k3u, k3v = self._rhs(h2, u2, v2)
        hn = (1.0 / 3.0) * h + (2.0 / 3.0) * (h2 + dt * k3h)
        un = (1.0 / 3.0) * u + (2.0 / 3.0) * (u2 + dt * k3u)
        vn = (1.0 / 3.0) * v + (2.0 / 3.0) * (v2 + dt * k3v)
        return SWState(h=hn, u=un, v=vn)

    def max_dt(self, state: SWState, safety: float = 0.3) -> float:
        """CFL-stable step ``safety · min(dx,dy)/(√(gH) + |u|_max)`` (s).

        The gravity-wave speed ``√(gH)`` plus the peak flow speed is the fastest signal;
        ``safety ≈ 0.3`` is a conservative RK3 C-grid factor. ``f·dt`` is automatically
        ``≪ 1`` at this step (mid-latitude ``f ~ 1e−4``), so the explicit Coriolis is stable
        without a semi-implicit treatment.
        """
        speed = self.gravity_wave_speed + float(np.max(np.sqrt(state.u ** 2 + state.v ** 2)))
        return safety * min(self.grid.dx, self.grid.dy) / speed

    def solve(self, state: SWState, t_end: float, dt: Optional[float] = None,
              safety: float = 0.3) -> SWState:
        """March from ``0`` to ``t_end`` (s). With ``dt=None`` uses the CFL step (fixed at the
        initial state's speed); otherwise the given ``dt`` (last step trimmed to land on ``t_end``).

        ``dt=None`` assumes a **non-accelerating** flow (the fixed step is set once from the initial
        speed); for a flow that spins up well past its initial speed (e.g. a jet forced from rest),
        pass an explicit ``dt`` or step in segments, or the per-step CFL guard will raise mid-run.
        """
        if t_end < 0.0:
            raise ValueError("t_end must be non-negative")
        if dt is None:
            dt = self.max_dt(state, safety=safety)
        s = state
        t = 0.0
        while t < t_end - 1e-12 * max(1.0, t_end):
            h = min(dt, t_end - t)
            s = self.step(s, h)
            t += h
        return s
