"""N-layer rotating shallow-water — the **baroclinic** extension of :mod:`engines.fluid`.

The single-layer :class:`~engines.fluid.shallowwater.ShallowWater` is *dry*: it has no
available potential energy and no vertical shear, so it **categorically cannot** be
baroclinically unstable. This module stacks **N free-surface SW layers** of slightly
different density (Phillips 1954) — coupled *only* through the hydrostatic **Montgomery
pressure** — which is the minimal honest vertical structure that supports real storms
(rung 3 of the GCM climb, §5). The **interface displacement is the dynamical
temperature/buoyancy**, so heat transport becomes *intrinsic to the dynamics* (the leap past
the rung-1 *passive* tracer).

Design (the CONTRACT's promised leading-axis seam, realized as a **sibling** engine):

* **Leading layer axis.** :class:`LayeredState` stacks the same plain C-grid fields
  ``(h, u, v)`` with a leading layer axis — shape ``(nl, ny, nx)``. The per-layer
  finite-difference stencils are *identical* to the single-layer engine (vectorized over the
  layer axis); the **only** inter-layer coupling is the Montgomery pressure stack.
* **Single-layer reduction is byte-identical.** ``LayeredShallowWater`` with ``nl=1`` (and no
  background) reproduces ``ShallowWater`` **bit-for-bit** (the by-construction rung-0 reduction
  that has held every rung); the test asserts ``np.array_equal``. The single-layer engine is
  left *untouched* — this is a sibling, so existing consumers are unaffected.
* **Optional thermal-wind background (default off).** A doubly-periodic domain cannot carry a
  uniform meridional gradient *as a field* (it would seam at the wrap — that is physics, not a
  bug). So an unstable basic state enters as **constant background coefficients**: per-layer
  mean zonal flow ``U_k`` (a Doppler term ``−U_k ∂_x``) and the thermal-wind thickness gradient
  ``G_k`` (a baroclinic term ``−G_k v'`` in continuity), with the **prognostic fields being the
  perturbations**. ``background=None`` → the plain nonlinear engine (and the bit-for-bit
  reduction). ``background`` set → the spike's linear terms **plus** the engine's eddy-eddy
  nonlinearity: at small amplitude the most-unstable mode grows at the analytic linear rate
  (Phase A, anchored on :mod:`engines.fluid.stability`); at large amplitude it saturates
  (Phase B — the standard "fixed mean shear + prognostic eddies" two-layer turbulence setup).

Scope at this rung (Phase A — the linear growth rate; the named, unbuilt rest):

* **No dissipation / no flux limiter** — inherited from the energy-conserving single-layer
  engine. Phase A runs are small-amplitude and short (pre-saturation), so this is fine; the
  saturated turbulent Phase-B run needs the named-but-unbuilt **hyperviscosity** (a turbulent
  run cascades enstrophy to the grid). The RHS is structured to leave room for that operator.
* **f-plane background.** The background-balanced mode is consistent only on the f-plane
  (``β=0``), matching the f-plane stability operator; ``β`` may be set for the dry/wave dynamics
  but the growth-rate anchor uses ``β=0`` (a β-capable PV-gradient treatment with a *finite*
  critical shear is the named within-rung extension).
* **No passive tracer** — at this rung the interface displacement *is* the temperature, so the
  rung-1 scalar tracer is not carried (a named later extension).
* **External-mode CFL.** The free surface carries a fast barotropic gravity wave ``√(g·H_tot)``
  that sets the explicit step while the slow baroclinic mode is all that matters (the named cost;
  the rigid-lid elliptic-solve fork is the named within-rung upgrade if it bites).

Units: SI throughout, as the single-layer engine. Sources (extending ``[[shallow-water-source]]``):
Vallis 2017 *AOFD* / Cushman-Roisin & Beckers (two-layer SW); Phillips 1954 / Eady 1949.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np

from .shallowwater import Grid2D
from .stability import TwoLayerStability


# --------------------------------------------------------------------------- #
# Finite-difference shifts on the LAST two axes (so they vectorize over a leading
# layer axis AND reduce, for a 2-D array, to the single-layer engine's stencils).
# --------------------------------------------------------------------------- #
def _xp(a: np.ndarray) -> np.ndarray:   # east  (i+1)
    return np.roll(a, -1, axis=-1)


def _xm(a: np.ndarray) -> np.ndarray:   # west  (i-1)
    return np.roll(a, 1, axis=-1)


def _yp(a: np.ndarray) -> np.ndarray:   # north (j+1)
    return np.roll(a, -1, axis=-2)


def _ym(a: np.ndarray) -> np.ndarray:   # south (j-1)
    return np.roll(a, 1, axis=-2)


# --------------------------------------------------------------------------- #
# State — the leading-axis realization of the stacked-field boundary (ADR 0001)
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class LayeredState:
    """N-layer shallow-water prognostic state: stacked C-grid fields with a leading layer axis.

    ``h`` (layer thickness, cell centers), ``u`` (E–W faces), ``v`` (N–S faces) — each of shape
    ``(nl, ny, nx)``, layer 0 the **top**. Same stable array boundary as :class:`SWState`, with
    the leading layer axis (ADR 0001 / the GCM-climb seam). ``h`` is the **total** thickness
    ``H_k + h'`` (so the divergence and pressure see the full column); with a background, ``u, v``
    are the **perturbation** velocities (the mean ``U_k`` rides in the engine's Doppler coefficient,
    not the field — see the module docstring).
    """

    h: np.ndarray
    u: np.ndarray
    v: np.ndarray

    def copy(self) -> "LayeredState":
        return LayeredState(
            h=np.array(self.h, dtype=float),
            u=np.array(self.u, dtype=float),
            v=np.array(self.v, dtype=float),
        )


# --------------------------------------------------------------------------- #
# Background — the thermal-wind basic state, injected as constant coefficients
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class ThermalWindBackground:
    """The unstable basic state as constant background coefficients (default-off feature).

    ``U`` — per-layer uniform mean zonal flow (length ``nl``, m/s). ``G`` — the per-layer
    thermal-wind thickness gradient ``d(H̄_k)/dy`` (length ``nl``, dimensionless), derived from
    ``U`` so the layers are in geostrophic/thermal-wind balance. Built via
    :meth:`LayeredShallowWater.thermal_wind` (which fills ``G`` from ``U``). The fields stay
    perturbations; these inject the Doppler ``−U_k ∂_x`` and baroclinic ``−G_k v'`` terms.
    """

    U: np.ndarray
    G: np.ndarray


# --------------------------------------------------------------------------- #
# Solver
# --------------------------------------------------------------------------- #
class LayeredShallowWater:
    """N-layer rotating shallow-water solver (free-surface, Montgomery-coupled); see module docstring.

    Parameters
    ----------
    grid : Grid2D
        The uniform doubly-periodic C-grid (use :func:`engines.fluid.uniform_grid`).
    g : float
        Gravity (m/s²) — the surface/barotropic restoring (``g'₀ = g``).
    mean_depth : sequence of float
        Per-layer rest thicknesses ``H_k`` (m), top → bottom; ``len == nl``.
    reduced_gravities : sequence of float
        Reduced gravities ``g'_k = g·Δρ/ρ`` at the ``nl−1`` *internal* interfaces (m/s²),
        top → bottom; ``len == nl − 1``. A single value for the two-layer model.
    f0, beta : float
        Coriolis ``f = f₀ + β(y − y_ref)``. The background mode is f-plane (``β=0``); see scope.
    y_ref : float, optional
        Reference latitude for ``f``; defaults to the domain centre.
    bottom : ndarray, optional
        Bottom topography ``h_b`` at cell centers (m); ``None`` is a flat bottom.
    background : ThermalWindBackground, optional
        The basic-state coefficients; ``None`` (default) runs the plain nonlinear dynamics (and
        the bit-for-bit single-layer reduction). Build with :meth:`thermal_wind`.
    """

    def __init__(
        self,
        grid: Grid2D,
        g: float,
        mean_depth: Sequence[float],
        reduced_gravities: Sequence[float],
        f0: float = 0.0,
        beta: float = 0.0,
        y_ref: Optional[float] = None,
        bottom: Optional[np.ndarray] = None,
        background: Optional[ThermalWindBackground] = None,
    ) -> None:
        if g <= 0.0:
            raise ValueError(f"g must be positive, got {g}")
        H = np.asarray(mean_depth, dtype=float)
        if H.ndim != 1 or H.size < 1:
            raise ValueError("mean_depth must be a 1-D sequence of per-layer thicknesses")
        if np.any(H <= 0.0):
            raise ValueError(f"every layer thickness must be positive, got {H}")
        gp = np.asarray(reduced_gravities, dtype=float)
        if gp.size != H.size - 1:
            raise ValueError(
                f"reduced_gravities must have nl-1={H.size - 1} entries (internal interfaces), "
                f"got {gp.size}"
            )
        if H.size > 1 and np.any(gp <= 0.0):
            raise ValueError(f"internal reduced gravities must be positive (stable stratification), got {gp}")
        self.grid = grid
        self.g = float(g)
        self.H = H
        self.nl = int(H.size)
        # Interface reduced gravities g'_k with g'₀ = g (surface): length nl, indexed by layer.
        self._gp = np.concatenate(([self.g], gp))
        self.f0 = float(f0)
        self.beta = float(beta)
        self.y_ref = float(0.5 * grid.Ly if y_ref is None else y_ref)
        if bottom is None:
            self.h_b = np.zeros((grid.ny, grid.nx))
        else:
            self.h_b = np.asarray(bottom, dtype=float)
            if self.h_b.shape != (grid.ny, grid.nx):
                raise ValueError(f"bottom must have shape {(grid.ny, grid.nx)}, got {self.h_b.shape}")
        # Coriolis at corners (where ζ/q live), shape (ny, nx) — broadcasts over the layer axis.
        y_corner = grid.y_corners()
        self.f_corner = (self.f0 + self.beta * (y_corner - self.y_ref))[:, None] * np.ones((1, grid.nx))
        self.background = background
        if background is not None:
            if self.beta != 0.0:
                # The background is balanced only on the f-plane: once f varies in y the omitted
                # mean Coriolis −f(y)·U_k is unbalanced (a spurious drift, not the β-effect). The
                # β-capable PV-gradient treatment is the named within-rung extension — enforce the
                # documented constraint rather than silently drift (cf. the enforced CFL guard).
                raise ValueError(
                    "a thermal-wind background is consistent only on the f-plane (beta=0); "
                    "got beta != 0 (a β-capable basic state is the named within-rung extension)"
                )
            U = np.asarray(background.U, dtype=float)
            G = np.asarray(background.G, dtype=float)
            if U.shape != (self.nl,) or G.shape != (self.nl,):
                raise ValueError(f"background U and G must each have shape ({self.nl},)")
            self._U = U.reshape(self.nl, 1, 1)
            self._G = G.reshape(self.nl, 1, 1)

    # -- derived physical scales -------------------------------------------- #
    @property
    def external_gravity_wave_speed(self) -> float:
        """The fast **barotropic** gravity-wave speed ``√(g·H_tot)`` (m/s) — sets the explicit CFL."""
        return float(np.sqrt(self.g * np.sum(self.H)))

    def internal_deformation_radius(self) -> float:
        """The two-layer internal deformation radius ``√(g'·H_e)/|f₀|`` (m); two-layer only."""
        if self.nl != 2:
            raise ValueError("internal_deformation_radius is defined here for the two-layer model")
        He = self.H[0] * self.H[1] / (self.H[0] + self.H[1])
        return float(np.sqrt(self._gp[1] * He) / abs(self.f0))

    def stability(self) -> TwoLayerStability:
        """The matching :class:`TwoLayerStability` analytic operator (two-layer only) — the anchor."""
        if self.nl != 2:
            raise ValueError("the linear-stability anchor is the two-layer operator (nl == 2)")
        return TwoLayerStability(f0=self.f0, g=self.g, gp=float(self._gp[1]), H1=self.H[0], H2=self.H[1])

    def thermal_wind(self, mean_flow: Sequence[float]) -> ThermalWindBackground:
        """Build the background from per-layer mean zonal flow ``U_k`` via thermal-wind balance.

        Geostrophy gives each interface slope from the shear across it; the per-layer thickness
        gradient ``G_k = d(H̄_k)/dy`` follows. The surface slope is ``−f₀U₀/g``; the slope of the
        interface below layer ``k−1`` is ``f₀(U_{k−1}−U_k)/g'_k``; the bottom is flat. For the
        two-layer model these reduce to :meth:`TwoLayerStability.basic_state_gradients`, so the
        engine injects exactly the ``G_k`` the analytic growth rate assumes.
        """
        U = np.asarray(mean_flow, dtype=float)
        if U.shape != (self.nl,):
            raise ValueError(f"mean_flow must have shape ({self.nl},), got {U.shape}")
        # Interface-height slopes s[k] = d(η_k)/dy, k = 0..nl (η_0 = surface, η_nl = bottom = flat).
        s = np.zeros(self.nl + 1)
        s[0] = -self.f0 * U[0] / self.g
        for k in range(1, self.nl):
            s[k] = self.f0 * (U[k - 1] - U[k]) / self._gp[k]
        # s[nl] = 0 (flat bottom). Layer-thickness gradient G_k = s[k] − s[k+1].
        G = s[:-1] - s[1:]
        return ThermalWindBackground(U=U, G=G)

    # -- the Montgomery pressure stack (the ONLY inter-layer coupling) ------- #
    def _montgomery(self, h: np.ndarray) -> np.ndarray:
        """Per-layer Montgomery potential ``M_k`` (m²/s²), shape ``(nl, ny, nx)``.

        ``M₀ = g·η_top``; ``M_k = M_{k−1} + g'_k·η_{k−1}`` where ``η_j`` is the height of
        interface ``j`` above the reference (``h_b`` plus the thicknesses below it). The momentum
        pressure-gradient force is ``−∇M_k``. For ``nl=1`` this is ``g·(h_b + h)`` — exactly the
        single-layer Bernoulli pressure (bit-for-bit by IEEE commutativity).
        """
        # S[k] = Σ_{j≥k} h_j (column thickness from layer k down to the bottom).
        S = np.cumsum(h[::-1], axis=0)[::-1]
        M = np.empty_like(h)
        M[0] = self.g * (self.h_b + S[0])               # g·η_top
        for k in range(1, self.nl):
            M[k] = M[k - 1] + self._gp[k] * (self.h_b + S[k])
        return M

    # -- the right-hand side (method of lines) ------------------------------ #
    def _rhs(self, h: np.ndarray, u: np.ndarray, v: np.ndarray):
        """Tendencies ``(dh, du, dv)`` of the N-layer vector-invariant system, shape ``(nl, ny, nx)``.

        Per layer the stencils are *identical* to the single-layer engine (vectorized over the
        leading axis); layers couple only through :meth:`_montgomery`. With a background, the two
        constant-coefficient terms (Doppler ``−U_k ∂_x``, baroclinic ``−G_k v'``) are added — and
        *only* then, so ``background=None`` is byte-identical to the single-layer RHS.
        """
        dx, dy = self.grid.dx, self.grid.dy

        # Mass fluxes at velocity faces.
        h_u = 0.5 * (h + _xm(h))
        h_v = 0.5 * (h + _ym(h))
        U = h_u * u
        V = h_v * v

        # Continuity (flux form), per layer.
        dh = -((_xp(U) - U) / dx + (_yp(V) - V) / dy)

        # Potential vorticity at corners, per layer.
        zeta = (v - _xm(v)) / dx - (u - _ym(u)) / dy
        h_corner = 0.25 * (h + _xm(h) + _ym(h) + _xm(_ym(h)))
        q = (self.f_corner + zeta) / h_corner

        # Bernoulli: Montgomery pressure (the inter-layer coupling) + kinetic energy.
        u2c = 0.5 * (u ** 2 + _xp(u) ** 2)
        v2c = 0.5 * (v ** 2 + _yp(v) ** 2)
        B = self._montgomery(h) + 0.5 * (u2c + v2c)

        # Coriolis + vorticity flux, symmetric corner-PV averaging (same as the single-layer engine).
        cor_u = 0.25 * (q * (V + _xm(V)) + _yp(q) * (_yp(V) + _xm(_yp(V))))
        cor_v = -0.25 * (q * (U + _ym(U)) + _xp(q) * (_xp(U) + _ym(_xp(U))))

        du = cor_u - (B - _xm(B)) / dx
        dv = cor_v - (B - _ym(B)) / dy

        # Background coefficients (added only if set → bit-for-bit reduction when None).
        if self.background is not None:
            du = du - self._U * (_xp(u) - _xm(u)) / (2 * dx)        # Doppler advection by the mean
            dv = dv - self._U * (_xp(v) - _xm(v)) / (2 * dx)
            cv = 0.5 * (v + _yp(v))                                  # v at N–S faces → centers
            dh = dh - self._U * (_xp(h) - _xm(h)) / (2 * dx) - self._G * cv
        return dh, du, dv

    # -- time stepping (SSP-RK3) -------------------------------------------- #
    def step(self, state: LayeredState, dt: float) -> LayeredState:
        """Advance ``state`` by one SSP-RK3 step ``dt`` (s); returns a new state (no mutation).

        Raises on non-positive ``dt`` and on ``dt`` above the external-mode CFL limit. With
        ``background=None`` and ``nl=1`` this is **bit-for-bit** identical to
        :meth:`ShallowWater.step`.
        """
        if dt <= 0.0:
            raise ValueError("dt must be positive")
        cfl_limit = self.max_dt(state, safety=1.0)
        if dt > cfl_limit:
            raise ValueError(
                f"dt={dt:g} exceeds the external-mode CFL stability limit {cfl_limit:g} "
                "(= min(dx,dy)/(√(g·H_tot)+|U|+|u'|)); reduce dt (the recommended step is max_dt())"
            )
        h, u, v = state.h, state.u, state.v

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
        return LayeredState(h=hn, u=un, v=vn)

    def max_dt(self, state: LayeredState, safety: float = 0.3) -> float:
        """CFL-stable step ``safety · min(dx,dy)/(√(g·H_tot) + |U|_max + |u'|_max)`` (s).

        The fast **external** gravity wave ``√(g·H_tot)`` is the binding signal (the named
        cost); the mean-flow Doppler speed and the peak perturbation speed are added.
        """
        speed = self.external_gravity_wave_speed
        speed += float(np.max(np.sqrt(state.u ** 2 + state.v ** 2)))
        if self.background is not None:
            speed += float(np.max(np.abs(self._U)))
        return safety * min(self.grid.dx, self.grid.dy) / speed

    def solve(self, state: LayeredState, t_end: float, dt: Optional[float] = None,
              safety: float = 0.3) -> LayeredState:
        """March from ``0`` to ``t_end`` (s). ``dt=None`` uses the CFL step fixed at the initial state."""
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

    # -- diagnostics -------------------------------------------------------- #
    def layer_mass(self, state: LayeredState) -> np.ndarray:
        """Per-layer mass ``∫h_k dA`` (length ``nl``) — each conserved to machine precision (flux form)."""
        return np.sum(state.h, axis=(-2, -1)) * self.grid.cell_area

    def perturbation_energy(self, state: LayeredState) -> float:
        """Quadratic perturbation 'energy' ``Σ_k ∫(u_k² + v_k² + h'_k²) dA`` with ``h'_k = h_k − H_k``.

        The growth diagnostic for the linear baroclinic mode (it grows as ``e^{2σt}``); ``σ`` is
        read as half the slope of ``ln(energy)`` in the clean exponential window. (Not the physical
        total energy — the layered KE+APE invariant — which a background run does **not** conserve:
        baroclinic instability extracts APE from the mean, so the perturbation energy *grows*. That
        growth is the signal.)
        """
        hp = state.h - self.H.reshape(self.nl, 1, 1)
        dens = state.u ** 2 + state.v ** 2 + hp ** 2
        return float(np.sum(dens) * self.grid.cell_area)
