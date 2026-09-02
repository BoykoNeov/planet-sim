"""The **complete equilibrium diagram** of the ice-albedo EBM — every branch, both folds, and the
**small-ice-cap instability** (rung 0+, the Snowball's second bifurcation).

Phase 1 traced the Snowball hysteresis by a **parameter-continuation sweep** (:func:`planet.albedo.snowball_hysteresis`):
relax to equilibrium at each solar constant, warm-started from the last. A sweep can only ever sit on
a **stable** climate, so it sees the two jumps (the freeze, the re-melt) but never the *unstable*
branch that separates them, nor the exact location of the folds it falls off. And its every
equilibrium rides the Strang-split relaxation, whose fixed point carries an **O(Δt) splitting bias in
the profile shape** (the transport substep is backward-Euler) — the reason Phase 1 had to shrink the
present-day step to ``n_tau = 0.01`` to land the ice line near 73°.

This module solves the **inverse problem** instead, and it is *exact* (dt-free) for the discrete model:

    Don't ask "where does the ice line sit for this sun?" — ask "**for an ice line here, what sun
    holds it in equilibrium?**"

With the ice line ``x_s`` *prescribed*, the step-function albedo is a *known* field, the EBM is
**linear**, and one tridiagonal solve (the engine's own transport operator, :meth:`~planet.ebm.EnergyBalanceModel._transport_tridiag`)
gives the profile per unit solar constant, ``u(x)``. The ice-line condition ``T(x_s) = T_f`` with
``T = S₀·u − A/B`` then closes in one line:

    S₀(x_s) = (T_f + A/B) / u(x_s).

Sweeping ``x_s`` from the equator to the pole traces **every equilibrium with an ice line** — stable
*and* unstable — as one curve ``S₀(x_s)``; the ice-free planet (``x_s = 1``) and the Snowball
(``x_s = 0``) branches cap its two ends. This is the classic Budyko/Sellers/North *S-curve*
(North 1975, Fig. 3), here read off the *same* discrete operator the marcher uses.

The slope-stability theorem (Cahalan & North 1979) — stability without marching
--------------------------------------------------------------------------------
For this class of model an equilibrium is **stable if and only if ``dS₀/dx_s > 0``**: on a branch
where a *brighter* sun means a *smaller* cap (the intuitive direction) the climate is stable; where
the curve bends back (a smaller cap needs a *dimmer* sun) it is unstable. So the *shape of the curve
is the stability analysis*: its local **minimum** is the fold the dimming planet falls off (the
**Snowball catastrophe**, Phase 1's freeze), and its local **maximum** near the pole is a *second*
fold — the **small-ice-cap instability (SICI, North 1984)**: caps smaller than a critical angular
radius ``θ_c`` cannot be held by any sun, so a brightening planet's cap does not shrink smoothly to
nothing — it **vanishes in a jump** once it reaches ``θ_c``. Between the two folds lies the whole
**finite-cap window**: the only band of solar constants in which a planet with a polar ice cap can
exist at all. Present-day Earth's cap sits *inside and near the top of* that window.

Validation triad (plan §3)
--------------------------
* **Tight (analytic).** The finite-volume curve converges to North's **Legendre-mode** solution
  (:func:`legendre_equilibrium_curve` — the even-mode expansion ``T = Σ Tₙ Pₙ``, the transport
  diagonal ``−n(n+1)D`` in each mode, the albedo step projected by exact Gauss–Legendre quadrature)
  at the engine's spatial order. The reduction: at a *given* sun the curve's stable finite-cap crossing
  is the **dt → 0 limit** of Phase 1's relaxation (``present_day_climate``) — the exact equilibrium
  the O(Δt) fixed point converges onto as the step shrinks (which also *quantifies* the bias).
* **Tight (structural).** The slope-stability theorem is *checked, not assumed*: a nonlinear
  relaxation seeded from a curve equilibrium **returns** to it on a stable segment and **departs to
  another branch** on an unstable one. And the Phase-1 continuation sweep — an independent method —
  freezes, thaws and grows its first cap **within one sweep step of the folds/branch-ends** read here.
* **Conservation.** Every equilibrium on the curve balances net TOA exactly (a linear solve's residual).
* **Loose (the payoff, calibrated).** Where the folds sit in S₀ and the critical cap radius ``θ_c`` ride
  the cited ``A, B, D, α, T_f`` ([[ebm-radiation-source]]); the *existence* of the second fold, the
  narrowness of the finite-cap window, and present-day's proximity to the SICI are the banked findings.
  ``θ_c`` sits near **10° for weak transport** and **grows with ``D``** once ``D ≳ 0.4`` (a stronger
  transport washes out a small cap's local cooling — ~14° at ``D = 1``), the finite-cap window narrows
  with it, and for strong enough ``D`` (≈ 1.4 here) the two folds **merge** and the finite-cap branch
  disappears — North's result that an efficiently-mixed planet is *either* ice-free *or* a Snowball.

Units — SI, climlab-conventional: S₀ in W m⁻², T in °C, ``x = sin φ`` on [0, 1], latitudes in degrees,
cap radius ``θ = 90° − φ_ice``. Sources: North 1975 (*J. Atmos. Sci.* 32, the S-curve); Cahalan & North
1979 (*J. Atmos. Sci.* 36, 1178 — the slope-stability theorem); North 1984 (*J. Atmos. Sci.* 41, 3390 — the
small ice cap instability); Budyko 1969 / Sellers 1969.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field, replace

import numpy as np
from scipy.linalg import solve_banded

from .albedo import EBMParams
from .ebm import EnergyBalanceModel, insolation, legendre_P2, ice_line_latitude


# --------------------------------------------------------------------------- #
# Result types — plain arrays (the loose-coupling currency).
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Fold:
    """A turning point of ``S₀(x_s)`` — where a branch of equilibria ends and the climate jumps.

    ``kind="min"`` is a local *minimum* of ``S₀`` (a dimming planet on the stable branch above it falls
    off here — the Snowball catastrophe); ``kind="max"`` is a local *maximum* (a brightening planet's cap
    cannot shrink past it — the small-ice-cap instability). ``x_ice``/``latitude_deg`` locate the fold
    (parabola-refined between grid faces); ``S0`` is the sun at which it occurs (W m⁻²).
    """

    x_ice: float
    latitude_deg: float
    S0: float
    kind: str

    @property
    def cap_radius_deg(self) -> float:
        """Angular radius of the cap at the fold, ``90° − φ_ice`` — for a ``"max"`` fold, the critical ``θ_c``."""
        return 90.0 - self.latitude_deg


@dataclass(frozen=True)
class Equilibrium:
    """One equilibrium climate at a given sun: its ice line and its stability.

    ``kind`` ∈ ``{"ice-free", "finite-cap", "snowball"}``; ``stable`` from the slope-stability theorem
    (the two capping branches are always stable where they exist).
    """

    S0: float
    x_ice: float
    latitude_deg: float
    stable: bool
    kind: str


@dataclass(frozen=True)
class EquilibriumCurve:
    """The complete equilibrium set ``S₀(x_s)`` of the ice-albedo EBM — the S-curve, all branches.

    ``x_ice`` the prescribed ice-line positions (the grid *faces*, ``sin φ`` from 0 = a Snowball to
    1 = ice-free); ``S0`` the solar constant holding each in equilibrium (W m⁻²); ``stable`` the
    slope-stability verdict ``dS₀/dx_s > 0``; ``T`` the equilibrium profiles ``[n_faces, n_cells]``
    (°C); ``global_mean_T`` their area means. Plain arrays; :attr:`folds` and the branch-end thresholds
    are derived on demand.
    """

    x_ice: np.ndarray
    S0: np.ndarray
    stable: np.ndarray
    T: np.ndarray
    global_mean_T: np.ndarray
    x_cells: np.ndarray
    params: EBMParams = field(default_factory=EBMParams)

    # -- coordinates --------------------------------------------------------- #
    def latitude_deg(self) -> np.ndarray:
        """Ice-line latitudes ``φ = asin(x_s)`` (degrees, 0 = Snowball → 90 = ice-free)."""
        return np.degrees(np.arcsin(np.clip(self.x_ice, 0.0, 1.0)))

    # -- the branch ends (the two always-stable capping branches) ------------ #
    @property
    def ice_free_threshold_S0(self) -> float:
        """``S₀(x_s = 1)``: the ice-free planet exists (pole at/above ``T_f``) for suns **brighter** than this."""
        return float(self.S0[-1])

    @property
    def snowball_threshold_S0(self) -> float:
        """``S₀(x_s = 0)``: the Snowball exists (equator at/below ``T_f``) for suns **dimmer** than this — the re-melt."""
        return float(self.S0[0])

    # -- folds --------------------------------------------------------------- #
    @property
    def folds(self) -> tuple[Fold, ...]:
        """Every turning point of ``S₀(x_s)``, parabola-refined, ordered from the equator to the pole."""
        S0, xs = self.S0, self.x_ice
        out = []
        for i in range(1, S0.size - 1):
            y0, y1, y2 = S0[i - 1], S0[i], S0[i + 1]
            if (y1 < y0 and y1 <= y2) or (y1 > y0 and y1 >= y2):        # a strict local extremum
                denom = y0 - 2.0 * y1 + y2
                shift = 0.5 * (y0 - y2) / denom if denom != 0.0 else 0.0
                h = xs[i] - xs[i - 1]
                x_f = float(xs[i] + shift * h)
                S0_f = float(y1 - 0.25 * (y0 - y2) * shift)
                kind = "min" if y1 < y0 else "max"
                out.append(Fold(x_f, math.degrees(math.asin(min(1.0, max(0.0, x_f)))), S0_f, kind))
        return tuple(out)

    @property
    def snowball_fold(self) -> Fold | None:
        """The global **minimum** of the curve — the fold a dimming planet freezes over at (the Snowball catastrophe)."""
        mins = [f for f in self.folds if f.kind == "min"]
        return min(mins, key=lambda f: f.S0) if mins else None

    @property
    def small_ice_cap_fold(self) -> Fold | None:
        """The polar **maximum** of the curve — the small-ice-cap instability; its ``cap_radius_deg`` is ``θ_c``."""
        maxes = [f for f in self.folds if f.kind == "max" and f.x_ice > 0.5]
        return max(maxes, key=lambda f: f.x_ice) if maxes else None

    @property
    def finite_cap_window(self) -> tuple[float, float] | None:
        """``(S₀_min, S₀_max)`` — the band of suns admitting a *stable* polar cap (Snowball fold → SICI fold); ``None`` if no such branch."""
        lo, hi = self.snowball_fold, self.small_ice_cap_fold
        if lo is None or hi is None or hi.S0 <= lo.S0:
            return None
        return (lo.S0, hi.S0)

    # -- all equilibria at one sun ------------------------------------------- #
    def equilibria_at(self, S0: float) -> list[Equilibrium]:
        """Every equilibrium climate at solar constant ``S0``: the crossings of the curve + the capping branches."""
        out = []
        if S0 >= self.ice_free_threshold_S0:
            out.append(Equilibrium(S0, 1.0, 90.0, True, "ice-free"))
        if S0 <= self.snowball_threshold_S0:
            out.append(Equilibrium(S0, 0.0, 0.0, True, "snowball"))
        d = self.S0 - S0
        for i in range(d.size - 1):
            if d[i] == 0.0 or d[i] * d[i + 1] < 0.0:
                t = d[i] / (d[i] - d[i + 1]) if d[i] != d[i + 1] else 0.0
                x_c = float(self.x_ice[i] + t * (self.x_ice[i + 1] - self.x_ice[i]))
                if x_c <= 0.0 or x_c >= 1.0:
                    continue
                stable = bool(self.S0[i + 1] > self.S0[i])                 # slope-stability theorem
                out.append(Equilibrium(S0, x_c, math.degrees(math.asin(x_c)), stable, "finite-cap"))
        return out

    def stable_finite_cap_at(self, S0: float) -> Equilibrium | None:
        """The stable finite-cap equilibrium at ``S0`` (Earth's branch), or ``None`` outside the window."""
        cands = [e for e in self.equilibria_at(S0) if e.kind == "finite-cap" and e.stable]
        return max(cands, key=lambda e: e.x_ice) if cands else None


# --------------------------------------------------------------------------- #
# The inverse solve — one linear system per prescribed ice line.
# --------------------------------------------------------------------------- #
def _ice_free_albedo(x: np.ndarray, p: EBMParams) -> np.ndarray:
    return p.a0 + p.a2 * legendre_P2(x)


def equilibrium_curve(params: EBMParams | None = None) -> EquilibriumCurve:
    """Trace ``S₀(x_s)`` over every grid face — the complete, dt-free equilibrium diagram.

    For each face ``k`` the cells poleward of it carry the ice albedo ``a_ice`` and those equatorward
    the ice-free ``a₀ + a₂P₂``; the linear steady EBM ``(L_T − B)u = −s(x)(1 − α)`` (``s`` the insolation
    per unit S₀) is solved by the engine-pinned tridiagonal, and the ice-line condition ``T(x_s) = T_f``
    gives ``S₀ = (T_f + A/B)/u(x_s)`` with ``u`` read at the face as the mean of its two neighbouring
    cells (2nd order; the end faces take the adjacent cell, exact to O(Δx²) by the no-flux symmetry).
    Honors ``params.face`` (harmonic / exact) and ``params.n_cells``. The stability of every point is
    the slope-stability verdict ``dS₀/dx_s > 0`` (central differences; one-sided at the ends).
    """
    p = params if params is not None else EBMParams()
    model = EnergyBalanceModel(A=p.A, B=p.B, D=p.D, T_freeze=p.T_freeze,
                               water_depth=p.water_depth, n_cells=p.n_cells, face=p.face)
    x = model.x
    faces = np.asarray(model.grid.edges, dtype=float)
    n = x.size
    sub, diag, sup = model._transport_tridiag()
    ab = np.zeros((3, n))
    ab[0, 1:] = sup[:-1]
    ab[1, :] = diag - p.B
    ab[2, :-1] = sub[1:]
    s_unit = insolation(x, 1.0, p.s2)                       # insolation per unit S₀
    ice_free = _ice_free_albedo(x, p)
    Tf_shift = p.T_freeze + p.A / p.B

    S0 = np.empty(n + 1)
    T = np.empty((n + 1, n))
    Tbar = np.empty(n + 1)
    for k in range(n + 1):
        alb = ice_free.copy()
        alb[k:] = p.ai                                       # ice poleward of face k
        u = solve_banded((1, 1), ab, -s_unit * (1.0 - alb))
        u_face = u[0] if k == 0 else (u[-1] if k == n else 0.5 * (u[k - 1] + u[k]))
        S0[k] = Tf_shift / u_face
        T[k] = S0[k] * u - p.A / p.B
        Tbar[k] = float(np.mean(T[k]))                       # equal-area cells ⟹ mean = area mean
    slope = np.gradient(S0, faces)
    stable = slope > 0.0
    return EquilibriumCurve(faces, S0, stable, T, Tbar, x, p)


# --------------------------------------------------------------------------- #
# North's Legendre-mode solution — the independent analytic anchor.
# --------------------------------------------------------------------------- #
def legendre_equilibrium_curve(x_ice, params: EBMParams | None = None, n_modes: int = 200) -> np.ndarray:
    """``S₀(x_s)`` from the **even-Legendre-mode** expansion (North 1975, generalized to N modes).

    On the hemisphere with no-flux ends the eigenfunctions of the transport are the even Legendre
    polynomials, ``d/dx[(1−x²)dPₙ/dx] = −n(n+1)Pₙ``, so each mode of ``T = Σₙ Tₙ Pₙ`` decouples:
    ``Tₙ·[n(n+1)D + B] = S₀·Hₙ(x_s) − A·δₙ₀`` with ``Hₙ(x_s) = (2n+1)∫₀¹ s(x)(1 − α(x; x_s))Pₙ(x)dx``
    (``s`` the insolation per unit S₀). The ice-line condition ``T(x_s) = T_f`` then gives

        S₀(x_s) = (T_f + A/B) / Σₙ Hₙ(x_s)Pₙ(x_s)/[n(n+1)D + B].

    The albedo step makes ``Hₙ`` an integral of a *piecewise* polynomial, evaluated **exactly** by
    Gauss–Legendre quadrature on each side of ``x_s`` (the integrand is a polynomial of degree
    ``n + 4`` on each piece). North's original is the ``n_modes = 2`` truncation (``P₀, P₂``); with
    ``n_modes`` large this is the continuous model's exact curve, the finite-volume
    :func:`equilibrium_curve` must converge to it at the engine's spatial order. Returns ``S₀`` at each
    ``x_ice`` (W m⁻²). Uniform-``D`` only (a callable ``D(x)`` is not mode-diagonal).
    """
    p = params if params is not None else EBMParams()
    if callable(p.D):
        raise ValueError("the Legendre-mode solution needs a uniform scalar D")
    x_ice = np.atleast_1d(np.asarray(x_ice, dtype=float))
    n_max = 2 * (int(n_modes) - 1)                          # even modes 0, 2, …, n_max
    modes = np.arange(0, n_max + 1, 2)
    # Gauss–Legendre nodes for a polynomial of degree ≤ n_max + 4 on each piece.
    n_gl = (n_max + 4) // 2 + 2
    g_nodes, g_weights = np.polynomial.legendre.leggauss(n_gl)
    denom = modes * (modes + 1) * p.D + p.B
    Tf_shift = p.T_freeze + p.A / p.B

    def H_all(xs: float) -> np.ndarray:
        """``Hₙ(x_s)`` for every even mode — the two-piece exact quadrature."""
        H = np.zeros(modes.size)
        for (lo, hi, alb_ice) in ((0.0, xs, False), (xs, 1.0, True)):
            if hi <= lo:
                continue
            xq = 0.5 * (hi - lo) * g_nodes + 0.5 * (hi + lo)
            wq = 0.5 * (hi - lo) * g_weights
            alb = np.full(xq.shape, p.ai) if alb_ice else _ice_free_albedo(xq, p)
            f = insolation(xq, 1.0, p.s2) * (1.0 - alb)
            P = np.polynomial.legendre.legvander(xq, n_max)[:, modes]    # [n_q, n_modes]
            H += (2 * modes + 1) * ((wq * f) @ P)
        return H

    out = np.empty(x_ice.size)
    for i, xs in enumerate(x_ice):
        P_at = np.polynomial.legendre.legvander(np.array([xs]), n_max)[0, modes]
        out[i] = Tf_shift / float(np.sum(H_all(float(xs)) * P_at / denom))
    return out


# --------------------------------------------------------------------------- #
# Stability by marching — the theorem checked against the nonlinear relaxation.
# --------------------------------------------------------------------------- #
def relax_from_curve(curve: EquilibriumCurve, index: int, perturb_K: float,
                     n_tau: float = 0.02, tol: float = 1e-8, max_iter: int = 60000):
    """Seed the nonlinear (ice-feedback) relaxation from the curve's ``index``-th equilibrium, nudged by ``perturb_K``.

    Returns the relaxed :class:`~planet.ebm.ClimateState` at that equilibrium's own ``S₀``. On a
    *stable* segment the relaxation stays (the ice line returns to within the O(Δt) fixed-point
    bias); on an *unstable* segment it departs to a neighbouring stable branch — how the
    slope-stability theorem is checked, not assumed. A gentle ``n_tau`` keeps the splitting bias small.
    """
    p = replace(curve.params, S0=float(curve.S0[index]))
    model = p.model()
    T0 = curve.T[index] + float(perturb_K)
    return model.equilibrate(p.absorbed_fn(), T0, n_tau=n_tau, tol=tol, max_iter=max_iter)


# --------------------------------------------------------------------------- #
# The critical cap radius as a function of transport — the D-sweep.
# --------------------------------------------------------------------------- #
def critical_cap_sweep(D_values, params: EBMParams | None = None, n_cells: int = 720):
    """``θ_c(D)`` — the small-ice-cap critical radius (degrees) and the finite-cap window, per transport ``D``.

    Returns ``(theta_c, window_lo, window_hi)`` arrays over ``D_values``; ``NaN`` where the finite-cap
    branch does not exist (the folds have merged — too efficient a transport leaves only the ice-free
    and Snowball climates). The direction banked: ``θ_c`` ≈ 10° at weak transport, growing with ``D``
    past ``D ≈ 0.4`` while the window narrows to nothing. Runs on a **fine grid**
    (``n_cells``, default 720): the uniform-``x`` grid is coarsest in *latitude* near the pole (~2° per
    cell at 80° for 180 cells), exactly where this fold sits, and the parabola-refined fold converges
    only to ~1° there; 720 cells brings it to ~0.1°.
    """
    p = params if params is not None else EBMParams()
    theta, lo, hi = [], [], []
    for D in np.asarray(D_values, dtype=float):
        c = equilibrium_curve(replace(p, D=float(D), n_cells=int(n_cells)))
        w = c.finite_cap_window
        f = c.small_ice_cap_fold
        theta.append(f.cap_radius_deg if (f is not None and w is not None) else np.nan)
        lo.append(w[0] if w is not None else np.nan)
        hi.append(w[1] if w is not None else np.nan)
    return np.array(theta), np.array(lo), np.array(hi)


def relaxation_bias_sweep(n_taus, params: EBMParams | None = None, **present_kw):
    """Ice line of Phase 1's relaxed present-day climate vs its step ``n_tau`` — the O(Δt) bias, quantified.

    Returns the relaxed ice-line latitude (degrees) for each ``n_tau``; read against the curve's
    exact :meth:`EquilibriumCurve.stable_finite_cap_at` at the same sun, the sequence converges onto it
    as the step shrinks — the reduction that ties the marcher to the exact diagram.
    """
    from .albedo import present_day_climate
    p = params if params is not None else EBMParams()
    return np.array([present_day_climate(p, n_tau=float(nt), **present_kw).ice_line_lat
                     for nt in np.asarray(n_taus, dtype=float)])
