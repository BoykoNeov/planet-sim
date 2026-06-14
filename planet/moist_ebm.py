"""The MSE-diffusing moist energy-balance model — emergent polar amplification (Planet rung 2.5).

Rung 2 (:mod:`planet.moist`) added moisture as a **pure diagnostic** that never touched the
temperature equation, so the validated Phase-1 climate stayed bit-for-bit fixed. Rung 2.5 is the
named step *up* from there: a moist EBM in which **``T`` itself responds** to moisture transport, so
the headline is an **emergent climate response** the diagnostic could not produce — **polar-amplified
warming** (Hwang & Frierson 2010; Flannery 1984; Siler, Roe & Armour 2018). It is a **separate model
that runs alongside rung-0**, not a replacement: the dry EBM stays the default everywhere; this is the
opt-in moist sibling, exactly as circulation-informed precip and the energy-constrained rate were opt-ins.

The mechanism — moisture makes the diffusivity grow where it is warm
-------------------------------------------------------------------
A moist atmosphere diffuses **moist static energy** ``m = c_p·T + L·q`` (sensible + latent) down its
own gradient, with ``q = RH·q_sat(T)`` at fixed relative humidity. In temperature-equivalent units
(dividing ``m`` by ``c_p``) the transported quantity is ``T + (L/c_p)·q``, so by the chain rule the
moist down-gradient flux is, term-by-term, a **temperature** diffusion with a *moisture-amplified
effective diffusivity*

    D_eff(T) = D_s · (1 + (L/c_p)·RH·dq_sat/dT)  ≡  D_s · (1 + β(T)),

where ``D_s`` is the sensible (dry) coefficient and ``β(T) = (L/c_p)·RH·dq_sat/dT``
(:func:`moisture_amplification`) is the latent amplification. Because ``q_sat`` rises steeply with
temperature (Clausius–Clapeyron, ~7 %/K), ``β`` is **large in the warm tropics and ≈ 0 at the cold
pole** — for Earth defaults ``D_eff`` runs ~1.3 at the equator down to ~0.35 at the pole. As the
climate warms, the tropical ``β`` grows fastest, the tropics export *more* heat poleward, and **the
poles warm more than the tropics**: polar amplification, emergent from moisture alone — **no
ice-albedo feedback, no change in ``D_s``**. The dry EBM (constant ``D``) cannot do this: under a
uniform forcing it warms **exactly uniformly** (the clean null, :func:`polar_amplification`).

The headline is redistribution, not extra heat (the framing, the tight anchor)
------------------------------------------------------------------------------
With constant albedo the global-mean response to a uniform OLR forcing ``ΔA`` is **pinned**:
``⟨δT⟩ = ΔA/B`` to machine precision, for *any* transport ``D`` (the diffusion conserves ``∫T dx``, so
transport cannot change the global mean). **Moisture redistributes that fixed ``⟨δT⟩`` poleward** — it
does not add warming at the pole; it moves warming from the equator to the pole around a pinned mean.
That global-mean identity is asserted **tight** (conservation), and the polar amplification is reported
as the *shape* of the redistribution around it.

Reusing the spine — D_eff is a state-dependent coefficient, frozen each step like α(T)
--------------------------------------------------------------------------------------
This is the diffusion spine's **fifth** reuse. The moist climate is **one nonlinear relaxation**, not
nested solves: each Strang substep freezes ``D_eff`` at the current ``T`` and rebuilds the conservative
transport operator ``∂/∂x[(1−x²)·D_eff(T)·∂T/∂x]`` (the coefficient **inside** the divergence, so the
engine's no-flux invariant still gives machine-exact energy conservation) — the **identical idiom** the
ice-albedo ``α(T)`` already uses in :meth:`planet.ebm.EnergyBalanceModel.equilibrate` (a state-dependent
coefficient re-frozen each substep). The radiation half-steps are the same analytic Strang sandwich.
The relaxation is **self-contained** — it does **not** modify :mod:`planet.ebm` (the validated rung-0
hot path is untouched), the ~one radiation helper duplicated being the correct price of that discipline.
``face="harmonic"`` is pinned: the reduction below is bit-for-bit because the per-step cells
``(D_s/C)(1−x²)`` match the dry model's once-built harmonic cells; ``face="exact"`` would pre-distort
them and silently break it.

**The headline magnitude, however, comes from the dt-FREE solve** (:func:`moist_steady_direct`), not the
relaxation. The same frozen-``D_eff`` coefficient, dropped into a **direct linear steady solve** ``(L_T −
B·I)T = A − S(1−α)`` and Picard-iterated, converges with **no time-stepping** — so with **no
operator-splitting error**. This is the spine's reuse of :meth:`planet.ebm.EnergyBalanceModel.steady_linear`
(the nonlinear generalisation of it, exactly as rung-4's :mod:`planet.radiative_ebm` is its Newton
generalisation). It matters because the Strang relaxation's steady state carries an **O(Δt) shape bias**
(the backward-Euler transport split against the exact radiation half-step): the global mean stays exact,
but the *shape* — and hence the polar-amplification *ratio* — drifts with the step, the default
``n_tau = 0.5`` suppressing the endpoint factor to ≈ 1.5 against the dt-free ≈ 2.05. The relaxation is kept
for the bit-for-bit rung-0 reduction and for animations (it has the transient); the magnitude is the
direct solve's.

The recalibration — re-derive D_s so we do not double-count latent transport (the wall)
---------------------------------------------------------------------------------------
Rung-0's ``D = 0.555`` is an **effective** diffusivity calibrated to the observed present-day gradient;
it *already* contains the real atmosphere's latent transport lumped into one number. Diffusing latent
heat **explicitly** on top of that would double-count it, so the moist EBM re-derives a **smaller
sensible** ``D_s`` (:func:`recalibrate_sensible_D`) such that the moist present-day climate reproduces
the dry present-day **equator-to-pole contrast** (≈ 0.28 vs 0.555 for Earth + RH 0.8). The global mean
``⟨T⟩`` is *automatically* equal in both (energy balance fixes it from ``A, B, ᾱ`` independent of ``D``),
so the contrast is the natural single recalibration scalar.

* **The trade (named, not a win).** A *single scalar* ``D_s`` cannot reproduce the whole dry ``T(x)``:
  the moist temperature profile matches the dry **mean and contrast** but differs in **higher-moment
  shape** (the matched-contrast moist profile is flatter in the interior with the curvature concentrated
  toward the edges). This is the honest cost of the recalibration — the same "trade, not a win" every
  rung banks.
* **The recalibration target is a modeling choice (named).** Matching endpoint-contrast leans on the
  two most polar cells (where the harmonic-face pole floor lives — :mod:`planet.ebm`); but the polar
  amplification **factor is invariant** to the choice (matching the P₂ amplitude ``T₂`` instead moves it
  < 2 %), because PA is set by the *shape* of ``β(T)``, not by the overall scaling.

Validation triad (plan §3) — re-classed for honesty
----------------------------------------------------
* **Tight (analytical / structural / conservation).** (a) ``β(T)`` is the **exact** analytic derivative
  of the Clausius–Clapeyron ``q_sat`` (the Whittaker-partition precedent — an exact function, not a
  fit). (b) The transport operator is the EBM's own ``∂/∂x[(1−x²)∂/∂x]`` with ``D_eff`` **inside** — it
  reproduces the **P₂ Legendre eigenvalue** ``−6`` and conserves ``∫T dx`` to machine precision. (c) The
  **global-mean warming is pinned**: ``⟨δT⟩ = ΔA/B`` to machine precision for both models. (d) The
  **frozen-D_eff attribution null**: freezing ``D_eff`` at its present profile and warming gives
  **exactly uniform** ``δT`` (PA = 1), proving the emergent PA is the ``dD_eff/dT`` feedback, not the
  recalibrated ``D``-shape.
* **Real-but-loose physics (the unlock).** The **polar amplification** itself: the poles warm more than
  the tropics, emergent from moisture transport alone. The **direction** (PA > 1, robustly across RH) is
  banked; the **magnitude** is loose — and the *metric matters*, so it is reported two ways from the
  **dt-free** :func:`moist_steady_direct` climate: the **single-endpoint** ratio ``δT(pole)/δT(equator)``
  ≈ **2.05** (Earth, RH 0.8 — the most generous, and the polar cell sits on the harmonic-face bias), and
  the **area-band** ratio ``mean(δT|φ≥60°)/mean(δT|φ≤30°)`` ≈ **1.80**; "~2.0" names the endpoint ratio
  specifically. *This is the model's **converged** number, not a **validated** one*: the magnitude is set
  by the formulation (the inside-divergence MSE diffusion, the contrast recalibration), so whether
  moisture transport alone "should" give ≈1.5 or ≈2.0 is a property of the model — and the observed ~2–3×
  additionally needs the ice-albedo + lapse-rate feedbacks held out of scope here.
  **(Earlier banked as ≈1.5 / ≈1.4 — that was the Strang relaxation's O(Δt) shape bias at the default
  ``n_tau = 0.5``, which suppresses the amplification; the dt-free direct solve removes it. The
  *direction* was never in doubt — and the relaxation's global mean was exact throughout.)**
* **Plumbing (by construction).** RH = 0 **and** ``D_s = D_TRANSPORT`` reduce **both** paths to the dry
  rung-0 solve **bit-for-bit** (``β ≡ 0`` ⟹ the per-step operator is the dry one): the relaxation
  (:func:`moist_equilibrium`) to :meth:`planet.ebm.EnergyBalanceModel.equilibrate`, and the dt-free direct
  solve (:func:`moist_steady_direct`) to :meth:`planet.ebm.EnergyBalanceModel.steady_linear`.

Named scope edges
-----------------
* **Only ``D`` is re-opened; ``A`` and ``B`` are held fixed.** The plan's "rung 2.5 re-opens the
  ``(A, B, D)`` calibration" is refined here: ``D_s`` is re-derived (the transport double-count), ``A``
  is the forcing knob (uniform ``ΔA`` = the CO₂ proxy), and **``B`` stays fixed** — re-deriving the
  water-vapour content of the linear OLR ``A + B·T`` is a **local-radiation** problem (a real
  vertically-resolved radiative transfer), the named **rung-4** wall, not opened here.
* **Fixed RH; constant albedo for the clean experiment.** RH is fixed (the diffusive-moist-EBM closure);
  the polar-amplification experiment holds the albedo **constant** (ice OFF) to isolate the *moisture*
  mechanism (the Hwang–Frierson clean experiment) — the ice-albedo and lapse-rate feedbacks that would
  push the factor to the observed ~2–3× are deliberately out of scope.
* **Forcing is uniform ``ΔA``, not ``ΔS₀``.** ``S₀`` is equator-weighted by ``(1 + s₂P₂)``, which would
  impose a tropical structure that *fights* polar amplification; a uniform OLR offset is the clean knob.
* **Separate from the rung-2 ``P − E`` budget.** This is a ``T``-response model; it does not touch
  :mod:`planet.moist`'s moisture-convergence diagnostic (those tests stay green).

Units — SI, climlab-conventional (W m⁻², °C, x = sin φ dimensionless)
---------------------------------------------------------------------
``D_s``/``D_eff`` in W m⁻² K⁻¹; ``A``, ``ΔA`` in W m⁻²; ``B`` in W m⁻² K⁻¹; ``q``, ``β`` dimensionless;
``T`` in °C (converted to K inside the C–C functions); ``x = sin φ`` on [0, 1]; latitudes in degrees.
See [[moist-ebm-source]], [[ebm-radiation-source]], [[precip-parameterization-source]], [[planet-plan]] §10.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np
from scipy.linalg import solve_banded
from scipy.optimize import brentq

from engines.diffusion import Diffusion1D, Neumann, uniform_grid

from . import moist
from .albedo import EBMParams
from .ebm import (
    ClimateState, CW_WATER, RHO_WATER, ice_line_latitude, insolation, legendre_P2,
)
from .transport import CP_AIR

# --------------------------------------------------------------------------- #
# Pinned constant. L/c_p converts a specific humidity to the temperature-equivalent latent energy that
# moist static energy m = c_p·T + L·q carries (≈ 2490 K). L_VAPOR and CP_AIR are the *already-pinned*
# rung-2 / rung-1 constants ([[moist-ebm-source]], [[shallow-water-source]]) — reused, not re-declared.
# --------------------------------------------------------------------------- #
L_OVER_CP = moist.L_VAPOR / CP_AIR        # K — latent energy per unit specific humidity, in temp units


# --------------------------------------------------------------------------- #
# The moisture amplification of the diffusivity — the exact analytic C–C derivative (the TIGHT leg).
# --------------------------------------------------------------------------- #
def dqsat_dT(T_celsius):
    """Exact derivative ``dq_sat/dT`` (per K) of the Clausius–Clapeyron saturation humidity.

    The analytic derivative of :func:`planet.moist.saturation_specific_humidity`
    ``q_sat = ε·e_sat/(p − (1−ε)e_sat)`` via the chain rule with ``de_sat/dT = e_sat·L_v/(R_v·T_K²)``:

        dq_sat/dT = [ε·p / (p − (1−ε)e_sat)²] · e_sat · L_v/(R_v·T_K²).

    An **exact** function (the Whittaker-partition precedent), not a finite-difference fit — pinned in
    the triad against both a numerical derivative and the ~7 %/K C–C log-slope. ``T`` in °C (→ K inside).
    """
    T_K = np.asarray(T_celsius, dtype=float) + moist.T0_KELVIN
    e_sat = moist.saturation_vapor_pressure(T_celsius)
    denom = moist.P_SURFACE - (1.0 - moist.EPSILON) * e_sat
    de_dT = e_sat * moist.L_VAPOR / (moist.R_VAPOR * T_K ** 2)
    dq_de = moist.EPSILON * moist.P_SURFACE / denom ** 2
    return dq_de * de_dT


def moisture_amplification(T_celsius, RH: float = moist.RH_DEFAULT):
    """The latent diffusivity amplification ``β(T) = (L/c_p)·RH·dq_sat/dT`` (dimensionless).

    The factor by which fixed-RH moisture amplifies the down-gradient temperature diffusivity (the
    moist static energy ``T + (L/c_p)q`` versus the dry ``T``). Large where it is warm (``dq_sat/dT``
    is steep) and ≈ 0 where it is cold — for Earth defaults ``β`` runs ~3 at a 30 °C equator down to
    ~0.1 at a −20 °C pole. The whole moist-EBM physics rides on ``β`` varying across the planet.
    """
    return L_OVER_CP * float(RH) * dqsat_dT(T_celsius)


def effective_diffusivity(T_celsius, D_s: float, RH: float = moist.RH_DEFAULT):
    """The moisture-amplified effective transport coefficient ``D_eff = D_s·(1 + β(T))`` (W m⁻² K⁻¹).

    The dry sensible diffusivity ``D_s`` amplified by the latent factor :func:`moisture_amplification`.
    This is the coefficient the moist relaxation places **inside** the spherical transport operator
    ``∂/∂x[(1−x²)·D_eff·∂T/∂x]`` (so the divergence stays conservative). Reduces to the dry ``D_s`` as
    ``RH → 0``.
    """
    return float(D_s) * (1.0 + moisture_amplification(T_celsius, RH))


# --------------------------------------------------------------------------- #
# The constant-albedo forcing — the clean experiment that isolates the moisture mechanism.
# --------------------------------------------------------------------------- #
def constant_albedo_absorbed(params: EBMParams, A: Optional[float] = None
                             ) -> Callable[[np.ndarray, np.ndarray], np.ndarray]:
    """The ice-OFF absorbed shortwave ``S(x)·(1 − [a₀ + a₂P₂])`` — albedo independent of ``T``.

    The Hwang–Frierson clean experiment: holding the albedo at its ice-free poleward-brightening value
    (no ``T < Tf`` ice jump) so the emergent polar amplification is attributable to **moisture alone**,
    not to ice retreat. ``A`` is accepted for call-site symmetry (the OLR offset lives in the relaxation,
    not the shortwave) and ignored here.
    """
    def absorbed(x, T):
        ice_free = params.a0 + params.a2 * legendre_P2(x)
        return insolation(x, params.S0, params.s2) * (1.0 - ice_free)
    return absorbed


# --------------------------------------------------------------------------- #
# The moist EBM — one nonlinear relaxation, D_eff(T) frozen each Strang step (like α(T)).
# --------------------------------------------------------------------------- #
def moist_equilibrium(params: EBMParams, D_s: float, absorbed_fn,
                      RH: float = moist.RH_DEFAULT, T_init=None, A: Optional[float] = None,
                      n_tau: float = 0.5, tol: float = 1e-9, max_iter: int = 20000) -> ClimateState:
    """Relax the MSE-diffusing moist EBM to equilibrium; return its :class:`~planet.ebm.ClimateState`.

    One nonlinear Strang-split relaxation in which the transport coefficient is the moisture-amplified
    ``D_eff(T) = D_s·(1 + β(T))`` **re-frozen each substep** (the ice-albedo ``α(T)`` idiom): each step
    is a radiation half-step, a full implicit transport step through the conservative operator
    ``∂/∂x[(1−x²)·D_eff(T)·∂T/∂x]`` (rebuilt from the current ``T`` — ``D_eff`` inside the divergence),
    and a second radiation half-step. Self-contained — it does **not** touch :mod:`planet.ebm`.

    **Carries an O(Δt) shape bias — not the headline path.** The exact radiation half-step is split
    against a *backward-Euler* transport step, so the converged steady state's **shape** carries a
    first-order operator-splitting error (the global mean stays exact): at the default ``n_tau = 0.5`` the
    bias *suppresses* the emergent polar amplification (endpoint ≈ 1.5 vs the dt-free ≈ 2.05), shrinking as
    ``n_tau → 0``. Use this path for the bit-for-bit rung-0 reduction and for animations (it produces the
    intermediate transient); use the **dt-free** :func:`moist_steady_direct` for the headline magnitude.

    Parameters
    ----------
    params : EBMParams
        The parameter bundle (``B``, ``water_depth``, ``n_cells``, albedo, ``S0`` …). ``face`` is
        **pinned to "harmonic"** internally (the bit-for-bit reduction relies on it).
    D_s : float
        The dry **sensible** transport coefficient (W m⁻² K⁻¹); the moist ``D_eff`` amplifies it. Use
        :func:`recalibrate_sensible_D` to obtain the value matching the present-day climate.
    absorbed_fn : callable ``(x, T) -> ndarray``
        The absorbed shortwave (W m⁻²) — :func:`constant_albedo_absorbed` for the clean experiment.
    RH : float
        Fixed relative humidity (default 0.8). ``RH = 0`` ⟹ ``β ≡ 0`` ⟹ the dry rung-0 relaxation.
    T_init : array | float | None
        Initial field (°C); ``None`` → an Earth-like warm-equator/cold-pole start.
    A : float | None
        OLR offset (W m⁻²); ``None`` → ``params.A``. The uniform ``ΔA`` forcing changes this.
    n_tau, tol, max_iter
        Strang step (in ``τ_rad``), convergence tolerance (°C), iteration cap — the rung-0 defaults.
    """
    grid = uniform_grid(1.0, params.n_cells)
    x = grid.centers
    C = RHO_WATER * CW_WATER * params.water_depth
    B = params.B
    A_olr = params.A if A is None else float(A)
    dt = n_tau * (C / B)
    decay = math.exp(-0.5 * dt * B / C)
    if T_init is None:
        T_init = 30.0 + (-30.0 - 30.0) * x                         # warm equator → cold pole (Earth-like)

    def rad_half(T):
        T_eq = (absorbed_fn(x, T) - A_olr) / B
        return T_eq + (T - T_eq) * decay

    T = np.array(np.broadcast_to(np.asarray(T_init, dtype=float), x.shape), dtype=float)
    converged, it = False, 0
    for it in range(1, max_iter + 1):
        T_old = T
        T = rad_half(T)
        # D_eff(T) frozen at this substep's T, placed INSIDE the divergence: cells (D_eff/C)(1−x²).
        Dcells = (effective_diffusivity(T, D_s, RH) / C) * (1.0 - x ** 2)
        T = Diffusion1D(grid, Dcells, Neumann(0.0), Neumann(0.0)).step(T, dt)
        T = rad_half(T)
        if np.max(np.abs(T - T_old)) < tol:
            converged = True
            break
    Tbar = float(np.mean(T))                                       # L = 1 ⟹ area mean (∫₀¹ T dx)
    return ClimateState(
        x=x, T=T, global_mean_T=Tbar,
        ice_line_lat=ice_line_latitude(x, T, params.T_freeze),    # the Tf-isotherm diagnostic (albedo fixed)
        net_toa=float(np.mean(absorbed_fn(x, T)) - A_olr - B * Tbar),
        converged=converged, iterations=it,
    )


# --------------------------------------------------------------------------- #
# The dt-FREE moist steady solve — Picard on the frozen-D_eff linear solve (the headline path).
#
# The moist steady state is nonlinear *only* through D_eff(T), so freezing D_eff and solving the linear
# (L_T[D_eff] − B·I)T = A − S(1−α) in one tridiagonal solve, then re-freezing and repeating, converges
# (Picard) onto the true nonlinear steady state with **no time-stepping** — hence **no operator-splitting
# dt bias**. This is the nonlinear generalisation of :meth:`planet.ebm.EnergyBalanceModel.steady_linear`
# (exactly as rung-4's :mod:`planet.radiative_ebm` is its Newton generalisation), and it is the path the
# headline :func:`polar_amplification` uses. The Strang relaxation :func:`moist_equilibrium` is kept as the
# general/animatable path, but its steady state carries an **O(Δt) shape bias** (the backward-Euler
# transport substep, split against the exact radiation half-step): at the default ``n_tau = 0.5`` that bias
# *suppresses* the emergent polar amplification (endpoint ratio ≈ 1.5 vs the dt-free ≈ 2.05), the global
# mean staying exact throughout. The linear assembly below mirrors ``steady_linear``'s (harmonic faces) and
# is pinned **bit-for-bit** against it by the RH = 0 reduction test.
# --------------------------------------------------------------------------- #
def _steady_tridiag_solve(x: np.ndarray, dx: np.ndarray, dxc: np.ndarray,
                          Dcoef: np.ndarray, B: float, rhs: np.ndarray) -> np.ndarray:
    """One direct linear steady solve ``(L_T − B·I)T = rhs`` with cell coefficient ``Dcoef = D_eff·(1−x²)``.

    Assembles the conservative transport operator ``L_T = d/dx[D_eff(1−x²) d/dx]`` with **harmonic-mean
    faces** (the engine's idiom), subtracts the radiative damping ``B``, and solves the tridiagonal system
    in one banded solve. Identical in form to :meth:`planet.ebm.EnergyBalanceModel.steady_linear` — pinned
    bit-for-bit by the RH = 0 reduction (constant ``D_eff`` ⟹ the dry linear EBM).
    """
    Dface = 2.0 * Dcoef[:-1] * Dcoef[1:] / (Dcoef[:-1] + Dcoef[1:])      # harmonic-mean interior faces
    Tt = Dface / dxc                                                     # interior-face transmissibilities
    n = x.size
    sub = np.zeros(n); diag = np.zeros(n); sup = np.zeros(n)
    sup[:-1] += Tt / dx[:-1]; diag[:-1] += -Tt / dx[:-1]
    sub[1:] += Tt / dx[1:];   diag[1:] += -Tt / dx[1:]
    diag = diag - B                                                      # (L_T − B·I): diagonally dominant
    ab = np.zeros((3, n))
    ab[0, 1:] = sup[:-1]; ab[1, :] = diag; ab[2, :-1] = sub[1:]
    return solve_banded((1, 1), ab, rhs)


def moist_steady_direct(params: EBMParams, D_s: float, absorbed_fn,
                        RH: float = moist.RH_DEFAULT, A: Optional[float] = None, T_init=None,
                        tol: float = 1e-11, max_iter: int = 500) -> ClimateState:
    """The **dt-free** moist steady climate by Picard iteration on the frozen-``D_eff`` linear solve.

    Each Picard step freezes ``D_eff(T_k) = D_s·(1 + β(T_k))``, assembles ``(L_T[D_eff] − B·I)`` with
    ``D_eff`` **inside** the conservative divergence, and solves ``(L_T − B·I)T_{k+1} = A − S(1−α)`` in one
    tridiagonal solve — repeating until ``max|T_{k+1} − T_k| < tol``. Because there is **no time-stepping**,
    there is **no operator-splitting error**: the fixed point is the true nonlinear steady state (= the
    ``Δt → 0`` limit of :func:`moist_equilibrium`), reached in ~20 iterations. This is the
    **headline path** — :func:`polar_amplification` uses it. The nonlinear generalisation of
    :meth:`planet.ebm.EnergyBalanceModel.steady_linear`; valid only for a **constant-albedo** (ice-OFF)
    absorbed field (raises otherwise), exactly as ``steady_linear`` does. ``RH = 0`` ⟹ constant ``D_eff``
    ⟹ one effective solve ⟹ the dry linear EBM, bit-for-bit ``steady_linear`` (the plumbing reduction).
    """
    grid = uniform_grid(1.0, params.n_cells)
    x = grid.centers
    dx = grid.widths
    dxc = np.diff(grid.centers)
    B = params.B
    A_olr = params.A if A is None else float(A)
    a_cold = absorbed_fn(x, np.full(x.shape, -100.0))
    a_warm = absorbed_fn(x, np.full(x.shape, 100.0))
    if not np.allclose(a_cold, a_warm):
        raise ValueError("moist_steady_direct requires a state-independent (constant-albedo) absorbed "
                         "field; the ice feedback must go through moist_equilibrium (the relaxation)")
    absorbed = np.asarray(a_cold, dtype=float)
    rhs = A_olr - absorbed
    if T_init is None:
        T_init = 30.0 + (-30.0 - 30.0) * x                            # warm equator → cold pole (Earth-like)
    T = np.array(np.broadcast_to(np.asarray(T_init, dtype=float), x.shape), dtype=float)
    converged, it = False, 0
    for it in range(1, max_iter + 1):
        Dcoef = effective_diffusivity(T, D_s, RH) * (1.0 - x ** 2)     # D_eff(T) inside the divergence
        T_new = _steady_tridiag_solve(x, dx, dxc, Dcoef, B, rhs)
        if np.max(np.abs(T_new - T)) < tol:
            T = T_new
            converged = True
            break
        T = T_new
    Tbar = float(np.mean(T))                                          # L = 1 ⟹ area mean (∫₀¹ T dx)
    return ClimateState(
        x=x, T=T, global_mean_T=Tbar,
        ice_line_lat=ice_line_latitude(x, T, params.T_freeze),
        net_toa=float(np.mean(absorbed_fn(x, T)) - A_olr - B * Tbar),
        converged=converged, iterations=it,
    )


def equator_pole_contrast(state: ClimateState) -> float:
    """Equator-to-pole temperature contrast ``T(0°) − T(90°)`` (°C) — the recalibration target."""
    return float(state.T[0] - state.T[-1])


def P2_amplitude(state: ClimateState) -> float:
    """North P₂ amplitude ``T₂ = ⟨T·P₂⟩/⟨P₂²⟩`` (°C) — the smooth-shape alternative recalibration scalar.

    Area means ``⟨·⟩ = ∫₀¹·dx`` with ``⟨P₂²⟩ = 1/5``. Matching ``T₂`` instead of the endpoint contrast
    leaves the polar-amplification factor essentially unchanged (PA is set by the shape of ``β(T)``) —
    the named "target is a modeling choice" invariance.
    """
    P2 = legendre_P2(state.x)
    return float(np.mean(state.T * P2) / np.mean(P2 ** 2))


# --------------------------------------------------------------------------- #
# The recalibration — re-derive D_s so the moist present climate matches the dry one (the wall).
# --------------------------------------------------------------------------- #
def recalibrate_sensible_D(params: Optional[EBMParams] = None, RH: float = moist.RH_DEFAULT,
                           target: str = "contrast", D_lo: float = 0.05, D_hi: float = 1.0) -> float:
    """The sensible ``D_s`` whose moist present-day climate matches the dry rung-0 ``target`` (W m⁻² K⁻¹).

    Re-derives the dry diffusivity so explicit latent transport does **not** double-count the latent heat
    already lumped into rung-0's effective ``D = 0.555`` (the named wall). Roots ``f(D_s) = (moist target)
    − (dry target)`` by bisection, both climates run with **constant albedo** at the present ``A``. The
    ``target`` is the equator-to-pole ``"contrast"`` (default; ``"T2"`` matches the P₂ amplitude — the
    factor is invariant to the choice). For Earth defaults + RH 0.8 this gives ``D_s ≈ 0.28`` (< 0.555,
    because the moisture-amplified ``D_eff`` transports more, so less dry ``D`` is needed to match).
    """
    if params is None:
        params = EBMParams()
    measure = {"contrast": equator_pole_contrast, "T2": P2_amplitude}[target]
    absorbed = constant_albedo_absorbed(params)
    dry_target = measure(moist_steady_direct(params, params.D, absorbed, RH=0.0))   # dt-free dry reference

    def gap(D_s: float) -> float:
        return measure(moist_steady_direct(params, D_s, absorbed, RH=RH)) - dry_target

    return float(brentq(gap, D_lo, D_hi, xtol=1e-4))


# --------------------------------------------------------------------------- #
# The headline experiment — polar amplification under a uniform ΔA, dry-uniform null beside it.
# --------------------------------------------------------------------------- #
# Band edges for the area-averaged polar-amplification metric (the less-generous, less-bias-exposed
# companion to the single-endpoint ratio): mean warming poleward of 60° over mean warming equatorward
# of 30°. The endpoint ratio reads off the two extreme cells (the warmest/coldest, and the polar cell
# sits on the harmonic-face bias ebm.py flags — the bias largely cancels in a δT *difference*, but the
# band average is reported alongside so the most-generous number is not mistaken for *the* number.
POLAR_BAND_DEG = 60.0
TROPICAL_BAND_DEG = 30.0


def _band_pa(phi: np.ndarray, dT: np.ndarray) -> float:
    """Band-averaged polar amplification ``mean(δT | φ≥60°) / mean(δT | φ≤30°)``."""
    return float(np.mean(dT[phi >= POLAR_BAND_DEG]) / np.mean(dT[phi <= TROPICAL_BAND_DEG]))


@dataclass(frozen=True)
class PolarAmplification:
    """Emergent polar amplification: dry-vs-moist warming under a uniform OLR forcing ``ΔA`` (°C).

    ``phi`` latitudes (deg); ``D_s`` the recalibrated sensible diffusivity; ``RH``/``dA`` the closure
    and forcing. ``*_present``/``*_warm`` are the **dt-free** (:func:`moist_steady_direct`) equilibrium
    climates before/after the forcing for the ``dry`` (constant ``D``) and ``moist`` (``D_eff(T)``)
    models; ``delta_T_dry``/``delta_T_moist`` the warming profiles. ``mean_delta_T`` is the **pinned**
    global-mean warming ``ΔA/B`` (both models share it to machine precision — the conservation anchor).
    The polar amplification is reported **two ways** (name the metric): ``pa_moist`` / ``pa_dry`` are the
    **single-endpoint** ratios ``δT(pole)/δT(equator)`` (the headline ~2.05 for moist; ≈ 1 the dry null),
    and ``pa_moist_band`` / ``pa_dry_band`` are the **area-band** ratios ``mean(δT|φ≥60°)/mean(δT|φ≤30°)``
    (the less-generous, less-bias-exposed companion, ~1.80 for moist) — both honest "polar amplification";
    the band average is carried so the most-generous endpoint number is not mistaken for *the* number.
    ``polar_excess`` is the moist pole warming **minus** the pinned mean (the redistributed excess). Plain
    arrays.
    """

    phi: np.ndarray
    D_s: float
    RH: float
    dA: float
    dry_present: ClimateState
    dry_warm: ClimateState
    moist_present: ClimateState
    moist_warm: ClimateState
    delta_T_dry: np.ndarray
    delta_T_moist: np.ndarray
    mean_delta_T: float
    pa_dry: float
    pa_moist: float
    pa_dry_band: float
    pa_moist_band: float
    polar_excess: float


def polar_amplification(params: Optional[EBMParams] = None, RH: float = moist.RH_DEFAULT,
                        dA: float = 10.0, D_s: Optional[float] = None) -> PolarAmplification:
    """Warm dry and moist EBMs by a uniform ``ΔA`` and measure the emergent polar amplification.

    The rung-2.5 headline. Builds the constant-albedo present-day climate for the **dry** EBM
    (``D = params.D``) and the **moist** EBM (``D_eff(T)`` at the recalibrated ``D_s``) — all four climates
    by the **dt-free** :func:`moist_steady_direct` (no operator-splitting bias) — warms both by a uniform
    OLR reduction ``A → A − ΔA`` (the CO₂ proxy), and returns the :class:`PolarAmplification`. The dry
    model warms **exactly uniformly** (``δT = ΔA/B`` — transport of a uniform field is zero), so its
    ``pa_dry ≈ 1`` is the clean null; the moist model **redistributes** that same pinned ``⟨δT⟩``
    poleward (``pa_moist`` ~ 2.05 for Earth defaults). ``D_s`` defaults to
    :func:`recalibrate_sensible_D` (so the moist present climate matches the dry contrast); pass a value
    to skip the recalibration.
    """
    if params is None:
        params = EBMParams()
    if D_s is None:
        D_s = recalibrate_sensible_D(params, RH)
    absorbed = constant_albedo_absorbed(params)
    A_warm = params.A - float(dA)

    dry_present = moist_steady_direct(params, params.D, absorbed, RH=0.0)
    dry_warm = moist_steady_direct(params, params.D, absorbed, RH=0.0, T_init=dry_present.T, A=A_warm)
    moist_present = moist_steady_direct(params, D_s, absorbed, RH=RH)
    moist_warm = moist_steady_direct(params, D_s, absorbed, RH=RH, T_init=moist_present.T, A=A_warm)

    dTd = dry_warm.T - dry_present.T
    dTm = moist_warm.T - moist_present.T
    phi = dry_present.latitude_deg()
    return PolarAmplification(
        phi=phi, D_s=float(D_s), RH=float(RH), dA=float(dA),
        dry_present=dry_present, dry_warm=dry_warm,
        moist_present=moist_present, moist_warm=moist_warm,
        delta_T_dry=dTd, delta_T_moist=dTm,
        mean_delta_T=float(np.mean(dTm)),
        pa_dry=float(dTd[-1] / dTd[0]),
        pa_moist=float(dTm[-1] / dTm[0]),
        pa_dry_band=_band_pa(phi, dTd),
        pa_moist_band=_band_pa(phi, dTm),
        polar_excess=float(dTm[-1] - np.mean(dTm)),
    )
