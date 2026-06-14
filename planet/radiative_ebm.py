"""The gray-OLR energy-balance model — emergent latitudinal radiative structure (Planet rung 4 completion).

Rung 4 (:mod:`planet.radiation`) made the outgoing longwave **emergent** from a gray radiative–convective
column, but only as a **single global-mean column**: it derived where climlab's slope ``B = 2`` comes from
(Planck − water-vapour + the omitted lapse-rate feedback) without ever letting that radiation *drive* a
climate. This module is the named rung-4 completion: it wires the emergent ``OLR(Ts, τ)`` **per latitude**
into the energy-balance spine, so the planet's temperature is set by *real radiation* rather than the
prescribed linear ``A + B·T``. It is a **separate sibling that runs alongside rung-0** — :mod:`planet.ebm`
is untouched and stays the default everywhere (the :mod:`planet.moist_ebm` / :mod:`planet.sphere_ebm` /
:mod:`planet.baroclinic_qg` discipline). Built spike-first; the headline **sign** was decided by the spike,
not assumed (see below).

The headline — the OLR slope is NOT a global constant, and it makes warming TROPICAL
-----------------------------------------------------------------------------------
Rung-0 assumes one number ``B = 2`` everywhere. The emergent ``OLR(Ts)`` is **nonlinear**, so its *local*
slope ``B_loc(Ts) = dOLR/dTs`` (:func:`local_radiative_slope`) **varies across the planet** — and that is
the whole story. Two effects compete in ``B_loc``:

* **Water-vapour feedback** is stronger where it is warm (Clausius–Clapeyron), which pulls ``B_loc`` *down*
  at the equator → favours **tropical** amplification.
* **The Planck term** ``4σT³`` is larger where it is warm, which pulls ``B_loc`` *up* at the equator →
  favours **polar** amplification.

Which wins is genuinely uncertain a priori, so it was **measured, not assumed** (the advisor's discriminator:
"whichever latitude has the smallest ``B_loc`` warms most"). The water-vapour feedback wins decisively: at
the climlab-matched loading ``B_loc`` runs ~**1.1 at a 30 °C equator** up to ~**2.4 in the mid-latitudes**
and back to ~**1.8 at a −30 °C pole**, so the *equator* has the smallest local damping. Under a uniform
forcing the warming concentrates there — **tropical amplification** (endpoint ``δT(pole)/δT(equator) ≈ 0.7``,
:class:`TropicalAmplification`), the **mirror image** of rung-2.5's moisture-*transport* polar amplification
(:mod:`planet.moist_ebm`, dt-free ~1.8–2.05). The two rungs make a clean "two mechanisms pull opposite ways" pair:
the water-vapour *radiative* feedback alone favours the tropics, while it is *transport* (and the
lapse-rate and ice feedbacks held out of scope) that make Earth's poles amplify.

**But the sign itself rides the water-vapour loading — it is not robust the way rung-2.5's is.** Turn the
loading down and the Planck term wins: the bare ``B_noWV`` *rises* with temperature, so a dry,
Planck-dominated planet damps **most** at the warm equator and is **polar** (the spike confirms ``amp > 1``
below a crossover loading ≈ 0.15). Earth's calibrated loading (≈ 0.35) sits well into the tropical regime,
but a drier planet flips back. So *both* the sign and the magnitude ride the loading (the wall) — whereas
rung-2.5's polar direction is genuinely robust to its RH wall (its ``β`` is always larger in the tropics).
The "mirror" is real *at Earth's loading*; it is not a symmetry of equally-robust mechanisms.

Why a Newton solve, not the Strang relaxation spine (a numerical finding)
------------------------------------------------------------------------
Rung-0 and rung-2.5 reach equilibrium by **Strang operator splitting**, whose radiation half-step is the
*analytically exact* solution of the linear ``−B·T`` relaxation. That makes the relaxed steady state's
**global mean** dt-free (exact), but **not its shape**: split against a *backward-Euler* transport step the
*shape* carries an **O(Δt)** splitting bias — which is exactly why both rungs read their headline off a
**direct** solve (:meth:`planet.ebm.EnergyBalanceModel.steady_linear`, and rung-2.5's Picard
:func:`planet.moist_ebm.moist_steady_direct`), never the relaxation. With a **nonlinear** ``OLR(Ts)`` even
that mean-exactness is lost: a frozen-slope linearised half-step leaves an operator-splitting error that
**does not vanish at equilibrium** — not even in the mean (the spike found the relaxed steady state drifting
with the step size, ``⟨T⟩ ≈ 18.9 → 18.0`` as ``n_τ → 0``, converging onto the Newton answer), and
near the warm-equator runaway edge (next paragraph) the local half-step goes outright unstable. So the
radiative EBM is solved by **coupled Newton iteration** on the steady residual

    R(T) = L_T·T + S(x)(1−α) − OLR(T) = 0,        J = L_T − diag(B_loc(T)),

i.e. the **nonlinear generalisation of rung-0's "direct" linear solve** (:meth:`planet.ebm.EnergyBalanceModel.steady_linear`)
— the *same* engine-pinned transport tridiagonal ``L_T`` (:meth:`~planet.ebm.EnergyBalanceModel._transport_tridiag`,
so the transport cannot drift from the engine), with the linear ``−B·I`` replaced by the local Jacobian
``−diag(B_loc)``. The transport coupling stabilises the equator, so Newton converges (~6 iterations) even
where the *local* column would run away. With constant albedo the radiation is the only nonlinearity — the
clean Hwang–Frierson-style experiment that isolates it (the moist-EBM precedent).

The runaway edge, now made local (a finding the per-latitude wire exposes)
--------------------------------------------------------------------------
A single global-mean column hides a hot trap: at the warm end ``B_loc`` collapses toward **zero** (the
water-vapour feedback approaching the Komabayashi–Ingersoll runaway, the named rung-4 edge). At the
**default** water-vapour loading (:data:`planet.radiation.WATER_VAPOUR_FRACTION` = 0.5, the value that
makes the global-mean ``B ≈ 1.33``) the equatorial column is already *past* its local runaway (``B_loc`` < 0
for a > ~32 °C surface), so a per-latitude wire has **no stable local equilibrium at the equator** — the
relaxation diverges there. This module therefore runs at the **climlab-matched loading**
(:func:`climlab_matched_column`: the water-vapour fraction for which the present-day global-mean slope is
exactly ``B = 2``, ≈ 0.35), where ``B_loc`` stays positive everywhere and a stable climate exists. That the
per-latitude wire surfaces a local runaway the global column did not is itself the rung-4 "linearisation
breaks far from present" edge, now concrete and latitude-resolved.

Validation triad (plan §3) — what is asserted tight vs loose
------------------------------------------------------------
* **Reduction (bit-for-bit / structural).** A **linear** ``olr_fn`` (:func:`linear_olr_fn`, ``A + B·T``)
  makes the residual affine, so one Newton step reproduces :meth:`planet.ebm.EnergyBalanceModel.steady_linear`
  to machine precision — the radiative EBM *is* rung-0's direct solve with the radiation swapped, and the
  whole departure below is attributable to the OLR **curvature** alone.
* **Conservation (tight).** At Newton convergence the residual is ~machine, and since the transport conserves
  ``∫T dx`` (``⟨L_T·T⟩ = 0``) the global energy balance ``⟨S(1−α)⟩ = ⟨OLR(T)⟩`` (net-TOA = 0) holds to that
  tolerance. **Unlike rung-2.5 the global-mean warming is *not* pinned at ``ΔA/B``**: ``OLR`` is concave, so
  by Jensen ``⟨OLR(T)⟩ ≠ OLR(⟨T⟩)`` and the water-vapour feedback amplifies the *mean* response too
  (``⟨δT⟩ > ΔA/B_tan``) — measured, not asserted (the advisor's catch; the moist-EBM "redistribution around a
  pinned mean" sentence deliberately does **not** transfer).
* **The unlock (real but loose).** **Tropical** amplification (``δT(pole)/δT(equator) < 1``) **at the
  Earth-calibrated loading** — *both the sign and the magnitude ride the water-vapour loading* (the wall): a
  drier, Planck-dominated planet is **polar** (``amp > 1`` below a crossover loading ≈ 0.15), so this is a
  "tropical at Earth's loading" result, not a loading-independent one. The *mean-state* contrast is
  **essentially unchanged** from rung-0 at the same ``D`` (the loading-matched *average* slope is ``≈ 2``, so
  the equator-to-pole spread matches) — the latitudinal signal lives in the *warming response*, not the
  present climate. What the present state *does* show is a **uniform warm shift** (``⟨T⟩`` ~2 °C above
  rung-0): ``OLR`` is concave, so by Jensen the planet must run warmer on average to emit the same mean.
  *(Part of the ~0.68 magnitude is runaway-proximity — the warmed equator reaches ~39 °C where ``B_loc`` ≈
  0.5, still positive/stable but close to the local runaway. Compared throughout to rung-0's dt-free*
  ``steady_linear`` *reference, not its relaxation default, whose shape carries an O(Δt) splitting error.)*

Named scope edges
-----------------
* **Only the radiation is re-opened.** ``D`` is held at rung-0's value: the present-day contrast is already
  ≈ rung-0's at this loading (above), so there is nothing to recalibrate ``D`` *for*. ``D`` sets the
  amplification *magnitude*, not its *sign*, so recalibrating it (the rung-2.5 move, done there for a
  *different* reason — latent double-counting) is an **available knob, not exercised here**.
* **Constant albedo; fixed lapse rate; clear-sky.** The ice-albedo feedback is held off (the clean
  experiment), and the gray column's fixed Γ supplies no lapse-rate feedback and no clouds — the named
  rung-4 within-rung upgrades (a moist-adiabatic Γ, band physics, clouds) are not opened here.
* **Forcing is a uniform ``ΔA`` (an OLR offset), not ``ΔS₀``** — the clean knob (``S₀`` carries a tropical
  structure that would confound the latitudinal signal), as in rung-2.5.

Units — SI, climlab-conventional (W m⁻², °C, x = sin φ)
-------------------------------------------------------
``OLR``/``S``/``ΔA`` in W m⁻²; ``B_loc``/``D`` in W m⁻² K⁻¹; ``T`` in **°C** (converted to K inside the
gray column, whose ``σT⁴`` is a kelvin law); ``x = sin φ`` on [0, 1]; latitudes in degrees.
See [[planet-rung4-radiation]], [[moist-ebm-source]], [[ebm-radiation-source]], [[planet-plan]] §10.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np
from scipy.linalg import solve_banded
from scipy.optimize import brentq

from . import radiation as rad
from .albedo import EBMParams
from .ebm import (
    A_OLR, B_OLR, D_TRANSPORT, WATER_DEPTH,
    ClimateState, EnergyBalanceModel, ice_line_latitude,
)
from .moist_ebm import constant_albedo_absorbed

# An ``olr_fn`` maps a temperature field (°C) to ``(OLR, B_loc)`` — the outgoing longwave (W m⁻²) and its
# local slope dOLR/dTs (W m⁻² K⁻¹), both per cell. It is the single injection point that makes the EBM
# radiation-agnostic: the linear one recovers rung-0, the gray one carries the emergent nonlinearity.
OlrFn = Callable[[np.ndarray], "tuple[np.ndarray, np.ndarray]"]

PRESENT_GLOBAL_B = 2.0          # W m⁻² K⁻¹ — the climlab present-day slope the loading is matched to
PRESENT_GLOBAL_TS = rad.PRESENT_SURFACE_T   # K — the surface temperature the slope is matched at (288 K)


# --------------------------------------------------------------------------- #
# The emergent OLR and its local slope — looped over radiation.py's scalar column (single source of truth).
# --------------------------------------------------------------------------- #
def column_olr(column: rad.GrayRadiationColumn, T_celsius, water_vapour: bool = True) -> np.ndarray:
    """Emergent ``OLR`` (W m⁻²) for each surface temperature in ``T_celsius`` (°C), via the gray column.

    A plain loop over :meth:`planet.radiation.GrayRadiationColumn.outgoing_longwave` (the validated rung-4
    forward model) — *one source of truth* for the radiation, not a re-implementation. Each latitude is its
    own column at its local surface temperature, so the water-vapour loading scales with that latitude's
    warmth (Clausius–Clapeyron), which is exactly what gives ``OLR`` its latitude-varying slope.
    """
    Ts_K = np.atleast_1d(np.asarray(T_celsius, dtype=float)) + 273.15
    return np.array([column.outgoing_longwave(float(t), water_vapour=water_vapour) for t in Ts_K])


def local_radiative_slope(column: rad.GrayRadiationColumn, T_celsius, dT: float = 1.0,
                          water_vapour: bool = True) -> np.ndarray:
    """The **local** OLR slope ``B_loc(Ts) = dOLR/dTs`` (W m⁻² K⁻¹) per cell — the discriminator.

    Central difference of :func:`column_olr` about each cell's temperature. This is the curve that decides
    the headline sign: the latitude with the **smallest** ``B_loc`` warms most under a uniform forcing.
    The water-vapour feedback makes it small at the warm equator (→ tropical amplification); ``water_vapour
    =False`` returns the bare Planck slope (which alone would do the opposite).
    """
    T = np.asarray(T_celsius, dtype=float)
    hi = column_olr(column, T + dT, water_vapour=water_vapour)
    lo = column_olr(column, T - dT, water_vapour=water_vapour)
    return (hi - lo) / (2.0 * dT)


def gray_olr_fn(column: rad.GrayRadiationColumn, forcing: float = 0.0, dT: float = 1.0) -> OlrFn:
    """Build the gray ``olr_fn(T) -> (OLR, B_loc)`` for the Newton solve, with a uniform ``forcing`` offset.

    ``forcing`` (W m⁻²) is subtracted from ``OLR`` — a uniform reduction in outgoing longwave is the clean
    warming knob (the CO₂ proxy, ``ΔA`` in rung-2.5's language), applied to the radiation so the solver
    stays a pure root-find. ``dT`` is the half-width of the slope's central difference.
    """
    def f(T_celsius: np.ndarray):
        return column_olr(column, T_celsius) - forcing, local_radiative_slope(column, T_celsius, dT=dT)
    return f


def linear_olr_fn(A: float = A_OLR, B: float = B_OLR, forcing: float = 0.0) -> OlrFn:
    """Build the linear ``olr_fn(T) -> (A + B·T − forcing, B)`` — the rung-0 OLR, for the reduction.

    Feeding this to :class:`RadiativeEBM` makes the steady residual affine, so the Newton solve reproduces
    :meth:`planet.ebm.EnergyBalanceModel.steady_linear` bit-for-bit — the structural reduction that pins the
    radiative EBM to rung-0 and attributes the whole gray departure to the OLR curvature.
    """
    def f(T_celsius: np.ndarray):
        T = np.asarray(T_celsius, dtype=float)
        return A + B * T - forcing, np.full_like(T, float(B))
    return f


def climlab_matched_column(target_B: float = PRESENT_GLOBAL_B, Ts: float = PRESENT_GLOBAL_TS,
                           lo: float = 0.05, hi: float = 0.95, **column_kw) -> rad.GrayRadiationColumn:
    """The gray column whose **present-day global-mean slope is** ``target_B`` (≈ climlab's 2) — the loading.

    Bisects the water-vapour fraction (:data:`planet.radiation.WATER_VAPOUR_FRACTION`, the wall) so the
    single-column slope at ``Ts`` (288 K) equals ``target_B``; the optical depth is then re-calibrated to the
    present greenhouse by :func:`planet.radiation.calibrate_column`. This is the principled operating point
    for the *wire* (distinct from rung-4's default 0.5 loading): it makes the global-mean reduce to rung-0's
    ``B = 2`` **and** keeps ``B_loc`` positive everywhere (sub-runaway), so a stable per-latitude climate
    exists. Returns ≈ 0.35 for Earth defaults.
    """
    def gap(wvf: float) -> float:
        col = rad.calibrate_column(wv_fraction=wvf, **column_kw)
        return float(local_radiative_slope(col, [Ts - 273.15])[0]) - target_B

    wvf = float(brentq(gap, lo, hi, xtol=1e-4))
    return rad.calibrate_column(wv_fraction=wvf, **column_kw)


# --------------------------------------------------------------------------- #
# The radiative EBM — coupled Newton steady solve (nonlinear generalisation of ebm.steady_linear).
# --------------------------------------------------------------------------- #
class RadiativeEBM:
    """A latitudinal EBM whose outgoing longwave is an injected ``olr_fn`` — solved by coupled Newton.

    Reuses :class:`planet.ebm.EnergyBalanceModel` purely for its **engine-pinned transport** (the
    geometry ``x = sin φ``, the heat capacity ``C``, and the tridiagonal operator ``L_T``); the radiation is
    replaced by ``olr_fn(T) -> (OLR, B_loc)``. With a linear ``olr_fn`` it *is* rung-0's direct solve; with
    the gray ``olr_fn`` it carries the emergent nonlinearity. The albedo enters through the same
    ``absorbed_fn(x, T)`` callable as rung-0 (held constant for the clean experiment), so the only
    nonlinearity in the residual is the radiation.
    """

    def __init__(self, olr_fn: OlrFn, D: float = D_TRANSPORT, n_cells: int = 180,
                 water_depth: float = WATER_DEPTH, face: str = "harmonic"):
        self.olr_fn = olr_fn
        # Build rung-0's machinery for the transport ONLY (A, B here are unused — radiation comes from
        # olr_fn). _transport_tridiag is the engine-pinned operator, so transport cannot drift.
        self._ebm = EnergyBalanceModel(D=D, n_cells=n_cells, water_depth=water_depth, face=face)
        self.x = self._ebm.x
        self.C = self._ebm.C
        self._sub, self._diag, self._sup = self._ebm._transport_tridiag()

    def _LT_apply(self, T: np.ndarray) -> np.ndarray:
        """Apply the transport operator ``L_T·T`` (W m⁻²) — the tridiagonal mat-vec (conserves ``∫T dx``)."""
        out = self._diag * T
        out[:-1] += self._sup[:-1] * T[1:]
        out[1:] += self._sub[1:] * T[:-1]
        return out

    def equilibrium(self, absorbed_fn, T_init=None, tol: float = 1e-10, max_iter: int = 200,
                    clip: Optional[float] = 10.0) -> ClimateState:
        """Newton-solve the steady residual ``L_T·T + S(1−α) − OLR(T) = 0``; return its :class:`ClimateState`.

        Each iteration linearises the radiation (``OLR ≈ OLR(T) + B_loc·δT``) and solves the banded system
        ``(L_T − diag(B_loc))·δT = −R`` for the update, optionally clipped to ``±clip`` °C for robustness in
        the steep-water-vapour equatorial regime. Converges when ``max|δT| < tol``. The transport coupling
        stabilises the equator, so this succeeds where an operator-split relaxation would not.

        Parameters
        ----------
        absorbed_fn : callable ``(x, T) -> ndarray``
            The absorbed shortwave ``S(x)(1−α)`` (W m⁻²) — :func:`planet.moist_ebm.constant_albedo_absorbed`
            for the clean (ice-free) experiment.
        T_init : array | float | None
            Initial field (°C); ``None`` → a uniform 15 °C (Newton is robust from it).
        tol, max_iter, clip
            Newton tolerance (°C), iteration cap, and the per-step clip (``None`` to disable).
        """
        n = self.x.size
        T = (np.full(n, 15.0) if T_init is None
             else np.array(np.broadcast_to(np.asarray(T_init, dtype=float), self.x.shape), dtype=float))
        converged, it = False, 0
        for it in range(1, max_iter + 1):
            OLR, B_loc = self.olr_fn(T)
            absorbed = np.asarray(absorbed_fn(self.x, T), dtype=float)
            R = self._LT_apply(T) + absorbed - OLR
            diag = self._diag - B_loc                       # Jacobian J = L_T − diag(B_loc)
            ab = np.zeros((3, n))
            ab[0, 1:] = self._sup[:-1]
            ab[1, :] = diag
            ab[2, :-1] = self._sub[1:]
            dT = solve_banded((1, 1), ab, -R)
            if clip is not None:
                dT = np.clip(dT, -clip, clip)
            T = T + dT
            if np.max(np.abs(dT)) < tol:
                converged = True
                break
        OLR, _ = self.olr_fn(T)
        Tbar = self._ebm.global_mean(T)
        net_toa = float(np.mean(absorbed_fn(self.x, T)) - np.mean(OLR))
        return ClimateState(
            x=self.x, T=T, global_mean_T=Tbar,
            ice_line_lat=ice_line_latitude(self.x, T, self._ebm.T_freeze),
            net_toa=net_toa, converged=converged, iterations=it,
        )


# --------------------------------------------------------------------------- #
# The headline experiment — tropical amplification, the present-tangent uniform-B null beside it.
# --------------------------------------------------------------------------- #
POLAR_BAND_DEG = 60.0
TROPICAL_BAND_DEG = 30.0


def _band_amp(phi: np.ndarray, dT: np.ndarray) -> float:
    """Band-averaged amplification ``mean(δT | φ≥60°) / mean(δT | φ≤30°)`` (< 1 ⟹ tropical)."""
    return float(np.mean(dT[phi >= POLAR_BAND_DEG]) / np.mean(dT[phi <= TROPICAL_BAND_DEG]))


@dataclass(frozen=True)
class TropicalAmplification:
    """Emergent tropical amplification: gray-vs-uniform-B warming under a uniform OLR forcing (°C).

    ``phi`` latitudes (deg); ``D``/``forcing`` the transport and the uniform ``ΔA``. ``B_loc_present`` is the
    present-day local-slope profile (the discriminator: smallest at the equator). ``gray_*``/``null_*`` are
    the present/warmed climates for the **gray** (nonlinear ``OLR``) and the **null** (its present-day
    *tangent*, a uniform ``B``); ``delta_T_*`` the warming profiles. The null warms ~uniformly
    (``amp_null ≈ 1``, the clean baseline) while the gray model concentrates warming in the tropics. The
    amplification is reported two ways: ``amp_*`` the **endpoint** ``δT(pole)/δT(equator)`` (the headline,
    ~0.7 for gray) and ``amp_*_band`` the **area-band** companion. ``mean_delta_T_gray`` is the *measured*
    global-mean warming and ``dA_over_B_tan`` the naive ``ΔA/B_tan`` it is **not** pinned to (Jensen + the
    water-vapour feedback amplify the mean). Plain arrays/scalars.
    """

    phi: np.ndarray
    D: float
    forcing: float
    B_loc_present: np.ndarray
    gray_present: ClimateState
    gray_warm: ClimateState
    null_present: ClimateState
    null_warm: ClimateState
    delta_T_gray: np.ndarray
    delta_T_null: np.ndarray
    mean_delta_T_gray: float
    mean_delta_T_null: float
    dA_over_B_tan: float
    amp_gray: float
    amp_null: float
    amp_gray_band: float
    amp_null_band: float


def tropical_amplification(column: Optional[rad.GrayRadiationColumn] = None,
                           params: Optional[EBMParams] = None, forcing: float = 10.0,
                           D: Optional[float] = None) -> TropicalAmplification:
    """Warm the gray and a present-tangent uniform-B EBM by a uniform ``ΔA`` and measure the amplification.

    The rung-4-completion headline. Builds the constant-albedo present-day climate with the **gray**
    ``OLR(Ts)`` and warms it by a uniform OLR reduction ``forcing`` (the CO₂ proxy); the **null** is the
    present-day *tangent* OLR — a uniform ``B = ⟨B_loc⟩`` with offset through the present mean state — warmed
    by the same forcing (it warms exactly uniformly, ``amp_null ≈ 1``, the clean baseline). The gray model
    concentrates the warming in the tropics (``amp_gray < 1``), the mirror of rung-2.5's polar amplification.
    ``column`` defaults to :func:`climlab_matched_column` (global-mean ``B = 2``); ``D`` to rung-0's value.
    """
    if column is None:
        column = climlab_matched_column()
    if params is None:
        params = EBMParams()
    if D is None:
        D = params.D
    absorbed = constant_albedo_absorbed(params)

    gray = RadiativeEBM(gray_olr_fn(column), D=D, n_cells=params.n_cells,
                        water_depth=params.water_depth)
    gray_present = gray.equilibrium(absorbed)
    gray_warm = RadiativeEBM(gray_olr_fn(column, forcing=forcing), D=D, n_cells=params.n_cells,
                             water_depth=params.water_depth).equilibrium(absorbed, T_init=gray_present.T)

    # The null: the present-day tangent — a uniform B (area-mean local slope) with the offset set so the
    # tangent line passes through the present mean state, then warmed by the same forcing.
    B_loc_present = local_radiative_slope(column, gray_present.T)
    B_tan = float(np.mean(B_loc_present))
    A_tan = float(np.mean(column_olr(column, gray_present.T)) - B_tan * np.mean(gray_present.T))
    null = RadiativeEBM(linear_olr_fn(A_tan, B_tan), D=D, n_cells=params.n_cells,
                        water_depth=params.water_depth)
    null_present = null.equilibrium(absorbed, T_init=gray_present.T)
    null_warm = RadiativeEBM(linear_olr_fn(A_tan, B_tan, forcing=forcing), D=D, n_cells=params.n_cells,
                             water_depth=params.water_depth).equilibrium(absorbed, T_init=null_present.T)

    dTg = gray_warm.T - gray_present.T
    dTn = null_warm.T - null_present.T
    phi = gray_present.latitude_deg()
    return TropicalAmplification(
        phi=phi, D=float(D), forcing=float(forcing), B_loc_present=B_loc_present,
        gray_present=gray_present, gray_warm=gray_warm,
        null_present=null_present, null_warm=null_warm,
        delta_T_gray=dTg, delta_T_null=dTn,
        mean_delta_T_gray=float(np.mean(dTg)), mean_delta_T_null=float(np.mean(dTn)),
        dA_over_B_tan=float(forcing / B_tan),
        amp_gray=float(dTg[-1] / dTg[0]), amp_null=float(dTn[-1] / dTn[0]),
        amp_gray_band=_band_amp(phi, dTg), amp_null_band=_band_amp(phi, dTn),
    )
