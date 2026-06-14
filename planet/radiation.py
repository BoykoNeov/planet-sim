"""Gray radiative transfer — where the prescribed OLR ``A + B·T`` comes from (Planet rung 4).

Every rung from 2.5 on names the same wall: **``B`` (the OLR slope) is held fixed** — the linear
outgoing-longwave parameterization ``OLR = A + B·T`` (:mod:`planet.ebm`) is *prescribed* from the
climlab/North constants, so the rung-0 EBM cannot tell you **where ``B = 2`` comes from**, cannot
produce a **water-vapour feedback**, and cannot **compute a CO₂ forcing**. Rung 4 makes that OLR
**emergent**: it computes the longwave flux to space from a **gray radiative–convective column**
(Schwarzschild two-stream radiative transfer over an optical depth set by the greenhouse-gas amount),
and reads ``A`` and ``B`` off the result. It is a **separate column model that runs alongside rung-0**,
not a replacement — :mod:`planet.ebm` is untouched and stays the default everywhere (the
:mod:`planet.moist_ebm` / :mod:`planet.sphere_ebm` / :mod:`planet.baroclinic_qg` sibling discipline).

The gray two-stream equations and their radiative-equilibrium solution (DERIVED, not recalled)
----------------------------------------------------------------------------------------------
With a **flux optical depth** ``τ`` measured downward from the top of atmosphere (the two-stream
diffusivity factor folded into ``τ``), the gray Schwarzschild equations for the upward/downward
hemispheric longwave fluxes ``F↑, F↓`` against the blackbody source ``σT⁴`` are

    dF↑/dτ = F↑ − σT⁴,        dF↓/dτ = σT⁴ − F↓.

Let ``Φ = F↑ − F↓`` (net up) and ``Σ = F↑ + F↓``. Then ``dΦ/dτ = Σ − 2σT⁴`` and ``dΣ/dτ = Φ``.
**Radiative equilibrium** is zero flux divergence, ``dΦ/dτ = 0``, so the net flux ``Φ`` is constant and
equals the outgoing longwave ``F`` at every level, and ``Σ = 2σT⁴``. Integrating ``dΣ/dτ = F`` with the
top boundary ``F↓(0) = 0`` (no longwave from space, so ``Σ(0) = F``) gives ``Σ(τ) = F(1 + τ)`` and hence

    σT⁴(τ) = ½ F (1 + τ),        F = OLR = σTₑ⁴   (set by solar balance),

with the **skin temperature** ``T(0) = Tₑ / 2^¼`` and, from ``F↑(τ_s) = (Σ+Φ)/2 = σTg⁴``, the
**ground temperature** carrying the surface–air discontinuity

    σTg⁴ = ½ F (2 + τ_s).

These three closed forms are the **tight anchor** — derived here from the equations (not a recalled
coefficient, the rung-3 ``K²=2F`` lesson), and reproduced by the numerical two-stream solver
(:func:`solve_gray_equilibrium`) to ~2nd order in the layer thickness with ``OLR`` machine-exact.

The emergent OLR and the decomposition ``B ≈ 4σTₑ³ − (water-vapour feedback)`` (the unlock)
-------------------------------------------------------------------------------------------
For a **radiative–convective** column (convection holds a lapse-rate floor ``Γ``; the atmosphere is no
longer in pure radiative equilibrium) the outgoing longwave for a given surface temperature ``Ts`` and
optical-depth profile ``τ(p)`` is the transmission-weighted emission (:func:`outgoing_longwave`)

    OLR = σTs⁴ e^(−τ_s) + ∫₀^{τ_s} σT(τ)⁴ e^(−τ) dτ.

Calibrating the present-day total optical depth so ``OLR(288 K) = 239 W m⁻²`` fixes the 33 K greenhouse
**by construction** — so the emergent OLR and the linear ``A + B·T`` agree in *value* at present, and the
only open quantity is the **slope** ``B = dOLR/dTs`` (``A`` and ``B`` are linked through the operating
point, ``A = 239 − B·T̄``; we recover the forced point and the slope is the finding). That slope
**decomposes**:

* **Clear-sky, no water-vapour feedback** (optical depth held fixed): ``B ≈ 4σTₑ³ ≈ 3.8 W m⁻² K⁻¹`` — the
  Planck slope of the ``Tₑ ≈ 255 K`` emission level (the bare-surface ``4σTs³ ≈ 5.4`` cut down by the
  greenhouse). Lands at the Soden & Held (2006) Planck order ``|λ₀| ≈ 3.1–3.2`` and well **above**
  climlab's prescribed ``2``.
* **With water vapour** (optical depth rises with ``Ts`` through Clausius–Clapeyron, reusing
  :func:`planet.moist.saturation_specific_humidity`): the emission level lifts to colder air, ``OLR``
  rises less per kelvin, and the slope **drops by ≈ 2 W m⁻² K⁻¹** — the positive feedback the fixed-``B``
  rung-0 structurally cannot represent. That drop lands at the Soden & Held water-vapour feedback order
  ``λ_wv ≈ 1.8`` (:data:`SH_WATER_VAPOUR`) — **order-validated, not tuned** — leaving a gray net
  ``B ≈ 1.3`` ≈ the clear-sky Planck+water-vapour combination ``3.2 − 1.8``.

So **climlab's ``B = 2`` ≈ Planck slope − water-vapour feedback + the lapse-rate feedback the gray column
omits**: the gray net (1.3) sits *below* climlab's 2 by ~ the Soden & Held lapse-rate feedback
``|λ_LR| ≈ 0.84`` (:data:`SH_LAPSE_RATE`) — which a **fixed** convective lapse rate (a uniform profile
shift under warming) **cannot produce** (climlab's obs-tuned ``B`` folds it, and clouds, in). Every term is
pinned to an independent feedback estimate, so the decomposition is **non-circular** (the rung-2.5
frozen-``D_eff`` and ITCZ closed-form attribution flavour, stronger than a bare "trade"). **NB the ``0.84``
here is the *imported* global-mean Soden & Held value — a touchstone for what the fixed-``Γ`` default
omits; turning the feedback *emergent* (``moist_adiabat=True``, :meth:`feedback_kernel`) gives the
**tropical** ``≈ 1.5`` and **overshoots** 2 (with-WV ``B ≈ 3.1``), not lands at it — see the lapse-rate
scope edge below and :mod:`planet.demo_lapse_rate`.**

Triad (plan §3) — what is asserted tight vs loose
-------------------------------------------------
* **Analytical limit (tight).** The numerical two-stream solver reproduces the derived gray-RE profile,
  skin and ground temperatures to ~2nd order in layer thickness, with ``OLR = σTₑ⁴`` machine-exact
  (energy conservation); and the no-feedback slope sits near the ``4σTₑ³`` Planck touchstone.
* **The unlock (real but loose).** The emergent ``B`` decomposition (Planck slope ≈ 3.4, a water-vapour
  feedback ≈ 2 dropping the net to ≈ 1.3) — **order-validated** against Soden & Held (2006) (Planck 3.2,
  ``λ_wv`` 1.8, the residual-to-climlab ≈ the *imported* global-mean lapse-rate ``λ_LR`` 0.84), **not
  tuned**. Loose: the exact magnitudes ride on the water-vapour optical-depth loading
  (:data:`WATER_VAPOUR_FRACTION`), the **wall**. (The *emergent* ``λ_LR`` from a moist adiabat is the
  **tropical** ``≈ 1.5`` and overshoots — :meth:`feedback_kernel`, the lapse-rate scope edge below.)
* **Reduction / plumbing.** Near present-day the emergent ``OLR(Ts)`` is **locally affine**, so the
  rung-0 ``A + B·T`` is its tangent line (:func:`linearized_olr`) — they share the operating point by
  calibration and rung-4 *derives* the slope rung-0 *prescribes*.

The named wall + scope edges
----------------------------
The **wall** is the **gray (band-independent) absorption** assumption together with the prescribed
**optical-depth ↔ greenhouse-gas** mapping (``τ_s`` and its water-vapour / CO₂ split are *calibrated*
to the present greenhouse, not derived from line-by-line spectroscopy — the cited-closure status of
``R_ATM_SLOPE`` at rung 2 and ``HADLEY_STRENGTH`` at the Hadley fix). Consequences, each named:

* **CO₂ forcing is saturating, not logarithmic — now FIXED by an opt-in spectral band.** A *gray*
  band gives a concave ``OLR(τ)`` (:meth:`GrayRadiationColumn.co2_forcing`): the right *sign and
  saturating shape*, locally doubling-like near present, but **not** the observed logarithmic law
  (Myhre ``5.35·ln C/C₀``) and at an **unrealistic magnitude** (a whole-band perturbation). The
  logarithmic law needs **spectral band wings** — built as :class:`SpectralCO2Band` (the within-rung
  slice): resolving the CO₂ band into bins with an **exponential** wing makes the per-doubling forcing
  **constant** (the Myhre log law) where the gray band saturates. The *functional form* is the win;
  the magnitude rides the band parameters (calibrated to order — the same wall). The gray
  :meth:`GrayRadiationColumn.co2_forcing` stays as the honest baseline (the band is a separate,
  opt-in construct).
* **Clear-sky only** (no clouds — the dominant real-world feedback uncertainty, out of scope).
* **Lapse-rate feedback — now BUILT, opt-in (and it OVERTURNED the scoped magnitude).** The default
  column's convective lapse rate is **fixed**, so warming is a uniform profile shift and there is *no*
  lapse-rate feedback — which is why the *default* gray net ``B ≈ 1.3`` sits below climlab's 2. Turning
  on a **moist adiabat** (``GrayRadiationColumn(moist_adiabat=True)``; :func:`moist_adiabat_temperature`)
  makes the feedback **emergent**: the adiabat flattens as it warms, so the surface warming amplifies in
  the upper troposphere and ``OLR(Ts)`` steepens. :meth:`feedback_kernel` measures it (Soden & Held
  kernel split) and the **sign and kind are banked** (``λ_LR > 0``; kernel closes to ~1e-3; resolution-
  converged). **But the §12 scoping guess "supplies ``λ_LR ≈ 0.84``, closing the gap to 2" was
  OVERTURNED:** the emergent value is **``≈ +1.5``**, so the moist-adiabat column **overshoots** — its
  with-water-vapour ``B ≈ 3.1`` sits *above* climlab's 2, not at it. **Reconciliation** (so the ``0.84``
  above does not read as a contradiction): the Soden & Held ``0.84`` is the **global-mean** lapse-rate
  feedback — a touchstone for what the *fixed*-``Γ`` default omits — whereas a single global moist-adiabat
  column captures only the **tropical** branch (the deep tropics are moist-adiabatic; the extratropics are
  not), and the tropical local feedback is ~``1.0–1.5`` *before* the model's own loadings push it higher.
  The magnitude is therefore **loose for two named reasons**: (a) the single column applies the tropical
  mechanism everywhere, missing the extratropical (bottom-heavy-warming) branch that pulls the global mean
  down to ``0.84``; and (b) it rides the prescribed vertical ``τ`` shape + :data:`WATER_VAPOUR_FRACTION`
  (the wall), which set where the emission level sits — the same loading the column's null is not perfectly
  clean about (fixed ``Γ`` already shows a small ``≈ −0.25`` tropopause-migration residual). Note the
  clean water-vapour/lapse-rate *separation* in :meth:`feedback_kernel` is partly a model artifact: here
  ``τ_wv`` tracks the **surface** ``Ts``, not the profile, so the upper-troposphere moisture–temperature
  coupling that links the two feedbacks in reality is absent. The remaining within-rung upgrades (the
  per-latitude wire :mod:`planet.radiative_ebm`, a moist-adiabat *with* latitudinal structure) are what
  would recover the extratropical branch and the global mean.
* **Single column.** The headline (emergent ``B``, water-vapour feedback, forcing) lives in the column.
  Wiring ``OLR(Ts, τ)`` *per latitude* into :mod:`planet.ebm` — real radiation *driving* the climate as
  an opt-in sibling EBM — is the natural rung-4 completion, **left to a user call** (not foreclosed). It
  is deferred not because the feedback is wrong (at the climlab-matched loading the global-mean ``B`` is 2
  exactly) but because the emergent ``OLR(Ts)`` is **nonlinear**: the per-latitude slope differs (cold
  pole vs warm equator), so the wire re-opens the meridional profile — a *feature* (emergent latitudinal
  radiative structure) as much as a risk.
* **Linearization breaks far from present** — a steep enough water-vapour loading drives a
  Komabayashi–Ingersoll runaway (the hot analogue of the snowball); ``B`` is not linearized across it.

Units — SI, with the climlab °C convention at the EBM seam
----------------------------------------------------------
Fluxes (``OLR``, ``F``) in **W m⁻²**; temperatures **in kelvin** inside the radiation (``σT⁴`` is a
kelvin law) and converted to **°C** only at :func:`linearized_olr` (the ``A + B·T`` seam, where
:mod:`planet.ebm` works in °C); optical depth ``τ`` dimensionless; pressure in **Pa**.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from planet.moist import (
    EPSILON,
    L_VAPOR,
    R_VAPOR,
    saturation_specific_humidity,
    saturation_vapor_pressure,
)

# --------------------------------------------------------------------------- #
# Pinned physical + calibration constants.
# Stefan–Boltzmann and the dry-air gas/gravity constants are CODATA/standard; the present-day
# operating point (Ts, OLR) is the observed global mean (Trenberth–Fasullo–Kiehl 2009: OLR ≈ 239
# W m⁻², Ts ≈ 288 K → Te ≈ 255 K, the 33 K greenhouse). The optical-depth split is the calibrated
# WALL — cited to *order*, not derived (Pierrehumbert PoPC §4; Goody & Yung *Atmospheric Radiation*).
# --------------------------------------------------------------------------- #
STEFAN_BOLTZMANN = 5.670374419e-8   # W m⁻² K⁻⁴ (CODATA)
PRESENT_SURFACE_T = 288.0           # K — observed global-mean surface temperature
PRESENT_OLR = 239.0                 # W m⁻² — observed global-mean outgoing longwave (= absorbed solar)
LAPSE_RATE = 6.5e-3                 # K m⁻¹ — convective-adjustment tropospheric lapse rate
STRATOSPHERE_T = 200.0              # K — cold-trap floor (no convective control above the tropopause)
P_SURFACE = 1.0e5                   # Pa — reference surface pressure
R_DRY_AIR = 287.0                   # J kg⁻¹ K⁻¹
C_P_AIR = 1004.0                    # J kg⁻¹ K⁻¹ — specific heat of dry air at constant pressure
GRAVITY = 9.81                      # m s⁻²
WATER_VAPOUR_FRACTION = 0.5         # — the WALL: fraction of the present optical depth that is
#                                     temperature-dependent water vapour (vs well-mixed CO₂). Sets the
#                                     water-vapour feedback magnitude; calibrated to order, not derived.

# Soden & Held (2006), "An Assessment of Climate Feedbacks in Coupled Ocean–Atmosphere Models",
# J. Climate 19:3354 — the multi-model-mean feedback parameters the emergent slopes are order-validated
# against (the [[…-source]] discipline, pinned at build, not carried from memory). Used as touchstones,
# NOT fitted to: the gray no-WV slope ≈ SH_PLANCK, the WV feedback ≈ SH_WATER_VAPOUR, and the gap from the
# gray net B to climlab's 2 ≈ SH_LAPSE_RATE (the feedback a fixed lapse rate cannot produce).
SH_PLANCK = 3.2                     # W m⁻² K⁻¹ — |λ₀| (their −3.1…−3.2): the Planck radiative damping
SH_WATER_VAPOUR = 1.8              # W m⁻² K⁻¹ — λ_wv: the clear-sky water-vapour feedback (largest positive)
SH_LAPSE_RATE = 0.84               # W m⁻² K⁻¹ — |λ_LR|: the GLOBAL-MEAN lapse-rate feedback gray's fixed Γ
#                                    OMITS. The emergent moist-adiabat column (feedback_kernel) OVERSHOOTS
#                                    it (~1.5, the tropical branch) — global-mean touchstone, not a target.


def emission_temperature(olr: float = PRESENT_OLR) -> float:
    """Emission temperature ``Tₑ = (OLR/σ)^¼`` (K) — the blackbody temperature of the planet to space."""
    return (olr / STEFAN_BOLTZMANN) ** 0.25


# --------------------------------------------------------------------------- #
# The tight anchor — the analytic gray radiative-equilibrium solution (derived in the module docstring).
# --------------------------------------------------------------------------- #
def gray_equilibrium_temperature(tau: np.ndarray | float, Te: float) -> np.ndarray:
    """Gray-RE air temperature ``T(τ) = [½Tₑ⁴(1+τ)]^¼`` (K) — the analytic profile.

    From ``σT⁴(τ) = ½σTₑ⁴(1+τ)``: the atmosphere warms downward from the skin value as the optical
    depth grows. ``τ`` is the flux optical depth measured from the top of atmosphere.
    """
    tau = np.asarray(tau, dtype=float)
    return (0.5 * Te ** 4 * (1.0 + tau)) ** 0.25


def skin_temperature(Te: float) -> float:
    """Skin temperature ``T(0) = Tₑ/2^¼`` (K) — the coldest, optically-thin top of the gray atmosphere."""
    return Te / 2.0 ** 0.25


def ground_temperature(tau_s: float, Te: float) -> float:
    """Gray-RE ground temperature ``Tg = [½Tₑ⁴(2+τ_s)]^¼`` (K) — carries the surface–air discontinuity.

    The ``2 + τ_s`` (vs the air's ``1 + τ_s`` just above) is the radiative-equilibrium surface jump: the
    ground is warmer than the air in contact with it. This coefficient is exactly where a recalled gray
    model goes wrong — derived in the module docstring, and confirmed by :func:`solve_gray_equilibrium`
    converging to it.
    """
    return (0.5 * Te ** 4 * (2.0 + tau_s)) ** 0.25


# --------------------------------------------------------------------------- #
# The moist adiabat — the temperature-dependent lapse rate that supplies the lapse-rate feedback.
# Derived by its LIMITS, not a recalled coefficient (the rung-3 K²=2F lesson): r_s→0 gives the dry
# adiabat g/c_p ≈ 9.8 K/km; a warm-moist column flattens toward ~4 K/km (the tropics). Wallace & Hobbs
# / Holton; the saturated (pseudo-)adiabatic lapse rate.
# --------------------------------------------------------------------------- #
def saturated_mixing_ratio(T_kelvin, p_pa):
    """Saturation mixing ratio ``r_s = ε·e_s/(p − e_s)`` (kg/kg) at the **local** pressure ``p``.

    Built from the Clausius–Clapeyron :func:`planet.moist.saturation_vapor_pressure` at the level's
    *own* pressure — not :func:`planet.moist.saturation_specific_humidity`, which bakes in a single
    reference pressure and would be wrong away from the surface.
    """
    e_s = saturation_vapor_pressure(np.asarray(T_kelvin, dtype=float) - 273.15)
    return EPSILON * e_s / np.maximum(p_pa - e_s, 1.0)


def moist_adiabatic_lapse_rate(T_kelvin, p_pa):
    """Saturated moist-adiabatic lapse rate ``Γ_m`` (K m⁻¹) at temperature ``T`` and pressure ``p``.

        Γ_m = g·(1 + L·r_s/(R_d·T)) / (c_p + L²·r_s/(R_v·T²))

    The dry limit ``r_s → 0`` recovers ``g/c_p ≈ 9.8 K/km``; latent-heat release in a warm, moist column
    (large ``r_s``) flattens it toward ``~4 K/km``. This temperature dependence — a flatter lapse rate
    when warmer — is the whole mechanism: warming the surface warms the upper troposphere *more*, which
    is the (negative) lapse-rate feedback a *fixed* ``Γ`` cannot produce.
    """
    r_s = saturated_mixing_ratio(T_kelvin, p_pa)
    num = 1.0 + L_VAPOR * r_s / (R_DRY_AIR * T_kelvin)
    den = C_P_AIR + L_VAPOR ** 2 * r_s / (R_VAPOR * T_kelvin ** 2)
    return GRAVITY * num / den


def moist_adiabat_temperature(Ts: float, p: np.ndarray, strat_T: float) -> np.ndarray:
    """Temperature profile ``T(p)`` (K) integrated up a moist adiabat from the surface, capped at ``strat_T``.

    ``p`` is TOA→surface (ascending pressure; index ``-1`` is the surface). Heights use the same fixed
    scale height as :meth:`GrayRadiationColumn._profile`, so ``moist_adiabat=False`` (a constant ``Γ``)
    and ``moist_adiabat=True`` (this) are the *same column* differing only in the lapse rate. Integrated
    upward with a predictor–corrector (Heun) step in height; converged in ``n_levels`` (see the spike).
    """
    scale_height = R_DRY_AIR * emission_temperature() / GRAVITY
    z = scale_height * np.log(P_SURFACE / p)            # z[-1] = 0 at the surface, rising toward the TOA
    T = np.empty_like(p, dtype=float)
    T[-1] = Ts
    for i in range(len(p) - 2, -1, -1):
        dz = z[i] - z[i + 1]                            # > 0 (going up)
        g1 = moist_adiabatic_lapse_rate(T[i + 1], p[i + 1])
        T_pred = T[i + 1] - g1 * dz
        g2 = moist_adiabatic_lapse_rate(max(T_pred, strat_T), p[i])
        T[i] = T[i + 1] - 0.5 * (g1 + g2) * dz
    return np.maximum(T, strat_T)


@dataclass(frozen=True)
class LapseRateFeedback:
    """Kernel split of the OLR slope ``B = dOLR/dTs`` into Planck + lapse-rate + water-vapour terms.

    The clean (Soden & Held) decomposition on **one** column: ``planck`` is the response to a *uniform*
    profile warming (τ fixed), ``lapse_rate`` is the response to the profile's *departure* from uniform
    warming (τ fixed; ``> 0`` and the headline for a moist adiabat, ``≈ 0`` for fixed ``Γ``), and
    ``water_vapour`` is the response to ``τ(Ts)`` alone (profile fixed; ``< 0``). They sum to ``total``
    to first order — :attr:`closure_residual` is the (small) second-order remainder.
    """

    total: float
    planck: float
    lapse_rate: float
    water_vapour: float

    @property
    def closure_residual(self) -> float:
        """``total − (planck + lapse_rate + water_vapour)`` — the second-order closure check (≈ 0)."""
        return self.total - (self.planck + self.lapse_rate + self.water_vapour)


def solve_gray_equilibrium(tau_s: float, Te: float, n_layers: int,
                           tol: float = 1e-13, max_iter: int = 50000):
    """Numerical two-stream radiative-equilibrium solve — the independent check on the analytic anchor.

    ``n_layers`` isothermal gray slabs of equal optical thickness ``Δτ = τ_s/n`` (emissivity
    ``ε = 1 − e^(−Δτ)``) are relaxed to radiative equilibrium: each layer's source is updated to the
    mean of the longwave entering its top and bottom (``F↓ + F↑ = 2σT⁴``), and the ground balances the
    absorbed solar ``σTₑ⁴`` plus the back-radiation. Uses *no* analytic input — it integrates the flux
    equations — so reproducing :func:`gray_equilibrium_temperature` / :func:`ground_temperature` is a
    genuine validation (it converges at ~2nd order in ``Δτ``; ``OLR = σTₑ⁴`` is machine-exact).

    Returns ``(tau_mid, T_layers, Tg, OLR)`` — layer-midpoint optical depths and temperatures (K), the
    ground temperature (K), and the top-of-atmosphere outgoing longwave (W m⁻²).
    """
    if n_layers < 1:
        raise ValueError(f"n_layers must be >= 1, got {n_layers}")
    dtau = tau_s / n_layers
    eps = 1.0 - math.exp(-dtau)
    trans = 1.0 - eps
    asr = STEFAN_BOLTZMANN * Te ** 4                 # absorbed solar = σTₑ⁴, all to the surface
    source = np.full(n_layers, 0.5 * asr)            # layer sources σT⁴, init at the skin value
    ground = asr
    f_down = np.zeros(n_layers)                      # F↓ at the bottom interface of each layer
    f_up = np.zeros(n_layers)                        # F↑ at the top interface of each layer
    for _ in range(max_iter):
        prev = 0.0                                   # F↓(top of atmosphere) = 0
        for i in range(n_layers):
            f_down[i] = trans * prev + eps * source[i]
            prev = f_down[i]
        prev = ground                                # F↑ entering the lowest layer from the surface
        for i in range(n_layers - 1, -1, -1):
            f_up[i] = trans * prev + eps * source[i]
            prev = f_up[i]
        down_into_top = np.concatenate(([0.0], f_down[:-1]))
        up_into_bottom = np.concatenate((f_up[1:], [ground]))
        new_source = 0.5 * (down_into_top + up_into_bottom)
        new_ground = asr + f_down[-1]
        delta = max(np.max(np.abs(new_source - source)), abs(new_ground - ground))
        source, ground = new_source, new_ground
        if delta < tol:
            break
    tau_mid = (np.arange(n_layers) + 0.5) * dtau
    T_layers = (source / STEFAN_BOLTZMANN) ** 0.25
    Tg = (ground / STEFAN_BOLTZMANN) ** 0.25
    return tau_mid, T_layers, Tg, float(f_up[0])


# --------------------------------------------------------------------------- #
# The RCE forward model — the emergent OLR(Ts, τ) the unlock is read off (a calibrated column).
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class GrayRadiationColumn:
    """A calibrated gray radiative–convective column: emergent ``OLR(Ts)``, its slope, and CO₂ forcing.

    The convective temperature profile is ``T(p) = max(Ts − Γ·z(p), T_strat)`` with ``z = H·ln(p_s/p)``;
    the optical depth from the top to pressure ``p`` is a **well-mixed CO₂** part ``∝ p/p_s`` plus a
    **bottom-heavy water-vapour** part ``∝ (p/p_s)²`` that scales with ``q_sat(Ts)`` (Clausius–Clapeyron)
    when ``water_vapour`` is on. ``total_tau`` is calibrated by :func:`calibrate_column` so
    ``OLR(PRESENT_SURFACE_T) = PRESENT_OLR``; ``wv_fraction`` splits it (the wall). ``co2_factor``
    multiplies only the CO₂ part (the forcing knob).

    ``moist_adiabat`` (default off → bit-for-bit the fixed-``Γ`` column above) swaps the constant ``Γ``
    for a temperature-dependent moist adiabat (:func:`moist_adiabat_temperature`). That is the **only**
    knob that lets the column produce a **lapse-rate feedback** (:meth:`feedback_kernel`): a moist adiabat
    flattens as it warms, so warming the surface warms the upper troposphere *more*, steepening
    ``OLR(Ts)``. Recalibrate ``total_tau`` for the new profile (``calibrate_column(moist_adiabat=True)``).
    """

    total_tau: float
    wv_fraction: float = WATER_VAPOUR_FRACTION
    lapse_rate: float = LAPSE_RATE
    strat_T: float = STRATOSPHERE_T
    n_levels: int = 200
    moist_adiabat: bool = False

    def _profile(self, Ts: float):
        """Pressure levels (TOA→surface, Pa) and the convective temperature profile ``T(p)`` (K).

        Fixed ``Γ`` by default; a moist adiabat when ``moist_adiabat`` is on (same heights, same column —
        only the lapse rate differs, so the off path is bit-for-bit the rung-4-core column).
        """
        p = np.linspace(0.02 * P_SURFACE, P_SURFACE, self.n_levels)
        if self.moist_adiabat:
            return p, moist_adiabat_temperature(Ts, p, self.strat_T)
        scale_height = R_DRY_AIR * emission_temperature() / GRAVITY
        z = scale_height * np.log(P_SURFACE / p)
        return p, np.maximum(Ts - self.lapse_rate * z, self.strat_T)

    def _optical_depth(self, Ts: float, co2_factor: float, water_vapour: bool) -> np.ndarray:
        """Cumulative flux optical depth from the TOA to each level (well-mixed CO₂ + C–C water vapour)."""
        p, _ = self._profile(Ts)
        tau_co2 = (1.0 - self.wv_fraction) * self.total_tau * co2_factor * (p / P_SURFACE)
        wv_scale = 1.0
        if water_vapour:
            q = saturation_specific_humidity(Ts - 273.15)
            q0 = saturation_specific_humidity(PRESENT_SURFACE_T - 273.15)
            wv_scale = q / q0
        tau_wv = self.wv_fraction * self.total_tau * wv_scale * (p / P_SURFACE) ** 2
        return tau_co2 + tau_wv

    def _olr_from(self, T: np.ndarray, tau: np.ndarray) -> float:
        """OLR (W m⁻²) from an *explicit* profile ``T(p)`` (K) and cumulative optical depth ``tau``.

        ``OLR = σT_s⁴ e^(−τ_s) + Σ_layers ε_k σT_k⁴ e^(−τ_top,k)``. Factored out of
        :meth:`outgoing_longwave` so :meth:`feedback_kernel` can evaluate the OLR for hand-built
        (profile, τ) pairs — a uniformly-warmed profile at fixed τ, etc. (``T[-1]`` is the surface, which
        equals ``Ts`` for both the fixed-``Γ`` and moist-adiabat profiles, so this is bit-for-bit the old
        body for the default column).
        """
        dtau = np.diff(tau)
        eps = 1.0 - np.exp(-dtau)
        T_layer = 0.5 * (T[:-1] + T[1:])
        trans_to_top = np.exp(-tau[:-1])
        surface = STEFAN_BOLTZMANN * float(T[-1]) ** 4 * math.exp(-float(tau[-1]))
        emitted = np.sum(eps * STEFAN_BOLTZMANN * T_layer ** 4 * trans_to_top)
        return surface + float(emitted)

    def outgoing_longwave(self, Ts: float, co2_factor: float = 1.0, water_vapour: bool = True) -> float:
        """Top-of-atmosphere outgoing longwave ``OLR`` (W m⁻²) by transmission-weighted emission.

        ``OLR = σTs⁴ e^(−τ_s) + Σ_layers ε_k σT_k⁴ e^(−τ_top,k)`` — the surface seen through the whole
        column plus each layer's emission attenuated by the optical depth above it. At ``Ts =
        PRESENT_SURFACE_T`` the water-vapour scaling is unity, so the calibrated present-day ``OLR`` is
        the same with or without ``water_vapour`` (the calibration is feedback-independent).
        """
        p, T = self._profile(Ts)
        tau = self._optical_depth(Ts, co2_factor, water_vapour)
        return self._olr_from(T, tau)

    def feedback_kernel(self, dT: float = 1.0) -> LapseRateFeedback:
        """Split ``B = dOLR/dTs`` into Planck + lapse-rate + water-vapour terms (Soden & Held kernels).

        On **one** (recalibrated) column, the warming response is decomposed by holding pieces fixed
        (advisor: cleaner than a two-column difference, which conflates the lapse-rate term with the
        shift in the Planck base):

        * **Planck** — warm the base profile *uniformly* (``T(p) → T(p) + dT``), ``τ`` frozen at present.
        * **Lapse rate** — the actual profile's *departure* from that uniform warming, ``τ`` frozen.
          ``> 0`` for a moist adiabat (upper-troposphere amplification ⇒ steeper ``OLR``), the headline
          feedback; ``≈ 0`` for fixed ``Γ`` (the null — only a small tropopause-migration residual).
        * **Water vapour** — the ``τ(Ts)`` change alone, base profile frozen. ``< 0``.

        They sum to the full slope to first order (:attr:`LapseRateFeedback.closure_residual` ≈ 0).
        """
        Ts0 = PRESENT_SURFACE_T
        _, T0 = self._profile(Ts0)
        tau0 = self._optical_depth(Ts0, 1.0, True)

        total = (self.outgoing_longwave(Ts0 + dT) - self.outgoing_longwave(Ts0 - dT)) / (2.0 * dT)
        planck = (self._olr_from(T0 + dT, tau0) - self._olr_from(T0 - dT, tau0)) / (2.0 * dT)

        def _olr_profile(Ts):                       # profile responds, τ frozen at present
            _, T = self._profile(Ts)
            return self._olr_from(T, tau0)
        profile = (_olr_profile(Ts0 + dT) - _olr_profile(Ts0 - dT)) / (2.0 * dT)
        lapse = profile - planck

        def _olr_tau(Ts):                           # τ responds, base profile frozen
            return self._olr_from(T0, self._optical_depth(Ts, 1.0, True))
        water_vapour = (_olr_tau(Ts0 + dT) - _olr_tau(Ts0 - dT)) / (2.0 * dT)

        return LapseRateFeedback(total=total, planck=planck, lapse_rate=lapse, water_vapour=water_vapour)

    def feedback_slope(self, water_vapour: bool = True, dT: float = 1.0, Ts: float = PRESENT_SURFACE_T) -> float:
        """Emergent OLR slope ``B = dOLR/dTs`` (W m⁻² K⁻¹) by central difference about ``Ts``.

        ``water_vapour=False`` holds the optical depth fixed → the Planck slope ``≈ 4σTₑ³`` (above
        climlab's 2); ``water_vapour=True`` lets ``τ`` rise with ``Ts`` → the slope drops through ~2
        (the positive feedback). The difference is the water-vapour feedback (:func:`feedback_decomposition`).
        """
        hi = self.outgoing_longwave(Ts + dT, water_vapour=water_vapour)
        lo = self.outgoing_longwave(Ts - dT, water_vapour=water_vapour)
        return (hi - lo) / (2.0 * dT)

    def feedback_decomposition(self):
        """``(B_no_wv, B_with_wv, water_vapour_feedback)`` — the rung's headline, ``B ≈ Planck − WV``.

        ``water_vapour_feedback = B_no_wv − B_with_wv > 0`` is the amount the water-vapour feedback
        subtracts from the bare Planck slope to land near climlab's prescribed ``B``. Direction banked;
        the magnitude rides on :data:`WATER_VAPOUR_FRACTION` (the wall).
        """
        b_dry = self.feedback_slope(water_vapour=False)
        b_moist = self.feedback_slope(water_vapour=True)
        return b_dry, b_moist, b_dry - b_moist

    def co2_forcing(self, co2_factor: float = 2.0, Ts: float = PRESENT_SURFACE_T) -> float:
        """Radiative forcing ``ΔF = OLR(1×) − OLR(co2_factor×)`` (W m⁻²) at fixed ``Ts`` (water vapour fixed).

        The instantaneous top-of-atmosphere forcing from scaling the CO₂ optical depth. **Gray and
        saturating**: ``ΔF`` per doubling rises then falls (concave ``OLR(τ)``), not the constant-per-
        doubling logarithmic law and at an unrealistic whole-band magnitude — the named band-physics edge.
        """
        base = self.outgoing_longwave(Ts, co2_factor=1.0, water_vapour=False)
        forced = self.outgoing_longwave(Ts, co2_factor=co2_factor, water_vapour=False)
        return base - forced

    def linearized_olr(self, dT: float = 5.0):
        """``(A, B)`` of the tangent ``OLR ≈ A + B·T`` (°C) at present — the reduction to rung-0's form.

        The emergent ``OLR(Ts)`` is locally affine near present, so rung-0's linear OLR is its tangent
        line: ``B`` is the present slope (with water vapour) and ``A = OLR(present) − B·T̄`` in the
        climlab **°C** convention (``T̄`` the present surface temperature in °C). Shares the operating
        point with climlab by calibration; the recovered ``B`` is the rung's finding, not ``2`` by
        construction. ``dT`` (K) is the half-width of the central difference.
        """
        B = self.feedback_slope(water_vapour=True, dT=dT)
        T_bar_c = PRESENT_SURFACE_T - 273.15
        A = self.outgoing_longwave(PRESENT_SURFACE_T) - B * T_bar_c
        return A, B


def calibrate_column(wv_fraction: float = WATER_VAPOUR_FRACTION, **kw) -> GrayRadiationColumn:
    """Build the present-day column: bisect ``total_tau`` so ``OLR(PRESENT_SURFACE_T) = PRESENT_OLR``.

    The single calibration that fixes the 33 K greenhouse; with it the emergent OLR passes through the
    observed present-day operating point by construction (so the *value* matches rung-0 and only the
    *slope* is open). Water vapour is unity at the reference ``Ts``, so the calibration is feedback-free.
    """
    lo, hi = 0.05, 40.0
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        col = GrayRadiationColumn(total_tau=mid, wv_fraction=wv_fraction, **kw)
        if col.outgoing_longwave(PRESENT_SURFACE_T, water_vapour=False) > PRESENT_OLR:
            lo = mid                                 # too transparent → more optical depth
        else:
            hi = mid
    return GrayRadiationColumn(total_tau=0.5 * (lo + hi), wv_fraction=wv_fraction, **kw)


# --------------------------------------------------------------------------- #
# The spectral-band log law — why the CO₂ forcing is LOGARITHMIC, not the gray SATURATING band.
# A within-rung upgrade: it replaces the gray column's band-independent absorption (one optical
# depth for the whole Planck spectrum) with a band-RESOLVED CO₂ absorption whose strength falls
# off EXPONENTIALLY in the wings, and shows the per-doubling forcing flatten to a constant.
# --------------------------------------------------------------------------- #
# Planck radiation constants (CODATA) — for the spectral Planck flux πB_ν(T), whose integral over
# wavenumber is σT⁴ (so summing per-bin emissions over the whole spectrum recovers the gray σT⁴,
# the reduction-to-gray anchor below).
PLANCK_H = 6.62607015e-34        # J s
SPEED_OF_LIGHT = 2.99792458e8    # m s⁻¹
BOLTZMANN_K = 1.380649e-23       # J K⁻¹

# The CO₂ 15-µm band, parameterized to ORDER (the wall — the same cited-closure status as the gray
# τ↔greenhouse mapping; these are not line-by-line spectroscopy). The band centre and the present
# band-centre optical depth set where the saturated core sits; the WING SCALE is the one ingredient
# that makes the forcing logarithmic (an exponential wing ⇒ the emission level spreads as l·ln C).
CO2_BAND_CENTRE_CM = 667.0       # cm⁻¹ — the ν₂ bending-mode band centre (15 µm)
CO2_BAND_HALF_WIDTH_CM = 217.0   # cm⁻¹ — modelled band half-extent (450–884 cm⁻¹); the far wing cutoff
CO2_BAND_WING_CM = 8.0           # cm⁻¹ — exponential wing decay scale l: k(ν) ∝ e^(−|ν−ν₀|/l). THE
#                                  ingredient; sets the forcing magnitude (~2l·π[B(Ts)−B_strat]).
CO2_BAND_CENTRE_TAU = 1000.0     # — band-centre surface optical depth at present CO₂ (deeply
#                                  saturated core, so the forcing comes only from the moving wings)

# Myhre et al. (1998) simplified expression ΔF = α·ln(C/C₀) — the observed logarithmic CO₂ forcing
# the emergent band law is order-validated against (a touchstone, NOT fitted to; the [[…-source]]
# discipline). α ≈ 5.35 W m⁻² ⇒ ≈ 3.71 W m⁻² per doubling.
MYHRE_COEFFICIENT = 5.35                          # W m⁻² — ΔF = MYHRE_COEFFICIENT·ln(C/C₀)
MYHRE_PER_DOUBLING = MYHRE_COEFFICIENT * math.log(2.0)   # ≈ 3.71 W m⁻² per CO₂ doubling


def planck_flux_per_wavenumber(wavenumber_per_m, T: float) -> np.ndarray:
    """Spectral flux ``πB_ν(T)`` per unit wavenumber (W m⁻² / m⁻¹); wavenumber ``ν`` in m⁻¹.

        πB_ν(T) = π · 2hc²ν³ / (e^{hcν/k_BT} − 1)

    Integrated over all wavenumbers this equals ``σT⁴`` (the Stefan–Boltzmann law) — so summing the
    per-bin emission ``πB_ν·Δν`` over a grid spanning the whole spectrum reproduces the gray ``σT⁴``
    source, which is exactly the reduction that ties :class:`SpectralCO2Band` back to the gray column.
    """
    n = np.asarray(wavenumber_per_m, dtype=float)
    x = PLANCK_H * SPEED_OF_LIGHT * n / (BOLTZMANN_K * T)
    return np.pi * 2.0 * PLANCK_H * SPEED_OF_LIGHT ** 2 * n ** 3 / np.expm1(x)


def _transmission_emission(tau: np.ndarray, surface_source: float, layer_source: np.ndarray) -> float:
    """One band's outgoing longwave: surface seen through the column + each layer's emission.

        OLR = S_surf·e^(−τ_s) + Σ_k ε_k·S_k·e^(−τ_top,k),   ε_k = 1 − e^(−Δτ_k)

    The same transmission-weighted-emission kernel as :meth:`GrayRadiationColumn._olr_from`, but with
    an explicit per-level *source* (``σT⁴`` for the gray whole-spectrum case, ``πB_ν·Δν`` for a
    spectral bin) rather than ``σT⁴`` hard-wired. Written independently of ``_olr_from`` — feeding it
    the gray ``σT⁴`` source reproduces ``_olr_from`` to machine precision, the cross-implementation
    check that anchors the band machinery to the gray column (:class:`SpectralCO2Band` reduction test).
    """
    dtau = np.diff(tau)
    eps = 1.0 - np.exp(-dtau)
    trans_to_top = np.exp(-tau[:-1])
    surface = surface_source * math.exp(-float(tau[-1]))
    emitted = np.sum(eps * layer_source * trans_to_top)
    return surface + float(emitted)


@dataclass(frozen=True)
class SpectralCO2Band:
    """Band-resolved CO₂ forcing: an exponential-wing absorption band over the column → the log law.

    The gray column treats CO₂ as a single optical depth for the *whole* Planck spectrum, so adding
    CO₂ pushes the *entire* emission to the cold upper atmosphere and the forcing **saturates** (a
    concave ``OLR(τ)``; per doubling ``ΔF`` 48→…→20 W m⁻², decreasing — :meth:`GrayRadiationColumn.
    co2_forcing`). Real CO₂ absorbs in a **band** whose strength falls off in the wings; the emission
    to space comes from the level where the band-resolved optical depth ``≈ 1``. This model resolves
    that band into ``n_bins`` spectral bins, each a gray sub-problem solved with the *same*
    transmission-weighted emission kernel (:func:`_transmission_emission`) over the column's fixed-``Γ``
    temperature profile and the well-mixed ``(p/p_s)`` CO₂ vertical shape, with the spectral Planck
    function ``πB_ν`` (:func:`planck_flux_per_wavenumber`) as the per-bin source.

    The band-centre optical depth is deeply saturated (``band_centre_tau ≫ 1``), so the forcing comes
    only from the **wings**. With an exponential wing ``k(ν) = k_c·e^(−|ν−ν₀|/l)`` the frequency at
    which ``τ(ν) = 1`` moves outward by ``l·Δ(ln C)`` per change in CO₂, so the spectral width that
    newly saturates is **constant per doubling** → ``ΔF`` is **constant per doubling** = the Myhre
    logarithmic law ``ΔF = α·ln(C/C₀)`` (:data:`MYHRE_COEFFICIENT`).

    **Scope — forcing only, CO₂ only.** This is the slice that fixes the *saturating-vs-logarithmic*
    edge, which lives entirely in the CO₂ forcing; the ``B = Planck − water-vapour + lapse-rate``
    slope decomposition is independent and stays on :class:`GrayRadiationColumn`. No water vapour here.

    **Triad.**

    * **Reduction (the independent anchor).** Collapsing the spectral resolution recovers the gray
      column: the band kernel with a single whole-spectrum bin and the ``σT⁴`` source reproduces
      :meth:`GrayRadiationColumn._olr_from` to machine precision (:func:`_transmission_emission`
      written independently of it), and a **uniform** ``k`` (no wings) makes :meth:`co2_forcing`
      *saturate* like gray — the exponential wing is the whole ingredient.
    * **The unlock (real but loose).** With exponential wings the per-doubling ``ΔF`` is **constant**
      across the realistic ``0.5×–8×`` range and lands in the **Myhre band** (:data:`MYHRE_PER_DOUBLING`
      ≈ 3.7 W m⁻²), versus gray's decreasing 48→20. *Loose:* the magnitude rides the band parameters
      (``l``, ``band_centre_tau``, the band half-width) — calibrated to **order**, not line-by-line
      (the **wall**); "CO₂ wings are ≈ exponential over the relevant range" is itself an empirical input.
    * **Derivation (consistency).** The cold-to-space (``τ = 1``) limit gives the slope
      ``dF/d ln C ≈ 2l·π[B_ν(Ts) − B_ν(T_strat)]`` (:meth:`log_law_coefficient`) — it matches the
      sharp emission-level estimate to ~1%; the column's finite-layer emission raises the realized
      coefficient ~20–30% above it. A derivation/consistency leg (both assume exponential wings), *not*
      an independent anchor.

    **The range is bounded (named edges).** The log law holds only in the middle regime — band centre
    saturated **and** wings not yet exhausted. Below ``C ≈ 1/band_centre_tau`` the band centre itself
    un-saturates and the forcing is linear/√ (grows per doubling, not constant); above the point where
    the active wing reaches the finite band edge the wings run out and it saturates again. Both edges
    sit far outside ``0.5×–8×`` for the present parameters.
    """

    column: GrayRadiationColumn = field(default_factory=lambda: GrayRadiationColumn(total_tau=4.0))
    band_centre_cm: float = CO2_BAND_CENTRE_CM
    half_width_cm: float = CO2_BAND_HALF_WIDTH_CM
    wing_scale_cm: float = CO2_BAND_WING_CM
    band_centre_tau: float = CO2_BAND_CENTRE_TAU
    n_bins: int = 300
    uniform: bool = False               # — flatten the wings (the null): the forcing then saturates

    def _grid(self):
        """``(centres_m, dnu_m, strengths)`` — bin centres (m⁻¹), widths (m⁻¹), and band-centre τ per bin.

        ``strengths[i] = band_centre_tau·e^(−|ν_i−ν₀|/l)`` (or flat if :attr:`uniform`) — the surface
        optical depth in bin ``i`` at present CO₂; the exponential wing is the log-law ingredient.
        """
        lo = self.band_centre_cm - self.half_width_cm
        hi = self.band_centre_cm + self.half_width_cm
        edges_cm = np.linspace(lo, hi, self.n_bins + 1)
        centres_cm = 0.5 * (edges_cm[:-1] + edges_cm[1:])
        dnu_m = np.diff(edges_cm) * 100.0                    # cm⁻¹ → m⁻¹
        if self.uniform:
            strengths = np.full_like(centres_cm, self.band_centre_tau)
        else:
            strengths = self.band_centre_tau * np.exp(-np.abs(centres_cm - self.band_centre_cm)
                                                      / self.wing_scale_cm)
        return centres_cm * 100.0, dnu_m, strengths

    def band_olr(self, Ts: float, co2_factor: float = 1.0) -> float:
        """Outgoing longwave from the CO₂ band alone (W m⁻²) at surface temperature ``Ts``.

        Sums each spectral bin's transmission-weighted emission over the column's fixed-``Γ`` profile,
        with ``τ_i(p) = co2_factor·strength_i·(p/p_s)`` and the spectral Planck source ``πB_ν·Δν``. The
        CO₂-transparent window (outside the band) is co₂-independent and omitted — it cancels in the
        forcing, the only quantity this model is built to deliver.
        """
        p, T = self.column._profile(Ts)
        shape = p / P_SURFACE
        T_layer = 0.5 * (T[:-1] + T[1:])
        centres_m, dnu_m, strengths = self._grid()
        olr = 0.0
        for nu, dn, k in zip(centres_m, dnu_m, strengths):
            tau = co2_factor * k * shape
            surface_source = planck_flux_per_wavenumber(nu, float(T[-1])) * dn
            layer_source = planck_flux_per_wavenumber(nu, T_layer) * dn
            olr += _transmission_emission(tau, surface_source, layer_source)
        return olr

    def co2_forcing(self, co2_factor: float = 2.0, Ts: float = PRESENT_SURFACE_T) -> float:
        """Radiative forcing ``ΔF = OLR_band(1×) − OLR_band(co2_factor×)`` (W m⁻²) at fixed ``Ts``.

        Same sign and convention as :meth:`GrayRadiationColumn.co2_forcing` (positive forcing = reduced
        OLR, fixed ``Ts``, no water-vapour feedback), so the two are an apples-to-apples comparison —
        gray saturates, this is logarithmic.
        """
        return self.band_olr(Ts, 1.0) - self.band_olr(Ts, co2_factor)

    def forcing_per_doubling(self, factors=(0.5, 1, 2, 4, 8), Ts: float = PRESENT_SURFACE_T):
        """``ΔF`` for each successive doubling across ``factors`` (W m⁻²) — the log-law signature.

        Returns the array of consecutive forcings ``OLR(f_{i})·... → F(f_{i+1}) − F(f_i)`` per
        ``log₂(f_{i+1}/f_i)``; **constant** (the Myhre log law) in the flat middle, versus gray's
        decreasing sequence. ``factors`` must be ascending and spaced by doublings for the labels to read.
        """
        F = np.array([self.co2_forcing(f, Ts=Ts) for f in factors])
        factors = np.asarray(factors, dtype=float)
        n_doublings = np.log2(factors[1:] / factors[:-1])
        return np.diff(F) / n_doublings

    def log_law_coefficient(self, Ts: float = PRESENT_SURFACE_T) -> float:
        """Analytic ``dF/d ln C ≈ 2l·π[B_ν(Ts) − B_ν(T_strat)]`` at band centre (W m⁻², the τ=1 limit).

        The cold-to-space derivation: an exponential wing's ``τ = 1`` level spreads by ``l·d ln C`` per
        wing, exposing a band of width ``2l·d ln C`` whose emission drops from the surface to the cold
        ``T_strat`` — so the forcing slope is ``2l·π[B_ν(Ts) − B_ν(T_strat)]``. Matches the sharp
        emission-level estimate to ~1%; the column's finite-layer emission realizes ~20–30% more. A
        consistency/derivation check (it assumes the same exponential wing), not an independent anchor.
        """
        nu0_m = self.band_centre_cm * 100.0
        contrast = (planck_flux_per_wavenumber(nu0_m, Ts)
                    - planck_flux_per_wavenumber(nu0_m, self.column.strat_T))
        return 2.0 * (self.wing_scale_cm * 100.0) * contrast
