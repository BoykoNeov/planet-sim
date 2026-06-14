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
frozen-``D_eff`` and ITCZ closed-form attribution flavour, stronger than a bare "trade").

Triad (plan §3) — what is asserted tight vs loose
-------------------------------------------------
* **Analytical limit (tight).** The numerical two-stream solver reproduces the derived gray-RE profile,
  skin and ground temperatures to ~2nd order in layer thickness, with ``OLR = σTₑ⁴`` machine-exact
  (energy conservation); and the no-feedback slope sits near the ``4σTₑ³`` Planck touchstone.
* **The unlock (real but loose).** The emergent ``B`` decomposition (Planck slope ≈ 3.4, a water-vapour
  feedback ≈ 2 dropping the net to ≈ 1.3) — **order-validated** against Soden & Held (2006) (Planck 3.2,
  ``λ_wv`` 1.8, the residual-to-climlab ≈ the lapse-rate ``λ_LR`` 0.84), **not tuned**. Loose: the exact
  magnitudes ride on the water-vapour optical-depth loading (:data:`WATER_VAPOUR_FRACTION`), the **wall**.
* **Reduction / plumbing.** Near present-day the emergent ``OLR(Ts)`` is **locally affine**, so the
  rung-0 ``A + B·T`` is its tangent line (:func:`linearized_olr`) — they share the operating point by
  calibration and rung-4 *derives* the slope rung-0 *prescribes*.

The named wall + scope edges
----------------------------
The **wall** is the **gray (band-independent) absorption** assumption together with the prescribed
**optical-depth ↔ greenhouse-gas** mapping (``τ_s`` and its water-vapour / CO₂ split are *calibrated*
to the present greenhouse, not derived from line-by-line spectroscopy — the cited-closure status of
``R_ATM_SLOPE`` at rung 2 and ``HADLEY_STRENGTH`` at the Hadley fix). Consequences, each named:

* **CO₂ forcing is saturating, not logarithmic.** A *gray* band gives a concave ``OLR(τ)``
  (:func:`co2_forcing`): the right *sign and saturating shape*, locally doubling-like near present, but
  **not** the observed logarithmic law (Myhre ``5.35·ln C/C₀``) and at an **unrealistic magnitude**
  (a whole-band perturbation) — the logarithmic law needs **spectral band wings**, the named within-rung
  band upgrade.
* **Clear-sky only** (no clouds — the dominant real-world feedback uncertainty, out of scope).
* **No lapse-rate feedback.** The convective lapse rate is **fixed**, so warming is a uniform profile
  shift — which produces *no* lapse-rate feedback. That is exactly why the gray net ``B ≈ 1.3`` sits
  below climlab's 2 by ``≈ SH_LAPSE_RATE``; a temperature-dependent (e.g. moist-adiabatic) lapse rate is
  the named within-rung upgrade that would supply it.
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
from dataclasses import dataclass

import numpy as np

from planet.moist import saturation_specific_humidity

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
SH_LAPSE_RATE = 0.84               # W m⁻² K⁻¹ — |λ_LR|: the lapse-rate feedback gray's fixed Γ OMITS


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
    """

    total_tau: float
    wv_fraction: float = WATER_VAPOUR_FRACTION
    lapse_rate: float = LAPSE_RATE
    strat_T: float = STRATOSPHERE_T
    n_levels: int = 200

    def _profile(self, Ts: float):
        """Pressure levels (TOA→surface, Pa) and the convective temperature profile ``T(p)`` (K)."""
        p = np.linspace(0.02 * P_SURFACE, P_SURFACE, self.n_levels)
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

    def outgoing_longwave(self, Ts: float, co2_factor: float = 1.0, water_vapour: bool = True) -> float:
        """Top-of-atmosphere outgoing longwave ``OLR`` (W m⁻²) by transmission-weighted emission.

        ``OLR = σTs⁴ e^(−τ_s) + Σ_layers ε_k σT_k⁴ e^(−τ_top,k)`` — the surface seen through the whole
        column plus each layer's emission attenuated by the optical depth above it. At ``Ts =
        PRESENT_SURFACE_T`` the water-vapour scaling is unity, so the calibrated present-day ``OLR`` is
        the same with or without ``water_vapour`` (the calibration is feedback-independent).
        """
        p, T = self._profile(Ts)
        tau = self._optical_depth(Ts, co2_factor, water_vapour)
        dtau = np.diff(tau)
        eps = 1.0 - np.exp(-dtau)
        T_layer = 0.5 * (T[:-1] + T[1:])
        trans_to_top = np.exp(-tau[:-1])
        surface = STEFAN_BOLTZMANN * Ts ** 4 * math.exp(-float(tau[-1]))
        emitted = np.sum(eps * STEFAN_BOLTZMANN * T_layer ** 4 * trans_to_top)
        return surface + float(emitted)

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
