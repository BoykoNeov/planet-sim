"""The column moist-EBM diagnostic — emergent moisture budget + the energy-constrained precip rate
(Planet rung 2, Phase A).

Rung 1 closed the EBM ⇄ circulation loop and wired the storm-track *position* to the emergent jet
(:mod:`planet.transport`, :mod:`planet.eddy_flux`, :mod:`planet.circ_precip`), but precipitation
stayed a **prescribed** kinematic field (:mod:`planet.precip`) with no water variable. Rung 2 adds the
first **moisture** physics — as a **column moist budget**, the fork the scoping spike settled (a
*resolved* storm-track precip pattern needs the vertical = ascent→condensation, which the single dry
layer lacks → that is **rung 3**, confirmed by ``outputs/rung2_moisture_convergence_spike.py``). This
module is that budget: a diagnostic moisture field ``q(φ)=RH·q_sat(T)`` over the rung-0 temperature,
a **down-gradient latent transport that reuses rung-1's eddy diffusivity** (no new free parameter),
and the precipitation diagnosed from the atmospheric water budget. It is a **fourth reuse of the
diffusion spine's structure** (the EBM's spherical ``(1−x²)`` transport operator, now applied to
moisture) and a **pure diagnostic** — it does **not** enter the temperature equation, so the Phase-1
climate and its triad stay green. See [[precip-parameterization-source]], [[ebm-radiation-source]],
[[planet-plan]] §10.

Two deliverables, deliberately split (the advisor's reframe)
------------------------------------------------------------
The bankable physics separates cleanly into two pieces with different status; fusing them into a
single emergent precipitation *field* would force an unphysical evaporation pattern (see "Why P−E,
not P" below), so they are kept apart:

1. **The RATE (the headline unlock — robust, opt-in).** The global-mean precipitation rate.
   :mod:`planet.precip` scales its whole pattern by a Clausius–Clapeyron **7 %/K** amplitude — the
   *moisture-capacity* rate, which it *already* flags in its own scope-edge as **not** the global
   precip rate. Global-mean precipitation is instead **energy-constrained**: the atmosphere can only
   shed so much latent heating, ``L⟨P⟩ ≈ R_atm − SH`` (atmospheric radiative cooling minus the
   sensible-heat flux), giving a much slower **≈ 2–3 %/K** (:func:`energy_constrained_rate`). This is
   the gap ``precip.py`` named. The energy-constrained amplitude (:func:`energy_constrained_factor`,
   :func:`energy_constrained_precip_field`) is **opt-in** and **= 1 at the present reference**, so the
   present-day map is unchanged and rung-0 ``precip.py`` (the 7 %/K) stays the default — exactly as
   circulation-informed precip stayed opt-in.

2. **The emergent moisture budget P−E (the trade — a diagnostic, not the default).** The atmospheric
   water budget in steady state is ``∂/∂t⟨water vapour⟩ = E − P − ∇·(moisture transport) = 0``, so

       P − E = −∇·(moisture transport) = the **moisture convergence**.

   With the down-gradient closure (below) this is a pure function of the rung-0 temperature
   (:func:`moisture_convergence`, :func:`moisture_budget`). It carries the *emergent* spatial content
   — and a named, honest finding (the "extratropical-only" trade below).

The latent transport reuses rung-1's eddy κ — and L cancels
-----------------------------------------------------------
The **same** midlatitude eddies that stir heat stir moisture, so the moisture transport reuses
rung-1's eddy diffusivity through the ``κ→D`` bridge (:mod:`planet.transport`) — **no new free
``D_q``**. The atmosphere diffuses the latent component of moist static energy ``L·q`` with the same
operator the EBM uses for sensible heat ``c_p·T`` (the diffusive moist EBM — Flannery 1984; Hwang &
Frierson 2010; Siler, Roe & Armour 2018). The κ→D bridge ``D = C_atm·κ/a²`` with the **atmospheric**
column heat capacity ``C_atm = c_p·p_s/g`` then gives, term-by-term, the **latent-energy** convergence

    L·(P − E) = (L / c_p)·D·∂/∂x[(1 − x²) ∂q/∂x]      (W m⁻²),

so the **moisture-mass** convergence (dividing by ``L``) is

    P − E = (D / c_p)·∂/∂x[(1 − x²) ∂q/∂x]            (kg m⁻² s⁻¹  →  cm/yr),

— the latent heat ``L`` **cancels**, leaving the EBM's *own* transport coefficient ``D`` and the
atmospheric ``c_p``: the moisture diffusivity is rung-1's eddy κ, reused, not a new knob. The
divergence is built in **conservative face-flux form** on the EBM ``x = sin φ`` grid with insulated
(Neumann 0) poles — the *identical* spherical operator the EBM transport reuses from
:mod:`engines.diffusion` — so ``∫(P − E) dx = 0`` to **machine precision** (the global ``∫E = ∫P``
plumbing leg falls straight out of moisture-mass conservation).

The energy-constrained rate and its sub-grid wall (R_atm)
---------------------------------------------------------
The rate ``L⟨P⟩ = R_atm − SH`` is closed by the **atmospheric radiative-cooling** response
``R_atm(T̄)`` — the **first sub-grid closure of the staircase, the named rung-2 wall**. With no
vertically-resolved radiation here, ``dR_atm/dT̄`` is **prescribed** (``R_ATM_SLOPE ≈ 2 W m⁻² K⁻¹``,
Held & Soden 2006; Allen & Ingram 2002) — **not** derived, and explicitly **not** ``B_OLR`` (also
2 W m⁻² K⁻¹ but a *different* quantity: the top-of-atmosphere longwave feedback, not the
atmospheric-column cooling sensitivity). Normalising by the present global precipitation
``⟨P⟩₀ ≈ 100 cm/yr`` gives the fractional rate ``R_ATM_SLOPE / (L⟨P⟩₀) ≈ 2.5 %/K``. The energy budget
is **linear** in T̄ (``R_atm`` linear), so the amplitude factor is **linear**, not the exponential of
the C–C moisture-capacity form — the honest functional difference, not merely a smaller exponent. The
rate's honesty lives entirely in the cited closure ``R_ATM_SLOPE``.

Why P−E, not a full emergent P (the advisor's catch)
----------------------------------------------------
A full emergent ``P`` would need an evaporation field ``E(φ)`` (``P = E + (P − E)``), and **no honest
zonal ``E`` keeps ``P ≥ 0``**: down-gradient diffusion *exports* ~2 m/yr from the moist equator (the
ITCZ-backwards term below), so a uniform ``E`` drives ``P`` negative there, and an ``E ∝ q_sat``
"rescue" only does so with an absurd ~6 m/yr of equatorial evaporation (real evaporation peaks in the
subtropics, not the equator). So this module reports **P − E** (the robust, sign-meaningful moisture
convergence) and the **global rate** ``⟨P⟩`` separately — it does **not** synthesise a full ``P``
field. Rung-0 ``precip.py`` remains the precipitation *map*.

Validation triad (plan §3) — re-classed for honesty
----------------------------------------------------
* **Tight (analytical/structural).** ``q_sat`` is the **exact Clausius–Clapeyron function** (the
  Whittaker-partition precedent: an exact, testable function, not a fit) — monotone, ~7 %/K local
  slope. The conservative moisture operator reproduces the **P₂ Legendre eigenvalue**
  (``∂/∂x[(1−x²)∂P₂/∂x] = −6 P₂``) to grid order — i.e. it *is* the EBM transport operator (the same
  anchor :mod:`planet.transport` uses for the channel geometry).
* **Real-but-loose physics (the unlock).** The energy-constrained ``≈ 2–3 %/K`` global rate — a
  *cited-closure* result (it follows from the prescribed ``R_ATM_SLOPE``, not from first principles),
  asserted in a loose band and contrasted with the faster C–C 7 %/K.
* **Consistency / plumbing (named as such).** Global ``∫(P − E) dx = 0`` (machine-exact, from the
  conservative flux form — moisture-mass conservation ⟹ ``∫E = ∫P``); the reduction to rung 0 as
  ``q → 0`` (a vanishing moisture layer ⟹ ``P − E → 0`` — **by construction**, the
  :func:`planet.transport.two_way_pass` plumbing honesty class, not an independent test); and
  ``energy_constrained_factor = 1`` at the present reference (rung-0 amplitude recovered).
* **Benchmark (loose) — the named "extratropical-only" trade.** Down-gradient moisture diffusion gets
  the **extratropics** right (midlatitude and polar ``P > E`` convergence) but gets the **deep tropics
  backwards** (it *exports* moisture from the moist equator — the real ITCZ is *up-gradient* Hadley
  convergence, a mean-circulation feature out of scope) **and** mislocates the subtropical evaporative
  belt (the steep equator–pole contrast hyper-peaks C–C ``q`` at the equator, pushing the
  moisture-flux maximum equatorward of the canonical subtropics, so the subtropics come out as weak
  ``P > E`` rather than the observed evaporative ``E > P``). So Phase A banks the moist **energetics +
  the extratropical budget + the global rate**, *not* a wholesale-better precip map — the same
  "trade, not a win" the staircase keeps finding. The benchmark therefore asserts **only** the
  equatorial export and the extratropical convergence (never subtropical ``E > P``).

Named scope edges
-----------------
* **ITCZ / Hadley deferred** — the up-gradient tropical moisture convergence is mean circulation, not
  eddy diffusion (a rung-3+ feature). The prescribed ``precip.py`` ITCZ is better in the deep tropics.
* **No resolved storm-track precip pattern** — that needs the vertical (ascent→condensation) = rung 3
  (spike-confirmed); rung 2 banks the *column* budget, not the *resolved* pattern.
* **Fixed RH; ``R_atm`` prescribed** — the sub-grid wall. A fuller moist EBM that diffuses *moist
  static energy* ``m = c_pT + L·q`` so ``T`` itself responds (emergent polar amplification) is the
  named **rung 2.5** extension — deferred because it re-opens the Phase-1 ``(A, B, D)`` calibration
  (rung-0's ``D = 0.555`` is an *effective* diffusivity already absorbing latent transport; explicit
  MSE diffusion double-counts the latent heat implicit in the linear OLR ``A + B·T``).
* **Deep-polar ``q → 0`` noise.** Poleward of ~79° (where ``q`` and the area weight ``1 − x²`` both
  approach 0) the convergence is a tiny, sign-flipping residual — physically negligible and not part of
  the banked structure; the benchmark asserts the extratropical convergence as a *mean poleward of 40°*,
  not pointwise into the polar cap.

Units — SI internally; P, P−E reported in **cm/yr** (the Whittaker / ``precip.py`` axis)
----------------------------------------------------------------------------------------
``q`` is specific humidity (kg/kg, dimensionless); ``q_sat`` from C–C; ``T`` in °C (the climlab/EBM
convention; converted to K inside ``q_sat``); fluxes in W m⁻² and mass fluxes in kg m⁻² s⁻¹ internally,
reported as **cm of water per year**; ``x = sin φ`` on [0, 1]; latitudes in degrees.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import precip
from .ebm import ClimateState, D_TRANSPORT, legendre_P2
from .transport import CP_AIR

# --------------------------------------------------------------------------- #
# Pinned moisture / thermodynamic constants.
# Clausius–Clapeyron and the saturation reference: Hartmann, *Global Physical Climatology* §1.5;
# Bohren & Albrecht, *Atmospheric Thermodynamics*. The energy-constrained precip rate: Held & Soden
# 2006; Allen & Ingram 2002 (extends [[precip-parameterization-source]], which already names the gap).
# The diffusive moist-EBM (fixed RH, latent diffusion): Flannery 1984; Hwang & Frierson 2010;
# Siler, Roe & Armour 2018. Cited and pinned at build — NOT carried from memory.
# --------------------------------------------------------------------------- #
L_VAPOR = 2.5e6            # J kg⁻¹    — latent heat of vaporization (~0–15 °C; Hartmann)
R_VAPOR = 461.5           # J kg⁻¹ K⁻¹ — gas constant for water vapour
E_SAT_0 = 611.2           # Pa        — saturation vapour pressure at T0 (0 °C)
T0_KELVIN = 273.15        # K         — Celsius/Kelvin offset and the C–C reference temperature
EPSILON = 0.622           # —         — R_dry / R_vapor (molar-mass ratio); q_sat = ε·e_sat/p
P_SURFACE = 1.013e5       # Pa        — global-mean surface pressure (== transport.P_SURFACE)
RH_DEFAULT = 0.8          # —         — fixed relative humidity (observed near-surface ~0.7–0.8)

# The energy-constrained-rate closure — the rung-2 sub-grid WALL (named, prescribed, NOT derived).
R_ATM_SLOPE = 2.0         # W m⁻² K⁻¹ — atmospheric radiative-cooling sensitivity dR_atm/dT̄ (Held &
                          #             Soden 2006; Allen & Ingram 2002). NOT B_OLR (= 2 W m⁻² K⁻¹ too,
                          #             but the *TOA* longwave feedback — a different quantity).
GLOBAL_PRECIP_REF_CMYR = 100.0   # cm/yr — present global-mean precipitation ⟨P⟩₀ (~1 m/yr; Hartmann)

SECONDS_PER_YEAR = 3.1557e7
# 1 kg m⁻² s⁻¹ of water = 1 mm s⁻¹ depth (ρ_w = 10³ kg m⁻³); cm/yr = mm/s · (s/yr) / (10 mm/cm).
_KGM2S_TO_CMYR = SECONDS_PER_YEAR / 10.0


# --------------------------------------------------------------------------- #
# Clausius–Clapeyron saturation humidity — the TIGHT leg (an exact, testable function).
# --------------------------------------------------------------------------- #
def saturation_vapor_pressure(T_celsius):
    """Saturation vapour pressure ``e_sat(T)`` (Pa) — the integrated Clausius–Clapeyron relation.

    ``e_sat(T) = e₀·exp[(L_v/R_v)(1/T₀ − 1/T)]`` with ``e₀ = 611.2 Pa`` at ``T₀ = 273.15 K`` (the
    constant-``L_v`` integration of ``d ln e_sat/dT = L_v/(R_v T²)``). ``T`` in °C (converted to K).
    The exact thermodynamic function moisture rides on — not a fit (the Whittaker-partition precedent).
    """
    T = np.asarray(T_celsius, dtype=float) + T0_KELVIN
    return E_SAT_0 * np.exp((L_VAPOR / R_VAPOR) * (1.0 / T0_KELVIN - 1.0 / T))


def saturation_specific_humidity(T_celsius):
    """Saturation specific humidity ``q_sat(T)`` (kg/kg) — ``ε·e_sat / (p − (1−ε)·e_sat)``.

    The mass of water vapour per unit mass of moist air at saturation, from the C–C
    :func:`saturation_vapor_pressure` and the surface pressure. ~7 %/K locally near present-day T
    (the moisture-capacity rate :mod:`planet.precip` cites) — verified in the triad as the exact
    function's log-slope, distinct from the energy-constrained *precip* rate below.
    """
    e_sat = saturation_vapor_pressure(T_celsius)
    return EPSILON * e_sat / (P_SURFACE - (1.0 - EPSILON) * e_sat)


def specific_humidity(T_celsius, RH: float = RH_DEFAULT):
    """Diagnostic specific humidity ``q(φ) = RH·q_sat(T)`` (kg/kg) — fixed relative humidity.

    Moisture is slaved to the rung-0 temperature at a **fixed** relative humidity (the diffusive
    moist-EBM closure — Hwang & Frierson 2010): warmer air holds proportionally more vapour. ``RH``
    is dimensionless on (0, 1]; the default is the observed near-surface ~0.8.
    """
    return float(RH) * saturation_specific_humidity(T_celsius)


# --------------------------------------------------------------------------- #
# The energy-constrained global precipitation rate — the headline unlock (opt-in).
# --------------------------------------------------------------------------- #
def energy_constrained_rate(slope: float = R_ATM_SLOPE,
                            mean_precip_cmyr: float = GLOBAL_PRECIP_REF_CMYR) -> float:
    """The energy-constrained fractional precipitation rate ``d ln⟨P⟩/dT̄`` (per K) — ≈ 0.025.

    From ``L⟨P⟩ = R_atm − SH`` with a prescribed atmospheric-cooling sensitivity ``slope = dR_atm/dT̄``:
    the fractional rate is ``slope / (L·⟨P⟩₀)``. With ``slope ≈ 2 W m⁻² K⁻¹`` (Held & Soden 2006) and
    ``⟨P⟩₀ ≈ 100 cm/yr`` this is ≈ **2.5 %/K** — about a third of the Clausius–Clapeyron 7 %/K
    moisture-capacity rate (the atmosphere's *energy* budget, not its moisture *capacity*, limits the
    global hydrological cycle). A **cited-closure** result: it inherits the honesty of the prescribed
    ``slope`` (the named sub-grid wall), it is **not** derived from first principles here.
    """
    LP0 = L_VAPOR * (float(mean_precip_cmyr) / _KGM2S_TO_CMYR)        # present L⟨P⟩ (W m⁻²)
    return float(slope) / LP0


def energy_constrained_factor(global_mean_T: float, ref_T: float = precip.PRECIP_REF_TEMP_C,
                              slope: float = R_ATM_SLOPE,
                              mean_precip_cmyr: float = GLOBAL_PRECIP_REF_CMYR) -> float:
    """Energy-constrained global moisture-amplitude factor ``= 1 + rate·(T̄ − T_ref)`` (floored at 0).

    The opt-in replacement for :func:`planet.precip.clausius_clapeyron_factor` (the C–C 7 %/K). The
    energy budget is **linear** in T̄, so this is a **linear** factor (not the C–C exponential) with
    the energy-constrained :func:`energy_constrained_rate` slope. ``= 1`` exactly at the present
    reference (so the present-day map is unchanged; rung-0 ``precip.py`` stays the default), and
    floored at 0 to stay physical in a deep-Snowball cooling. The headline rung-2 unlock.
    """
    rate = energy_constrained_rate(slope, mean_precip_cmyr)
    return float(max(0.0, 1.0 + rate * (float(global_mean_T) - ref_T)))


def energy_constrained_precip_field(state: ClimateState, RH: float = RH_DEFAULT,
                                    slope: float = R_ATM_SLOPE) -> np.ndarray:
    """Opt-in precipitation ``P(φ)`` (cm/yr): the rung-0 *pattern* × the **energy-constrained** amplitude.

    Like :func:`planet.precip.precip_field`, but the global moisture amplitude is the
    energy-constrained :func:`energy_constrained_factor` (≈ 2–3 %/K) instead of the C–C 7 %/K. Keeps
    the circulation-set pattern (so ``P ≥ 0`` is free — see the module docstring on why a *new*
    emergent ``P`` field is not built) and reduces to ``precip.precip_field`` **bit-for-bit** at the
    present reference temperature. The opt-in rung-2 amplitude; rung-0 ``precip.py`` remains the default.
    ``RH`` is unused here (the pattern is prescribed) and kept only for signature parallelism with
    :func:`moisture_budget`.
    """
    pattern = precip.precip_pattern(state.latitude_deg())
    return pattern * energy_constrained_factor(state.global_mean_T, slope=slope)


# --------------------------------------------------------------------------- #
# The emergent moisture budget P − E — the diagnostic (the trade).
# --------------------------------------------------------------------------- #
def _spherical_flux_divergence(field: np.ndarray, x: np.ndarray) -> np.ndarray:
    """The EBM spherical transport operator ``∂/∂x[(1 − x²) ∂field/∂x]`` in conservative flux form.

    Built as differences of interior **face** fluxes ``(1 − x_face²)·Δfield/Δx`` on the uniform
    ``x = sin φ`` grid, with **insulated (Neumann 0)** ends — the *identical* operator the EBM transport
    reuses from :mod:`engines.diffusion` (so this is the diffusion spine's *structure*, a fourth time).
    Because the boundary fluxes are zero (equator by symmetry, pole where ``1 − x² → 0``), the discrete
    operator sums to **exactly zero** over the grid (``Σ div·Δx = 0``) — the machine-exact ``∫E = ∫P``
    plumbing leg. Reproduces the Legendre eigenvalue ``∂/∂x[(1−x²)∂P₂/∂x] = −6 P₂`` to grid order (the
    tight anchor that it *is* the EBM operator).
    """
    x = np.asarray(x, dtype=float)
    f = np.asarray(field, dtype=float)
    dx = x[1] - x[0]                                   # uniform EBM grid (uniform_grid(1.0, n))
    x_face = 0.5 * (x[:-1] + x[1:])                    # interior cell faces
    g = 1.0 - x_face ** 2                              # area weight (1 − x²) at faces (→ 0 at the pole)
    flux = g * (f[1:] - f[:-1]) / dx                   # interior-face down-gradient flux × (1 − x²)
    div = np.zeros_like(f)
    div[1:-1] = (flux[1:] - flux[:-1]) / dx            # interior cells: difference of bounding faces
    div[0] = flux[0] / dx                              # equator: lower face is Neumann 0 (symmetry)
    div[-1] = -flux[-1] / dx                           # pole: upper face is Neumann 0 (no flux)
    return div


def moisture_convergence(state: ClimateState, RH: float = RH_DEFAULT) -> np.ndarray:
    """The moisture convergence ``P − E`` (cm/yr) — down-gradient latent transport off the rung-0 climate.

    ``P − E = (D / c_p)·∂/∂x[(1 − x²) ∂q/∂x]`` with ``q = RH·q_sat(T)`` (see the module docstring for
    the L-cancellation: the latent heat drops out, leaving the EBM's *own* transport coefficient ``D``
    and the atmospheric ``c_p`` — rung-1's eddy diffusivity reused, no new ``D_q``). Positive where the
    column **gains** moisture (``P > E`` — midlatitude/polar convergence), negative where it **loses**
    (``E > P``). Built in conservative flux form, so ``∫(P − E) dx = 0`` to machine precision.

    **The honest reading (the trade):** the **extratropics** are right (poleward ``P > E``), but the
    **deep equator is backwards** (down-gradient diffusion *exports* moisture from the moist equator;
    the real ITCZ is up-gradient Hadley convergence) and the subtropical evaporative belt is mislocated
    — so this is an emergent *budget*, not a better precip *map* (see the module triad).
    """
    q = specific_humidity(state.T, RH)
    conv_mass = (D_TRANSPORT / CP_AIR) * _spherical_flux_divergence(q, state.x)   # kg m⁻² s⁻¹
    return conv_mass * _KGM2S_TO_CMYR                                             # cm/yr


@dataclass(frozen=True)
class MoistureBudget:
    """The rung-2 column moist-budget diagnostic for one rung-0 climate (plain arrays — loose coupling).

    ``phi`` latitudes (deg); ``q`` the diagnostic specific humidity (kg/kg, ``RH·q_sat(T)``);
    ``p_minus_e`` the moisture convergence ``P − E`` (cm/yr, :func:`moisture_convergence`);
    ``mean_precip`` the energy-constrained global-mean ``⟨P⟩`` (cm/yr, :func:`energy_constrained_factor`
    × the present reference); ``energy_rate`` the fractional ``d ln⟨P⟩/dT̄`` (per K, the unlock);
    ``net_p_minus_e`` the area integral ``∫(P − E) dx`` (cm/yr, ~0 to machine precision — the
    ``∫E = ∫P`` plumbing leg). ``equatorial_export`` (``P − E`` at the equator, < 0 — the
    ITCZ-backwards term) and ``extratropical_convergence`` (mean ``P − E`` poleward of 40°, > 0) are
    the two banked benchmark numbers (the named extratropical-only trade).
    """

    phi: np.ndarray
    q: np.ndarray
    p_minus_e: np.ndarray
    mean_precip: float
    energy_rate: float
    net_p_minus_e: float
    equatorial_export: float
    extratropical_convergence: float


def moisture_budget(state: ClimateState, RH: float = RH_DEFAULT,
                    slope: float = R_ATM_SLOPE) -> MoistureBudget:
    """Build the :class:`MoistureBudget` — the moisture field, the ``P − E`` convergence, and the rate.

    Reads the rung-0 :class:`~planet.ebm.ClimateState` (its temperature and grid) and returns the full
    diagnostic: the fixed-RH moisture field, the down-gradient ``P − E`` convergence
    (:func:`moisture_convergence`), the energy-constrained global mean and rate, and the two banked
    benchmark numbers (equatorial export, extratropical convergence). A **pure diagnostic** — it does
    not modify ``state`` or the climate, so the Phase-1 triad is untouched.
    """
    phi = state.latitude_deg()
    q = specific_humidity(state.T, RH)
    pme = moisture_convergence(state, RH)
    # Area integral ∫(P − E) dx on the equal-area (uniform-Δx) grid is the area MEAN — the EBM's own
    # `total` convention (rectangle rule), under which the conservative flux divergence sums to *exactly*
    # zero (Σ div = 0). np.trapezoid would break the telescoping and leave an O(boundary) quadrature
    # residual, so the machine-exact ∫E = ∫P plumbing leg uses the mean.
    net = float(np.mean(pme))
    extra_mask = phi >= 40.0
    return MoistureBudget(
        phi=phi, q=q, p_minus_e=pme,
        mean_precip=float(GLOBAL_PRECIP_REF_CMYR * energy_constrained_factor(state.global_mean_T, slope=slope)),
        energy_rate=energy_constrained_rate(slope),
        net_p_minus_e=net,
        equatorial_export=float(pme[0]),
        extratropical_convergence=float(np.mean(pme[extra_mask])),
    )
