"""Ice-albedo feedback and the Snowball-Earth hysteresis (Planet Phase 1, the banked payoff).

This module supplies the **nonlinearity** that turns the linear energy-balance machinery of
:mod:`planet.ebm` into the project's dramatic early win. The albedo jumps from its
ice-free value to a high ice value wherever the surface freezes (``T < Tf``); fed through the
EBM as the absorbed-shortwave forcing, that threshold makes the local radiation step a
*nonlinear* relaxation — which gives the system **two stable climates for one sun** over a
range of solar constants. Ramp the sun down and a temperate planet's ice line creeps
equatorward, then past a threshold the feedback runs away and the whole planet **freezes in a
catastrophic jump** (Snowball Earth); brighten it back up and it **stays frozen far past where
it froze**, because a white planet reflects too much to thaw — a wide **hysteresis loop**. One
knob (S₀), a bifurcation, the cheapest most teachable payoff in the capstone.

The hysteresis is a **parameter-continuation sweep**, not a single solve (plan, the crux): at
each S₀ we relax to equilibrium *warm-started from the previous S₀'s equilibrium*, so the
solver tracks whichever branch it is already on — exactly how the two branches and the jumps
between them are traced out.

The albedo model (climlab `StepFunctionAlbedo`, [[ebm-radiation-source]])
-------------------------------------------------------------------------
``α(x, T) = a₀ + a₂·P₂(x)``  where the surface is unfrozen (a smoothly poleward-brightening
ice-free albedo), and ``α = a_ice`` wherever ``T < Tf``. The absorbed shortwave the EBM
integrates is then ``S(x)·(1 − α(x, T))`` — the **one place** the ice nonlinearity enters; the
transport and OLR stay linear. Constants are the climlab defaults
(``a₀=0.30, a₂=0.078, a_ice=0.62, Tf=−10 °C``), pinned in :mod:`~planet.ebm`.

Validation (plan §3) — what is asserted tight vs loose
------------------------------------------------------
*Validated tight (structural/qualitative):* the **existence and structure of the hysteresis**
— two stable branches over a band of S₀, a catastrophic (discontinuous) freezing jump, and a
re-melt at a *higher* S₀ than the freeze (the loop has positive width). Nothing but the
feedback produces this; it is emergent, so it is asserted firmly. *Calibrated/flagged (loose):*
the **exact threshold S₀ values and the loop width** depend on the cited albedo/radiation
constants, so they are asserted only in loose bands (and benchmarked against climlab in
:mod:`~planet.climate_reference`). The present-day ice line (~70°) is the loose
benchmark this module's :func:`present_day_climate` reproduces.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np

from .ebm import (
    EnergyBalanceModel, ClimateState, insolation, legendre_P2,
    S0_EARTH, S2_INSOLATION, A_OLR, B_OLR, D_TRANSPORT, T_FREEZE,
    ALBEDO_A0, ALBEDO_A2, ALBEDO_ICE, WATER_DEPTH,
)


def planetary_albedo(x, T, T_freeze: float = T_FREEZE,
                     a0: float = ALBEDO_A0, a2: float = ALBEDO_A2, ai: float = ALBEDO_ICE):
    """Step-function planetary albedo ``α(x, T)`` (climlab `StepFunctionAlbedo`).

    Ice-free ``α = a₀ + a₂·P₂(x)`` (poleward-brightening), or the ice value ``a_ice`` wherever
    ``T < Tf``. Returns an array matching the broadcast of ``x`` and ``T``. This is the **only
    nonlinear, state-dependent** term in the model — the source of the multiple equilibria.
    """
    x = np.asarray(x, dtype=float)
    T = np.asarray(T, dtype=float)
    ice_free = a0 + a2 * legendre_P2(x)
    return np.where(T < T_freeze, ai, ice_free)


def absorbed_shortwave(x, T, S0: float = S0_EARTH, s2: float = S2_INSOLATION,
                       T_freeze: float = T_FREEZE, a0: float = ALBEDO_A0,
                       a2: float = ALBEDO_A2, ai: float = ALBEDO_ICE):
    """Absorbed shortwave ``S(x)·(1 − α(x, T))`` (W m⁻²) — the EBM's radiation source.

    The incident insolation times the coalbedo, with the ice-albedo feedback in ``α``. This is
    the ``absorbed_fn(x, T)`` :meth:`EnergyBalanceModel.equilibrate` integrates; freezing it at
    each substep's ``T`` (the engine does) is what makes the local step a nonlinear relaxation.
    """
    return insolation(x, S0, s2) * (1.0 - planetary_albedo(x, T, T_freeze, a0, a2, ai))


@dataclass(frozen=True)
class EBMParams:
    """The full planet parameter bundle (climlab Earth defaults) — forcing + machinery.

    Bundles the radiation/albedo forcing (``S0, s2, Tf, a0, a2, ai``) with the EBM machinery
    constants (``A, B, D, water_depth, n_cells``) so a sweep can vary one knob (``S0``) with
    :func:`dataclasses.replace` while holding the rest fixed. :meth:`model` builds the
    (S₀-independent) :class:`~planet.ebm.EnergyBalanceModel`; :meth:`absorbed_fn`
    builds the S₀-/albedo-dependent forcing callable.
    """

    S0: float = S0_EARTH
    s2: float = S2_INSOLATION
    A: float = A_OLR
    B: float = B_OLR
    D: float = D_TRANSPORT
    T_freeze: float = T_FREEZE
    a0: float = ALBEDO_A0
    a2: float = ALBEDO_A2
    ai: float = ALBEDO_ICE
    water_depth: float = WATER_DEPTH
    n_cells: int = 180
    face: str = "harmonic"          # engine face-diffusivity mode: "harmonic" (default) | "exact"

    def model(self) -> EnergyBalanceModel:
        """Build the S₀-independent EBM machinery (transport + OLR + the ice-line isotherm)."""
        return EnergyBalanceModel(A=self.A, B=self.B, D=self.D, T_freeze=self.T_freeze,
                                  water_depth=self.water_depth, n_cells=self.n_cells, face=self.face)

    def absorbed_fn(self):
        """The ``absorbed_fn(x, T)`` closure: ``S(x)(1−α(x,T))`` at this parameter set's S₀/albedo."""
        return lambda x, T: absorbed_shortwave(
            x, T, S0=self.S0, s2=self.s2, T_freeze=self.T_freeze,
            a0=self.a0, a2=self.a2, ai=self.ai)


def present_day_climate(params: EBMParams | None = None, ic_equator: float = 30.0,
                        ic_pole: float = -30.0, **eq_kw) -> ClimateState:
    """The present-day equilibrium climate — a temperate planet with a polar ice cap (ice line ~70°).

    At the present solar constant the EBM is **bistable**: both an ice-free planet *and* a planet
    with a finite polar cap are stable (and, dimmer, a Snowball). Earth sits on the **finite-cap
    branch** — so this relaxes from an **Earth-like initial condition** (warm equator, frozen pole,
    ``T_init = ic_equator + (ic_pole − ic_equator)·x``) that lands on it, giving the ice line near
    70° (the loose climlab benchmark). A *warm uniform* start would instead settle on the (equally
    valid) ice-free branch — that the choice of start picks the branch is the very bistability the
    Snowball loop traces, here visible at present insolation. The Phase-1 climate Phase 2's biomes
    will consume.
    """
    if params is None:
        params = EBMParams()
    model = params.model()
    T_init = ic_equator + (ic_pole - ic_equator) * model.x      # warm equator → cold (iced) pole
    return model.equilibrate(params.absorbed_fn(), T_init, **eq_kw)


@dataclass(frozen=True)
class HysteresisLoop:
    """The Snowball hysteresis loop — the down- and up-sweep branches and the jumps between them.

    ``S0_down``/``Tbar_down``/``iceline_down`` trace the sun **dimming** (high → low S₀,
    warm-started from a temperate planet → the catastrophic freeze); ``*_up`` trace it
    **brightening** back (low → high, warm-started from the frozen planet → the late re-melt).
    Global-mean temperature (°C) and ice-line latitude (degrees) are recorded against the solar
    constant (W m⁻²). The branches differ over the bistable band — the loop. Plain arrays.
    """

    S0_down: np.ndarray
    Tbar_down: np.ndarray
    iceline_down: np.ndarray
    S0_up: np.ndarray
    Tbar_up: np.ndarray
    iceline_up: np.ndarray
    params: EBMParams = field(default_factory=EBMParams)

    @staticmethod
    def _jump_S0(S0: np.ndarray, Tbar: np.ndarray) -> float:
        """Midpoint S₀ of the adjacent pair with the largest |ΔT̄| — the catastrophic transition."""
        i = int(np.argmax(np.abs(np.diff(Tbar))))
        return 0.5 * (S0[i] + S0[i + 1])

    @property
    def freeze_S0(self) -> float:
        """The solar constant at which the dimming planet freezes over (the down-sweep jump)."""
        return self._jump_S0(self.S0_down, self.Tbar_down)

    @property
    def melt_S0(self) -> float:
        """The solar constant at which the frozen planet thaws (the up-sweep jump)."""
        return self._jump_S0(self.S0_up, self.Tbar_up)

    @property
    def hysteresis_width(self) -> float:
        """Loop width ``melt_S0 − freeze_S0`` (W m⁻²) — positive: it re-melts later than it froze."""
        return self.melt_S0 - self.freeze_S0


def snowball_hysteresis(params: EBMParams | None = None, S0_min: float = 1000.0,
                        S0_max: float = 1900.0, n_steps: int = 60, warm_start_T: float = 40.0,
                        **eq_kw) -> HysteresisLoop:
    """Trace the Snowball hysteresis loop by a parameter-continuation sweep in S₀.

    Sweeps the solar constant **down** from ``S0_max`` to ``S0_min`` (warm-started from an
    ice-free hot planet → the freeze), then **up** again to ``S0_max`` (warm-started from the
    frozen end state → the re-melt). Each S₀ relaxes to equilibrium from the previous one
    (continuation), so the solver stays on its current branch — and the two passes trace the two
    stable climates and the catastrophic jumps between them. Returns a :class:`HysteresisLoop`.
    """
    if params is None:
        params = EBMParams()
    model = params.model()
    S0_grid = np.linspace(S0_min, S0_max, n_steps)

    def sweep(S0_sequence, T0):
        T = T0
        Tbar, iceline = [], []
        for S0 in S0_sequence:
            st = model.equilibrate(replace(params, S0=S0).absorbed_fn(), T, **eq_kw)
            Tbar.append(st.global_mean_T)
            iceline.append(st.ice_line_lat)
            T = st.T                                     # warm-start the next S₀ (continuation)
        return np.array(Tbar), np.array(iceline), T

    S0_down = S0_grid[::-1].copy()
    Tbar_down, iceline_down, T_frozen = sweep(S0_down, warm_start_T)
    Tbar_up, iceline_up, _ = sweep(S0_grid, T_frozen)
    return HysteresisLoop(
        S0_down=S0_down, Tbar_down=Tbar_down, iceline_down=iceline_down,
        S0_up=S0_grid.copy(), Tbar_up=Tbar_up, iceline_up=iceline_up, params=params,
    )
