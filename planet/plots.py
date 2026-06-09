"""Planet-local plot helpers — the render layer (Planet Phase 1; ADR 0002).

The **viz floor**: static matplotlib figures that *consume* the plain arrays
:mod:`projects.planet.ebm` / :mod:`projects.planet.albedo` produce. Per ADR 0002 this layer is
strictly downstream of correctness — a figure draws already-validated numbers, it is never
evidence of validity (the §3 triad tests do that). It is the only place in planet that imports a
plotting library; the compute modules stay headless so the test suite never needs matplotlib.

The headline view is the **mechanism** one ADR 0002 §5 calls for: the Snowball **hysteresis traced
by the continuation sweep**, so a learner sees *why a dimming sun jumps to a frozen planet, and why
it stays frozen on the way back* — the bistability *seen*, not stated — beside the present-day
``T(φ)`` profile with its ice line. These helpers start project-local; the 2-D field/heatmap and
time-animation primitives are the ADR-0002 candidates a third reuse would promote to a shared
``viz/`` by rule-of-three (ARCHITECTURE.md §6).

Requires the optional ``viz`` extra (``pip install -e .[viz]``).
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from .albedo import HysteresisLoop
from .ebm import ClimateState, T_FREEZE, S0_EARTH

# Stable colours: cool blue = the dimming/cooling (down) branch, warm red = the brightening (up) branch.
COOLING_COLOR = "#2f6fb0"       # the sun dims → ice advances → the freeze jump
WARMING_COLOR = "#d4711f"       # the sun brightens → the late re-melt jump
PRESENT_COLOR = "#27795b"       # present-day Earth (the finite-cap branch)
ICE_COLOR = "#9ec9e2"           # the frozen / ice-covered region shading
ICELINE_COLOR = "#c0392b"       # the ice line (T = T_freeze isotherm)
FREEZE_COLOR = "#7f8c8d"        # the freeze-isotherm reference line


def _latitude_deg(x: np.ndarray) -> np.ndarray:
    """Cell-centre area coordinates x = sin φ → latitude φ in degrees (equator 0° → pole 90°)."""
    return np.degrees(np.arcsin(np.clip(np.asarray(x, dtype=float), 0.0, 1.0)))


def hysteresis_axes(ax, loop: HysteresisLoop, present: ClimateState | None = None,
                    present_S0: float = S0_EARTH) -> None:
    """Draw the global-mean-T hysteresis loop (T̄ vs S₀) with the freeze/melt jumps marked."""
    ax.plot(loop.S0_down, loop.Tbar_down, "-o", ms=3, color=COOLING_COLOR,
            label="dimming sun → freeze")
    ax.plot(loop.S0_up, loop.Tbar_up, "-o", ms=3, color=WARMING_COLOR,
            label="brightening sun → melt")
    ax.axvline(loop.freeze_S0, ls="--", lw=1, color=COOLING_COLOR, alpha=0.7)
    ax.axvline(loop.melt_S0, ls="--", lw=1, color=WARMING_COLOR, alpha=0.7)
    ax.annotate("Snowball\njump", (loop.freeze_S0, loop.Tbar_down.min()),
                textcoords="offset points", xytext=(8, 14), fontsize=8, color=COOLING_COLOR)
    ax.annotate("re-melt", (loop.melt_S0, loop.Tbar_up.max()),
                textcoords="offset points", xytext=(-38, -4), fontsize=8, color=WARMING_COLOR)
    if present is not None:
        # Present-day Earth is the finite-cap branch — a DISTINCT, slightly colder stable state than
        # the near-ice-free branch the dimming sweep traces at the same S₀ (today's sun admits both).
        ax.plot([present_S0], [present.global_mean_T], "*", ms=14, color=PRESENT_COLOR,
                label="present-day (finite-cap branch)", zorder=5)
    ax.axhline(0.0, color="#cccccc", lw=0.6, zorder=0)
    ax.set_xlabel("solar constant  S₀  (W m⁻²)")
    ax.set_ylabel("global-mean temperature  T̄  (°C)")
    ax.set_title("Snowball hysteresis: two climates for one sun")
    ax.legend(fontsize=8, loc="upper left")


def iceline_axes(ax, loop: HysteresisLoop) -> None:
    """Draw the ice-line latitude vs S₀ for both branches (the loop in ice-line space)."""
    ax.plot(loop.S0_down, loop.iceline_down, "-o", ms=3, color=COOLING_COLOR)
    ax.plot(loop.S0_up, loop.iceline_up, "-o", ms=3, color=WARMING_COLOR)
    ax.set_xlabel("solar constant  S₀  (W m⁻²)")
    ax.set_ylabel("ice-line latitude  (°)")
    ax.set_ylim(-3, 93)
    ax.set_title("Ice line: from pole (ice-free) to equator (Snowball)")


def profile_axes(ax, state: ClimateState) -> None:
    """Draw the present-day equilibrium T(φ) profile, shading the ice cap, marking the ice line."""
    lat = _latitude_deg(state.x)
    ax.plot(lat, state.T, "-", color=PRESENT_COLOR, lw=2, label="T(φ) equilibrium")
    ax.axhline(T_FREEZE, ls=":", color=FREEZE_COLOR, lw=1.2, label=f"freeze isotherm ({T_FREEZE:.0f} °C)")
    # shade the ice cap (poleward of the ice line)
    iced = state.T < T_FREEZE
    if np.any(iced):
        ax.fill_between(lat, state.T, T_FREEZE, where=iced, color=ICE_COLOR, alpha=0.7, label="ice cap")
    if 0.0 < state.ice_line_lat < 90.0:
        ax.axvline(state.ice_line_lat, ls="--", color=ICELINE_COLOR, lw=1.2)
        ax.annotate(f"ice line\n{state.ice_line_lat:.0f}°", (state.ice_line_lat, state.T.max()),
                    textcoords="offset points", xytext=(-30, -8), fontsize=8, color=ICELINE_COLOR)
    ax.set_xlabel("latitude  φ  (°)")
    ax.set_ylabel("temperature  T  (°C)")
    ax.set_xlim(0, 90)
    ax.set_title("Present-day climate (the finite-cap branch)")
    ax.legend(fontsize=8, loc="lower left")


def snowball_figure(loop: HysteresisLoop, present: ClimateState, present_S0: float = S0_EARTH):
    """The banked Phase-1 artifact: the hysteresis loop + ice-line loop + present-day T(φ) profile.

    Left (spanning): the T̄-vs-S₀ hysteresis loop with the two branches, the catastrophic jumps,
    and present-day Earth marked. Right: the ice-line-vs-S₀ loop (top) and the present-day
    temperature profile with its ice cap (bottom). *Knob in, frozen-or-temperate planet out.*
    """
    fig, axd = plt.subplot_mosaic(
        [["loop", "profile"], ["loop", "iceline"]], figsize=(13, 7), constrained_layout=True,
    )
    hysteresis_axes(axd["loop"], loop, present, present_S0)
    profile_axes(axd["profile"], present)
    iceline_axes(axd["iceline"], loop)
    fig.suptitle("Planet Phase 1 — the latitudinal EBM & the Snowball bifurcation", fontsize=13)
    return fig
