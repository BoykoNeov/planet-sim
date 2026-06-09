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
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch

from .albedo import HysteresisLoop
from .biomes import Biome, BIOME_COLORS, BIOME_NAMES, classify_field
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


# --------------------------------------------------------------------------- #
# Phase 2 — the biome map (the consequence payoff) + the mechanism (T, P) view.
# --------------------------------------------------------------------------- #
PRECIP_COLOR = "#2f6fb0"        # the precipitation profile (cool blue = water)


def _biome_cmap_norm():
    """A discrete colormap + boundary norm over the :class:`Biome` integer codes (for the map/diagram)."""
    codes = sorted(int(b) for b in Biome)
    cmap = ListedColormap([BIOME_COLORS[c] for c in codes])
    norm = BoundaryNorm([codes[0] - 0.5] + [c + 0.5 for c in codes], cmap.N)
    return cmap, norm


def _biome_legend(ax, present_codes) -> None:
    """Draw a biome legend on ``ax`` for the biomes present (ordered by code → polar → tropical)."""
    present = sorted(set(int(c) for c in np.asarray(present_codes).ravel()))
    handles = [Patch(facecolor=BIOME_COLORS[c], edgecolor="#444444", label=BIOME_NAMES[c]) for c in present]
    ax.legend(handles=handles, fontsize=8, loc="center", frameon=False, ncols=1, title="Biomes")
    ax.axis("off")


def biome_map_axes(ax, state: ClimateState, codes: np.ndarray) -> None:
    """Paint the latitude-banded biome map of the planet (the showcase artifact).

    The v1 planet is **zonal-mean** (plan §9.3): the biome code at each latitude is broadcast across
    longitude, so the globe paints **latitude bands** — honestly, not a premature 2-D field. The
    hemisphere is mirrored about the equator (the annual-mean climate is hemispherically symmetric),
    giving a full pole-to-pole band map. Knob in → these bands migrate.
    """
    lat = state.latitude_deg()
    full_lat = np.concatenate([-lat[::-1], lat])               # mirror to the southern hemisphere
    full_codes = np.concatenate([codes[::-1], codes]).astype(float)
    lon = np.linspace(-180.0, 180.0, 64)
    grid = np.repeat(full_codes[:, None], lon.size, axis=1)
    cmap, norm = _biome_cmap_norm()
    ax.pcolormesh(lon, full_lat, grid, cmap=cmap, norm=norm, shading="auto")
    ax.set_xlabel("longitude (zonal-mean → bands)")
    ax.set_ylabel("latitude  φ  (°)")
    ax.set_yticks([-90, -60, -30, 0, 30, 60, 90])
    ax.set_xticks([-180, -90, 0, 90, 180])
    ax.set_title("Biome bands of the planet")


def whittaker_axes(ax, state: ClimateState, precip_cm: np.ndarray) -> None:
    """Draw the (T, P) plane shaded by biome, with the planet's climate trajectory overlaid.

    The **mechanism view** (ADR 0002 §5): the background is the classifier's partition of the
    temperature–precipitation plane (this module's own Whittaker diagram), and the planet's
    equator→pole ``(T(φ), P(φ))`` curve is drawn through it — so a learner *sees* why the warm, dry
    subtropics land in the desert region (the curve dips through it near 25–30°) and why the wet warm
    equator is rain forest. Points are coloured by latitude.
    """
    Tg = np.linspace(-20.0, 32.0, 240)
    Pg = np.linspace(0.0, 400.0, 240)
    TT, PP = np.meshgrid(Tg, Pg)
    cmap, norm = _biome_cmap_norm()
    ax.pcolormesh(Tg, Pg, classify_field(TT, PP).astype(float), cmap=cmap, norm=norm,
                  shading="auto", alpha=0.85)
    lat = state.latitude_deg()
    sc = ax.scatter(state.T, precip_cm, c=lat, cmap="inferno", s=18, edgecolor="k",
                    linewidth=0.3, zorder=5)
    ax.plot(state.T, precip_cm, "-", color="#222222", lw=0.8, alpha=0.6, zorder=4)
    cb = ax.figure.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label("latitude φ (°)", fontsize=8)
    # mark the subtropical-desert dip = the driest WARM point (NOT the globally-driest pole),
    # the local precip minimum between the wet ITCZ and the midlatitude storm tracks.
    warm = np.asarray(state.T) > 12.0
    warm_idx = np.where(warm)[0]
    dry = int(warm_idx[np.argmin(np.asarray(precip_cm)[warm_idx])]) if warm_idx.size else int(np.argmin(precip_cm))
    ax.annotate("subtropics\n(~25–30°): desert", (state.T[dry], precip_cm[dry]),
                textcoords="offset points", xytext=(10, 18), fontsize=8, color="#7a5b1e",
                arrowprops=dict(arrowstyle="->", color="#7a5b1e", lw=0.8))
    ax.set_xlabel("mean annual temperature  T  (°C)")
    ax.set_ylabel("annual precipitation  P  (cm/yr)")
    ax.set_xlim(-20, 32)
    ax.set_ylim(0, 400)
    ax.set_title("Whittaker (T, P) plane — the planet's climate trajectory")


def biome_profile_axes(ax, state: ClimateState, precip_cm: np.ndarray) -> None:
    """Draw the T(φ) and P(φ) profiles vs latitude on twin axes (the two classifier inputs)."""
    lat = state.latitude_deg()
    ax.plot(lat, state.T, "-", color=PRESENT_COLOR, lw=2, label="T(φ)  (°C)")
    ax.axhline(T_FREEZE, ls=":", color=FREEZE_COLOR, lw=1.0)
    ax.set_xlabel("latitude  φ  (°)")
    ax.set_ylabel("temperature  T  (°C)", color=PRESENT_COLOR)
    ax.set_xlim(0, 90)
    axp = ax.twinx()
    axp.plot(lat, precip_cm, "-", color=PRECIP_COLOR, lw=2, label="P(φ)  (cm/yr)")
    axp.set_ylabel("precipitation  P  (cm/yr)", color=PRECIP_COLOR)
    axp.set_ylim(0, max(260.0, float(np.max(precip_cm)) * 1.1))
    ax.set_title("The two inputs: temperature & precipitation")


def biomes_figure(state: ClimateState, precip_cm: np.ndarray, codes: np.ndarray):
    """The banked Phase-2 artifact: the biome-band map + the Whittaker mechanism view + the inputs.

    Left (spanning): the **biome-band map** of the planet (the showcase). Top-right: the **Whittaker
    (T, P) plane** shaded by biome with the planet's climate trajectory drawn through it (the mechanism
    — *why* deserts sit at 30°). Bottom-right: the ``T(φ)`` and ``P(φ)`` profiles that feed it. Plus a
    biome legend. *Knob in, bands of life out.*
    """
    fig, axd = plt.subplot_mosaic(
        [["map", "whittaker"], ["map", "profile"], ["legend", "profile"]],
        figsize=(14, 8), constrained_layout=True,
        gridspec_kw={"height_ratios": [1.0, 0.85, 0.35], "width_ratios": [1.0, 1.15]},
    )
    biome_map_axes(axd["map"], state, codes)
    whittaker_axes(axd["whittaker"], state, precip_cm)
    biome_profile_axes(axd["profile"], state, precip_cm)
    _biome_legend(axd["legend"], codes)
    fig.suptitle("Planet Phase 2 — climate → habitability: the biome map", fontsize=13)
    return fig


# --------------------------------------------------------------------------- #
# Phase 3 — the shallow-water engine: geostrophic adjustment + a westward Rossby wave.
# --------------------------------------------------------------------------- #
def adjustment_axes(ax_init, ax_bal, adj) -> None:
    """Two panels: the initial unbalanced height bump, and the balanced remnant over L_R.

    A circle of radius ``L_R`` is drawn on the balanced panel to show the adjusted vortex has
    spread to the deformation radius — the geostrophic-adjustment scale made visible.
    """
    xk, yk = adj.x / 1e3, adj.y / 1e3                       # km
    draw = 100.0 * (1.0 - adj.eta_balanced.max() / adj.eta_init.max())
    # each panel on its OWN scale: the balanced remnant (~0.1× the bump) would vanish on a shared
    # one — the title carries the drawdown, the panels carry the *shape* (the L_R-scale vortex).
    panels = (
        (ax_init, adj.eta_init, "initial height anomaly (unbalanced)"),
        (ax_bal, adj.eta_balanced, f"balanced remnant ({draw:.0f}% radiated away)"),
    )
    for ax, field, title in panels:
        vmax = float(np.max(np.abs(field))) or 1.0
        im = ax.pcolormesh(xk, yk, field, cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="auto")
        ax.set_aspect("equal")
        ax.set_xlabel("x (km)"); ax.set_ylabel("y (km)")
        ax.set_title(title, fontsize=9)
        ax.figure.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="η (m)")
    cx, cy = adj.x[adj.x.size // 2] / 1e3, adj.y[adj.y.size // 2] / 1e3
    th = np.linspace(0, 2 * np.pi, 100)
    ax_bal.plot(cx + (adj.L_R / 1e3) * np.cos(th), cy + (adj.L_R / 1e3) * np.sin(th),
                "--", color="#222222", lw=1.2, label=f"L_R ≈ {adj.L_R/1e3:.0f} km")
    ax_bal.legend(fontsize=8, loc="upper right")


def rossby_axes(ax, ros) -> None:
    """A zonal cross-section of η through mid-latitude at successive times — the crest moves WEST.

    The clearest read of the dispersion: the same wave crest shifts toward −x (westward) as time
    advances, at the measured phase speed (printed against the analytic value)."""
    jmid = ros.y.size // 2
    xk = ros.x / 1e3
    n = len(ros.snapshots)
    for i, (snap, t) in enumerate(zip(ros.snapshots, ros.snapshot_times)):
        shade = 0.15 + 0.75 * i / max(1, n - 1)
        ax.plot(xk, snap[jmid, :], color=(0.1, 0.2, 0.55, shade),
                label=f"t = {t/86400:.1f} d")
    ax.set_xlabel("x (km)")
    ax.set_ylabel("η at mid-latitude (m)")
    ax.set_title(f"Rossby wave drifts WESTWARD  (c = {ros.c_measured:.1f} m/s, "
                 f"analytic {ros.c_analytic:.1f})", fontsize=9)
    ax.annotate("", xy=(0.12, 0.92), xytext=(0.42, 0.92), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="->", color="#c0392b", lw=1.6))
    ax.text(0.27, 0.95, "west", transform=ax.transAxes, ha="center", fontsize=8, color="#c0392b")
    ax.legend(fontsize=7, loc="lower right", ncols=2)


def conservation_axes(ax, adj) -> None:
    """The conservation diagnostics over the adjustment run: mass / energy / enstrophy drift.

    Mass holds at machine precision (flux form); energy and potential enstrophy hold to small
    bounded drifts — the honest 'invariants stay flat' panel the contract promises."""
    td = adj.times / 86400.0
    ax.plot(td, np.abs(adj.mass) + 1e-18, color="#27795b", label="mass |Δ| (machine-exact)")
    ax.plot(td, np.abs(adj.energy) + 1e-18, color="#2f6fb0", label="energy |Δ|")
    ax.plot(td, np.abs(adj.enstrophy) + 1e-18, color="#d4711f", label="enstrophy |Δ|")
    ax.set_yscale("log")
    ax.set_xlabel("time (days)")
    ax.set_ylabel("relative drift")
    ax.set_title("Conservation diagnostics hold flat", fontsize=9)
    ax.legend(fontsize=8, loc="center right")


def shallowwater_figure(adj, ros):
    """The banked Phase-3 artifact: geostrophic adjustment + a westward Rossby wave + conservation.

    Top row: the geostrophic adjustment — an unbalanced height bump (left) radiates gravity waves and
    settles to a balanced vortex of scale ``L_R`` (centre). Bottom-left: the conservation diagnostics
    (mass machine-exact; energy/enstrophy bounded). Right (spanning): a westward-propagating Rossby
    wave. *The rotating fluid engine, exercised with planetary numbers.*
    """
    fig, axd = plt.subplot_mosaic(
        [["init", "balanced", "rossby"], ["conserve", "conserve", "rossby"]],
        figsize=(16, 8.5), constrained_layout=True,
        gridspec_kw={"width_ratios": [1.0, 1.0, 1.2]},
    )
    adjustment_axes(axd["init"], axd["balanced"], adj)
    conservation_axes(axd["conserve"], adj)
    rossby_axes(axd["rossby"], ros)
    fig.suptitle("Planet Phase 3 — the rotating shallow-water engine (engines/fluid)", fontsize=13)
    return fig
