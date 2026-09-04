"""Planet-local plot helpers — the render layer (Planet Phase 1; ADR 0002).

The **viz floor**: static matplotlib figures that *consume* the plain arrays
:mod:`planet.ebm` / :mod:`planet.albedo` produce. Per ADR 0002 this layer is
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
from .ebm import ClimateState, T_FREEZE, S0_EARTH, ALBEDO_A0

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
# Rung 5A.2 — the orographic rain shadow placed on the sphere (the first step off the zonal mean).
# --------------------------------------------------------------------------- #
def orographic_scene_figure(scene):
    """The Rung-5A.2 regional scene: a rain shadow behind a mountain range under the westerly jet.

    Four panels on the scene's fine regional patch (:class:`planet.orographic_scene.OrographicScene`):
    the **terrain** with the prescribed cross-mountain wind arrow; the **orographic bonus** (the S&B
    windward rain + lee shadow); the **biome map** the enhancement re-classifies (the payoff — the
    mountain finally *changes the map*); and a mid-latitude **cross-section** of the zonal-mean baseline
    vs the enhanced total, over the elevation, that shows windward-wet / lee-dry directly. This is the
    honest deliverable: a *2-D precipitation* driven by geography, over a temperature climate that is
    still zonal-mean (:mod:`planet.orographic_scene`). Requires the ``viz`` extra.
    """
    lon, lat = scene.lon_deg, scene.lat_deg
    mid = scene.precip_cm.shape[0] // 2
    u_dir = "→ E (westerly)" if scene.wind_direction_deg == 270.0 else "← W (easterly)"

    fig, axd = plt.subplot_mosaic(
        [["terrain", "bonus"], ["biome", "section"]], figsize=(13, 9), constrained_layout=True,
    )

    # Terrain + the prescribed wind.
    ax = axd["terrain"]
    im = ax.pcolormesh(lon, lat, scene.elevation_m, cmap="terrain", shading="auto")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="elevation (m)")
    ax.set_title(f"terrain + prescribed wind ({scene.wind_speed:.0f} m/s {u_dir})", fontsize=9)
    ax.annotate("", xy=(0.72, 0.5), xytext=(0.28, 0.5), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color="#c0392b", lw=2.2))

    # The orographic bonus (windward rain + lee shadow).
    ax = axd["bonus"]
    im = ax.pcolormesh(lon, lat, scene.orographic_precip_cm, cmap="GnBu", shading="auto")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="orographic bonus (cm/yr)")
    ax.set_title("Smith–Barstad orographic precip — windward wet, lee shadow", fontsize=9)

    # The biome map the enhancement produces (the payoff), with baseline-changed cells outlined.
    ax = axd["biome"]
    cmap, norm = _biome_cmap_norm()
    ax.pcolormesh(lon, lat, scene.biome_codes, cmap=cmap, norm=norm, shading="auto")
    changed = scene.biome_codes != scene.baseline_biome_codes
    ax.contour(lon, lat, changed.astype(float), levels=[0.5], colors="#c0392b", linewidths=1.2)
    ax.set_title(f"biome map (payoff: {100 * scene.biome_changed_fraction:.0f}% of cells re-classified)",
                 fontsize=9)
    _biome_legend_inset(ax, scene.biome_codes)

    # Cross-section: baseline vs total precip, over the elevation. With Rung-5A.3 depletion on, the total
    # dips BELOW the baseline in the lee — the desert is the shaded gap; without it the lee sits at baseline.
    depleted = bool(np.any(scene.depletion_factor < 1.0 - 1e-9))
    ax = axd["section"]
    base_row, tot_row = scene.baseline_precip_cm[mid, :], scene.precip_cm[mid, :]
    ax.plot(lon, base_row, color="#7f8c8d", lw=1.6, ls="--", label="zonal-mean baseline")
    ax.plot(lon, tot_row, color=PRECIP_COLOR, lw=2.0, label="total")
    if depleted:
        ax.fill_between(lon, tot_row, base_row, where=tot_row < base_row, color="#c0392b", alpha=0.25,
                        label="lee desert (below baseline)")
    ax.set_ylabel("precipitation (cm/yr)")
    ax.set_title("cross-section — windward rain, lee " + ("desert (depleted)" if depleted else "shadow"),
                 fontsize=9)
    ax.legend(fontsize=8, loc="upper left")
    axe = ax.twinx()
    axe.fill_between(lon, 0, scene.elevation_m[mid, :], color="#c8b28a", alpha=0.45, zorder=0)
    axe.set_ylabel("elevation (m)", color="#8a6d3b")

    for key in ("terrain", "bonus", "biome"):
        axd[key].set_xlabel("longitude (°)"); axd[key].set_ylabel("latitude (°)")
    axd["section"].set_xlabel("longitude (°)")
    rung = "5A.3 — the rain-shadow desert" if depleted else "5A.2 — the orographic rain shadow"
    fig.suptitle(f"Planet Rung {rung}, placed on the sphere", fontsize=13)
    return fig


def _biome_legend_inset(ax, codes) -> None:
    """A compact in-axes biome legend for the regional biome panel (only the biomes present)."""
    present = sorted(set(int(c) for c in np.asarray(codes).ravel()))
    handles = [Patch(facecolor=BIOME_COLORS[c], edgecolor="#444444", label=BIOME_NAMES[c]) for c in present]
    ax.legend(handles=handles, fontsize=6.5, loc="lower right", frameon=True, framealpha=0.85, ncols=1)


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


# --------------------------------------------------------------------------- #
# Phase 4 — the one-way coupler: the EBM gradient forces an emergent geostrophic jet.
# --------------------------------------------------------------------------- #
WESTERLY_COLOR = "#b3361f"      # the emergent westerly jet
EASTERLY_COLOR = "#2f6fb0"      # the flanking easterly return (the periodic-channel consequence)
GEOSTROPHIC_COLOR = "#222222"   # the geostrophic estimate −(g/f)∂h/∂y


def jet_profile_axes(ax, jet) -> None:
    """The emergent zonal-wind profile ``u(φ)`` vs the geostrophic estimate — geostrophic balance made visible.

    The modelled zonal-mean wind (a westerly jet flanked by the periodic channel's compensating easterly
    return) lies on top of ``−(g/f)∂h/∂y`` in the jet core: the jet *is* in geostrophic balance with the
    forced height field.
    The jet maximum (the westerly) is marked, and the EBM gradient maximum is flagged — the jet sits
    *there*, not at the channel centre (the emergence).
    """
    phi = jet.phi
    ax.axhline(0.0, color="#999999", lw=0.8)
    ax.fill_between(phi, jet.u_profile, 0.0, where=jet.u_profile > 0, color=WESTERLY_COLOR, alpha=0.18)
    ax.fill_between(phi, jet.u_profile, 0.0, where=jet.u_profile < 0, color=EASTERLY_COLOR, alpha=0.18)
    ax.plot(phi, jet.u_profile, color=WESTERLY_COLOR, lw=2.0, label="modelled u (zonal mean)")
    ax.plot(phi, jet.u_geostrophic, "--", color=GEOSTROPHIC_COLOR, lw=1.3,
            label="geostrophic  −(g/f)∂h/∂y")
    ax.axvline(jet.jet_lat, color=WESTERLY_COLOR, lw=1.0, ls=":")
    ax.plot([jet.jet_lat], [jet.jet_speed], "o", color=WESTERLY_COLOR, ms=6)
    ax.annotate(f"westerly jet\n{jet.jet_speed:.0f} m/s @ {jet.jet_lat:.0f}°",
                xy=(jet.jet_lat, jet.jet_speed), xytext=(0.62, 0.82), textcoords="axes fraction",
                fontsize=8, color=WESTERLY_COLOR,
                arrowprops=dict(arrowstyle="->", color=WESTERLY_COLOR, lw=1.0))
    ax.axvline(jet.gradient_peak_lat, color="#7a7a7a", lw=1.0, ls="--")
    ax.text(jet.gradient_peak_lat, ax.get_ylim()[0], " EBM ∂T/∂φ max", rotation=90,
            va="bottom", ha="right", fontsize=7, color="#7a7a7a")
    ax.set_xlabel("latitude (°)")
    ax.set_ylabel("zonal wind u (m/s)   +east")
    ax.set_title(f"Emergent jet in geostrophic balance (core residual {100*jet.core_balance_residual:.1f}%)",
                 fontsize=9)
    ax.legend(fontsize=8, loc="lower left")


def coupling_chain_axes(ax, jet, state) -> None:
    """The forcing chain: the EBM temperature gradient sets the target height field (warm → high).

    Left axis: the EBM temperature at the channel latitudes (the climate gradient that does the
    forcing). Right axis: the windowed, zero-mean target height anomaly ``η_target`` the EBM hands the
    flow — high where warm. *Climate gradient in, height field out* — the coupler's one job.
    """
    phi = jet.phi
    T_chan = np.interp(phi, state.latitude_deg(), state.T)
    ax.plot(phi, T_chan, color="#c0392b", lw=1.8, label="EBM T(φ)")
    ax.set_xlabel("latitude (°)")
    ax.set_ylabel("EBM temperature (°C)", color="#c0392b")
    ax.tick_params(axis="y", labelcolor="#c0392b")
    ax.set_title("The forcing chain: warm climate → high target", fontsize=9)
    ax2 = ax.twinx()
    ax2.axhline(0.0, color="#999999", lw=0.6)
    ax2.plot(phi, jet.eta_target, color="#1f6f4a", lw=1.8, label="η_target (forced)")
    ax2.set_ylabel("target height anomaly η (m)", color="#1f6f4a")
    ax2.tick_params(axis="y", labelcolor="#1f6f4a")
    lines = ax.get_lines() + ax2.get_lines()[1:]
    ax.legend(lines, [l.get_label() for l in lines], fontsize=8, loc="upper right")


def jet_field_axes(ax, jet) -> None:
    """The 2-D jet on the β-plane channel: the zonal wind painted with flow arrows over it.

    The forcing is zonally symmetric, so the emergent jet is a coherent zonal band — red (eastward
    westerly) at midlatitudes, blue (westward easterly) on the flanks — the circulation the
    interactive map registers as a ``vector_overlay`` (Phase 4)."""
    xk = jet.x / 1e3
    lat = jet.phi
    vmax = float(np.max(np.abs(jet.u2d))) or 1.0
    im = ax.pcolormesh(xk, lat, jet.u2d, cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="auto")
    si, sj = max(1, jet.u2d.shape[1] // 12), max(1, jet.u2d.shape[0] // 16)
    ax.quiver(xk[::si], lat[::sj], jet.u2d[::sj, ::si], jet.v2d[::sj, ::si],
              color="#222222", scale=400, width=0.003, alpha=0.7)
    ax.set_xlabel("x (km)")
    ax.set_ylabel("latitude (°)")
    ax.set_title("The emergent circulation (zonal wind + flow)", fontsize=9)
    ax.figure.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="u (m/s)  +east")


def coupler_conservation_axes(ax, jet) -> None:
    """Conservation: mass machine-exact while forced; on RELEASE the bare engine conserves all invariants.

    The honest forced–dissipative story (the reframed conservation leg). During forcing, mass holds at
    machine precision but energy and potential enstrophy are *not* conserved — the forcing–drag balance
    is what selects the steady jet. When the forcing is switched **off** (shaded), the bare engine
    conserves mass / energy / enstrophy (its Phase-3 guarantees) — and the jet persists.
    """
    tf = jet.times / 86400.0
    t0 = tf[-1] if tf.size else 0.0
    tr = t0 + jet.times_release / 86400.0
    eps = 1e-18
    ax.plot(tf, np.abs(jet.mass) + eps, color="#27795b", lw=1.4, label="mass |Δ| (machine-exact)")
    ax.plot(tf, np.abs(jet.energy) + eps, color="#2f6fb0", lw=1.4, label="energy |Δ| (forced: not conserved)")
    ax.plot(tf, np.abs(jet.enstrophy) + eps, color="#d4711f", lw=1.4, label="enstrophy |Δ|")
    ax.plot(tr, np.abs(jet.mass_release) + eps, color="#27795b", lw=1.4)
    ax.plot(tr, np.abs(jet.energy_release) + eps, color="#2f6fb0", lw=1.4)
    ax.plot(tr, np.abs(jet.enstrophy_release) + eps, color="#d4711f", lw=1.4)
    if tr.size:
        ax.axvspan(t0, tr[-1], color="#cccccc", alpha=0.25)
        ax.text(0.5 * (t0 + tr[-1]), ax.get_ylim()[1] if False else 1e-1, "forcing OFF\n(release)",
                ha="center", va="top", fontsize=7, color="#555555")
    ax.axvline(t0, color="#555555", lw=0.8, ls="--")
    ax.set_yscale("log")
    ax.set_xlabel("time (days)")
    ax.set_ylabel("relative drift")
    ax.set_title("Conservation: mass forced-exact; release re-confirms the engine", fontsize=9)
    ax.legend(fontsize=7, loc="lower right")


def coupler_figure(jet, state):
    """The banked Phase-4 artifact: the EBM gradient forces an emergent, geostrophically-balanced jet.

    Top-left: the emergent zonal-wind profile on top of the geostrophic estimate (balance made visible).
    Top-right: the forcing chain (warm EBM → high target height). Bottom-left: the 2-D jet on the
    channel. Bottom-right: the conservation diagnostics — mass forced-exact, then the release test
    re-confirming the engine's invariants. *Climate in, circulation out — the two engines coupled.*
    """
    fig, axd = plt.subplot_mosaic(
        [["jet", "chain"], ["field", "conserve"]],
        figsize=(15, 9), constrained_layout=True,
    )
    jet_profile_axes(axd["jet"], jet)
    coupling_chain_axes(axd["chain"], jet, state)
    jet_field_axes(axd["field"], jet)
    coupler_conservation_axes(axd["conserve"], jet)
    fig.suptitle("Planet Phase 4 — one-way EBM → shallow-water coupler: the emergent jet", fontsize=13)
    return fig


# --------------------------------------------------------------------------- #
# §9.1 — the exoplanet knobs: stellar spectrum → ice albedo, planet size → transport.
# --------------------------------------------------------------------------- #
SUN_COLOR = "#d98a1f"           # the Sun-like (G-type) host star
MDWARF_COLOR = "#7b3fb0"        # the redder, cooler M-dwarf host star


def stellar_hysteresis_axes(ax, sun_loop: HysteresisLoop, mdwarf_loop: HysteresisLoop,
                            mdwarf_label: str) -> None:
    """Overlay the Snowball hysteresis loops of a Sun-like vs an M-dwarf planet — redder = harder to snowball.

    Two T̄-vs-S₀ loops (solid = dimming/freeze branch, dashed = brightening/re-melt branch), one per
    host star, with each star's freeze threshold marked. The M-dwarf's weaker ice-albedo feedback gives
    a **narrower loop shifted to a lower freeze threshold** — it must be dimmed *further* to snowball,
    and re-melts sooner. The mechanism made visible (ADR 0002 §5).
    """
    for loop, color, lbl in ((sun_loop, SUN_COLOR, "Sun (G2V)"), (mdwarf_loop, MDWARF_COLOR, mdwarf_label)):
        ax.plot(loop.S0_down, loop.Tbar_down, "-", color=color, lw=1.8, label=f"{lbl}: dimming → freeze")
        ax.plot(loop.S0_up, loop.Tbar_up, "--", color=color, lw=1.4, label=f"{lbl}: brightening → melt")
        ax.axvline(loop.freeze_S0, ls=":", color=color, lw=1.0, alpha=0.8)
        ax.annotate(f"freeze\n{loop.freeze_S0:.0f}", (loop.freeze_S0, loop.Tbar_down.min()),
                    textcoords="offset points", xytext=(6, 12), fontsize=7, color=color)
    ax.axhline(0.0, color="#cccccc", lw=0.6, zorder=0)
    ax.set_xlabel("solar constant  S₀  (W m⁻²)")
    ax.set_ylabel("global-mean temperature  T̄  (°C)")
    ax.set_title("A redder star is harder to snowball\n(weaker ice-albedo feedback → narrower loop)", fontsize=10)
    ax.legend(fontsize=7, loc="lower right")


def stellar_albedo_axes(ax, stellar_ai: dict) -> None:
    """Bar the effective ice albedo by stellar type — the bright-ice albedo falls toward the redder stars.

    The knob's mechanism in one panel: cooler (redder) host stars emit more near-IR, where ice is dark,
    so the broadband ice albedo drops; the ice-free ocean/land albedo line shows the contrast weakening
    but never inverting (the bounded, *modest* effect).
    """
    names = list(stellar_ai)
    vals = [stellar_ai[n] for n in names]
    colors = plt.cm.YlOrRd_r(np.linspace(0.15, 0.8, len(names)))
    ax.bar(names, vals, color=colors, edgecolor="#444444", width=0.7)
    ax.axhline(ALBEDO_A0, ls="--", color="#2f6fb0", lw=1.2, label=f"ice-free ocean/land α₀ = {ALBEDO_A0:.2f}")
    ax.set_ylabel("effective ice albedo  a_ice")
    ax.set_ylim(0.0, 0.75)
    ax.tick_params(axis="x", labelrotation=35, labelsize=7)
    ax.set_title("Ice albedo weakens toward redder stars\n(never below the ocean albedo)", fontsize=10)
    ax.legend(fontsize=7, loc="upper left")


def size_profiles_axes(ax, size_states, sizes) -> None:
    """Overlay the relaxed T(φ) for several planet sizes — a bigger planet sharpens the equator–pole gradient.

    Bigger planet → weaker per-area transport (``D ∝ 1/size²``) → a steeper gradient and a colder pole,
    so the ice cap reaches further equatorward — while the global mean is (nearly) fixed (the 0-D mean
    is size-invariant; the relaxed mean drifts only through the ice-albedo feedback). The label carries
    each size's mean T̄ and ice-line latitude.
    """
    colors = plt.cm.viridis(np.linspace(0.15, 0.82, len(sizes)))
    for st, size, c in zip(size_states, sizes, colors):
        lat = st.latitude_deg()
        label = f"{size:g} R⊕  (T̄ {st.global_mean_T:.1f}°C, ice {st.ice_line_lat:.0f}°)"
        ax.plot(lat, st.T, "-", color=c, lw=2.0, label=label)
        if 0.0 < st.ice_line_lat < 90.0:
            ax.axvline(st.ice_line_lat, ls=":", color=c, lw=1.0, alpha=0.7)
    ax.axhline(T_FREEZE, ls="--", color=FREEZE_COLOR, lw=1.0, label=f"freeze isotherm ({T_FREEZE:.0f} °C)")
    ax.set_xlabel("latitude  φ  (°)")
    ax.set_ylabel("temperature  T  (°C)")
    ax.set_xlim(0, 90)
    ax.set_title("A bigger planet sharpens the gradient\n(transport-only: D ∝ 1/size²)", fontsize=10)
    ax.legend(fontsize=7, loc="lower left")


def exoplanet_figure(sun_loop: HysteresisLoop, mdwarf_loop: HysteresisLoop, mdwarf_label: str,
                     stellar_ai: dict, size_states, sizes):
    """The banked §9.1 artifact: the two exoplanet knobs — stellar spectrum and planet size.

    Top-left: the Snowball hysteresis loops of a Sun-like vs an M-dwarf planet (redder = harder to
    snowball). Top-right: the effective ice albedo by stellar type (the knob's mechanism, bounded above
    the ocean albedo). Bottom (spanning): the relaxed T(φ) for a range of planet sizes (a bigger planet
    sharpens the equator-to-pole gradient, transport-only). *Other-world knobs in, other-world climate out.*
    """
    fig, axd = plt.subplot_mosaic(
        [["stellar", "albedo"], ["size", "size"]],
        figsize=(14, 9), constrained_layout=True,
        gridspec_kw={"height_ratios": [1.0, 0.9]},
    )
    stellar_hysteresis_axes(axd["stellar"], sun_loop, mdwarf_loop, mdwarf_label)
    stellar_albedo_axes(axd["albedo"], stellar_ai)
    size_profiles_axes(axd["size"], size_states, sizes)
    fig.suptitle("Planet §9.1 — the exoplanet knobs: stellar spectrum & planet size", fontsize=13)
    return fig


# --------------------------------------------------------------------------- #
# §9.1 obliquity knob — axial tilt → the annual-mean-insolation P₂ coefficient s₂
# --------------------------------------------------------------------------- #
from .obliquity import OBLIQUITY_EARTH, OBLIQUITY_FAITHFUL_MAX   # noqa: E402 (grouped with its figure)

S2_CURVE_COLOR = "#6a4c93"      # the geometric s₂(ε) curve
ANCHOR_COLOR = "#c0392b"        # the exact analytic anchors (−5/8, Earth)


def obliquity_s2_axes(ax, eps_grid, s2_grid, s2_earth) -> None:
    """Plot the geometric s₂(ε) curve — how axial tilt grades the annual-mean insolation equator→pole.

    ``s₂`` rises from **exactly −5/8** at zero tilt (sun pinned at the equator) toward zero and then
    **positive** past the ≈55° critical obliquity, where the poles receive more annual sun than the
    equator (shaded). Earth's tilt is marked at ``s₂ ≈ −0.48`` (the independent climlab cross-check), and
    the 0–45° band the live slider drives is shaded — beyond it the EBM's single-P₂-mode insolation only
    approximates the strongly-flattened profile (the named scope edge).
    """
    ax.axhline(0.0, color="#cccccc", lw=0.8, zorder=0)
    ax.axvspan(0.0, OBLIQUITY_FAITHFUL_MAX, color="#f0ecf7", zorder=0,
               label=f"wired slider range (0–{OBLIQUITY_FAITHFUL_MAX:g}°)")
    cross = float(np.interp(0.0, s2_grid, eps_grid))          # where the gradient reverses (s₂ = 0)
    ax.axvspan(cross, eps_grid[-1], color="#fdebd0", zorder=0, label="poles warmer than equator (s₂ > 0)")
    ax.plot(eps_grid, s2_grid, "-", color=S2_CURVE_COLOR, lw=2.2, label="s₂(ε)  (annual-mean insolation)")
    # the exact analytic anchor at ε = 0 (−5/8) and the Earth cross-check point
    ax.plot(0.0, -0.625, "o", color=ANCHOR_COLOR, ms=7, zorder=5)
    ax.annotate("−5/8 exactly\n(no tilt)", (0.0, -0.625), textcoords="offset points",
                xytext=(10, 4), fontsize=8, color=ANCHOR_COLOR)
    ax.plot(OBLIQUITY_EARTH, s2_earth, "s", color=ANCHOR_COLOR, ms=7, zorder=5)
    ax.annotate(f"Earth {OBLIQUITY_EARTH:.2f}°\ns₂ ≈ {s2_earth:.2f}  (≈ climlab −0.48)",
                (OBLIQUITY_EARTH, s2_earth), textcoords="offset points", xytext=(10, -28), fontsize=8,
                color=ANCHOR_COLOR)
    ax.annotate(f"gradient reverses ≈ {cross:.0f}°", (cross, 0.0), textcoords="offset points",
                xytext=(-4, 10), fontsize=8, color="#a6611a", ha="right")
    ax.set_xlabel("obliquity  ε  (°)")
    ax.set_ylabel("insolation P₂ coefficient  s₂")
    ax.set_xlim(0, 90)
    ax.set_title("Axial tilt sets the insolation gradient\n(exact −5/8 at 0°; reverses past ≈55°)", fontsize=10)
    ax.legend(fontsize=7, loc="lower right")


def obliquity_climate_axes(ax, climate_states, tilts) -> None:
    """Overlay the relaxed T(φ) for several obliquities — more tilt flattens the planet, the ice retreats.

    Less tilt → a steeper gradient and a colder pole (a larger ice cap); more tilt spreads the sunlight
    poleward, warming the pole until the cap is gone. Each label carries the tilt, the global mean T̄ and
    the ice-line latitude (Earth's 23.44° at ~70°, the climlab benchmark).
    """
    colors = plt.cm.cividis(np.linspace(0.1, 0.85, len(tilts)))
    for st, tilt, c in zip(climate_states, tilts, colors):
        lat = st.latitude_deg()
        label = f"ε {tilt:g}°  (T̄ {st.global_mean_T:.1f}°C, ice {st.ice_line_lat:.0f}°)"
        ax.plot(lat, st.T, "-", color=c, lw=2.0, label=label)
        if 0.0 < st.ice_line_lat < 90.0:
            ax.axvline(st.ice_line_lat, ls=":", color=c, lw=1.0, alpha=0.7)
    ax.axhline(T_FREEZE, ls="--", color=FREEZE_COLOR, lw=1.0, label=f"freeze isotherm ({T_FREEZE:.0f} °C)")
    ax.set_xlabel("latitude  φ  (°)")
    ax.set_ylabel("temperature  T  (°C)")
    ax.set_xlim(0, 90)
    ax.set_title("More tilt flattens the planet\n(the pole warms, the ice cap retreats)", fontsize=10)
    ax.legend(fontsize=7, loc="lower left")


def obliquity_figure(eps_grid, s2_grid, s2_earth, climate_states, tilts):
    """The banked §9.1 artifact: the obliquity knob — axial tilt → the insolation gradient → the climate.

    Left: the geometric ``s₂(ε)`` curve (the mechanism — exact −5/8 at no tilt, the Earth cross-check,
    the ≈55° gradient reversal). Right: the relaxed ``T(φ)`` at a range of obliquities (the consequence —
    more tilt warms the pole and retreats the ice cap). *Axial tilt in, climate gradient out.*
    """
    fig, axd = plt.subplot_mosaic(
        [["s2", "climate"]], figsize=(14, 5.5), constrained_layout=True,
    )
    obliquity_s2_axes(axd["s2"], eps_grid, s2_grid, s2_earth)
    obliquity_climate_axes(axd["climate"], climate_states, tilts)
    fig.suptitle("Planet §9.1 — the obliquity knob: axial tilt sets the annual-mean-insolation gradient", fontsize=13)
    return fig


# --------------------------------------------------------------------------- #
# Rung A — the animated eddy life cycle (the program's FIRST time-animation primitive; §9.5).
# The MECHANISM, two panels: the tracer stirred by the released eddies (left) + the cumulative
# meridional transport (right) — throughput rages while net stays small, so eddy_flux's
# ~90%-reversible finding is made VISIBLE (the overclaim a bare stirring movie would commit).
# --------------------------------------------------------------------------- #
THROUGHPUT_COLOR = "#b3361f"    # the |F̄| throughput — the raging swirls
NET_COLOR = "#1f6f4a"           # the net ∫F̄ transport — the small down-gradient residual
WINDOW_COLOR = "#7a7a7a"        # the κ-diagnosis window onset
SATURATION_COLOR = "#6a4c93"    # the eddy-KE saturation time


def eddy_life_animation(eddy, *, interval: int = 120):
    """The banked rung-A artifact: the emergent eddy life cycle animated — the MECHANISM, two panels.

    Left — **the stirring**: the passive temperature tracer ``θ`` advected by the released,
    barotropically-unstable jet on the **midlatitude β-plane channel** (a longitude × latitude *band*,
    NOT a globe — the honest scope), with the **eddy velocity** ``(u−ū, v−v̄)`` overlaid as arrows — the
    swirls doing the stirring. Right — **the transport budget**: the cumulative meridional eddy heat
    flux, the **throughput** ``Σ∫|F̄|dt`` climbing steeply while the **net** ``Σ|∫F̄dt|`` (integrated per
    latitude first) stays a small fraction — so :mod:`planet.eddy_flux`'s headline finding, that the
    instantaneous flux is **~90 % reversible** (``irreversible_fraction ~0.1``), is made *visible*: the
    swirls rage, the net barely moves. Without this second panel a stirring movie would *contradict* the
    module's own finding (ADR 0002 §5: visualize the mechanism, not the output). A dashed line marks
    where κ is diagnosed from (the window onset), reconciling the full-release curve with the banked
    windowed number.

    Returns a :class:`matplotlib.animation.FuncAnimation`; save it with a Pillow writer (GIF, CI-safe)
    or ffmpeg (MP4). Requires the frames side-channel — raises :class:`ValueError` if ``eddy.frames``
    is ``None`` (recompute with ``eddy_life_cycle(..., n_frames=N)``).
    """
    fr = eddy.frames
    if fr is None:
        raise ValueError("eddy.frames is None — recompute with eddy_life_cycle(..., n_frames=N)")
    from matplotlib.animation import FuncAnimation

    irr = eddy.irreversible_fraction
    rev_pct = round(100 * (1 - irr))     # the share of the stirring that cancels out (≈90%)
    res_pct = round(100 * irr)           # the residual that survives as net κ transport (≈10%)

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 6.6), constrained_layout=True)
    # reserve the bottom band of the figure for the plain-language caption: constrained_layout lays the
    # two panels into the top ~83%, leaving the lower ~17% for fig.text — so the standalone GIF stands on
    # its own (a standing preference; planet memory: viz-prose-novice-intermediate).
    fig.get_layout_engine().set(rect=(0.0, 0.17, 1.0, 0.83))

    # -- left: the tracer field stirred by the eddies. Colour range AND quiver scale are fixed across
    #    frames so the eye reads the genuine GROWTH, not matplotlib's per-frame autoscale. -- #
    xk = fr.x / 1e3                                          # km
    lat = fr.phi                                             # deg
    tmin, tmax = float(fr.theta.min()), float(fr.theta.max())
    mesh = axL.pcolormesh(xk, lat, fr.theta[0], cmap="RdBu_r", vmin=tmin, vmax=tmax, shading="nearest")
    fig.colorbar(mesh, ax=axL, fraction=0.046, pad=0.04, label="tracer θ (°C)")
    si = max(1, fr.x.size // 16)
    sj = max(1, fr.phi.size // 16)
    up0 = fr.u[0] - fr.u[0].mean(axis=1, keepdims=True)
    vp0 = fr.v[0] - fr.v[0].mean(axis=1, keepdims=True)
    eddy_max = float(max(np.abs(fr.u - fr.u.mean(axis=2, keepdims=True)).max(),
                         np.abs(fr.v - fr.v.mean(axis=2, keepdims=True)).max(), 1e-6))
    q = axL.quiver(xk[::si], lat[::sj], up0[::sj, ::si], vp0[::sj, ::si], color="#222222", alpha=0.7,
                   scale_units="width", scale=eddy_max / 0.07, width=0.004)
    axL.set_xlabel("x (km)")
    axL.set_ylabel("latitude φ (°)")
    axL.set_title("temperature θ stirred by the eddies\n(one midlatitude band — a β-plane channel, not a globe)",
                  fontsize=9)
    tlabel = axL.text(0.02, 0.96, "", transform=axL.transAxes, fontsize=9, va="top", color="#222222",
                      bbox=dict(boxstyle="round", fc="white", ec="#cccccc", alpha=0.8))

    # -- right: the cumulative transport budget — throughput rages, net stays small -- #
    # plain-language legend (novice→intermediate): "back-and-forth stirring" / "net heat moved" carry the
    # meaning, the technical term (throughput / κ residual) rides parenthetically — the cryptic "swirls
    # rage" / "small residual" editorializing + the raw formulas drop to the caption (where each is glossed).
    axR.plot(fr.times, fr.thru_cum, color=THROUGHPUT_COLOR, lw=2.0,
             label="back-and-forth stirring (throughput)")
    axR.plot(fr.times, fr.net_cum, color=NET_COLOR, lw=2.0,
             label="net heat moved poleward (κ residual)")
    ytop = axR.get_ylim()[1]
    axR.axvline(fr.window_start, ls="--", color=WINDOW_COLOR, lw=1.2)
    axR.text(fr.window_start, ytop, " κ diagnosed →", rotation=90, va="top", ha="left",
             fontsize=7, color=WINDOW_COLOR)
    if eddy.saturation_period > fr.times[0]:
        axR.axvline(eddy.saturation_period, ls=":", color=SATURATION_COLOR, lw=1.2)
        axR.text(eddy.saturation_period, ytop, " eddy-KE peak", rotation=90, va="top", ha="right",
                 fontsize=7, color=SATURATION_COLOR)
    cursor = axR.axvline(fr.times[0], color="#222222", lw=1.0)
    pt_thru, = axR.plot([fr.times[0]], [fr.thru_cum[0]], "o", color=THROUGHPUT_COLOR, ms=6)
    pt_net, = axR.plot([fr.times[0]], [fr.net_cum[0]], "o", color=NET_COLOR, ms=6)
    axR.set_xlabel("release time (inertial periods)")
    axR.set_ylabel("cumulative heat transport  (interior-band sum, K·m)")
    axR.set_title(f"Total stirring races up, net transport stays small\n"
                  f"(only ~{res_pct}% survives as real poleward heat)", fontsize=9)
    axR.legend(fontsize=7, loc="upper left")

    # The one piece of NEW prose: a plain-language caption (novice→intermediate) so the standalone GIF
    # stands on its own — the same de-jargoning + formula-with-gloss the eddy-globe exemplar carries
    # (planet memory: viz-prose-novice-intermediate). matplotlib fig.text can't tint spans like the Plotly
    # globe, so the legend's curve colours carry the colour cue; the formulas use mathtext ($\bar F$)
    # rather than a combining-macron unicode char that fonts render poorly. Wrapped by hand to the figure
    # width; the %s are driven off the same `irreversible_fraction` as the panel so they agree.
    caption = (
        "How to read it.  Left: temperature (θ) stirred by swirling eddies on one midlatitude band — a flat "
        "β-plane channel, not a\n"
        f"whole globe. The eddies push heat poleward one moment and back the next, so most of it cancels: "
        f"about {rev_pct}% is reversible.\n"
        r"Right: total stirring (throughput, $\Sigma\!\int\!|\bar{F}|\,dt$ — the flux size summed over time "
        r"and latitude, ignoring direction)" "\n"
        r"races up, while net transport ($\Sigma|\!\int\!\bar{F}\,dt|$ — what's left once the poleward and "
        r"equatorward parts cancel) barely moves." "\n"
        f"Only that surviving ~{res_pct}% — the eddy diffusivity κ — is real poleward heat transport."
    )
    fig.text(0.5, 0.095, caption, ha="center", va="center", fontsize=10.5, color="#33373b")

    def update(k):
        mesh.set_array(fr.theta[k].ravel())
        upk = fr.u[k] - fr.u[k].mean(axis=1, keepdims=True)
        vpk = fr.v[k] - fr.v[k].mean(axis=1, keepdims=True)
        q.set_UVC(upk[::sj, ::si], vpk[::sj, ::si])
        tlabel.set_text(f"t = {fr.times[k]:.0f} periods")
        cursor.set_xdata([fr.times[k], fr.times[k]])
        pt_thru.set_data([fr.times[k]], [fr.thru_cum[k]])
        pt_net.set_data([fr.times[k]], [fr.net_cum[k]])
        return mesh, q, tlabel, cursor, pt_thru, pt_net

    return FuncAnimation(fig, update, frames=fr.times.size, interval=interval, blit=False)


# --------------------------------------------------------------------------- #
# Rung 5B.1 — the seasonal cycle & continentality.
# --------------------------------------------------------------------------- #
SEASON_LAND = "#c0622d"          # the small-C land tile — big, prompt seasonal swing
SEASON_OCEAN = "#2f6fb0"         # the large-C ocean tile — small, lagged swing
SEASON_SUN = "#e0b000"           # the insolation forcing


def seasonal_figure(result):
    """The Rung-5B.1 payoff: heat capacity woken, continentality cast (:mod:`planet.seasonal`).

    Three panels off the converged annual limit cycle. **The cycle** (a midlatitude band, default 45°N):
    the land tile's temperature swings hard and nearly in step with the sun, while the ocean tile at the
    *same latitude* barely moves and lags ~2 months — continentality, from the ``C`` contrast alone (the
    sun is drawn on a twin axis). **The march** — a latitude–month Hovmöller of the land-tile temperature,
    the seasons sweeping pole to pole, the two hemispheres half a year out of phase (the antisymmetry).
    **The law** — seasonal amplitude vs latitude for both tiles: the land curve towers over the ocean
    curve everywhere, growing toward the poles. Requires the ``viz`` extra.
    """
    m, clim = result.model, result.climate
    lat = clim.latitude_deg()
    months = clim.days / (365.25 / 12.0)
    i = m.nearest_index(result.band_lat_deg)

    fig, axd = plt.subplot_mosaic(
        [["cycle", "hov"], ["amp", "hov"]], figsize=(13, 8.5), constrained_layout=True,
    )

    # -- the cycle at one midlatitude band ------------------------------------ #
    ax = axd["cycle"]
    ax.plot(months, clim.T_land[i], color=SEASON_LAND, lw=2.2, label="land tile")
    ax.plot(months, clim.T_ocean[i], color=SEASON_OCEAN, lw=2.2, label="ocean tile")
    ax.set_ylabel("surface temperature (°C)")
    ax.set_title(f"the seasonal cycle at {abs(lat[i]):.0f}°{'N' if lat[i] >= 0 else 'S'} — "
                 f"land swings {clim.amplitude('land')[i]:.0f} K, ocean {clim.amplitude('ocean')[i]:.1f} K",
                 fontsize=9)
    axs = ax.twinx()
    axs.plot(months, m.insolation_series()[i], color=SEASON_SUN, lw=1.4, ls=":", label="insolation")
    axs.set_ylabel("insolation (W/m²)", color="#9a7b00")
    ax.legend(fontsize=8, loc="upper left"); axs.legend(fontsize=8, loc="upper right")

    # -- the pole-to-pole seasonal march (Hovmöller, land tile) --------------- #
    ax = axd["hov"]
    im = ax.pcolormesh(months, lat, clim.T_land, cmap="RdBu_r", shading="auto",
                       vmin=-np.max(np.abs(clim.T_land)), vmax=np.max(np.abs(clim.T_land)))
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="land-tile T (°C)")
    ax.axhline(0.0, color="k", lw=0.6, alpha=0.5)
    ax.set_title("the seasons sweep pole to pole (hemispheres ½ year out of phase)", fontsize=9)
    ax.set_ylabel("latitude (°)")

    # -- the amplitude law vs latitude ---------------------------------------- #
    ax = axd["amp"]
    ax.plot(lat, clim.amplitude("land"), color=SEASON_LAND, lw=2.2, label="land tile")
    ax.plot(lat, clim.amplitude("ocean"), color=SEASON_OCEAN, lw=2.2, label="ocean tile")
    ax.fill_between(lat, clim.amplitude("ocean"), clim.amplitude("land"),
                    color=SEASON_LAND, alpha=0.10)
    ax.set_xlabel("latitude (°)"); ax.set_ylabel("seasonal amplitude (K)")
    ax.set_title("continentality — the land tile towers over the ocean at every latitude", fontsize=9)
    ax.legend(fontsize=8, loc="upper center")

    for key in ("cycle", "hov"):
        axd[key].set_xlabel("month of year"); axd[key].set_xlim(0, 12)
        axd[key].set_xticks(range(0, 13, 2))

    fig.suptitle("Planet Rung 5B.1 — the seasonal cycle wakes heat capacity, casting continentality",
                 fontsize=13)
    return fig


# --------------------------------------------------------------------------- #
# Rung 5B.2 — the 2-D map: continentality resolved in longitude.
# --------------------------------------------------------------------------- #
def seasonal_map_figure(result):
    """The Rung-5B.2 payoff: continentality as a MAP (:mod:`planet.seasonal_map`, NMS83).

    Four panels off the converged 2-D annual limit cycle. **The seasonal-range map** (top-left): peak-to-
    peak seasonal range over the globe, coastlines overlaid — continental interiors blaze with tens of
    kelvin while the oceans stay maritime-cool, the interior/coast/ocean sample points marked. **The
    annual-mean map** (top-right): the *same* field averaged over the year — visibly **zonally flat**, the
    mask invisible (average over the year and ``C`` cancels; the land/sea signal is *all* seasonal
    amplitude). **The cross-section** (bottom-left): seasonal range vs longitude along the midlatitude
    band, land shaded — the range peaks deep in the continent and is drawn down over the sea, the coast
    moderated between. **The cycle** (bottom-right): the seasonal temperature cycle at the interior, coast,
    and ocean points — the interior swings hard and prompt, the ocean barely moves and lags. Requires the
    ``viz`` extra.
    """
    m, c = result.model, result.climate
    lat, lon = c.latitude_deg(), c.longitude_deg()
    months = c.days / (365.25 / 12.0)
    i = result.band_index()
    interior, coast, ocean = result.sample_columns()
    rng = c.seasonal_range()
    amean = c.annual_mean()
    mask = c.land_mask.astype(float)
    sample_colors = (SEASON_LAND, "#7a4bd0", SEASON_OCEAN)   # interior / coast / ocean

    fig, axd = plt.subplot_mosaic(
        [["range", "amean"], ["cross", "cycle"]], figsize=(13.5, 8.5), constrained_layout=True,
    )

    # -- the seasonal-range map (the headline) -------------------------------- #
    ax = axd["range"]
    im = ax.pcolormesh(lon, lat, rng, cmap="inferno", shading="auto")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="seasonal range (K)")
    ax.contour(lon, lat, mask, levels=[0.5], colors="w", linewidths=0.7, alpha=0.8)
    for k, col in zip((interior, coast, ocean), sample_colors):
        ax.plot(lon[k], lat[i], "o", mfc=col, mec="w", mew=1.2, ms=8)
    ax.set_title("seasonal range — continental interiors blaze, oceans stay maritime", fontsize=9)
    ax.set_ylabel("latitude (°)")

    # -- the annual-mean map (zonally flat — the NMS headline) ---------------- #
    ax = axd["amean"]
    vlim = np.max(np.abs(amean))
    im = ax.pcolormesh(lon, lat, amean, cmap="RdBu_r", vmin=-vlim, vmax=vlim, shading="auto")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="annual-mean T (°C)")
    ax.contour(lon, lat, mask, levels=[0.5], colors="k", linewidths=0.7, alpha=0.5)
    spread = float(np.max(amean.max(axis=1) - amean.min(axis=1)))
    ax.set_title(f"annual mean — zonally FLAT (east–west spread {spread:.2f} K): the mask is invisible",
                 fontsize=9)

    # -- the cross-section: range vs longitude at the band ------------------- #
    ax = axd["cross"]
    ax.plot(lon, rng[i], color="#333333", lw=2.0)
    ax.fill_between(lon, 0, rng[i].max() * 1.05, where=c.land_mask[i], color=SEASON_LAND,
                    alpha=0.12, step="mid", label="land")
    for k, col, name in zip((interior, coast, ocean), sample_colors, ("interior", "coast", "ocean")):
        ax.plot(lon[k], rng[i, k], "o", mfc=col, mec="k", mew=0.8, ms=9, label=name)
    ax.set_ylim(0, rng[i].max() * 1.05)
    ax.set_xlim(0, 360)
    ax.set_xlabel("longitude (°)"); ax.set_ylabel("seasonal range (K)")
    ax.set_title(f"at {abs(lat[i]):.0f}°{'N' if lat[i] >= 0 else 'S'} — the range peaks in the interior, "
                 f"drops over the sea", fontsize=9)
    ax.legend(fontsize=8, loc="upper right")

    # -- the seasonal cycle at the three points ------------------------------ #
    ax = axd["cycle"]
    for k, col, name in zip((interior, coast, ocean), sample_colors, ("interior", "coast", "ocean")):
        ax.plot(months, c.T[i, k], color=col, lw=2.0,
                label=f"{name} (range {rng[i, k]:.0f} K)")
    ax.set_xlim(0, 12); ax.set_xticks(range(0, 13, 2))
    ax.set_xlabel("month of year"); ax.set_ylabel("surface temperature (°C)")
    ax.set_title("the cycle — interior swings hard and prompt, ocean barely moves and lags", fontsize=9)
    ax.legend(fontsize=8, loc="upper left")

    fig.suptitle("Planet Rung 5B.2 — seasons meet continents: continentality becomes a map (NMS83)",
                 fontsize=13)
    return fig


# --------------------------------------------------------------------------- #
# Rung 5B.1+ — the seasonal ice-albedo feedback: a migrating ice edge, ice asymmetry, bistability.
# --------------------------------------------------------------------------- #
def seasonal_ice_figure(result):
    """The seasonal ice-albedo payoff: a migrating ice edge, continentality in ice (:mod:`planet.seasonal`).

    Four panels off the converged ice-albedo limit cycle. **The cycle** (top-left, a midlatitude band):
    the land- and ocean-tile temperatures over the year with the freezing isotherm drawn — the land tile
    dips below freezing for a wide slab of the year (shaded: seasonal land ice) while the ocean tile stays
    open. **The migrating ice edge** (top-right): a latitude–month Hovmöller of the land tile with the
    ``T_freeze`` contour — a white ice cap that grows in winter and retreats in summer, in both hemispheres
    half a year apart. **The ice asymmetry** (bottom-left): fraction of the year frozen vs latitude, land
    vs ocean — the land curve reaches far equatorward, the ocean barely freezes. **Bistability**
    (bottom-right): the warm-start (finite-ice) and cold-start (snowball) annual-mean profiles — one sun,
    two climates, the seed picks the branch (Phase 1's Snowball, inside the seasonal cycle). Requires the
    ``viz`` extra.
    """
    m, c = result.model, result.warm
    lat = c.latitude_deg()
    months = c.days / (365.25 / 12.0)
    i = result.band_index()
    Tf = T_FREEZE

    fig, axd = plt.subplot_mosaic(
        [["cycle", "hov"], ["frac", "bistab"]], figsize=(13.5, 8.5), constrained_layout=True,
    )

    # -- the cycle at one midlatitude band, with the freezing line + seasonal ice ---- #
    ax = axd["cycle"]
    ax.plot(months, c.T_land[i], color=SEASON_LAND, lw=2.2, label="land tile")
    ax.plot(months, c.T_ocean[i], color=SEASON_OCEAN, lw=2.2, label="ocean tile")
    ax.axhline(Tf, color="#4a4a4a", lw=1.0, ls="--", alpha=0.8)
    ax.fill_between(months, c.T_land[i], Tf, where=c.T_land[i] < Tf,
                    color=ICE_COLOR, alpha=0.8, label="land frozen")
    ax.annotate("freezing", (0.2, Tf), textcoords="offset points", xytext=(0, 3),
                fontsize=7, color="#4a4a4a")
    ax.set_ylabel("surface temperature (°C)")
    ax.set_title(f"the cycle at {abs(lat[i]):.0f}°{'N' if lat[i] >= 0 else 'S'} — land freezes "
                 f"{100 * c.ice_fraction('land')[i]:.0f}% of the year, ocean stays open", fontsize=9)
    ax.legend(fontsize=8, loc="lower right")

    # -- the migrating ice edge (Hovmöller, land tile, with the T_freeze contour) ---- #
    ax = axd["hov"]
    im = ax.pcolormesh(months, lat, c.T_land, cmap="RdBu_r", shading="auto",
                       vmin=-np.max(np.abs(c.T_land)), vmax=np.max(np.abs(c.T_land)))
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="land-tile T (°C)")
    ax.contour(months, lat, c.T_land, levels=[Tf], colors="w", linewidths=1.6)
    ax.axhline(0.0, color="k", lw=0.6, alpha=0.5)
    ax.set_title("the ice edge migrates — white cap grows in winter, retreats in summer", fontsize=9)
    ax.set_ylabel("latitude (°)")

    # -- the ice asymmetry: fraction of year frozen vs latitude ---------------------- #
    ax = axd["frac"]
    ax.plot(lat, 100 * c.ice_fraction("land"), color=SEASON_LAND, lw=2.2, label="land tile")
    ax.plot(lat, 100 * c.ice_fraction("ocean"), color=SEASON_OCEAN, lw=2.2, label="ocean tile")
    ax.fill_between(lat, 100 * c.ice_fraction("ocean"), 100 * c.ice_fraction("land"),
                    color=ICE_COLOR, alpha=0.25)
    ax.set_xlabel("latitude (°)"); ax.set_ylabel("% of year frozen")
    ax.set_ylim(-3, 103)
    ax.set_title("continentality as ice — land freezes far equatorward, ocean barely does", fontsize=9)
    ax.legend(fontsize=8, loc="upper center")

    # -- bistability: warm (finite-ice) vs cold (snowball) annual-mean profiles ------ #
    ax = axd["bistab"]
    ax.plot(lat, result.warm.annual_mean("mean"), color="#b5402e", lw=2.2,
            label=f"warm start (⟨T⟩ {result.warm.annual_mean('mean').mean():.0f} °C)")
    ax.plot(lat, result.snowball.annual_mean("mean"), color="#3a6ea5", lw=2.2,
            label=f"cold start (⟨T⟩ {result.snowball.annual_mean('mean').mean():.0f} °C)")
    ax.axhline(Tf, color="#4a4a4a", lw=1.0, ls="--", alpha=0.8)
    ax.set_xlabel("latitude (°)"); ax.set_ylabel("annual-mean T (°C)")
    ax.set_title("bistability — one sun, two climates (the seed picks the branch)", fontsize=9)
    ax.legend(fontsize=8, loc="center")

    for key in ("cycle", "hov"):
        axd[key].set_xlabel("month of year"); axd[key].set_xlim(0, 12)
        axd[key].set_xticks(range(0, 13, 2))

    fig.suptitle("Planet Rung 5B.1+ — seasonal ice-albedo: a migrating ice edge, continentality in ice",
                 fontsize=13)
    return fig


# --------------------------------------------------------------------------- #
# Rung 0+ — the complete equilibrium diagram: every branch, both folds, the small-ice-cap instability.
# --------------------------------------------------------------------------- #
STABLE_COLOR = "#1f4e79"         # stable branches (solid)
UNSTABLE_COLOR = "#b03a2e"       # unstable branches (dashed)
FOLD_COLOR = "#d4711f"           # the folds (turning points)


def bifurcation_figure(result):
    """The equilibrium diagram (:mod:`planet.bifurcation`): the S-curve, its folds, and what they cost Phase 1.

    Four panels. **The S-curve** (top-left): ice-line latitude vs solar constant for *every* equilibrium —
    stable branches solid, unstable dashed (the slope-stability theorem read off the curve), the ice-free
    and Snowball branches as the horizontal caps, the two folds marked, and the Phase-1 continuation sweep
    overlaid as dots: it rides the stable branches and jumps at exactly the folds. **The small-ice-cap
    zoom** (top-right): the narrow finite-cap window around today's sun — Earth's cap sits a couple of
    W/m² below the cliff where caps smaller than ``θ_c`` cannot exist. **θ_c vs transport** (bottom-left):
    a stronger transport needs a bigger cap, the window narrows, and past a critical ``D`` the finite-cap
    branch is gone. **The retired bias** (bottom-right): Phase 1's relaxed ice line vs its step, converging
    onto the exact curve. Requires the ``viz`` extra.
    """
    c, loop = result.curve, result.loop
    lat = c.latitude_deg()
    lo, hi = c.snowball_fold, c.small_ice_cap_fold
    pres = result.present

    fig, axd = plt.subplot_mosaic(
        [["scurve", "zoom"], ["theta", "bias"]], figsize=(13.5, 8.8), constrained_layout=True,
    )

    def draw_branches(ax, lw=2.2):
        # split the curve into runs of equal stability so solid/dashed segments join cleanly
        stable = c.stable
        start = 0
        for k in range(1, stable.size + 1):
            if k == stable.size or stable[k] != stable[start]:
                sl = slice(start, min(k + 1, stable.size))
                ax.plot(c.S0[sl], lat[sl], color=STABLE_COLOR if stable[start] else UNSTABLE_COLOR,
                        ls="-" if stable[start] else "--", lw=lw)
                start = k

    # -- the S-curve ------------------------------------------------------------ #
    ax = axd["scurve"]
    draw_branches(ax)
    ax.plot([c.ice_free_threshold_S0, loop.S0_up.max()], [90, 90], color=STABLE_COLOR, lw=2.2)
    ax.plot([loop.S0_up.min(), c.snowball_threshold_S0], [0, 0], color=STABLE_COLOR, lw=2.2)
    ax.plot(loop.S0_down, loop.iceline_down, "o", ms=3.2, color=COOLING_COLOR, alpha=0.8,
            label="Phase-1 sweep — dimming")
    ax.plot(loop.S0_up, loop.iceline_up, "o", ms=3.2, color=WARMING_COLOR, alpha=0.8,
            label="Phase-1 sweep — brightening")
    for f, name in ((lo, "Snowball fold"), (hi, "small-ice-cap fold")):
        ax.plot(f.S0, f.latitude_deg, "D", ms=8, color=FOLD_COLOR, mec="k", mew=0.8, zorder=5)
        ax.annotate(name, (f.S0, f.latitude_deg), textcoords="offset points",
                    xytext=(10, -14 if f is lo else 8), fontsize=8, color="#4a3000")
    ax.axvline(S0_EARTH, color=PRESENT_COLOR, lw=1.2, ls=":", alpha=0.9)
    ax.plot(S0_EARTH, pres.latitude_deg, "*", ms=13, color=PRESENT_COLOR, mec="k", mew=0.6, zorder=6,
            label=f"today: finite cap at {pres.latitude_deg:.0f}°")
    ax.plot([], [], color=STABLE_COLOR, lw=2.2, label="stable (dS₀/dx_s > 0)")
    ax.plot([], [], color=UNSTABLE_COLOR, lw=2.2, ls="--", label="unstable")
    ax.set_xlim(loop.S0_up.min(), loop.S0_up.max()); ax.set_ylim(-3, 93)
    ax.set_xlabel("solar constant S₀ (W m⁻²)"); ax.set_ylabel("ice-line latitude (°)")
    ax.set_title("every equilibrium — the S-curve; the sweep rides the stable branches and jumps at the folds",
                 fontsize=9)
    ax.legend(fontsize=7.5, loc="center right")

    # -- the small-ice-cap zoom ------------------------------------------------- #
    ax = axd["zoom"]
    draw_branches(ax, lw=2.4)
    ax.plot([c.ice_free_threshold_S0, hi.S0 + 25], [90, 90], color=STABLE_COLOR, lw=2.4)
    ax.axvspan(lo.S0, hi.S0, color=ICE_COLOR, alpha=0.25, label="finite-cap window")
    ax.plot(hi.S0, hi.latitude_deg, "D", ms=9, color=FOLD_COLOR, mec="k", mew=0.8, zorder=5)
    ax.axvline(S0_EARTH, color=PRESENT_COLOR, lw=1.2, ls=":")
    ax.plot(S0_EARTH, pres.latitude_deg, "*", ms=14, color=PRESENT_COLOR, mec="k", mew=0.6, zorder=6)
    ax.annotate(f"today: {hi.S0 - S0_EARTH:.1f} W/m² ({100 * (hi.S0 - S0_EARTH) / S0_EARTH:.2f} %)\n"
                f"below the cliff", (S0_EARTH, pres.latitude_deg), textcoords="offset points",
                xytext=(-120, -30), fontsize=8, color=PRESENT_COLOR,
                arrowprops=dict(arrowstyle="-", color=PRESENT_COLOR, lw=0.8))
    ax.annotate(f"θ_c = {hi.cap_radius_deg:.1f}°: no sun holds\na smaller cap", (hi.S0, hi.latitude_deg),
                textcoords="offset points", xytext=(-10, 18), fontsize=8, ha="right", color="#4a3000")
    ax.set_xlim(S0_EARTH - 40, hi.S0 + 20); ax.set_ylim(55, 92)
    ax.set_xlabel("solar constant S₀ (W m⁻²)"); ax.set_ylabel("ice-line latitude (°)")
    ax.set_title("the second cliff — a brightening cap shrinks to θ_c, then vanishes in a jump", fontsize=9)
    ax.legend(fontsize=8, loc="lower left")

    # -- θ_c vs D ----------------------------------------------------------------- #
    ax = axd["theta"]
    ok = ~np.isnan(result.theta_c)
    ax.plot(result.D_values[ok], result.theta_c[ok], "o-", color=FOLD_COLOR, lw=2.0, ms=6, label="critical cap radius θ_c")
    ax.set_xlabel("transport D (W m⁻² K⁻¹)"); ax.set_ylabel("θ_c (°)", color=FOLD_COLOR)
    ax2 = ax.twinx()
    width = result.window_hi - result.window_lo
    ax2.plot(result.D_values[ok], width[ok], "s--", color=STABLE_COLOR, lw=1.6, ms=5, label="finite-cap window width")
    ax2.set_ylabel("window width (W m⁻²)", color=STABLE_COLOR)
    if (~ok).any():
        Dgone = float(result.D_values[~ok].min())
        ax.axvspan(Dgone, result.D_values.max() * 1.02, color="#dddddd", alpha=0.6)
        ax.text(Dgone, 0.5 * np.nanmax(result.theta_c), "  no stable\n  finite cap", fontsize=8, va="center")
    ax.axvline(0.555, color=PRESENT_COLOR, lw=1.0, ls=":")
    ax.set_xlim(result.D_values.min() - 0.05, result.D_values.max() * 1.02)
    ax.set_title("θ_c ≈ 10° at weak transport, growing past D ≈ 0.4 — and past a critical D the branch is gone",
                 fontsize=9)
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper left")

    # -- the retired bias -------------------------------------------------------- #
    ax = axd["bias"]
    ax.semilogx(result.n_taus, result.relaxed_ice_line, "o-", color=COOLING_COLOR, lw=2.0, ms=6,
                label="Phase-1 relaxation (n_tau = Δt/τ_rad)")
    ax.axhline(pres.latitude_deg, color=PRESENT_COLOR, lw=2.0, label=f"exact curve: {pres.latitude_deg:.1f}°")
    ax.invert_xaxis()
    ax.set_xlabel("relaxation step Δt / τ_rad  (→ smaller)"); ax.set_ylabel("present-day ice line (°)")
    ax.set_title("the marcher's O(Δt) fixed-point bias, quantified and retired: no step, no bias", fontsize=9)
    ax.legend(fontsize=8, loc="lower right")

    fig.suptitle("Planet Rung 0+ — the complete equilibrium diagram: two folds, one narrow window for a polar cap",
                 fontsize=13)
    return fig


# --------------------------------------------------------------------------- #
# Rung 5B.3 — seasonal ice on the map: the snow map, and a mask the annual mean can see.
# --------------------------------------------------------------------------- #
_MONTHS = ("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")


def _ice_overlay(ax, lon, lat, frozen, mask):
    """Paint the frozen cells white on top of a temperature map, and draw the coastlines."""
    ice = np.ma.masked_where(~frozen, np.ones_like(frozen, dtype=float))
    ax.pcolormesh(lon, lat, ice, cmap=ListedColormap(["#f4f8fb"]), vmin=0, vmax=1, shading="auto")
    ax.contour(lon, lat, mask.astype(float), levels=[0.5], colors="k", linewidths=0.7, alpha=0.6)


def seasonal_ice_map_figure(result):
    """The Rung-5B.3 payoff: the seasonal ice MAP and the rectified annual mean (:mod:`planet.seasonal_map`).

    Four panels off the converged 2-D ice-albedo limit cycle. **January / July** (top): the temperature map
    with the frozen cells painted white — winter snow sprawls across the northern continents while the
    ocean at the same latitude stays open; six months later the snow is gone and the south wears it.
    **Fraction of the year frozen** (bottom-left): the seasonal-ice map — continental interiors freeze
    for months, the polar ocean is perennial, the midlatitude ocean never freezes. **The annual mean sees
    the mask** (bottom-right): the annual-mean anomaly from its zonal mean — exactly zero for the fixed-
    albedo 5B.2 map, now a cold imprint of every continent: winter snow reflects sun the ocean keeps
    absorbing, so the seasonal cycle is *rectified* into a colder land mean. Requires the ``viz`` extra.
    """
    c = result.climate
    lat, lon = c.latitude_deg(), c.longitude_deg()
    mask = c.land_mask
    steps = result.month_steps()
    vlim = float(np.max(np.abs(c.T)))
    i = result.band_index()
    interior, ocean = result.sample_columns()

    fig, axd = plt.subplot_mosaic(
        [["jan", "jul"], ["frac", "anom"]], figsize=(13.5, 8.6), constrained_layout=True,
    )
    for key, mi in (("jan", 0), ("jul", 6)):
        ax = axd[key]
        im = ax.pcolormesh(lon, lat, c.T[:, :, steps[mi]], cmap="RdBu_r", vmin=-vlim, vmax=vlim, shading="auto")
        _ice_overlay(ax, lon, lat, c.frozen(steps[mi]), mask)
        ax.set_title(f"{_MONTHS[mi]} — surface temperature, frozen cells painted white "
                     f"({100 * c.frozen(steps[mi]).mean():.0f}% of the globe iced)", fontsize=9)
        ax.set_ylabel("latitude (°)")
    fig.colorbar(im, ax=[axd["jan"], axd["jul"]], fraction=0.025, pad=0.02, label="surface T (°C)")

    ax = axd["frac"]
    frac = c.ice_fraction()
    im = ax.pcolormesh(lon, lat, 100 * frac, cmap="Blues", vmin=0, vmax=100, shading="auto")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="% of year frozen")
    ax.contour(lon, lat, mask.astype(float), levels=[0.5], colors="k", linewidths=0.7, alpha=0.6)
    for k, col in ((interior, SEASON_LAND), (ocean, SEASON_OCEAN)):
        ax.plot(lon[k], lat[i], "o", mfc=col, mec="w", mew=1.2, ms=8)
    ax.set_title(f"the seasonal-ice map — at {abs(lat[i]):.0f}°N the interior is frozen "
                 f"{100 * frac[i, interior]:.0f}% of the year, the ocean {100 * frac[i, ocean]:.0f}%", fontsize=9)
    ax.set_xlabel("longitude (°)"); ax.set_ylabel("latitude (°)")

    ax = axd["anom"]
    anom = c.zonal_anomaly()
    alim = float(np.max(np.abs(anom)))
    im = ax.pcolormesh(lon, lat, anom, cmap="RdBu_r", vmin=-alim, vmax=alim, shading="auto")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="annual-mean T − zonal mean (K)")
    ax.contour(lon, lat, mask.astype(float), levels=[0.5], colors="k", linewidths=0.7, alpha=0.6)
    ax.set_title(f"the annual mean now SEES the mask — interior {anom[i, interior]:+.1f} K, ocean "
                 f"{anom[i, ocean]:+.1f} K (was 0 with a fixed albedo)", fontsize=9)
    ax.set_xlabel("longitude (°)")

    fig.suptitle("Planet Rung 5B.3 — seasonal ice on the map: winter snow, sea ice, and a rectified annual mean",
                 fontsize=13)
    return fig


def seasonal_ice_map_animation(result, *, interval: int = 160):
    """The year, month by month: the temperature map with the snow/sea-ice cover painted white.

    A single panel animated over the 12 month-centres (colour range fixed across frames so the eye reads
    the real seasonal swing), the coastlines drawn, a month badge, and a running readout of the iced
    fraction of the globe. Returns a :class:`matplotlib.animation.FuncAnimation`; save with a Pillow writer
    (GIF, CI-safe). Requires the ``viz`` extra.
    """
    from matplotlib.animation import FuncAnimation

    c = result.climate
    lat, lon = c.latitude_deg(), c.longitude_deg()
    mask = c.land_mask
    steps = result.month_steps()
    vlim = float(np.max(np.abs(c.T)))
    ice_cmap = ListedColormap(["#f4f8fb"])

    fig, ax = plt.subplots(figsize=(9.5, 5.6), constrained_layout=True)
    fig.get_layout_engine().set(rect=(0.0, 0.14, 1.0, 0.86))
    mesh = ax.pcolormesh(lon, lat, c.T[:, :, steps[0]], cmap="RdBu_r", vmin=-vlim, vmax=vlim, shading="auto")
    fig.colorbar(mesh, ax=ax, fraction=0.035, pad=0.03, label="surface T (°C)")
    ice_mesh = ax.pcolormesh(lon, lat, np.ma.masked_where(~c.frozen(steps[0]), np.ones(mask.shape)),
                             cmap=ice_cmap, vmin=0, vmax=1, shading="auto")
    ax.contour(lon, lat, mask.astype(float), levels=[0.5], colors="k", linewidths=0.7, alpha=0.6)
    ax.set_xlabel("longitude (°)"); ax.set_ylabel("latitude (°)")
    badge = ax.text(0.015, 0.97, "", transform=ax.transAxes, fontsize=16, fontweight="bold", va="top",
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#888", alpha=0.9))
    readout = ax.text(0.985, 0.97, "", transform=ax.transAxes, fontsize=9, va="top", ha="right",
                      bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#888", alpha=0.9))
    ax.set_title("Rung 5B.3 — the seasonal ice map: frozen cells painted white, month by month", fontsize=10)
    fig.text(0.5, 0.045,
             "Every grid point freezes on its own temperature (the rung-0 step-function ice-albedo, re-frozen each\n"
             "half-step). Small-C continents grow winter snow the ocean at the same latitude never does; the polar\n"
             "ocean's sea ice lingers. Idealized blocky continents; a diffusive 2-D EBM, not a weather model.",
             ha="center", va="center", fontsize=8.5, color="#333333")

    def update(frame):
        s = steps[frame]
        mesh.set_array(c.T[:, :, s].ravel())
        frozen = c.frozen(s)
        ice_mesh.set_array(np.ma.masked_where(~frozen, np.ones(mask.shape)).ravel())
        badge.set_text(_MONTHS[frame])
        readout.set_text(f"iced: {100 * frozen.mean():.0f}% of globe · "
                         f"{100 * frozen[mask].mean():.0f}% of land · {100 * frozen[~mask].mean():.0f}% of ocean")
        return mesh, ice_mesh, badge, readout

    return FuncAnimation(fig, update, frames=len(steps), interval=interval, blit=False)


# --------------------------------------------------------------------------- #
# Rung 5B.4 — the seasonal small-ice-cap instability.
# --------------------------------------------------------------------------- #
def seasonal_sici_figure(result):
    """The seasons dissolving the small-ice-cap fold (:mod:`planet.seasonal_sici`) — and where it comes back.

    Four panels. **The forbidden band** (top-left): the annual-mean equilibrium diagram of this very model
    near the pole — ice-line latitude against solar constant, stable branches solid, unstable dashed — with
    the shaded band of cap sizes it says no sun can hold, and the *seasonal* sweep's perennial ice edge
    laid over it, walking straight through. **One cell at a time** (top-right): the same sweep read as the
    un-interpolated perennial ice-**cell** count, dimming and brightening on top of each other; a fold
    would show as a step of several cells and a visible gap between the two legs. **Bringing it back**
    (bottom-left): plant a cap of the critical size at the fold's sun and deepen the ocean mixed layer —
    the seasonal swing dies, and past some depth the planted cap survives where a warm start grows none,
    i.e. two climates at one sun. **The grid is not the answer** (bottom-right): the loop width at three
    resolutions against the parent's own loop — the verdict does not move as the polar cell shrinks.
    Requires the ``viz`` extra.
    """
    r = result
    parent, loop = r.parent, r.loop
    theta_c = r.theta_c

    fig, axd = plt.subplot_mosaic(
        [["band", "cells"], ["depth", "res"]], figsize=(13.5, 8.8), constrained_layout=True,
    )

    # -- the annual mean's forbidden band, and the seasonal sweep crossing it -------- #
    ax = axd["band"]
    lat = parent.latitude_deg()
    near_pole = lat >= 60.0
    stable = parent.stable
    labelled = set()                                        # each branch kind labelled once, not per run
    start = 0
    for k in range(1, stable.size + 1):                     # solid/dashed runs, as the rung-0 figure draws
        if k == stable.size or stable[k] != stable[start]:
            sl = slice(start, min(k + 1, stable.size))
            sel = near_pole[sl]
            if sel.any():
                kind = bool(stable[start])
                ax.plot(parent.S0[sl][sel], lat[sl][sel],
                        color=STABLE_COLOR if kind else UNSTABLE_COLOR,
                        ls="-" if kind else "--", lw=2.2,
                        label=None if kind in labelled else
                        ("annual mean — stable cap" if kind else "annual mean — unstable (no cap lives here)"))
                labelled.add(kind)
            start = k
    ax.axhspan(90.0 - theta_c, 90.0, color=UNSTABLE_COLOR, alpha=0.12)
    ax.annotate(f"caps smaller than θ_c = {theta_c:.1f}°\ncannot exist at equilibrium",
                (0.03, 0.93), xycoords="axes fraction", fontsize=8, va="top", color=UNSTABLE_COLOR)
    down, up = loop.down, loop.up
    for cont, colour, tag in ((down, COOLING_COLOR, "seasonal — dimming"),
                              (up, WARMING_COLOR, "seasonal — brightening")):
        has_ice = cont.perennial_cap_deg > 0.0
        ax.plot(cont.S0[has_ice], 90.0 - cont.perennial_cap_deg[has_ice], "o", ms=4.5,
                color=colour, alpha=0.85, label=tag)
    ax.set_xlabel("solar constant S₀ (W m⁻²)")
    ax.set_ylabel("edge of the year-round ice (° latitude)")
    ax.set_ylim(60, 91)
    ax.set_xlim(min(down.S0.min(), parent.ice_free_threshold_S0) - 12.0, down.S0.max() + 3.0)
    ax.set_title("the seasonal planet holds caps the annual mean forbids", fontsize=9)
    ax.legend(fontsize=7.5, loc="lower left")

    # -- the un-interpolated read: perennial cells, one at a time -------------------- #
    ax = axd["cells"]
    ax.step(down.S0, down.n_perennial_cells, where="mid", color=COOLING_COLOR, lw=2.2,
            label="dimming")
    ax.step(up.S0, up.n_perennial_cells, where="mid", color=WARMING_COLOR, lw=1.4, ls="--",
            label="brightening (retraces it)")
    ax.set_xlabel("solar constant S₀ (W m⁻²)")
    ax.set_ylabel("cells frozen all year (one hemisphere)")
    ax.set_title(f"the cap grows {down.max_cell_jump} cell"
                 f"{'s' if down.max_cell_jump != 1 else ''} per step — a fold would switch on "
                 f"{r.fold_cells} at once", fontsize=9)
    ax.annotate("the albedo feedback only ever flips whole cells,\n"
                "so this is the honest read — the interpolated\n"
                "cap radius could look smooth over a jump",
                (0.40, 0.97), xycoords="axes fraction", fontsize=7.5, va="top", color="#4a4a4a")
    ax.legend(fontsize=8, loc="lower left")

    # -- deepen the ocean and the bistability returns -------------------------------- #
    ax = axd["depth"]
    depths = np.array(r.depths, dtype=float)
    gaps = np.array([s.gap_deg for s in r.seed_points])
    swings = np.array([s.warm_polar_amplitude_K for s in r.seed_points])
    colours = [UNSTABLE_COLOR if s.bistable else STABLE_COLOR for s in r.seed_points]
    ax.bar(np.arange(depths.size), gaps, color=colours, alpha=0.85, width=0.55)
    ax.axhline(r.seed_points[0].threshold_deg, color="#4a4a4a", lw=1.0, ls="--")
    ax.annotate("detection threshold (½ polar cell)", (0.02, r.seed_points[0].threshold_deg),
                xycoords=("axes fraction", "data"), textcoords="offset points", xytext=(0, 3),
                fontsize=7.5, color="#4a4a4a")
    ax.set_xticks(np.arange(depths.size))
    ax.set_xticklabels([f"{d:.0f} m" for d in depths])
    ax.set_xlabel("ocean mixed-layer depth (the seasonal-swing knob)")
    ax.set_ylabel("difference between the two climates (°)")
    ax.set_title("damp the seasons and the fold returns — two climates at one sun", fontsize=9)
    twin = ax.twinx()
    swing_colour = "#7a5195"                                # distinct from either bar colour
    twin.plot(np.arange(depths.size), swings, "o-", color=swing_colour, lw=1.6, ms=5)
    twin.set_ylabel("seasonal swing at 80° (K)", color=swing_colour)
    twin.tick_params(axis="y", colors=swing_colour)
    twin.set_ylim(0, max(swings.max() * 1.35, 0.1))

    # -- the resolution check --------------------------------------------------------- #
    ax = axd["res"]
    idx = np.arange(len(r.res_n_cells))
    ax.bar(idx, r.res_loop_width, color=STABLE_COLOR, alpha=0.85, width=0.55,
           label="seasonal loop width")
    ax.axhline(r.parent_loop_width, color=UNSTABLE_COLOR, lw=1.8, ls="--",
               label=f"the annual mean's own loop ({r.parent_loop_width:.1f} W m⁻²)")
    for i, (j, fc) in enumerate(zip(r.res_max_jump, r.res_fold_cells)):
        ax.annotate(f"{j}/step (a fold: {fc})", (i, 0.06 * max(r.parent_loop_width, 1.0)),
                    ha="center", fontsize=7.5, color="#333333")
    ax.set_xticks(idx)
    ax.set_xticklabels([f"{n} cells\n(polar cell {pc:.1f}°)"
                        for n, pc in zip(r.res_n_cells, r.res_polar_cell)], fontsize=8)
    for i, w in enumerate(r.res_loop_width):
        ax.annotate(f"{w:.1f}", (i, w), textcoords="offset points", xytext=(0, 4),
                    ha="center", fontsize=9, fontweight="bold", color=STABLE_COLOR)
    ax.set_ylabel("hysteresis loop width (W m⁻²)")
    ax.set_ylim(0, max(r.parent_loop_width * 1.35, 1.0))
    ax.set_title("not a grid artifact — the verdict holds as the polar cell shrinks", fontsize=9)
    ax.legend(fontsize=7.5, loc="upper right")

    fig.suptitle("Planet Rung 5B.4 — the seasons dissolve the small-ice-cap instability "
                 "(and a deep ocean brings it back)", fontsize=13)
    return fig


# --------------------------------------------------------------------------- #
# Rung 5A.4 — the mountain is cold: elevation → temperature → alpine biomes.
# --------------------------------------------------------------------------- #
SEA_LEVEL_COLOR = "#c0863d"      # the uncorrected zonal-mean (sea-level) temperature
TERRAIN_T_COLOR = "#1f4e79"      # the terrain-cooled surface temperature (what gets classified)
CONSTANT_G_COLOR = "#1f4e79"     # the pinned 6.5 K/km constant (the default that survived)
MOIST_G_COLOR = "#b03a2e"        # the emergent moist adiabat (the diagnostic that did not)
OBSERVED_BAND_COLOR = "#1f6f4a"  # the observed tropical freezing-level band


def alpine_biomes_figure(scene, sea_level_scene, diagnostic):
    """The Rung-5A.4 story: the terrain cools its own air, and the crest turns alpine.

    Four panels. **transect** — a mid-patch cross-section of the sea-level zonal-mean temperature
    (flat) against the terrain-cooled surface temperature (which dives over the ridge), with the
    Whittaker cold thresholds drawn so the crossing is visible over the drawn terrain: this is the
    whole mechanism in one line. **before** / **after** — the biome map under the *same* orographic
    rainfall without and with the cooling, the changed cells outlined, so the alpine cap is
    attributable to the lapse rate alone. **rates** — the verdict: the effective rate the emergent
    moist adiabat realises versus the pinned 6.5 K/km constant across latitude, and (twin axis) the
    freezing level each implies against the observed tropical band — the benchmark the constant lands just
    below and the moist adiabat overshoots by ~45 %, an ordering that is why the constant stays the default.

    Consumes an :class:`~planet.orographic_scene.OrographicScene` built with ``lapse=True``, the same
    scene built without it, and a :class:`~planet.elevation_temperature.LapseDiagnostic`. Requires the
    ``viz`` extra.
    """
    from .elevation_temperature import OBSERVED_TROPICAL_FREEZING_LEVEL_M
    from .biomes import BOREAL_MAX_C, TUNDRA_MAX_C

    lon, lat = scene.lon_deg, scene.lat_deg
    mid = scene.temperature_C.shape[0] // 2

    fig, axd = plt.subplot_mosaic(
        [["transect", "rates"], ["before", "after"]], figsize=(13, 9), constrained_layout=True,
    )

    # 1. The mechanism: the temperature transect over the ridge, against the cold thresholds.
    ax = axd["transect"]
    ax.plot(lon, scene.sea_level_temperature_C[mid, :], color=SEA_LEVEL_COLOR, lw=1.8, ls="--",
            label="sea-level zonal mean (rung 0)")
    ax.plot(lon, scene.temperature_C[mid, :], color=TERRAIN_T_COLOR, lw=2.2,
            label="terrain-cooled surface (rung 5A.4)")
    ax.axhline(0.0, color=FREEZE_COLOR, lw=1.0, ls=":")
    ax.annotate("freezing", xy=(lon[2], 0.0), fontsize=7, color=FREEZE_COLOR, va="bottom")
    for thr, name in ((BOREAL_MAX_C, "boreal"), (TUNDRA_MAX_C, "tundra")):
        ax.axhline(thr, color=ICELINE_COLOR, lw=0.9, ls="-.", alpha=0.7)
        ax.annotate(f"{name} threshold", xy=(lon[2], thr), fontsize=7, color=ICELINE_COLOR, va="bottom")
    ax.set_ylabel("surface temperature (°C)")
    ax.set_title(f"the crest cools {scene.elevation_cooling_K.max():.0f} K and crosses the cold thresholds",
                 fontsize=9)
    ax.legend(fontsize=8, loc="lower left")
    axe = ax.twinx()
    axe.fill_between(lon, 0, scene.elevation_m[mid, :], color="#c8b28a", alpha=0.45, zorder=0)
    axe.set_ylabel("elevation (m)", color="#8a6d3b")

    # 2. The verdict: the two lapse rates, and the freezing-level benchmark that chose between them.
    ax = axd["rates"]
    d = diagnostic
    ax.plot(d.latitude_deg, 1e3 * d.gamma_moist, color=MOIST_G_COLOR, lw=2.0,
            label="emergent moist adiabat (effective)")
    ax.plot(d.latitude_deg, 1e3 * d.gamma_constant, color=CONSTANT_G_COLOR, lw=2.0, ls="--",
            label="pinned constant Γ (the default)")
    ax.set_xlabel("latitude (°)"); ax.set_ylabel("lapse rate over the ridge (K/km)")
    ax.set_title("the emergent rate CONFIRMS the constant at mid-latitudes, diverges elsewhere", fontsize=9)
    ax.legend(fontsize=8, loc="upper left")
    axf = ax.twinx()
    axf.plot(d.latitude_deg, d.freezing_constant_m / 1e3, color=CONSTANT_G_COLOR, lw=1.2, ls=":")
    axf.plot(d.latitude_deg, d.freezing_moist_m / 1e3, color=MOIST_G_COLOR, lw=1.2, ls=":")
    lo, hi = OBSERVED_TROPICAL_FREEZING_LEVEL_M
    axf.axhspan(lo / 1e3, hi / 1e3, color=OBSERVED_BAND_COLOR, alpha=0.13, zorder=0)
    axf.scatter([d.latitude_deg[0]] * 2, [d.freezing_constant_m[0] / 1e3, d.freezing_moist_m[0] / 1e3],
                s=26, color=[CONSTANT_G_COLOR, MOIST_G_COLOR], zorder=4)
    axf.annotate(f"observed tropical freezing level {lo / 1e3:.1f}–{hi / 1e3:.1f} km\n"
                 "the constant lands just below it; the moist adiabat overshoots ~45 %",
                 xy=(d.latitude_deg[-1], hi / 1e3), xytext=(0.28, 0.545), textcoords="axes fraction",
                 fontsize=7, color=OBSERVED_BAND_COLOR, ha="left", va="top")
    axf.set_ylabel("freezing level (km, dotted)", color="#555555")

    # 3+4. The payoff: the same rain, classified without and with the terrain cooling.
    cmap, norm = _biome_cmap_norm()
    for key, sc_, title in (("before", sea_level_scene, "same rain, NO terrain cooling (rung 5A.3)"),
                            ("after", scene, "with terrain cooling — the alpine cap (rung 5A.4)")):
        ax = axd[key]
        ax.pcolormesh(lon, lat, sc_.biome_codes, cmap=cmap, norm=norm, shading="auto")
        ax.set_xlabel("longitude (°)"); ax.set_ylabel("latitude (°)")
        ax.set_title(title, fontsize=9)
        _biome_legend_inset(ax, sc_.biome_codes)
    changed = scene.biome_codes != scene.sea_level_biome_codes
    axd["after"].contour(lon, lat, changed.astype(float), levels=[0.5], colors="#c0392b", linewidths=1.2)
    axd["after"].set_title(f"with terrain cooling — {100 * scene.alpine_fraction:.0f}% re-classified "
                           "by the lapse rate alone", fontsize=9)

    fig.suptitle("Planet Rung 5A.4 — the mountain is COLD, not only wet: elevation → temperature → alpine biomes",
                 fontsize=13)
    return fig
