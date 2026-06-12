"""Rung B — the emergent eddy life cycle on the **globe** (the Plotly-globe animation, §9.5).

Rung A (:mod:`planet.plots.eddy_life_animation`) animated the released eddy life cycle as a flat
two-panel matplotlib movie. Rung B lifts the *same* banked frames (:class:`planet.eddy_flux.EddyFrames`
— the diagnostic-pure ``n_frames`` side-channel) onto the **existing Plotly globe**
(:mod:`planet.planetmap`'s unit sphere): the tracer ``θ`` stirred by the released jet, painted on the
planet and animated with Plotly's native play/slider. *No new stack* — it reuses the ``[webviz]``
Plotly the biome map already depends on (plan §9.5).

The two honesty edges, carried geometrically — not just in a caption
------------------------------------------------------------------
The flow is a **doubly-periodic midlatitude β-plane band patch**, not a global field, and its
instantaneous flux is **~90 % reversible** (``irreversible_fraction ~0.1``). A naive globe — wrapping
the periodic channel around all 360° of longitude and/or mirroring it into both hemispheres — would
fabricate a *planet-wide circulation* the model never produced. So Rung B is deliberately constrained:

* **One honest band, at its true width.** The channel's ``x`` is a Cartesian (m) periodic coordinate;
  its physical width ``Lx`` subtends only ``Δlon = Lx/(a·cosφ_c) ≈ 55°`` of longitude — *not* a full
  wrap. The band is laid as a **bounded longitude sector** at the channel's real latitudes
  (:func:`_band_geometry`), a single patch on an otherwise-bare base sphere; it is **not** mirrored to
  the southern hemisphere (the eddy field is the project's only longitudinally-structured field, so the
  :func:`planet.planetmap.circulation_layer` two-hemisphere broadcast — valid only for a *zonal-mean*
  jet — does not transfer). Earth's radius ``a`` is recovered from the frames' own ``(y, φ)`` pair (the
  linear β-plane metric ``φ = φ_ref + deg((y−y_ref)/a)``), so no new constant or schema field is added.
* **The flux indicator stays.** Beside the globe, the same cumulative-transport panel Rung A carries —
  the **throughput** ``Σ∫|F̄|dt`` raging upward while the **net** ``Σ|∫F̄dt|`` stays a small fraction —
  so the swirls-rage-but-net-barely-moves finding is on screen, not just narrated. Without it a
  streaming-particle globe would silently overclaim "ocean currents carrying heat" (ADR 0002 §5:
  visualize the mechanism, not the output).

Layering (ADR 0004 / the headless discipline)
----------------------------------------------
This module is **plotly-free at import** — Plotly is imported *lazily* inside :func:`eddy_globe_figure`
/ :func:`save_eddy_globe_html`, exactly like :mod:`planet.planetmap`. It reuses
:func:`planet.planetmap._sphere_xyz` (the globe parametrization) and adds *reach, not correctness*: the
physics is sealed in :mod:`planet.eddy_flux`; this module's only test is an execution smoke-test plus a
geometry pin (the band is a bounded midlat sector, not a pole-to-pole / 360° global field).

Run headless (saves the HTML globe):

    python -m planet.demo_eddy_globe
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .planetmap import _sphere_xyz

# Colours mirrored from the Rung-A two-panel (plots.py) — kept LOCAL so this module pulls no matplotlib.
THROUGHPUT_COLOR = "#b3361f"    # the |F̄| throughput — the raging swirls
NET_COLOR = "#1f6f4a"           # the net ∫F̄ transport — the small down-gradient residual
WINDOW_COLOR = "#7a7a7a"        # the κ-diagnosis window onset
SATURATION_COLOR = "#6a4c93"    # the eddy-KE saturation time
BASE_SPHERE_COLOR = "#dce3ea"   # the bare planet the eddy band sits on
BAND_RAISE = 1.012              # the eddy band radius — just above the base sphere (avoid z-fighting)


def _band_geometry(fr):
    """Map the doubly-periodic β-plane channel onto a TRUE-WIDTH longitude sector (NOT a 360° wrap).

    Returns ``((x, y, z), lat2d, lon2d, lon1d)``: the raised unit-sphere mesh of the eddy band (each
    ``(ny, nx)``) plus the latitude/longitude meshes and the 1-D sector longitudes (deg). Latitude is
    the channel's *own* ``fr.phi`` (a single midlatitude band — never mirrored to the other hemisphere);
    longitude is the channel ``x`` rescaled by the spherical metric ``Δlon = Δx/(a·cosφ_c)`` and centred
    on ``lon = 0``, so the patch occupies only its honest angular width (~55°), not the whole planet.
    Earth's radius ``a`` is recovered from the frames' linear β-plane metric — no new constant needed.
    """
    a = (fr.y[1] - fr.y[0]) / np.radians(fr.phi[1] - fr.phi[0])          # Earth radius from φ = φ_ref + deg((y−y_ref)/a)
    phi_c = float(fr.phi.mean())
    lon1d = np.degrees((fr.x - fr.x.mean()) / (a * np.cos(np.radians(phi_c))))   # centred sector (deg)
    lon2d, lat2d = np.meshgrid(lon1d, fr.phi)                            # (ny, nx): rows = lat, cols = lon
    x, y, z = _sphere_xyz(lat2d, lon2d)
    return (BAND_RAISE * x, BAND_RAISE * y, BAND_RAISE * z), lat2d, lon2d, lon1d


def _base_sphere(go):
    """A bare, uniformly-shaded planet — the globe the one honest eddy band is a *patch* on."""
    lon = np.linspace(-180.0, 180.0, 121)
    lat = np.linspace(-90.0, 90.0, 61)
    LON, LAT = np.meshgrid(lon, lat)
    x, y, z = _sphere_xyz(LAT, LON)
    return go.Surface(
        x=x, y=y, z=z, surfacecolor=np.zeros_like(x),
        colorscale=[[0.0, BASE_SPHERE_COLOR], [1.0, BASE_SPHERE_COLOR]],
        showscale=False, hoverinfo="skip", name="planet",
        lighting=dict(ambient=0.85, diffuse=0.3, specular=0.05),
    )


def eddy_globe_figure(eddy, *, frame_ms: int = 120):
    """The banked Rung-B artifact: the emergent eddy life cycle on the globe — globe + flux indicator.

    Left — **the stirring, on the planet**: the passive tracer ``θ`` advected by the released,
    barotropically-unstable jet, painted as a **bounded midlatitude longitude band** (its true ~55°
    width, a single hemisphere — *not* a planet-wide flow) on an otherwise-bare globe, animated through
    the release. Right — **the transport budget** (the same panel Rung A carries): the cumulative
    meridional eddy heat flux, the **throughput** ``Σ∫|F̄|dt`` climbing steeply while the **net**
    ``Σ|∫F̄dt|`` stays a small fraction — so :mod:`planet.eddy_flux`'s headline finding, that the
    instantaneous flux is **~90 % reversible**, is *visible*: the band streams, the net barely moves.
    A moving cursor ties the two together; a dashed line marks where κ is diagnosed (the window onset).

    Returns a :class:`plotly.graph_objects.Figure` with native play/slider controls. Requires the
    optional ``[webviz]`` extra (Plotly) **and** the frames side-channel — raises :class:`ValueError`
    if ``eddy.frames`` is ``None`` (recompute with ``eddy_life_cycle(..., n_frames=N)``).
    """
    fr = eddy.frames
    if fr is None:
        raise ValueError("eddy.frames is None — recompute with eddy_life_cycle(..., n_frames=N)")
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    (bx, by, bz), lat2d, lon2d, _lon1d = _band_geometry(fr)
    tmin, tmax = float(fr.theta.min()), float(fr.theta.max())
    hover = np.vectorize(lambda la, lo: f"lat {la:.0f}°, lon {lo:.0f}°")(lat2d, lon2d)

    def band_surface(theta):
        return go.Surface(
            x=bx, y=by, z=bz, surfacecolor=theta, colorscale="RdBu_r", cmin=tmin, cmax=tmax,
            colorbar=dict(title="θ (°C)", len=0.55, x=0.02, xanchor="left", thickness=14),
            text=hover, hoverinfo="text", name="eddy θ",
            lighting=dict(ambient=0.9, diffuse=0.25, specular=0.05),
        )

    fig = make_subplots(rows=1, cols=2, column_widths=[0.6, 0.4],
                        specs=[[{"type": "scene"}, {"type": "xy"}]])

    # --- static + animated traces; record the indices the frames will update -------------------- #
    fig.add_trace(_base_sphere(go), row=1, col=1)                 # idx 0 — the bare planet (static)
    fig.add_trace(band_surface(fr.theta[0]), row=1, col=1)       # idx 1 — the eddy band (ANIMATED: surfacecolor)
    band_i = 1

    # The transport budget, normalized to a FRACTION of total throughput so the panel reads 0..1 — the
    # raw cumulative is ~3e7 K·m (a meaningless magnitude that collapses the axis); the RATIO is the
    # story: throughput climbs to 1.0 while net plateaus at the irreversible fraction (~0.1). κ is
    # diagnosed elsewhere, so this rescale is legibility-only (ADR 0002: the panel is "reach not correctness").
    total = float(fr.thru_cum[-1]) or 1.0
    thru_frac = fr.thru_cum / total
    net_frac = fr.net_cum / total
    ymax = 1.08
    fig.add_trace(go.Scatter(x=fr.times, y=thru_frac, mode="lines",
                             line=dict(color=THROUGHPUT_COLOR, width=2.5),
                             name="throughput  Σ∫|F̄|dt  (the swirls rage)"), row=1, col=2)
    fig.add_trace(go.Scatter(x=fr.times, y=net_frac, mode="lines",
                             line=dict(color=NET_COLOR, width=2.5),
                             name="net  Σ|∫F̄dt|  (the small residual)"), row=1, col=2)
    fig.add_trace(go.Scatter(x=[fr.window_start, fr.window_start], y=[0.0, ymax], mode="lines",
                             line=dict(color=WINDOW_COLOR, dash="dash", width=1.4),
                             showlegend=False, hoverinfo="skip"), row=1, col=2)
    if eddy.saturation_period > fr.times[0]:
        fig.add_trace(go.Scatter(x=[eddy.saturation_period, eddy.saturation_period], y=[0.0, ymax],
                                 mode="lines", line=dict(color=SATURATION_COLOR, dash="dot", width=1.4),
                                 showlegend=False, hoverinfo="skip"), row=1, col=2)
    cursor_i = len(fig.data)
    fig.add_trace(go.Scatter(x=[fr.times[0], fr.times[0]], y=[0.0, ymax], mode="lines",
                             line=dict(color="#222222", width=1.2),
                             showlegend=False, hoverinfo="skip"), row=1, col=2)
    thru_i = len(fig.data)
    fig.add_trace(go.Scatter(x=[fr.times[0]], y=[thru_frac[0]], mode="markers",
                             marker=dict(color=THROUGHPUT_COLOR, size=9),
                             showlegend=False, hoverinfo="skip"), row=1, col=2)
    net_i = len(fig.data)
    fig.add_trace(go.Scatter(x=[fr.times[0]], y=[net_frac[0]], mode="markers",
                             marker=dict(color=NET_COLOR, size=9),
                             showlegend=False, hoverinfo="skip"), row=1, col=2)

    # --- the frames: update ONLY the band's surfacecolor + the cursor/markers (x/y/z stay static) - #
    frames = []
    for k in range(fr.times.size):
        t = float(fr.times[k])
        frames.append(go.Frame(
            name=str(k),
            # surfacecolor-only keeps the HTML lean (x/y/z are static, merged from the base trace);
            # cmin/cmax are re-stated per frame so the fixed colour range can't autoscale on a merge.
            data=[go.Surface(surfacecolor=fr.theta[k], cmin=tmin, cmax=tmax),
                  go.Scatter(x=[t, t], y=[0.0, ymax]),
                  go.Scatter(x=[t], y=[thru_frac[k]]),
                  go.Scatter(x=[t], y=[net_frac[k]])],
            traces=[band_i, cursor_i, thru_i, net_i],
        ))
    fig.frames = frames

    # --- play / pause / slider (Plotly-native; redraw=True so the 3-D surface repaints) ---------- #
    play = dict(type="buttons", direction="left", showactive=False, x=0.0, xanchor="left",
                y=0.0, yanchor="top", pad=dict(t=8, r=8),
                buttons=[dict(label="▶ play", method="animate",
                              args=[None, dict(frame=dict(duration=frame_ms, redraw=True),
                                               fromcurrent=True, transition=dict(duration=0))]),
                         dict(label="❚❚ pause", method="animate",
                              args=[[None], dict(frame=dict(duration=0, redraw=False),
                                                 mode="immediate", transition=dict(duration=0))])])
    slider = dict(active=0, x=0.08, len=0.5, y=0.0, yanchor="top", pad=dict(t=40),
                  currentvalue=dict(prefix="release t = ", suffix=" periods", font=dict(size=12)),
                  steps=[dict(method="animate", label=f"{fr.times[k]:.0f}",
                              args=[[str(k)], dict(mode="immediate", transition=dict(duration=0),
                                                   frame=dict(duration=0, redraw=True))])
                         for k in range(fr.times.size)])

    irr = eddy.irreversible_fraction
    # showspikes=False kills the hover crosshair/projection lines that track the pointer (3-D scene spikes
    # default ON) — a standing preference: no pointer-following spike lines on any visualization.
    no_axis = dict(showbackground=False, showticklabels=False, showgrid=False, zeroline=False,
                   visible=False, showspikes=False)
    fig.update_layout(
        title=dict(
            text="Planet §9.5 — the emergent eddy life cycle, on the globe<br>"
                 "<sub>one midlatitude β-plane band (a ~55° patch — NOT a planet-wide flow); the "
                 f"instantaneous flux is ~{round(100 * (1 - irr)):.0f}% reversible — the band streams, "
                 "net transport is only the small κ residual</sub>",
            x=0.5, xanchor="center", font=dict(size=15)),
        width=1180, height=640, margin=dict(l=0, r=0, t=88, b=10),
        scene=dict(xaxis=no_axis, yaxis=no_axis, zaxis=no_axis, aspectmode="data",
                   camera=dict(eye=dict(x=1.5, y=0.5, z=0.95))),
        updatemenus=[play], sliders=[slider],
        # legend seated in the panel's empty top-left (throughput is still near zero there early on),
        # clear of the title block above and the rising curves to its right.
        legend=dict(x=0.635, y=0.80, xanchor="left", yanchor="top", font=dict(size=10),
                    bgcolor="rgba(255,255,255,0.7)", bordercolor="#cccccc", borderwidth=1),
    )
    fig.update_xaxes(title_text="release time (inertial periods)", showspikes=False,
                     range=[float(fr.times[0]), float(fr.times[-1])], row=1, col=2)
    fig.update_yaxes(title_text="cumulative transport (fraction of total throughput)", showspikes=False,
                     range=[0.0, ymax], row=1, col=2)
    # row/col resolves to the panel's real x/y axes — NOT a hardcoded "x2"/"y2", which would spawn a
    # phantom overlaid axis (the bug that produced the double x-axis + the collapsed 33.2M range).
    fig.add_annotation(x=fr.window_start, y=ymax, text="κ diagnosed",
                       textangle=-90, xanchor="left", yanchor="top", showarrow=False,
                       font=dict(size=9, color=WINDOW_COLOR), row=1, col=2)
    return fig


def save_eddy_globe_html(eddy, path, *, frame_ms: int = 120) -> Path:
    """Render the Rung-B globe and write a standalone, viewable HTML animation (Plotly only).

    The interactive analogue of the Rung-A GIF: a self-contained HTML globe with native play/slider,
    needing no server and no dependency beyond Plotly (CDN ``plotly.js``). Does **not** autoplay — the
    viewer presses play. Returns the written path.
    """
    fig = eddy_globe_figure(eddy, frame_ms=frame_ms)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(path), include_plotlyjs="cdn", auto_play=False)
    return path
