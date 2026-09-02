"""The seasonal-ice **globe** — the 2-D ice-albedo limit cycle on a rotatable sphere with a month slider (rung 5B.3).

The interactive twin of :func:`planet.plots.seasonal_ice_map_animation`: the same twelve month-centred
temperature maps, the same frozen cells painted white, but on the Plotly globe :mod:`planet.planetmap`
renders — rotate it, zoom it, drag the month slider or press play and watch the winter snow sweep across
the northern continents while the southern summer bares the land, then trade places six months later.

A pure **render** layer (ADR 0002): it consumes the plain arrays of a converged
:class:`~planet.seasonal_map.SeasonalMapClimate` and paints them; nothing here is evidence. Reuses the
globe geometry helpers of :mod:`planet.planetmap` (``_sphere_xyz``, ``_polecapped``) rather than a second
sphere parametrisation. The frozen cells are painted by pinning them to the **top** of the colour scale,
whose last stop is white — one ``surfacecolor`` array per frame, so the sphere geometry is built once and
each frame only swaps the colours (small, fast). Requires the optional ``[webviz]`` extra (Plotly).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .ebm import T_FREEZE
from .planetmap import _sphere_xyz, _polecapped

_MONTHS = ("January", "February", "March", "April", "May", "June",
           "July", "August", "September", "October", "November", "December")

# A diverging blue→white→red temperature scale whose LAST stop is ice-white: frozen cells are pinned to
# the top value (a sentinel above the warmest temperature), so they read as snow/sea ice on the same
# surface, no second trace needed. The temperature range ±vlim occupies the central ±_ICE_SENTINEL_FRACTION
# of the colour axis; the band above it is the ice colour.
_ICE_SENTINEL_FRACTION = 0.90
_TEMPERATURE_STOPS = [(-1.0, "#1c3f80"), (-0.5, "#5b9bd5"), (-0.15, "#cfe3f3"), (0.0, "#f7f7f7"),
                      (0.15, "#f9d0c4"), (0.5, "#e47c6c"), (1.0, "#b2182b")]


def _colorscale():
    f = 0.5 * _ICE_SENTINEL_FRACTION
    scale = [[0.5 + f * t, col] for t, col in _TEMPERATURE_STOPS]
    scale[0][0] = 0.0                                                   # the bottom pads with the coldest colour
    top = 0.5 + f
    return scale + [[top + 1e-4, "#f4f8fb"], [1.0, "#f4f8fb"]]


def _wrap_lon(lon_deg: np.ndarray, *fields: np.ndarray):
    """Close the sphere at the longitude seam: append the first column at ``lon + 360°`` to each field."""
    lon_ext = np.concatenate([lon_deg, [lon_deg[0] + 360.0]])
    return (lon_ext, *(np.hstack([np.asarray(f), np.asarray(f)[:, :1]]) for f in fields))


def month_steps(n_steps: int) -> np.ndarray:
    """Time-sample index of the middle of each of the 12 months (endpoint-excluded year)."""
    return ((np.arange(12) + 0.5) * n_steps / 12.0).astype(int)


def seasonal_ice_globe(climate, *, T_freeze: float = T_FREEZE, caption: str | None = None):
    """Build the month-slider globe for a converged 2-D ice-albedo limit cycle (a Plotly figure).

    ``climate`` is a :class:`~planet.seasonal_map.SeasonalMapClimate` (ideally from an ice-albedo march —
    a fixed-albedo cycle renders too, its sub-freezing cells painted white *diagnostically*). One
    ``go.Surface`` on the unit sphere; twelve frames swapping ``surfacecolor``; a slider + play/pause. The
    hover text names the month, the surface (land/sea), the temperature, and whether the cell is frozen.
    """
    import plotly.graph_objects as go

    lat_c = climate.latitude_deg()
    lon_c = climate.longitude_deg()
    steps = month_steps(climate.days.size)
    vlim = float(np.max(np.abs(climate.T)))
    # the colour axis: [-vlim, +vlim] for temperature, and a sentinel band above +vlim for "frozen"
    cmax = vlim / _ICE_SENTINEL_FRACTION
    cmin = -cmax
    lat, mask = _polecapped(lat_c, climate.land_mask.astype(float))
    lon, mask = _wrap_lon(lon_c, mask)                       # close the seam (cell-centred grid)
    LON, LAT = np.meshgrid(lon, lat)
    X, Y, Z = _sphere_xyz(LAT, LON)
    surface_name = np.where(mask > 0.5, "land", "sea")

    def frame_arrays(k: int):
        s = steps[k]
        T = climate.T[:, :, s]
        frozen = T < T_freeze
        colour = np.where(frozen, cmax, np.clip(T, -vlim, vlim))
        lat_p, colour_p, T_p, frozen_p = _polecapped(lat_c, colour, T, frozen.astype(float))
        _, colour_p, T_p, frozen_p = _wrap_lon(lon_c, colour_p, T_p, frozen_p)
        hover = np.empty(colour_p.shape, dtype=object)
        for i in range(colour_p.shape[0]):
            for j in range(colour_p.shape[1]):
                hover[i, j] = (f"{_MONTHS[k]} · lat {LAT[i, j]:.0f}°, lon {LON[i, j]:.0f}° · {surface_name[i, j]}"
                               f"<br>T = {T_p[i, j]:.1f} °C" + (" · frozen" if frozen_p[i, j] > 0.5 else ""))
        return colour_p, hover, float(frozen.mean())

    tick_step = 10.0 if vlim > 25.0 else 5.0
    ticks = [t for t in np.arange(-np.floor(vlim / tick_step) * tick_step, vlim + 1e-9, tick_step)]
    colorbar = dict(title="surface T (°C)", len=0.7, tickvals=ticks + [0.5 * (vlim + cmax)],
                    ticktext=[f"{t:.0f}" for t in ticks] + ["frozen"])
    contours = dict(x=dict(highlight=False), y=dict(highlight=False), z=dict(highlight=False))
    colour0, hover0, iced0 = frame_arrays(0)
    base = go.Surface(x=X, y=Y, z=Z, surfacecolor=colour0, text=hover0, hoverinfo="text",
                      colorscale=_colorscale(), cmin=cmin, cmax=cmax, colorbar=colorbar, contours=contours)

    frames, slider_steps = [], []
    for k in range(12):
        colour, hover, iced = frame_arrays(k)
        title = f"{_MONTHS[k]} — {100 * iced:.0f}% of the globe frozen"
        frames.append(go.Frame(name=_MONTHS[k][:3], data=[go.Surface(surfacecolor=colour, text=hover)],
                               layout=go.Layout(title=dict(text=f"Seasonal ice on the map — {title}"))))
        slider_steps.append(dict(method="animate", label=_MONTHS[k][:3],
                                 args=[[_MONTHS[k][:3]], dict(mode="immediate", frame=dict(duration=0, redraw=True),
                                                           transition=dict(duration=0))]))

    no_axis = dict(showbackground=False, showticklabels=False, showgrid=False, zeroline=False,
                   visible=False, showspikes=False)          # showspikes=False: the standing viz preference
    fig = go.Figure(data=[base], frames=frames)
    fig.update_layout(
        title=f"Seasonal ice on the map — {_MONTHS[0]} — {100 * iced0:.0f}% of the globe frozen",
        width=900, height=960, margin=dict(l=0, r=0, t=50, b=250),
        scene=dict(xaxis=no_axis, yaxis=no_axis, zaxis=no_axis, aspectmode="data",
                   camera=dict(eye=dict(x=1.4, y=1.4, z=0.9))),
        updatemenus=[dict(type="buttons", showactive=False, x=0.02, y=0.0, xanchor="left", yanchor="bottom",
                          buttons=[dict(label="▶ play", method="animate",
                                        args=[None, dict(frame=dict(duration=450, redraw=True), fromcurrent=True,
                                                         transition=dict(duration=0))]),
                                   dict(label="❚❚ pause", method="animate",
                                        args=[[None], dict(frame=dict(duration=0, redraw=False), mode="immediate",
                                                           transition=dict(duration=0))])])],
        sliders=[dict(active=0, currentvalue=dict(prefix="month: ", visible=True), pad=dict(t=10),
                      x=0.15, len=0.8, y=0.0, steps=slider_steps)],
    )
    if caption is None:
        caption = ("Rung 5B.3 — the 2-D seasonal energy-balance model with the ice-albedo feedback on every grid "
                   "point, marched to its annual limit cycle. Frozen cells (T below the freezing isotherm) are painted "
                   "white: winter snow sweeps over the small-heat-capacity continents while the ocean at the same "
                   "latitude stays open; the polar ocean's sea ice lingers all year. Drag the slider or press play. "
                   "Idealized blocky continents; a diffusive EBM, not a weather model — the physics is validated in "
                   "planet/tests, this globe only paints it.")
    from .planetmap import _wrap_html
    fig.add_annotation(xref="paper", yref="paper", x=0.5, y=-0.13, xanchor="center", yanchor="top",
                       showarrow=False, align="left", font=dict(size=13, color="#33373b"), text=_wrap_html(caption))
    return fig


def save_seasonal_ice_globe(climate, path, **kw) -> Path:
    """Render the month-slider globe and write a standalone HTML file (Plotly via CDN, no server)."""
    fig = seasonal_ice_globe(climate, **kw)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(path), include_plotlyjs="cdn")
    return path
