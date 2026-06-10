"""The deep-end interactive planet map — a layer registry painted by Plotly (Planet §9, ADR 0004).

Planet is the **first project in the program to reach the deep end** of ADR 0002 §4 — a payoff
that is inherently a *map*. Beyond the static matplotlib floor (:mod:`planet.plots`) and
the teaching notebook, this module is the **interactive planet map**: a rotatable globe whose
knob-sliders drive an instant recompute-and-remap. Its **first version is the biome map** — the
Phase-2 dramatic end-to-end win (:mod:`planet.demo_biomes`) made live.

It adds **reach, not correctness** (ADR 0002): every array it paints is produced by a model already
sealed behind its own validation triad (the EBM, the precip parameterization, the Whittaker
classifier). This module introduces **no new physics, no new constant** — it is a pure consumer of
:func:`planet.demo_biomes.compute`, the same way Steel's ``app.py`` is a pure re-composition
of its ``sweep`` harness. Its only *real* correctness property lives next door in
:mod:`planet.planet_spec` (the round-trip-identity of the exported state); the map's own
test is an **execution smoke-test** (ADR 0002 §2 / ADR 0004 #4).

The layer registry (ADR 0004 #1 — the structural heart)
-------------------------------------------------------
A view is an **ordered stack of self-contained layers**, each a :class:`Layer`
``(name, kind, data, units, style, z_order, inert)`` with ``kind`` one of :class:`LayerKind`
``{scalar_field, vector_overlay, annotation}``. The renderer (:func:`render`) is **generic over
kind**: a phase contributes to the map by *registering layers*, never by editing the renderer. v1
registers four climate layers — **temperature** & **precipitation** & **biome** (scalar fields) plus
the **ice line** (an annotation) — and the renderer paints all four with zero per-layer special
casing: that *is* the proof that "phases just register, the renderer is unchanged" (ADR 0004 #1).
Later phases register more — circulation streamlines / the jet axis (Phase 4, the
``vector_overlay`` kind), elevation & coastlines (the geography seam) — at the same seam.

``vector_overlay`` was a **declared-but-unpainted** kind through Phases 1–3 ("build the seam, not the
machinery" — ADR 0004): :func:`render` raised for it until there was a real flow to draw. **Phase 4
builds the machinery** (:func:`_vector_overlay_trace`): a Phase-4 coupled jet
(:func:`planet.coupler.couple_jet`) registers a ``circulation`` ``vector_overlay`` layer
(:func:`circulation_layer`), painted as flow arrows on the globe — the emergent jet *is* the consumer.
The jet is **not** part of the live-slider loop (it is a separate integration, the first compute too
heavy for the rung-0 instant remap; §9.2): a circulation view is *computed, then viewed*.

Three layers of code, by the ADR-0002 / Steel-``app.py`` discipline
-------------------------------------------------------------------
1. **The registry + the layer builders** (this module's top half) — plain functions over NumPy
   that turn a climate result into a :class:`PlanetView`. They import **neither Plotly nor
   ipywidgets**, so the module imports on a bare core install and the builders are unit-tested
   *always-green* (``tests/test_planetmap.py``), exactly like the ``sweep`` tests. This is also the
   half :mod:`planet.planet_spec` imports (the :class:`Layer` it serializes), so the
   interchange schema never pulls a render dependency.
2. **The renderer + the HTML artifact** (:func:`render`, :func:`save_html`) — Plotly imported
   **lazily** inside the function, smoke-tested with ``importorskip("plotly")`` (fast — no kernel,
   so **not** ``@slow``).
3. **:func:`interactive_map`** — the *only* place ``import ipywidgets`` lives, the live-slider loop
   for a notebook. It is the ``main()`` analogue: paper-thin, not unit-tested (the UI is reach).

The renderer-input seam: 2-D-ready now (ADR 0004 / plan §9.3)
------------------------------------------------------------
The renderer consumes a **2-D lat×lon field**. v1 is zonal-mean (the §3 scope edge), so each scalar
layer is the latitudinal profile **mirrored to the full globe** (the equator-to-pole EBM grid is one
hemisphere; the annual-mean climate is hemispherically symmetric) and **broadcast across longitude** —
so the globe paints honest **latitude bands**, not a premature 2-D field. This is exactly the shape a
longitudinally-varying field (the rung-5 2-D geography exit) will take, so that future change is *new
data*, not a renderer rewrite.

The interaction model is tier-dependent; the renderer is not (ADR 0004 #2)
--------------------------------------------------------------------------
Because the map only consumes arrays, it is invariant up the whole §5 GCM staircase — only the
*trigger* changes: at **rung 0** (here) compute is laptop-seconds, so a slider drives an instant
recompute-and-remap (:func:`interactive_map`); at the heavy upper rungs a live loop is impossible and
the trigger becomes *set parameters → launch a run → view the result*. The slider is a driver of
compute; the map is a consumer of arrays — only the former is tier-dependent.

The knobs (plan §9.1)
---------------------
The map wires the climate levers — **solar constant ``S₀``** (the Snowball lever), **CO₂** (a drop in
the OLR offset ``A`` — :mod:`planet.demo_biomes`'s warming knob), and **meridional transport
``D``** — plus the two **exoplanet** knobs (:mod:`planet.exoplanet`, now that each one's
source is pinned): the host **star's effective temperature ``T_star``** (a redder star lowers the ice
albedo → harder to snowball) and the **planet ``size``** in Earth radii (bigger → weaker transport → a
sharper gradient); and the **obliquity** (axial tilt → the insolation P₂ coefficient ``s₂`` —
:mod:`planet.obliquity`, now that its ``s₂(obliquity)`` relation is pinned and wired). All are
already-validated parameter knobs, so the default (Sun, Earth-size, Earth-tilt) recovers the present-day
map exactly.

The geography seam — preplanned, carried inert (ADR 0004 #4 / plan §9.3)
-----------------------------------------------------------------------
Every view carries an **elevation** scalar layer tagged ``inert=True`` — carried, displayed, and
round-tripped (via :mod:`planet.planet_spec`), but it **does not change the climate** at v1
(the consuming physics is a rung on the §5 staircase: a lapse-rate diagnostic at the cheap tier,
true 2-D orographic precip at rung 5). Default elevation is flat (zeros); an imported heightmap rides
the same layer. The honesty flag is named, not blurred.

Units (the §7 discipline, carried into the interchange)
-------------------------------------------------------
Temperature in **°C**, precipitation in **cm/yr**, biome as integer :class:`~planet.biomes.Biome`
codes, latitude/longitude in **degrees**, elevation in **m**. Each :class:`Layer` carries its unit
string so the exported state is self-describing.
"""
from __future__ import annotations

import sys
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path

import numpy as np

# --- run-as-script bootstrap: put the repo root on sys.path BEFORE the absolute imports below,
#     so `python -m planet.planetmap` (the HTML-artifact demo) and a notebook `%run`
#     both resolve `planet.*`. A no-op under pytest / `python -m`, where it is already there.
_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from planet import demo_biomes, exoplanet
from planet.albedo import EBMParams
from planet.biomes import Biome, BIOME_COLORS, BIOME_NAMES
from planet.ebm import A_OLR, D_TRANSPORT, S0_EARTH
from planet.exoplanet import T_SUN, exoplanet_params
from planet.obliquity import OBLIQUITY_EARTH, OBLIQUITY_FAITHFUL_MAX, OBLIQUITY_MIN, insolation_s2

N_LON = 73                     # default longitude samples (−180…180 inclusive) — bands, so coarse is fine
LIVE_N_TAU = 0.01              # relaxation step for a live recompute; matches demo_biomes' fine present-day
                               #   step — the steady state carries an O(Δt) operator-splitting bias, so we
                               #   do NOT coarsen it (the biome bands would shift; see ebm.py / planet memory)


# --------------------------------------------------------------------------- #
# 1. The layer-registry primitives (ADR 0004 #1) — NumPy only, no Plotly/ipywidgets.
# --------------------------------------------------------------------------- #
class LayerKind(str, Enum):
    """The three render primitives a view is built from (ADR 0002 §3 / ADR 0004 #1).

    ``SCALAR_FIELD`` — a 2-D lat×lon field painted as the globe surface (temperature, precip, biome,
    elevation). ``ANNOTATION`` — a line/marker overlay given as ``(k, 2)`` ``[lat, lon]`` points (the
    ice line). ``VECTOR_OVERLAY`` — a flow field (the emergent circulation), ``data`` a stacked
    ``(2, n_lat, n_lon)`` velocity ``[u, v]``, painted as flow arrows (:func:`render` /
    :func:`_vector_overlay_trace`) — **Phase 4** built this machinery; it was a declared-but-unpainted
    seam through Phases 1–3. A ``str``-valued enum so ``kind.value`` serializes straight into the
    planet-spec JSON.
    """

    SCALAR_FIELD = "scalar_field"
    VECTOR_OVERLAY = "vector_overlay"
    ANNOTATION = "annotation"


@dataclass(frozen=True, eq=False)
class Grid:
    """The map's grid geometry — full-globe latitudes and longitudes in **degrees**.

    ``lat`` runs pole-to-pole (−90 → +90: the hemisphere EBM grid mirrored about the equator);
    ``lon`` runs −180 → +180. ``eq=False`` because the fields are NumPy arrays (element-wise ``==``
    is ambiguous); array-aware equality lives in :mod:`planet.planet_spec` (the one place a
    *real* equality test — the round-trip invariant — is needed).
    """

    lat: np.ndarray
    lon: np.ndarray
    lat_units: str = "degrees_north"
    lon_units: str = "degrees_east"


@dataclass(frozen=True, eq=False)
class Layer:
    """One registry layer — the self-contained unit a phase registers (ADR 0004 #1).

    ``name`` is the key (unique within a view); ``kind`` a :class:`LayerKind`; ``data`` the array
    (a 2-D lat×lon field for ``SCALAR_FIELD``, ``(k, 2)`` ``[lat, lon]`` points for ``ANNOTATION``);
    ``units`` the self-describing unit string carried into the interchange; ``style`` a **JSON-safe**
    dict of render hints (colorscale, the ``categorical`` flag + biome colour/name maps); ``z_order``
    the paint order; ``inert`` marks a layer that is carried/displayed/round-tripped but **not consumed
    by the climate** (the geography seam, §9.3). ``eq=False`` for the same array reason as :class:`Grid`.
    """

    name: str
    kind: LayerKind
    data: np.ndarray
    units: str
    style: dict = field(default_factory=dict)
    z_order: int = 0
    inert: bool = False


@dataclass(frozen=True, eq=False)
class PlanetView:
    """An ordered stack of :class:`Layer` over a shared :class:`Grid` — the registry, and the export manifest.

    This single structure is what :func:`render` paints and what :mod:`planet.planet_spec`
    serializes (ADR 0004 #3: "the registry *is* the export manifest — one structure, not two"). The
    layers are stored in registration order; :meth:`ordered` returns them by ``z_order`` for painting.
    """

    grid: Grid
    layers: tuple[Layer, ...]

    def layer(self, name: str) -> Layer:
        """The layer registered under ``name`` (raises ``KeyError`` if absent)."""
        for ly in self.layers:
            if ly.name == name:
                return ly
        raise KeyError(f"no layer {name!r}; have {[ly.name for ly in self.layers]}")

    def ordered(self) -> list[Layer]:
        """The layers sorted by ``z_order`` (stable) — the paint order."""
        return sorted(self.layers, key=lambda ly: ly.z_order)

    def scalar_fields(self, include_inert: bool = False) -> list[Layer]:
        """The ``SCALAR_FIELD`` layers — the paintable globe surfaces (inert ones excluded by default)."""
        return [ly for ly in self.ordered()
                if ly.kind is LayerKind.SCALAR_FIELD and (include_inert or not ly.inert)]


def _mirror_profile(profile: np.ndarray) -> np.ndarray:
    """Mirror a hemisphere profile (equator→pole) to the full globe (south pole→north pole).

    The EBM grid is one hemisphere (``x = sin φ`` on [0, 1]); the annual-mean climate is
    hemispherically symmetric, so the southern hemisphere is the northern one reversed
    (``concat([profile[::-1], profile])`` — the same mirroring the static :mod:`~planet.plots`
    biome map uses). Pairs with :func:`_mirror_latitude`.
    """
    profile = np.asarray(profile)
    return np.concatenate([profile[::-1], profile])


def _mirror_latitude(lat_hemi_deg: np.ndarray) -> np.ndarray:
    """Mirror hemisphere latitudes [0…90] to the full globe [−90…90] (``concat([−lat[::-1], lat])``)."""
    lat_hemi_deg = np.asarray(lat_hemi_deg, dtype=float)
    return np.concatenate([-lat_hemi_deg[::-1], lat_hemi_deg])


def _broadcast_field(profile_hemi: np.ndarray, n_lon: int) -> np.ndarray:
    """A hemisphere profile → a full-globe 2-D lat×lon field (mirror, then repeat across longitude)."""
    full = _mirror_profile(profile_hemi)
    return np.repeat(full[:, None], n_lon, axis=1)


def _iceline_annotation(ice_line_lat: float, lon: np.ndarray) -> np.ndarray:
    """The ice-line annotation as ``(k, 2)`` ``[lat, lon]`` points — the two latitude circles at ±φ_ice.

    A NaN ``[nan, nan]`` row separates the northern and southern circles so the rendered line does not
    connect across the globe (Plotly breaks a 3-D line on NaN). A degenerate ice line — ice-free (90°)
    or Snowball (0°) — yields an **empty** ``(0, 2)`` annotation (no cap edge to draw). The NaN
    separator is *why* the round-trip equality uses ``equal_nan`` (ADR 0004 #4 / planet_spec).
    """
    lon = np.asarray(lon, dtype=float)
    if not (0.0 < float(ice_line_lat) < 90.0):
        return np.empty((0, 2), dtype=float)
    north = np.column_stack([np.full(lon.size, float(ice_line_lat)), lon])
    south = np.column_stack([np.full(lon.size, -float(ice_line_lat)), lon])
    gap = np.array([[np.nan, np.nan]])
    return np.vstack([north, gap, south])


def _biome_style() -> dict:
    """JSON-safe render hints for the categorical biome layer — the discrete colour/name maps.

    Keys are ``str(code)`` (JSON object keys must be strings); covers all nine biomes so the legend is
    stable as bands migrate. Pairs with the continuous-field ``{"colorscale": …}`` style.
    """
    return {
        "categorical": True,
        "colors": {str(int(b)): BIOME_COLORS[b] for b in Biome},
        "names": {str(int(b)): BIOME_NAMES[b] for b in Biome},
    }


def circulation_layer(jet, lat_full: np.ndarray, n_lon: int = N_LON) -> Layer:
    """A ``VECTOR_OVERLAY`` circulation layer from a Phase-4 coupled jet (the §9 Phase-4 registration).

    Maps the coupler's zonal-wind profile ``u(φ)`` (a midlatitude channel band, NH) onto the full-globe
    latitude grid — **mirrored to both hemispheres** (mid-latitude westerlies in each) and zero outside
    the channel band — and broadcasts it across longitude. The layer ``data`` is the stacked
    ``(2, n_lat, n_lon)`` velocity ``[u, v]`` (m/s; ``v ≡ 0`` for this zonally-symmetric jet, but the
    two-component shape is the general vector-field seam). The renderer paints it as flow arrows
    (:func:`render`). This is the layer Phase 4 *registers* — the deferred ``VECTOR_OVERLAY`` kind, now
    painted (build the seam, *then* the machinery).
    """
    abs_lat = np.abs(np.asarray(lat_full, dtype=float))
    u1d = np.interp(abs_lat, jet.phi, jet.u_profile, left=0.0, right=0.0)   # 0 outside the channel band
    u_field = np.repeat(u1d[:, None], n_lon, axis=1)
    data = np.stack([u_field, np.zeros_like(u_field)])                      # (2, n_lat, n_lon) = [u, v]
    return Layer("circulation", LayerKind.VECTOR_OVERLAY, data, "m/s",
                 style={"colorscale": "RdBu_r", "arrow_color": "#1a1a1a",
                        "label": f"zonal wind — jet {jet.jet_speed:.0f} m/s @ {jet.jet_lat:.0f}°"},
                 z_order=2)


def build_view(result, n_lon: int = N_LON, elevation: np.ndarray | None = None, jet=None) -> PlanetView:
    """Turn a Phase-2 climate result into the v1 :class:`PlanetView` — the biome-map layer stack.

    Consumes a :class:`~planet.demo_biomes.BiomeResult` (the validated climate → precip →
    biome composition) and registers the v1 layers: **temperature** & **precipitation** & **biome**
    (``SCALAR_FIELD`` — each the hemisphere profile mirrored + broadcast across longitude), the **ice
    line** (``ANNOTATION``), and an inert **elevation** layer (the geography seam, §9.3 — flat zeros by
    default, or an imported full-globe heightmap). With a Phase-4 ``jet`` it *also* registers the
    **circulation** (``VECTOR_OVERLAY``) — the seam the renderer paints with no restructuring (ADR 0004
    #1). No physics here — only re-shaping validated arrays into the renderer-input seam.

    Parameters
    ----------
    result : BiomeResult
        The climate state + precip field + biome codes (:func:`planet.demo_biomes.compute`).
    n_lon : int
        Longitude samples for the broadcast (bands, so coarse is fine).
    elevation : ndarray | None
        Optional full-globe ``(n_lat, n_lon)`` elevation field (m) for the inert geography layer;
        defaults to flat (zeros). Carried and round-tripped, never consumed by the climate at v1.
    jet : CoupledJet | None
        Optional Phase-4 emergent jet (:func:`planet.coupler.couple_jet`); if given, a
        ``circulation`` ``VECTOR_OVERLAY`` layer is registered. (Not part of the live-slider loop — the
        jet is a separate integration, the first compute too heavy for the rung-0 instant remap, §9.2.)
    """
    state = result.state
    lat_full = _mirror_latitude(state.latitude_deg())
    lon = np.linspace(-180.0, 180.0, n_lon)
    grid = Grid(lat=lat_full, lon=lon)

    n_lat = lat_full.size
    if elevation is None:
        elevation_field = np.zeros((n_lat, n_lon), dtype=float)
    else:
        elevation_field = np.asarray(elevation, dtype=float)
        if elevation_field.shape != (n_lat, n_lon):
            raise ValueError(f"elevation must be {(n_lat, n_lon)} (full-globe lat×lon), got {elevation_field.shape}")

    layers = [
        Layer("temperature", LayerKind.SCALAR_FIELD, _broadcast_field(state.T, n_lon), "°C",
              style={"colorscale": "RdBu_r"}, z_order=0),
        Layer("precipitation", LayerKind.SCALAR_FIELD, _broadcast_field(result.precip_cm, n_lon), "cm/yr",
              style={"colorscale": "GnBu"}, z_order=0),
        Layer("biome", LayerKind.SCALAR_FIELD, _broadcast_field(result.codes, n_lon).astype(int),
              "Whittaker biome code", style=_biome_style(), z_order=1),
        Layer("ice_line", LayerKind.ANNOTATION, _iceline_annotation(state.ice_line_lat, lon), "degrees",
              style={"color": "#c0392b", "label": "ice line (T = T_freeze)"}, z_order=3),
        Layer("elevation", LayerKind.SCALAR_FIELD, elevation_field, "m",
              style={"colorscale": "Earth"}, z_order=-1, inert=True),
    ]
    if jet is not None:
        layers.append(circulation_layer(jet, lat_full, n_lon))
    return PlanetView(grid=grid, layers=tuple(layers))


def climate_view(S0: float = S0_EARTH, A: float = A_OLR, D: float = D_TRANSPORT, *,
                 T_star: float = T_SUN, size: float = 1.0, obliquity_deg: float = OBLIQUITY_EARTH,
                 n_cells: int = 180, n_tau: float = LIVE_N_TAU, n_lon: int = N_LON,
                 elevation: np.ndarray | None = None) -> PlanetView:
    """The live recompute: knob values → a fresh equilibrium climate → the biome-map :class:`PlanetView`.

    Pure re-composition of :func:`planet.demo_biomes.compute` (the tested chain) at a knob
    setting — the rung-0 instant-remap trigger (ADR 0004 #2). The wired knobs are all **already
    validated**: the climate levers ``S₀`` (the Snowball lever), the OLR offset ``A`` (a drop ≈ a CO₂
    increase — :mod:`~planet.demo_biomes`'s warming knob), and the meridional transport ``D``;
    the **obliquity** ``obliquity_deg`` (axial tilt → the insolation P₂ coefficient ``s₂`` —
    :mod:`planet.obliquity`, now that its source is pinned; more tilt spreads sunlight poleward
    → a flatter planet, the ice cap retreats); plus the two **exoplanet** knobs
    (:mod:`planet.exoplanet`, §9.1) — the host **star's effective temperature ``T_star``** (a
    redder star lowers the ice albedo → harder to snowball) and the **planet ``size``** in Earth radii
    (bigger → weaker transport → a sharper gradient). Defaults (``T_star=T_SUN``, ``size=1``,
    ``obliquity_deg=OBLIQUITY_EARTH``) recover the Earth model exactly, so the present-day map is
    unchanged. Each call relaxes to equilibrium from the Earth-like initial condition
    (:func:`~planet.albedo.present_day_climate`'s finite-cap start), so it shows the *climate
    at these knobs* — **not** the hysteresis branch-tracking, which is the Snowball demo's continuation
    sweep (a different, path-dependent question).
    """
    base = EBMParams(S0=S0, A=A, D=D, s2=insolation_s2(obliquity_deg), n_cells=n_cells)
    params = exoplanet_params(T_star=T_star, size=size, base=base)
    result = demo_biomes.compute(params, n_tau=n_tau)
    return build_view(result, n_lon=n_lon, elevation=elevation)


# --------------------------------------------------------------------------- #
# 2. The Plotly renderer + the HTML artifact (Plotly imported lazily; importorskip-gated).
# --------------------------------------------------------------------------- #
def _sphere_xyz(lat_deg, lon_deg):
    """Unit-sphere Cartesian coordinates from latitude/longitude in degrees (the globe parametrization)."""
    lat = np.radians(np.asarray(lat_deg, dtype=float))
    lon = np.radians(np.asarray(lon_deg, dtype=float))
    return np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)


def _discrete_colorscale(layer: Layer):
    """A Plotly discrete colorscale + (cmin, cmax) + tick (value, name) pairs for a categorical layer."""
    colors = layer.style["colors"]
    codes = sorted(int(c) for c in colors)
    cmin, cmax = codes[0] - 0.5, codes[-1] + 0.5
    span = cmax - cmin
    scale = []
    for i, c in enumerate(codes):
        lo = (c - 0.5 - cmin) / span
        hi = (c + 0.5 - cmin) / span
        scale.append([lo, colors[str(c)]])
        scale.append([hi, colors[str(c)]])
    names = layer.style.get("names", {})
    ticks = [(float(c), names.get(str(c), str(c))) for c in codes]
    return scale, cmin, cmax, ticks


def _hovertext(grid: Grid, layer: Layer) -> np.ndarray:
    """Per-cell hover strings (lat°, lon°, value) — the default surface hover shows meaningless x/y/z."""
    LON, LAT = np.meshgrid(grid.lon, grid.lat)
    if layer.style.get("categorical"):
        names = layer.style.get("names", {})
        label = np.vectorize(lambda v: names.get(str(int(v)), str(int(v))))(layer.data)
        val = label
        fmt = lambda la, lo, v: f"lat {la:.0f}°, lon {lo:.0f}°<br>{layer.name}: {v}"
    else:
        val = layer.data
        fmt = lambda la, lo, v: f"lat {la:.0f}°, lon {lo:.0f}°<br>{layer.name}: {v:.1f} {layer.units}"
    return np.vectorize(fmt)(LAT, LON, val)


def _vector_overlay_trace(go, grid: Grid, layer: Layer):
    """A Plotly ``Cone`` trace painting a ``VECTOR_OVERLAY`` velocity field as flow arrows on the globe.

    ``layer.data`` is the stacked ``(2, n_lat, n_lon)`` horizontal velocity ``[u, v]`` (m/s). At each
    (sub-sampled) cell the eastward/northward components are rotated onto the **sphere-tangent** 3-D
    directions (east ``ê_λ``, north ``ê_φ``) and drawn as a cone (arrowhead) just above the surface, so
    a band of eastward cones at mid-latitudes *is* the emergent westerly jet. Cones are coloured by
    wind speed. This is the machinery the ``VECTOR_OVERLAY`` seam deferred to Phase 4 (ADR 0004 #1)."""
    u_field, v_field = np.asarray(layer.data[0], dtype=float), np.asarray(layer.data[1], dtype=float)
    lat = np.radians(grid.lat)
    lon = np.radians(grid.lon)
    # sub-sample so the globe is not crowded (≈16 lat × 12 lon arrows); skip negligible wind.
    sj = max(1, lat.size // 16)
    si = max(1, lon.size // 12)
    xs, ys, zs, us, vs, ws = [], [], [], [], [], []
    umax = float(np.max(np.abs(u_field))) or 1.0
    for j in range(0, lat.size, sj):
        for i in range(0, lon.size, si):
            u, v = u_field[j, i], v_field[j, i]
            if abs(u) + abs(v) < 0.05 * umax:
                continue
            la, lo = lat[j], lon[i]
            east = np.array([-np.sin(lo), np.cos(lo), 0.0])
            north = np.array([-np.sin(la) * np.cos(lo), -np.sin(la) * np.sin(lo), np.cos(la)])
            pos = 1.02 * np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])
            vec = u * east + v * north
            xs.append(pos[0]); ys.append(pos[1]); zs.append(pos[2])
            us.append(vec[0]); vs.append(vec[1]); ws.append(vec[2])
    return go.Cone(
        x=xs, y=ys, z=zs, u=us, v=vs, w=ws,
        colorscale=layer.style.get("colorscale", "RdBu_r"), showscale=False,
        sizemode="scaled", sizeref=0.8, anchor="tail",
        name=layer.style.get("label", layer.name), hoverinfo="name",
    )


def render(view: PlanetView, active: str = "biome"):
    """Paint a :class:`PlanetView` as an interactive Plotly globe — the generic, kind-dispatching renderer.

    The ``active`` ``SCALAR_FIELD`` layer becomes the globe surface (a ``go.Surface`` on the unit
    sphere, with a discrete colorscale for the categorical biome layer or a continuous one otherwise,
    plus custom lat/lon/value hover); every ``ANNOTATION`` layer is overlaid as a 3-D line. The
    renderer is **generic over** :class:`LayerKind` — it dispatches, it does not special-case a phase —
    so registering a new scalar or annotation layer needs no edit here (ADR 0004 #1).

    ``VECTOR_OVERLAY`` layers are painted as flow arrows (:func:`_vector_overlay_trace` — Plotly cones on
    the sphere-tangent plane) — the Phase-4 machinery the seam deferred through Phases 1–3. Requires the
    optional ``[webviz]`` extra (Plotly); raises ``ImportError`` without it (caught by callers).
    """
    import plotly.graph_objects as go

    grid = view.grid
    scalars = view.scalar_fields(include_inert=True)
    names = [ly.name for ly in scalars]
    if active not in names:
        raise KeyError(f"no scalar layer {active!r} to paint; have {names}")
    surface_layer = view.layer(active)

    LON, LAT = np.meshgrid(grid.lon, grid.lat)       # both (n_lat, n_lon); rows = lat, cols = lon
    X, Y, Z = _sphere_xyz(LAT, LON)
    hover = _hovertext(grid, surface_layer)
    if surface_layer.style.get("categorical"):
        scale, cmin, cmax, ticks = _discrete_colorscale(surface_layer)
        colorbar = dict(tickvals=[t[0] for t in ticks], ticktext=[t[1] for t in ticks],
                        title="biome", len=0.8)
        surf_kw = dict(colorscale=scale, cmin=cmin, cmax=cmax, colorbar=colorbar)
    else:
        surf_kw = dict(colorscale=surface_layer.style.get("colorscale", "Viridis"),
                       colorbar=dict(title=f"{surface_layer.name} ({surface_layer.units})", len=0.8))

    traces = [go.Surface(x=X, y=Y, z=Z, surfacecolor=surface_layer.data,
                         text=hover, hoverinfo="text", showscale=True, **surf_kw)]

    for ly in view.ordered():
        if ly.kind is LayerKind.SCALAR_FIELD:
            continue                                          # one scalar is the surface; others switch via `active`
        if ly.kind is LayerKind.ANNOTATION:
            if ly.data.size == 0:
                continue                                      # degenerate ice line (ice-free / Snowball): nothing to draw
            ax, ay, az = _sphere_xyz(ly.data[:, 0], ly.data[:, 1])
            traces.append(go.Scatter3d(
                x=ax, y=ay, z=az, mode="lines", name=ly.style.get("label", ly.name),
                line=dict(color=ly.style.get("color", "#000000"), width=5), hoverinfo="name"))
        elif ly.kind is LayerKind.VECTOR_OVERLAY:
            traces.append(_vector_overlay_trace(go, grid, ly))

    title = (f"Planet — {surface_layer.name} ({surface_layer.units})"
             if not surface_layer.style.get("categorical") else "Planet — biome map")
    no_axis = dict(showbackground=False, showticklabels=False, showgrid=False, zeroline=False, visible=False)
    fig = go.Figure(data=traces)
    fig.update_layout(
        title=title, width=820, height=720, margin=dict(l=0, r=0, t=40, b=0),
        scene=dict(xaxis=no_axis, yaxis=no_axis, zaxis=no_axis,
                   aspectmode="data", camera=dict(eye=dict(x=1.4, y=1.4, z=0.9))),
    )
    return fig


def save_html(view: PlanetView, path, active: str = "biome") -> Path:
    """Render ``view`` and write a standalone, viewable HTML globe — the banked artifact (Plotly only).

    The interactive analogue of the static ``demo_*.py`` PNGs: a self-contained HTML globe needing no
    server and no extra dependency beyond Plotly (no kaleido/PNG export). Returns the written path.
    """
    fig = render(view, active=active)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(path), include_plotlyjs="cdn")
    return path


# --------------------------------------------------------------------------- #
# 3. interactive_map() — the live-slider loop (the ONLY place ipywidgets is imported; not unit-tested).
# --------------------------------------------------------------------------- #
def interactive_map(n_tau: float = LIVE_N_TAU):
    """The in-notebook live map: knob-sliders → instant recompute-and-remap (ADR 0004 #2, rung 0).

    The ``main()`` analogue (Steel ``app.py`` discipline): the *only* place ``import ipywidgets`` lives,
    kept paper-thin — every value comes from the tested :func:`climate_view` / :func:`render` above, so
    the only statements that can raise are widget/Plotly calls. Not unit-tested (the live UI is reach,
    ADR 0002). Sliders: the climate levers ``S₀`` / CO₂ (→ ``A``) / transport ``D``, the two **exoplanet**
    knobs **star ``T_star``** & **planet ``size``** (§9.1 — now pinned and wired), the **obliquity**
    (axial tilt → ``s₂`` — :mod:`planet.obliquity`, now pinned and wired too), and a layer
    selector. Requires the ``[webviz]`` extra (Plotly + ipywidgets).
    """
    import ipywidgets as widgets
    import plotly.graph_objects as go
    from IPython.display import display

    s0 = widgets.FloatSlider(value=S0_EARTH, min=1000.0, max=1900.0, step=10.0,
                             description="S₀ (W/m²)", continuous_update=False, readout_format=".0f")
    a = widgets.FloatSlider(value=A_OLR, min=180.0, max=230.0, step=1.0,
                            description="A — CO₂↓", continuous_update=False, readout_format=".0f")
    d = widgets.FloatSlider(value=D_TRANSPORT, min=0.1, max=1.2, step=0.01,
                            description="transport D", continuous_update=False, readout_format=".2f")
    t_star = widgets.FloatSlider(value=T_SUN, min=exoplanet.STAR_TEFF_MIN, max=8000.0, step=50.0,
                                 description="star Teff (K)", continuous_update=False, readout_format=".0f",
                                 tooltip="host-star spectrum → ice albedo (redder = harder to snowball; §9.1)")
    size = widgets.FloatSlider(value=1.0, min=exoplanet.SIZE_MIN, max=2.0, step=0.05,
                               description="size (R⊕)", continuous_update=False, readout_format=".2f",
                               tooltip="planet size → transport only (bigger = sharper gradient; §9.1)")
    obliquity = widgets.FloatSlider(value=OBLIQUITY_EARTH, min=OBLIQUITY_MIN, max=OBLIQUITY_FAITHFUL_MAX, step=0.5,
                                    description="obliquity (°)", continuous_update=False, readout_format=".1f",
                                    tooltip="axial tilt → insolation gradient (more tilt = flatter, ice retreats; §9.1)")
    layer = widgets.Dropdown(options=["biome", "temperature", "precipitation", "elevation"],
                             value="biome", description="layer")

    fig = go.FigureWidget(render(climate_view(n_tau=n_tau)))

    def update(_=None):
        new = render(climate_view(S0=s0.value, A=a.value, D=d.value, T_star=t_star.value,
                                  size=size.value, obliquity_deg=obliquity.value, n_tau=n_tau),
                     active=layer.value)
        with fig.batch_update():
            fig.data = ()                                     # drop the old traces
            for tr in new.data:
                fig.add_trace(tr)
            fig.layout.title = new.layout.title

    for w in (s0, a, d, t_star, size, obliquity, layer):
        w.observe(update, names="value")

    controls = widgets.VBox([s0, a, d, t_star, size, obliquity, layer])
    display(widgets.HBox([controls, fig]))
    return fig


def main() -> None:
    """Headless demo: build the present-day biome-map globe and write the banked HTML artifact."""
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")              # °C, ₂, → on legacy codepages

    docs_html = _REPO_ROOT / "docs" / "figures" / "planet-map.html"
    output_html = _REPO_ROOT / "outputs" / "planet-map.html"
    view = climate_view()
    print("Planet interactive map (v1 — the biome map): present-day globe")
    print(f"  layers: {[ly.name for ly in view.layers]}")
    try:
        for target in (docs_html, output_html):
            save_html(view, target)
        print(f"  HTML globe saved → {docs_html.relative_to(_REPO_ROOT)}")
        print("  open it in a browser (rotate / zoom / hover); for the live sliders run "
              "interactive_map() in a notebook (pip install -e .[webviz]).")
    except ImportError:
        print("(plotly not installed — install the webviz extra to render the globe: "
              "pip install -e .[webviz])")


if __name__ == "__main__":
    main()
