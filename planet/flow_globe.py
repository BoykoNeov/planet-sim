"""Rung C — the **showcase**: a flow-on-a-globe **particle-streaming** renderer (§9.5).

Rung A animated the eddy life cycle as a flat two-panel movie; Rung B painted the same banked frames
as a *scalar field* on the Plotly globe — both **honest-by-construction** (the geometry cannot lie).
Rung C is the **immersive** view: particles streaming along the emergent flow on a *real*, rotatable
3-D planet — the NASA *Perpetual-Ocean* / Ventusky look the user is after. It is **reach / delivery,
not new teaching**, and it is the one renderer governed by the **honest-by-disclosure** carve-out
(ADR 0002 status note, 2026-06-12): it *may* render a view the model does not literally compute (a
continuously-streaming field that *implies* persistent currents, though the instantaneous eddy flux is
~90 % reversible) **provided a visible on-screen disclaimer documents the departure.**

A general-purpose renderer, not an eddy-specific one
----------------------------------------------------
The module is architected around a **renderer-agnostic data contract** (:class:`FlowField`), *not* the
eddy band: a lat×lon grid, per-cell ``(u, v)``, an optional ``scalar`` for colour, a **coverage** extent,
and an **honesty** disclaimer string. The eddy band is merely its *first* consumer
(:func:`flow_field_from_eddy`); the same contract would one day carry a full GCM/ESM wind-or-current
field. The contract commits to **nothing** about projection or particles — those are renderer-side —
which is why the renderer can be swapped (the originally-planned Canvas2D globe → this three.js sphere)
without touching the contract, §9.3, or the carve-out.

The two honesty edges, and how each is carried
-----------------------------------------------
(1) **One band, not a globe.** The flow is a doubly-periodic midlatitude β-plane patch (the same edge as
:func:`planet.eddy_globe`); its honest extent is carried *in the data* as :class:`Coverage`, and the
particles are seeded **only within it** — the rest of the planet is left bare. We do **not** fabricate a
global ``(u, v)`` from a 55° patch (that would be inventing data, not illustrating a richer model).
(2) **~90 % reversible.** Band-confinement keeps the coverage truthful but does **not** make the showcase
honest-by-construction: streaming particles still *imply* persistent currents the reversible flux does
not produce. So the disclaimer stays mandatory and carries the *"mostly sloshes / net is only the small κ
residual"* clause — it **is the entire license**, and is the one thing machine-checked
(``test_flow_globe.py``), on-screen and legible, even though the physics-fidelity of the render is not.

Vendoring (the §6 deliverable this renderer *owes*)
---------------------------------------------------
The scene/camera/WebGL framework is **three.js** (r137 UMD), **vendored inline** at
``planet/vendor/three.min.js`` and emitted verbatim into the HTML — so the artifact opens straight off
``file://`` with no network (the property a CDN ``<script src>`` would break). three.js' MIT licence is
attributed in the repo ``NOTICE`` and in the inlined ``@license`` banner that travels with the artifact.
Only three.js is vendored; the particle advection and the orbit camera here are **original**.

This module is **NumPy-only at import** (no plotly / matplotlib): it builds an HTML *string*. Run it via
the demo::

    python -m planet.demo_eddy_particles
"""
from __future__ import annotations

import html
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np

_VENDOR_DIR = Path(__file__).resolve().parent / "vendor"
_THREE_JS_PATH = _VENDOR_DIR / "three.min.js"

DEFAULT_N_PARTICLES = 6500       # streaming points — reads as flow without taxing the per-frame JS loop
_BAND_CROSSING_SECONDS = 6.0     # wall-clock seconds for the fastest particle to cross the band's width
DEFAULT_PARTICLE_SIZE = 0.035      # point size in sphere-radius units (the size slider's initial value)
DEFAULT_PARTICLE_OPACITY = 0.95    # overall particle opacity ceiling (the opacity slider's initial value)
DEFAULT_PARTICLE_SHARPNESS = 0.55  # sprite edge: 0 = soft bloom, 1 = hard disc (the sharpness slider's start)


# --------------------------------------------------------------------------------------------------- #
# The renderer-agnostic data contract — a generic vector-field-on-a-globe layer (the "one day a GCM"
# hook). It carries only what *any* such renderer needs; nothing about projection or particles.
# --------------------------------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Coverage:
    """The honest spatial extent a :class:`FlowField` actually covers — the band-vs-globe truth, in the data.

    Today the eddy field sets a bounded midlatitude sector (``is_global=False``); a future global model
    would set the full sphere (``is_global=True``). Particles are seeded only within this box, so the
    renderer never paints flow where the field has none.
    """

    lat_min: float
    lat_max: float
    lon_min: float
    lon_max: float
    is_global: bool = False


@dataclass(frozen=True)
class FlowField:
    """A vector field on a lat×lon globe patch — the Rung-C renderer's only input (§9.5).

    ``lat`` (ny,) / ``lon`` (nx,) the cell-centre coordinates (deg, monotone increasing); ``u``/``v``
    (ny, nx) the eastward/northward velocity at those centres (m/s); ``coverage`` the honest extent
    (above); ``honesty`` the on-screen disclaimer text — **the entire license** under the
    honest-by-disclosure carve-out. ``scalar`` (ny, nx) an optional field that colours the particles,
    ``scalar_label`` its name; ``radius_m`` the planet radius used by the lat/lon advection metric
    (``dλ/dt = u/(a cosφ)``, ``dφ/dt = v/a``).

    ``mask`` (ny, nx, bool) is the **per-cell validity mask** (§9.6 O1 — the first contract growth past
    R1): ``True`` = a valid data cell, ``False`` = no data there (land, or an unobserved cell of a real
    ocean product). :class:`Coverage` stays the *bounding box*; the mask is validity *inside* it — an
    ocean field's land is inside its box, so a box alone cannot carry it. Particles are seeded and kept
    **only in valid cells** (the R1 band-zeros honesty style: where-the-data-is is carried in the data,
    not a caption). ``None`` = all cells valid = the exact pre-mask behaviour (default-off discipline).
    """

    lat: np.ndarray
    lon: np.ndarray
    u: np.ndarray
    v: np.ndarray
    coverage: Coverage
    honesty: str
    scalar: Optional[np.ndarray] = None
    scalar_label: str = ""
    radius_m: float = 6.371e6
    mask: Optional[np.ndarray] = None


def flow_field_from_eddy(eddy) -> FlowField:
    """Map a released eddy life cycle (:class:`planet.eddy_flux.EddyFlux`) onto the generic contract.

    Uses the **saturated-flow** snapshot (the banked frame nearest ``eddy.saturation_period``) for a
    steady streaming look, and **reuses Rung B's band geometry** (:func:`planet.eddy_globe._band_geometry`
    / :func:`planet.eddy_globe._earth_radius`) so the two globes can never drift on *where* the band sits.
    The C-grid face velocities are collocated to cell centres exactly as
    :func:`planet.eddy_flux.eddy_heat_flux` does for ``v``. Raises :class:`ValueError` if the frames
    side-channel is absent (recompute with ``eddy_life_cycle(..., n_frames=N)``).
    """
    fr = eddy.frames
    if fr is None:
        raise ValueError("eddy.frames is None — recompute with eddy_life_cycle(..., n_frames=N)")
    from .eddy_globe import _band_geometry, _earth_radius

    (_xyz, _lat2d, _lon2d, lon1d) = _band_geometry(fr)
    a = float(_earth_radius(fr))
    k = int(np.argmin(np.abs(fr.times - eddy.saturation_period)))      # the saturated-flow frame
    # collocate the C-grid faces to cell centres (u on E–W faces → axis 1; v on N–S faces → axis 0).
    u_c = 0.5 * (fr.u[k] + np.roll(fr.u[k], -1, axis=1))
    v_c = 0.5 * (fr.v[k] + np.roll(fr.v[k], -1, axis=0))
    lat = np.asarray(fr.phi, dtype=float)
    lon = np.asarray(lon1d, dtype=float)
    theta = np.asarray(fr.theta[k], dtype=float)

    rev_pct = round(100.0 * (1.0 - eddy.irreversible_fraction))        # the share that just sloshes (~92%)
    res_pct = round(100.0 * eddy.irreversible_fraction)                # the surviving net κ transport (~8%)
    lon_span = float(lon.max() - lon.min())
    honesty = (
        f"This is ONE midlatitude band — a single β-plane latitude zone about {lon_span:.0f}° of "
        f"longitude wide, in one hemisphere — NOT a planet-wide flow. The rest of the globe is left "
        f"bare because the model resolves only this band. The streaming particles are illustrative: they "
        f"imply persistent currents, but the instantaneous eddy heat flux is about {rev_pct:.0f}% "
        f"reversible — it mostly sloshes back and forth, and the genuine net transport is only the "
        f"small remaining ~{res_pct:.0f}% (the emergent eddy diffusivity κ). “Currents carrying "
        f"heat” is an artistic reading the numbers do not validate."
    )
    coverage = Coverage(lat_min=float(lat.min()), lat_max=float(lat.max()),
                        lon_min=float(lon.min()), lon_max=float(lon.max()), is_global=False)
    return FlowField(lat=lat, lon=lon, u=u_c, v=v_c, coverage=coverage, honesty=honesty,
                     scalar=theta, scalar_label="θ (°C)", radius_m=a)


# --------------------------------------------------------------------------------------------------- #
# The renderer — emit a self-contained three.js HTML scene (data + three.js + app, all inlined).
# --------------------------------------------------------------------------------------------------- #
def _three_js() -> str:
    """The vendored three.js (r137 UMD) source, inlined verbatim — its ``@license`` banner intact."""
    return _THREE_JS_PATH.read_text(encoding="utf-8")


def _build_data(field: FlowField, n_particles: int, particle_size: float,
                particle_opacity: float, particle_sharpness: float,
                crossing_seconds: float = _BAND_CROSSING_SECONDS,
                sequential: bool = False) -> dict:
    """Pack a :class:`FlowField` into the JSON the in-browser renderer consumes (compact, rounded)."""
    lat = np.asarray(field.lat, dtype=float)
    lon = np.asarray(field.lon, dtype=float)
    u = np.asarray(field.u, dtype=float)
    v = np.asarray(field.v, dtype=float)
    scalar = None if field.scalar is None else np.asarray(field.scalar, dtype=float)
    mask = None if field.mask is None else np.asarray(field.mask, dtype=bool)

    center_lat = float(lat.mean())
    center_lon = float(lon.mean())
    span_lon = float(lon.max() - lon.min()) or 1.0
    umax = float(np.max(np.abs(u))) or 1.0
    # the field's peak *speed* — the normaliser for speed-weighted seeding (§9.6 O3c): particles respawn
    # with acceptance ∝ |u,v|/speed_max so the fast western-boundary currents visually dominate the way
    # they physically do. Independent of the `scalar` (which may be θ, not speed), so it is always honest.
    speed_max = float(np.max(np.hypot(u, v))) or 1.0
    # `accel` (model-seconds per wall-clock second) is auto-scaled so the *fastest* particle crosses the
    # field's longitude width in ~`crossing_seconds` — readable streaming regardless of the field's
    # actual speeds. In-browser step: Δλ_deg = deg(u/(a·cosφ)) · accel · dt_real (the honest metric × a
    # purely-visual time-acceleration; the showcase's physics-fidelity is relaxed, ADR 0002 note).
    # The default (`_BAND_CROSSING_SECONDS`) was tuned on the ~55° eddy band; a 360° global field wants a
    # proportionally longer crossing (else the whole ocean sprints), so the pace is a caller knob (§9.6 O2).
    dlamdt_max = math.degrees(umax / (field.radius_m * math.cos(math.radians(center_lat))))
    accel = (span_lon / crossing_seconds) / dlamdt_max if dlamdt_max > 0 else 1.0

    def flat(arr, nd):
        return [round(float(x), nd) for x in np.asarray(arr).ravel(order="C")]

    return {
        "lat": flat(lat, 4), "lon": flat(lon, 4),
        "u": flat(u, 3), "v": flat(v, 3),
        "mask": ([int(x) for x in mask.ravel(order="C")] if mask is not None else None),
        "scalar": (flat(scalar, 3) if scalar is not None else None),
        "scalar_min": round(float(scalar.min()), 3) if scalar is not None else 0.0,
        "scalar_max": round(float(scalar.max()), 3) if scalar is not None else 1.0,
        "scalar_label": field.scalar_label,
        "radius_m": field.radius_m,
        "accel": accel,
        "speed_max": round(speed_max, 4),           # normaliser for speed-weighted seeding (§9.6 O3c)
        "sequential": bool(sequential),             # True → the sequential speed ramp (else diverging RdBu_r)
        "coverage": {
            "lat_min": round(float(field.coverage.lat_min), 4),
            "lat_max": round(float(field.coverage.lat_max), 4),
            "lon_min": round(float(field.coverage.lon_min), 4),
            "lon_max": round(float(field.coverage.lon_max), 4),
            "is_global": bool(field.coverage.is_global),
        },
        "center_lat": round(center_lat, 4), "center_lon": round(center_lon, 4),
        "n_particles": int(n_particles),
        "particle_size": float(particle_size),
        "particle_opacity": float(particle_opacity),
        "particle_sharpness": float(particle_sharpness),
    }


_CSS = """\
:root { color-scheme: dark; }
* { box-sizing: border-box; }
html, body { margin: 0; height: 100%; }
body { font: 15px/1.5 system-ui, -apple-system, Segoe UI, Roboto, sans-serif;
       color: #e8ecf4; background: #0b1020; overflow: hidden; }
#stage { position: relative; width: 100vw; height: 100vh; }
#globe { display: block; width: 100%; height: 100%; cursor: grab; }
#globe:active { cursor: grabbing; }
.title { position: absolute; top: 14px; left: 0; right: 0; text-align: center;
         pointer-events: none; text-shadow: 0 1px 8px #000a; }
.title h1 { margin: 0; font-size: 1.2rem; font-weight: 650; letter-spacing: -.01em; }
.title p { margin: .25rem auto 0; max-width: 42rem; color: #c4cde2; font-size: .86rem; }
.hint { position: absolute; top: 16px; right: 16px; color: #8b95ad; font-size: .78rem;
        pointer-events: none; text-shadow: 0 1px 6px #000a; }
.disclaimer { position: absolute; left: 50%; bottom: 16px; transform: translateX(-50%);
              max-width: 58rem; width: calc(100% - 2rem); background: rgba(12, 18, 38, .9);
              border: 1px solid #36406a; border-radius: 10px; padding: .7rem 1rem;
              color: #d8e0f2; font-size: .84rem; line-height: 1.45; }
.disclaimer strong { color: #ffd166; }
.controls { position: absolute; top: 16px; left: 16px; display: flex; flex-direction: column;
            gap: .45rem; background: rgba(12, 18, 38, .82); border: 1px solid #2a3358;
            border-radius: 10px; padding: .6rem .75rem; font-size: .76rem; color: #c4cde2;
            text-shadow: 0 1px 6px #000a; }
.controls label { display: flex; align-items: center; gap: .6rem; justify-content: space-between; }
.controls input[type="range"] { width: 9rem; accent-color: #ffd166; }
"""

# The in-browser app. Plain JS (no f-string — the braces are JS). DATA + three.js are inlined ahead of
# this block. The particle advection — both the GPU ping-pong pass (the default) and the CPU step()
# fallback — and the spherical orbit camera are ORIGINAL; three.js supplies only the scene / camera /
# WebGL renderer / render targets / shader-material plumbing (vendored inline).
_APP_JS = r"""
const D = window.FLOW_DATA, T = window.THREE;
const lat = D.lat, lon = D.lon, ny = lat.length, nx = lon.length;
const U = D.u, V = D.v, S = D.scalar, smin = D.scalar_min, smax = D.scalar_max;
const MASK = D.mask;   // per-cell validity, 1 = data / 0 = land or unobserved (null = every cell valid)
const A = D.radius_m, ACCEL = D.accel, cov = D.coverage;
const DEG = 180 / Math.PI, RAD = Math.PI / 180;
const SEQ = !!D.sequential;            // §9.6 O3c: sequential speed ramp (ocean) vs diverging RdBu_r (θ)
const SPEEDMAX = D.speed_max || 1.0;   // peak |u,v| — the speed-weighted-seeding normaliser
const SEED_FLOOR = 0.08;               // min respawn acceptance so calm regions keep a light ambient fill
let density = 1.0;                      // the §9.5 density knob (fraction of particles drawn; 1 = all)

const stage = document.getElementById("stage");
const canvas = document.getElementById("globe");
const renderer = new T.WebGLRenderer({ canvas: canvas, antialias: true });
renderer.setPixelRatio(Math.min(2, window.devicePixelRatio || 1));
function size() { return [stage.clientWidth || 800, stage.clientHeight || 600]; }
let [W, H] = size();
renderer.setSize(W, H, false);

const scene = new T.Scene();
scene.background = new T.Color(0x0b1020);
const camera = new T.PerspectiveCamera(42, W / H, 0.01, 100);

// lat/lon (deg) -> unit-sphere xyz (y = north pole up). Internally consistent with the camera below.
function sph(latd, lond, r) {
  const p = latd * RAD, l = lond * RAD;
  return [r * Math.cos(p) * Math.cos(l), r * Math.sin(p), r * Math.cos(p) * Math.sin(l)];
}

// the bare planet: an OPAQUE sphere, so particles on the far side are correctly occluded (back-face).
// ===== (§9.6 O3a) the land/ocean base layer ======================================================== #
// When the field carries an O1 validity mask, drape the sphere with a two-tone land/ocean skin derived
// from that mask — coastlines are what make a stream read as *the Gulf Stream*. Alignment is
// honest-by-CONSTRUCTION: the base fragment shader inverts every surface point to (lat, lon) with the
// SAME mapping the particles use (`sph()`), then samples the SAME mask on the SAME bilinear + 0.5-threshold
// coastline rule as `validAt()` — so the coast under the particles can never drift from the coast under the
// base. No mask (the eddy band) → the plain solid sphere, exactly as before; a base-shader compile miss
// also degrades to the solid sphere (and says why), the same discipline as the advection path.
const BASE_VS = `precision highp float;
uniform mat4 modelViewMatrix, projectionMatrix;
attribute vec3 position; varying vec3 vPos;
void main() { vPos = position; gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0); }`;

const BASE_FS = `precision highp float;
uniform sampler2D uMask; uniform vec2 uLonRange, uLatRange;
uniform vec3 uOcean, uLand, uBare, uLight;
varying vec3 vPos;
const float DEG = 57.29577951308232;
void main() {
  vec3 n = normalize(vPos);
  float lat = asin(clamp(n.y, -1.0, 1.0)) * DEG;
  float lon = atan(n.z, n.x) * DEG;               // inverse of sph(): x = cos p cos l, z = cos p sin l
  float gx = (lon - uLonRange.x) / max(1e-6, uLonRange.y - uLonRange.x);
  float gy = (lat - uLatRange.x) / max(1e-6, uLatRange.y - uLatRange.x);
  vec3 col = uBare;                               // outside the data box → the plain planet, no invented land
  if (gx >= 0.0 && gx <= 1.0 && gy >= 0.0 && gy <= 1.0) {
    float m = texture2D(uMask, vec2(gx, gy)).r;   // linear-filtered mask; 0.5 = the coast (as validAt())
    col = m >= 0.5 ? uOcean : uLand;
  }
  float lambert = 0.42 + 0.58 * max(0.0, dot(n, normalize(uLight)));
  gl_FragColor = vec4(col * lambert, 1.0);
}`;

function solidBase() {                             // the pre-O3 look: the fallback when there is no mask
  const m = new T.Mesh(new T.SphereGeometry(1.0, 64, 48),
                       new T.MeshPhongMaterial({ color: 0x1f2c4d, shininess: 6, specular: 0x0c1326 }));
  scene.add(m); return m;
}
function buildBase() {
  if (!MASK) return solidBase();                  // no validity mask (the eddy band) → the plain sphere
  if (!compileOK(BASE_VS, BASE_FS).ok) {
    console.warn("[flow-globe] base-layer shader rejected — solid sphere"); return solidBase();
  }
  // the mask as a linear-filtered byte texture (row-major: lat_min→lat_max = data row 0→last, matching
  // _build_data's C-order ravel and DataTexture's un-flipped rows, so gy=0 reads lat_min).
  const md = new Uint8Array(nx * ny * 4);
  for (let k = 0; k < nx * ny; k++) {
    const val = MASK[k] ? 255 : 0; md[k*4] = val; md[k*4+1] = val; md[k*4+2] = val; md[k*4+3] = 255;
  }
  const maskTex = new T.DataTexture(md, nx, ny, T.RGBAFormat, T.UnsignedByteType);
  maskTex.minFilter = maskTex.magFilter = T.LinearFilter;
  maskTex.wrapS = maskTex.wrapT = T.ClampToEdgeWrapping; maskTex.needsUpdate = true;
  const mat = new T.RawShaderMaterial({
    uniforms: { uMask: { value: maskTex },
                uLonRange: { value: new T.Vector2(cov.lon_min, cov.lon_max) },
                uLatRange: { value: new T.Vector2(cov.lat_min, cov.lat_max) },
                uOcean: { value: new T.Color(0x0b2038) }, uLand: { value: new T.Color(0x2a2a20) },
                uBare: { value: new T.Color(0x1f2c4d) }, uLight: { value: new T.Vector3(2, 1.4, 2.2) } },
    vertexShader: BASE_VS, fragmentShader: BASE_FS });
  const m = new T.Mesh(new T.SphereGeometry(1.0, 96, 64), mat);   // a touch finer for crisper coasts
  scene.add(m); return m;
}
buildBase();
scene.add(new T.AmbientLight(0x8c9bc4, 0.95));
const sun = new T.DirectionalLight(0xffffff, 0.65); sun.position.set(2, 1.4, 2.2); scene.add(sun);

// a faint graticule (meridians every 30deg, parallels every 30deg) so the sphere reads as a planet and
// the eddy band reads as the small patch it is.
(function graticule() {
  const verts = [], R = 1.002;
  const push = (a, b) => verts.push(a[0], a[1], a[2], b[0], b[1], b[2]);
  for (let lo = -180; lo < 180; lo += 30) {
    let prev = null;
    for (let la = -90; la <= 90; la += 5) { const p = sph(la, lo, R); if (prev) push(prev, p); prev = p; }
  }
  for (let la = -60; la <= 60; la += 30) {
    let prev = null;
    for (let lo = -180; lo <= 180; lo += 5) { const p = sph(la, lo, R); if (prev) push(prev, p); prev = p; }
  }
  const g = new T.BufferGeometry(); g.setAttribute("position", new T.Float32BufferAttribute(verts, 3));
  scene.add(new T.LineSegments(g, new T.LineBasicMaterial({ color: 0x34406c, transparent: true, opacity: 0.45 })));
})();

// bilinear sample of a flat (ny x nx) field at (latd, lond), clamped to the grid.
function sample(F, latd, lond) {
  let fi = (latd - lat[0]) / (lat[ny - 1] - lat[0]) * (ny - 1);
  let fj = (lond - lon[0]) / (lon[nx - 1] - lon[0]) * (nx - 1);
  fi = Math.max(0, Math.min(ny - 1, fi)); fj = Math.max(0, Math.min(nx - 1, fj));
  const i0 = Math.floor(fi), j0 = Math.floor(fj);
  const i1 = Math.min(ny - 1, i0 + 1), j1 = Math.min(nx - 1, j0 + 1);
  const di = fi - i0, dj = fj - j0;
  const a = F[i0 * nx + j0], b = F[i0 * nx + j1], c = F[i1 * nx + j0], e = F[i1 * nx + j1];
  return (a * (1 - dj) + b * dj) * (1 - di) + (c * (1 - dj) + e * dj) * di;
}

// per-cell validity (the O1 mask): bilinearly sampled 0/1 against a 0.5 threshold — the SAME rule the
// GPU path applies to the linear-filtered mask channel, so the two paths agree at the coastline (to
// within half a source cell). No mask → everywhere valid → the exact pre-mask behaviour.
function validAt(latd, lond) { return !MASK || sample(MASK, latd, lond) >= 0.5; }

function mix3(a, b, s) { return [a[0] + (b[0] - a[0]) * s, a[1] + (b[1] - a[1]) * s, a[2] + (b[2] - a[2]) * s]; }
// Two ramps (§9.6 O3c). SEQ=false → RdBu_r, a DIVERGING blue→white→red for a signed θ field (Rung A/B).
// SEQ=true → a SEQUENTIAL blue→cyan→green→yellow speed ramp: a diverging map bleaches mid-speed and is
// semantically wrong for a 0→max field, so ocean speed gets its own monotone ramp (the demo opts in).
function cmap(t) {
  if (SEQ) {
    const s0 = [0.16, 0.34, 0.60], s1 = [0.24, 0.63, 0.75], s2 = [0.55, 0.86, 0.60], s3 = [0.98, 0.95, 0.55];
    if (t < 0.4) return mix3(s0, s1, t / 0.4);
    if (t < 0.75) return mix3(s1, s2, (t - 0.4) / 0.35);
    return mix3(s2, s3, (t - 0.75) / 0.25);
  }
  const lo = [0.23, 0.31, 0.75], mid = [0.96, 0.96, 0.96], hi = [0.79, 0.16, 0.18];
  return t < 0.5 ? mix3(lo, mid, t / 0.5) : mix3(mid, hi, (t - 0.5) / 0.5);
}

// --------------------------------------------------------------------------------------------------- #
// Particle advection — TWO implementations of the SAME lat/lon-metric integration: a GPU ping-pong pass
// (the default; state lives in a float texture, advanced entirely on the GPU) and a CPU step() loop (the
// fallback). We cannot run WebGL in CI, so the design rule is: a GPU failure must degrade to the working
// CPU globe — never a blank one — and the console must say which path ran and why. The active path's
// per-frame update is bound to `tick`; the appearance sliders to `applySize`/`applyOpacity`/`applySharp`.
// --------------------------------------------------------------------------------------------------- #
const N = D.n_particles;
let tick = null, applySize = null, applyOpacity = null, applySharp = null, applyDensity = null, onResizeHook = null;

// a deterministic, fixed-seed LCG → a reproducible *initial* state for either path (the motion is then
// stochastic, which the honest-by-disclosure carve-out permits — there is no byte-golden on this figure).
let seed = 1234567;
function rnd() { seed = (1103515245 * seed + 12345) % 2147483648; return seed / 2147483648; }

// ===== GPU ping-pong advection (the default) ======================================================= #
// Particle state is an RGBA32F texture, one texel per particle = (lon, lat, age, life). Each frame an
// off-screen fragment shader (UPDATE_FS) reads the current state texture, advects every particle by the
// same metric the CPU path uses, and writes the next state into a second target; the two targets
// ping-pong. A Points cloud then draws the particles, its vertex shader (DRAW_VS) reading each particle's
// position straight from the state texture. The velocity (+θ +mask) field rides along as a half-float
// DataTexture, RGBA = (u, v, θ, mask) — the O1 validity mask on the formerly-free 4th channel, so the
// respawn/advection logic rejects land texels with NO new texture. All shaders are GLSL1 (texture2D /
// attribute / varying) on
// RawShaderMaterial, so the source we hand three is exactly what we can validate against the live context.
const QUAD_VS = `precision highp float;
attribute vec3 position; attribute vec2 uv; varying vec2 vUv;
void main() { vUv = uv; gl_Position = vec4(position.xy, 0.0, 1.0); }`;

const INIT_FS = `precision highp float;
uniform sampler2D uSeed; varying vec2 vUv;
void main() { gl_FragColor = texture2D(uSeed, vUv); }`;

const UPDATE_FS = `precision highp float;
uniform sampler2D uState, uVel;
uniform float uDt, uAccel, uRadius, uRandom, uSpeedMax, uSeedFloor;
uniform vec2 uLonRange, uLatRange, uVelGrid;
varying vec2 vUv;
const float DEG = 57.29577951308232, RAD = 0.017453292519943295;
float hash(vec2 p) { return fract(sin(dot(p, vec2(12.9898, 78.233))) * 43758.5453); }
vec2 velUV(float lon, float lat) {
  float gx = clamp((lon - uLonRange.x) / max(1e-6, uLonRange.y - uLonRange.x), 0.0, 1.0);
  float gy = clamp((lat - uLatRange.x) / max(1e-6, uLatRange.y - uLatRange.x), 0.0, 1.0);
  return vec2((gx * (uVelGrid.x - 1.0) + 0.5) / uVelGrid.x,
              (gy * (uVelGrid.y - 1.0) + 0.5) / uVelGrid.y);   // texel-centre sampling, matches CPU sample()
}
void main() {
  vec4 st = texture2D(uState, vUv);
  float lon = st.x, lat = st.y, age = st.z, life = st.w;
  vec4 vel = texture2D(uVel, velUV(lon, lat));
  float cosp = max(0.05, cos(lat * RAD));
  lon += DEG * (vel.x / (uRadius * cosp)) * uAccel * uDt;   // dλ/dt = u/(a cosφ)
  lat += DEG * (vel.y / uRadius) * uAccel * uDt;            // dφ/dt = v/a
  age += uDt;
  bool gone = lon < uLonRange.x || lon > uLonRange.y || lat < uLatRange.x || lat > uLatRange.y
           || texture2D(uVel, velUV(lon, lat)).w < 0.5;     // drifted onto a masked (land) texel → recycle
  if (age > life || gone) {                                 // respawn inside coverage (never paints bare globe)
    lat = uLatRange.x + hash(vUv + vec2(uRandom, 0.123)) * (uLatRange.y - uLatRange.x);
    lon = uLonRange.x + hash(vUv + vec2(0.456, uRandom)) * (uLonRange.y - uLonRange.x);
    age = 0.0; life = 2.0 + hash(vUv * 1.7 + uRandom) * 4.0;
    // a draw that landed on a masked texel — OR one rejected by the speed weighting — is kept INVISIBLE
    // (age past life ⇒ zero fade in DRAW_VS) and re-rolled next frame with a fresh uRandom. A shader cannot
    // loop a rejection sample, so this invisible-retry idiom is the GPU form of BOTH "seed only valid cells"
    // (converges in ~2 frames even at OSCAR's 44% land) and "seed ∝ speed" (§9.6 O3c) — the two criteria
    // compose: land is rejected outright, then valid cells are accepted with probability |u,v|/speed_max
    // (floored so calm water still gets a light ambient fill), concentrating particles in the fast currents.
    vec4 rv = texture2D(uVel, velUV(lon, lat));
    if (rv.w < 0.5) age = life + 1.0;                       // masked (land) → invisible, re-roll
    else if (hash(vUv * 2.3 + vec2(uRandom, uRandom * 1.7))
             > max(uSeedFloor, length(rv.xy) / max(1e-6, uSpeedMax))) age = life + 1.0;   // speed-weighted reject
  }
  gl_FragColor = vec4(lon, lat, age, life);
}`;

const DRAW_VS = `precision highp float;
uniform mat4 modelViewMatrix, projectionMatrix;
uniform sampler2D uState, uVel;
uniform float uSize, uScale, uRadius, uHasScalar, uSeq, uDensity;
uniform vec2 uLonRange, uLatRange, uVelGrid, uScalarRange;
attribute vec3 position;                 // .xy = this particle's texel-centre UV into uState
attribute float aSeq;                    // this particle's index / N ∈ [0,1] — the density-knob cut
varying vec3 vColor; varying float vFade;
const float RAD = 0.017453292519943295;
vec2 velUV(float lon, float lat) {
  float gx = clamp((lon - uLonRange.x) / max(1e-6, uLonRange.y - uLonRange.x), 0.0, 1.0);
  float gy = clamp((lat - uLatRange.x) / max(1e-6, uLatRange.y - uLatRange.x), 0.0, 1.0);
  return vec2((gx * (uVelGrid.x - 1.0) + 0.5) / uVelGrid.x,
              (gy * (uVelGrid.y - 1.0) + 0.5) / uVelGrid.y);
}
vec3 cmap(float t, float seq) {           // seq=0 → RdBu_r (θ); seq=1 → sequential speed ramp (§9.6 O3c)
  if (seq > 0.5) {
    vec3 s0 = vec3(0.16, 0.34, 0.60), s1 = vec3(0.24, 0.63, 0.75), s2 = vec3(0.55, 0.86, 0.60), s3 = vec3(0.98, 0.95, 0.55);
    if (t < 0.4) return mix(s0, s1, t / 0.4);
    if (t < 0.75) return mix(s1, s2, (t - 0.4) / 0.35);
    return mix(s2, s3, (t - 0.75) / 0.25);
  }
  vec3 lo = vec3(0.23, 0.31, 0.75), mid = vec3(0.96, 0.96, 0.96), hi = vec3(0.79, 0.16, 0.18);
  return t < 0.5 ? mix(lo, mid, t / 0.5) : mix(mid, hi, (t - 0.5) / 0.5);
}
void main() {
  if (aSeq > uDensity) { gl_Position = vec4(2.0, 2.0, 2.0, 1.0); gl_PointSize = 0.0; vFade = 0.0; vColor = vec3(0.0); return; }
  vec4 st = texture2D(uState, position.xy);
  float lon = st.x, lat = st.y, age = st.z, life = st.w;
  float p = lat * RAD, l = lon * RAD;
  vec3 xyz = vec3(cos(p) * cos(l), sin(p), cos(p) * sin(l)) * 1.012;
  vec4 mv = modelViewMatrix * vec4(xyz, 1.0);
  gl_Position = projectionMatrix * mv;
  gl_PointSize = uSize * uScale / max(0.001, -mv.z);        // replicate three's size attenuation
  float t = 0.5;
  if (uHasScalar > 0.5) {
    float theta = texture2D(uVel, velUV(lon, lat)).z;
    t = clamp((theta - uScalarRange.x) / max(1e-6, uScalarRange.y - uScalarRange.x), 0.0, 1.0);
  }
  vColor = cmap(t, uSeq);
  vFade = min(1.0, age / 0.3) * min(1.0, (life - age) / 0.5);   // fade in from spawn, out toward death
}`;

const DRAW_FS = `precision highp float;
uniform float uOpacity, uSharp; varying vec3 vColor; varying float vFade;
void main() {
  float r = length(gl_PointCoord - 0.5) * 2.0;             // 0 at centre → 1 at the sprite edge
  if (r > 1.0) discard;                                    // round, not a square GL point
  float core = clamp(uSharp, 0.0, 0.97);                   // opaque out to the core radius, then linear to 0
  float a = 1.0 - clamp((r - core) / max(1e-4, 1.0 - core), 0.0, 1.0);
  float alpha = a * vFade * uOpacity;
  if (alpha < 0.01) discard;
  gl_FragColor = vec4(vColor, alpha);
}`;

function buildGPU() {
  // velocity (+θ +mask) as a linear-filtered half-float texture: RGBA = (u, v, θ, mask). Half precision
  // is ample for a few-m/s velocity, a colour channel, and a 0/1 validity flag (linear filtering blends
  // the flag near coasts; the 0.5 threshold in the shaders keeps the boundary to half a texel).
  const toHalf = T.DataUtils.toHalfFloat;
  const velData = new Uint16Array(nx * ny * 4);
  for (let k = 0; k < nx * ny; k++) {
    velData[k * 4] = toHalf(U[k]); velData[k * 4 + 1] = toHalf(V[k]);
    velData[k * 4 + 2] = toHalf(S ? S[k] : 0); velData[k * 4 + 3] = toHalf(MASK ? MASK[k] : 1);
  }
  const velTex = new T.DataTexture(velData, nx, ny, T.RGBAFormat, T.HalfFloatType);
  velTex.minFilter = velTex.magFilter = T.LinearFilter;
  velTex.wrapS = velTex.wrapT = T.ClampToEdgeWrapping; velTex.needsUpdate = true;

  // pack the particles into the smallest square float texture that holds them (one texel each).
  const texSize = Math.ceil(Math.sqrt(N)), M = texSize * texSize;
  const seedData = new Float32Array(M * 4);
  for (let i = 0; i < M; i++) {                           // same LCG draw order as spawn(): lat, lon, life
    const la = cov.lat_min + rnd() * (cov.lat_max - cov.lat_min);
    const lo = cov.lon_min + rnd() * (cov.lon_max - cov.lon_min);
    seedData[i * 4] = lo; seedData[i * 4 + 1] = la; seedData[i * 4 + 2] = 0; seedData[i * 4 + 3] = 2.0 + rnd() * 4.0;
  }
  const seedTex = new T.DataTexture(seedData, texSize, texSize, T.RGBAFormat, T.FloatType);
  seedTex.minFilter = seedTex.magFilter = T.NearestFilter; seedTex.needsUpdate = true;

  // the two ping-pong targets — FloatType so positions keep full precision, NearestFilter so a particle
  // reads exactly its own texel (LinearFilter would blend neighbouring particles into garbage motion).
  const rtOpts = { type: T.FloatType, format: T.RGBAFormat, minFilter: T.NearestFilter,
                   magFilter: T.NearestFilter, depthBuffer: false, stencilBuffer: false,
                   wrapS: T.ClampToEdgeWrapping, wrapT: T.ClampToEdgeWrapping };
  let rtCur = new T.WebGLRenderTarget(texSize, texSize, rtOpts);
  let rtNext = new T.WebGLRenderTarget(texSize, texSize, rtOpts);

  // an off-screen full-screen quad whose fragment shader advects every particle once per frame.
  const quadScene = new T.Scene(), quadCam = new T.Camera();
  const quad = new T.Mesh(new T.PlaneGeometry(2, 2), new T.RawShaderMaterial(
    { uniforms: { uSeed: { value: seedTex } }, vertexShader: QUAD_VS, fragmentShader: INIT_FS }));
  quad.frustumCulled = false;            // the shader writes clip coords directly — never cull the update pass
  quadScene.add(quad);
  const updateMat = new T.RawShaderMaterial({
    uniforms: { uState: { value: null }, uVel: { value: velTex }, uDt: { value: 0 }, uAccel: { value: ACCEL },
                uRadius: { value: A }, uRandom: { value: 0 }, uSpeedMax: { value: SPEEDMAX },
                uSeedFloor: { value: SEED_FLOOR },
                uLonRange: { value: new T.Vector2(cov.lon_min, cov.lon_max) },
                uLatRange: { value: new T.Vector2(cov.lat_min, cov.lat_max) },
                uVelGrid: { value: new T.Vector2(nx, ny) } },
    vertexShader: QUAD_VS, fragmentShader: UPDATE_FS });
  renderer.setRenderTarget(rtCur); renderer.render(quadScene, quadCam);   // seed both targets from the LCG state
  renderer.setRenderTarget(rtNext); renderer.render(quadScene, quadCam);
  renderer.setRenderTarget(null);
  // a DIAGNOSTIC (not a gate): read particle 0 straight back out of the float target. The feature gate +
  // compile check can't catch a driver that advertises EXT_color_buffer_float yet renders an INCOMPLETE
  // float target — there the state round-trips as zeros/NaN and the particles freeze at seed while the
  // console still says "GPU active". Logging the round-tripped state turns that into a paste-back tell.
  try {
    const probe = new Float32Array(4);
    renderer.readRenderTargetPixels(rtCur, 0, 0, 1, 1, probe);
    console.log("[flow-globe] GPU state round-trip — particle 0 (lon,lat,age,life) = " +
                probe[0].toFixed(2) + ", " + probe[1].toFixed(2) + ", " + probe[2].toFixed(2) + ", " +
                probe[3].toFixed(2) + "  (zeros/NaN here ⇒ the float target didn't round-trip ⇒ frozen)");
  } catch (e) { console.warn("[flow-globe] GPU state read-back unavailable:", e && e.message); }
  quad.material = updateMat;

  // the visible Points: each vertex's `position.xy` is its texel-centre UV; DRAW_VS reads the particle's
  // (lon, lat) from the state texture. Bounds live in the texture, so the CPU bounding sphere is
  // meaningless → disable frustum culling (else three culls the whole cloud at most camera angles).
  const refs = new Float32Array(M * 3), seqs = new Float32Array(M);
  for (let i = 0; i < M; i++) {
    refs[i * 3] = ((i % texSize) + 0.5) / texSize; refs[i * 3 + 1] = (Math.floor(i / texSize) + 0.5) / texSize;
    seqs[i] = (i + 0.5) / M;                              // stable [0,1] rank → the density-knob cut
  }
  const gpuGeo = new T.BufferGeometry();
  gpuGeo.setAttribute("position", new T.BufferAttribute(refs, 3));
  gpuGeo.setAttribute("aSeq", new T.BufferAttribute(seqs, 1));
  const drawMat = new T.RawShaderMaterial({
    uniforms: { uState: { value: rtCur.texture }, uVel: { value: velTex }, uSize: { value: D.particle_size },
                uScale: { value: 0.5 * (renderer.domElement.height || H) }, uOpacity: { value: D.particle_opacity },
                uSharp: { value: D.particle_sharpness }, uRadius: { value: A }, uHasScalar: { value: S ? 1 : 0 },
                uSeq: { value: SEQ ? 1 : 0 }, uDensity: { value: density },
                uLonRange: { value: new T.Vector2(cov.lon_min, cov.lon_max) },
                uLatRange: { value: new T.Vector2(cov.lat_min, cov.lat_max) },
                uVelGrid: { value: new T.Vector2(nx, ny) }, uScalarRange: { value: new T.Vector2(smin, smax) } },
    vertexShader: DRAW_VS, fragmentShader: DRAW_FS,
    transparent: true, depthTest: true, depthWrite: false, blending: T.NormalBlending });   // occlude far side
  const gpuPoints = new T.Points(gpuGeo, drawMat); gpuPoints.frustumCulled = false;
  scene.add(gpuPoints);

  tick = function (dt) {
    updateMat.uniforms.uDt.value = dt; updateMat.uniforms.uRandom.value = Math.random();
    updateMat.uniforms.uState.value = rtCur.texture;
    renderer.setRenderTarget(rtNext); renderer.render(quadScene, quadCam); renderer.setRenderTarget(null);
    const swap = rtCur; rtCur = rtNext; rtNext = swap;     // ping-pong
    drawMat.uniforms.uState.value = rtCur.texture;
  };
  applySize = (v) => { drawMat.uniforms.uSize.value = v; };
  applyOpacity = (v) => { drawMat.uniforms.uOpacity.value = v; };
  applySharp = (v) => { drawMat.uniforms.uSharp.value = v; };
  applyDensity = (v) => { density = v; drawMat.uniforms.uDensity.value = v; };   // the §9.5 density knob
  onResizeHook = () => { drawMat.uniforms.uScale.value = 0.5 * (renderer.domElement.height || H); };
}

// ===== CPU advection (the fallback) ================================================================ #
// The original, dependency-free path: integrate every particle in JS each frame and re-upload the
// position/colour buffers. Identical metric and look to the GPU path; this is what runs if WebGL2 /
// float-render targets / the shaders are unavailable, so the showcase is never a blank globe.
const pos = new Float32Array(N * 3), col = new Float32Array(N * 4);   // colour is RGBA — alpha carries the fade
const pLat = new Float32Array(N), pLon = new Float32Array(N), pAge = new Float32Array(N), pLife = new Float32Array(N);
function spawn(i) {
  for (let tries = 0; tries < 40; tries++) {       // reject land, THEN accept ∝ speed — the same two
    pLat[i] = cov.lat_min + rnd() * (cov.lat_max - cov.lat_min);   // criteria the GPU respawn composes
    pLon[i] = cov.lon_min + rnd() * (cov.lon_max - cov.lon_min);
    if (!validAt(pLat[i], pLon[i])) continue;      // masked (land) → re-roll (~2 tries at OSCAR's 44% land)
    const spd = Math.hypot(sample(U, pLat[i], pLon[i]), sample(V, pLat[i], pLon[i]));
    if (rnd() < Math.max(SEED_FLOOR, spd / SPEEDMAX)) break;   // speed-weighted accept → fast currents dominate
  }
  pAge[i] = 0; pLife[i] = 2.0 + rnd() * 4.0;       // seconds before respawn (staggered so trails don't blink together)
  // a near-all-masked field can exhaust the tries: keep the particle invisible (age past life ⇒ zero
  // alpha in tick) and re-roll next frame — the same idiom as the GPU respawn's masked-texel branch.
  if (!validAt(pLat[i], pLon[i])) pAge[i] = pLife[i] + 1;
}
// a ROUND particle sprite (a radial alpha falloff): a square GL point is the amateur tell. The sprite is
// white, so the per-vertex temperature colour survives; `sharp` ∈ [0,1] sets the opaque-core radius
// (0 = soft bloom, ~1 = near-hard disc). Cheap to regenerate (a 64² canvas), so the slider rebuilds it.
function particleSprite(sharp) {
  const s = 64, cvs = document.createElement("canvas"); cvs.width = cvs.height = s;
  const g = cvs.getContext("2d");
  const core = Math.max(0, Math.min(0.97, sharp));   // 0.97 cap keeps a 1-px anti-aliased rim (no jaggies)
  const grd = g.createRadialGradient(s / 2, s / 2, 0, s / 2, s / 2, s / 2);
  grd.addColorStop(0.0, "rgba(255,255,255,1)");
  grd.addColorStop(core, "rgba(255,255,255,1)");     // opaque out to `core`, then fall to transparent
  grd.addColorStop(1.0, "rgba(255,255,255,0)");
  g.fillStyle = grd; g.fillRect(0, 0, s, s);
  return new T.CanvasTexture(cvs);
}
function buildCPU() {
  const geo = new T.BufferGeometry();
  geo.setAttribute("position", new T.BufferAttribute(pos, 3));
  geo.setAttribute("color", new T.BufferAttribute(col, 4));            // RGBA: the 4th channel is the spawn/death fade
  const pmat = new T.PointsMaterial({ size: D.particle_size, map: particleSprite(D.particle_sharpness),
                                      vertexColors: true, transparent: true, opacity: D.particle_opacity,
                                      depthWrite: false });
  scene.add(new T.Points(geo, pmat));
  for (let i = 0; i < N; i++) spawn(i);
  tick = function (dt) {
    const nActive = Math.floor(N * density);            // the §5 density knob: draw the first `nActive`
    for (let i = 0; i < N; i++) {
      const la = pLat[i], lo = pLon[i];
      const uu = sample(U, la, lo), vv = sample(V, la, lo);
      const cosp = Math.max(0.05, Math.cos(la * RAD));
      pLon[i] += DEG * (uu / (A * cosp)) * ACCEL * dt;   // dλ/dt = u/(a cosφ)
      pLat[i] += DEG * (vv / A) * ACCEL * dt;            // dφ/dt = v/a
      pAge[i] += dt;
      const gone = pLat[i] < cov.lat_min || pLat[i] > cov.lat_max || pLon[i] < cov.lon_min || pLon[i] > cov.lon_max
                || !validAt(pLat[i], pLon[i]);     // drifted onto a masked (land) cell → recycle
      if (pAge[i] > pLife[i] || gone) spawn(i);
      const xyz = sph(pLat[i], pLon[i], 1.012);
      pos[i * 3] = xyz[0]; pos[i * 3 + 1] = xyz[1]; pos[i * 3 + 2] = xyz[2];
      let t = S ? (sample(S, pLat[i], pLon[i]) - smin) / ((smax - smin) || 1) : 0.5;
      t = Math.max(0, Math.min(1, t));
      const c = cmap(t);                                 // full-brightness colour, ALWAYS (never dimmed toward black)
      // the fade lives in ALPHA, not RGB: a fresh particle fades in from transparent, a dying one fully out.
      let alpha = Math.min(1, pAge[i] / 0.3) * Math.min(1, (pLife[i] - pAge[i]) / 0.5);
      if (i >= nActive) alpha = 0;                       // density knob hides the tail (advection still runs)
      col[i * 4] = c[0]; col[i * 4 + 1] = c[1]; col[i * 4 + 2] = c[2]; col[i * 4 + 3] = alpha;
    }
    geo.attributes.position.needsUpdate = true;
    geo.attributes.color.needsUpdate = true;
  };
  applySize = (v) => { pmat.size = v; };
  applyOpacity = (v) => { pmat.opacity = v; };
  applySharp = (v) => { const old = pmat.map; pmat.map = particleSprite(v); pmat.needsUpdate = true; if (old) old.dispose(); };
  applyDensity = (v) => { density = v; };               // the §5 density knob (CPU: hides the tail in tick)
}

// ===== pick the path: GPU if WebGL2 + float-render targets + the shaders compile; else CPU ========== #
// We cannot run WebGL here, so we validate the GLSL against the live context (compile + link) before
// trusting it — three logs but does NOT throw on a link failure — and degrade to CPU on any miss, naming
// the reason in the console so a blind play-through is debuggable from a paste-back, not a guess.
function compileOK(srcV, srcF) {
  const gl = renderer.getContext();
  function sh(type, src) {
    const s = gl.createShader(type); gl.shaderSource(s, src); gl.compileShader(s);
    const ok = gl.getShaderParameter(s, gl.COMPILE_STATUS);
    return { s, ok, log: ok ? "" : (gl.getShaderInfoLog(s) || "") };
  }
  const v = sh(gl.VERTEX_SHADER, srcV), f = sh(gl.FRAGMENT_SHADER, srcF);
  let ok = v.ok && f.ok, log = v.log + f.log;
  if (ok) {
    const p = gl.createProgram(); gl.attachShader(p, v.s); gl.attachShader(p, f.s); gl.linkProgram(p);
    ok = gl.getProgramParameter(p, gl.LINK_STATUS); if (!ok) log += (gl.getProgramInfoLog(p) || ""); gl.deleteProgram(p);
  }
  gl.deleteShader(v.s); gl.deleteShader(f.s); return { ok, log };
}
let useGPU = false, why = "";
if (!renderer.capabilities.isWebGL2) why = "WebGL2 unavailable";
else if (!renderer.getContext().getExtension("EXT_color_buffer_float")) why = "EXT_color_buffer_float unavailable";
else {
  const u = compileOK(QUAD_VS, UPDATE_FS), d = compileOK(DRAW_VS, DRAW_FS), z = compileOK(QUAD_VS, INIT_FS);
  if (u.ok && d.ok && z.ok) useGPU = true;
  else { why = "shader compile/link failed"; console.warn("[flow-globe] GPU shaders rejected:\n" + u.log + d.log + z.log); }
}
if (useGPU) { try { buildGPU(); } catch (e) { useGPU = false; why = "GPU init threw: " + (e && e.message); console.warn("[flow-globe]", e); } }
if (useGPU) console.log("[flow-globe] GPU ping-pong advection active");
else { buildCPU(); console.warn("[flow-globe] CPU advection fallback active — " + why); }

// --- hand-rolled spherical orbit camera (drag → azimuth/elevation, wheel → radius) --- #
let az = D.center_lon * RAD, el = D.center_lat * RAD, radius = 3.0;
const ELIM = 85 * RAD;
function placeCam() {
  const ce = Math.cos(el), se = Math.sin(el), ca = Math.cos(az), sa = Math.sin(az);
  camera.position.set(radius * ce * ca, radius * se, radius * ce * sa);
  camera.lookAt(0, 0, 0);
}
placeCam();
let drag = false, px = 0, py = 0;
canvas.addEventListener("pointerdown", (e) => { drag = true; px = e.clientX; py = e.clientY; });
window.addEventListener("pointerup", () => { drag = false; });
window.addEventListener("pointermove", (e) => {
  if (!drag) return;
  az -= (e.clientX - px) * 0.005; el += (e.clientY - py) * 0.005;
  el = Math.max(-ELIM, Math.min(ELIM, el));
  px = e.clientX; py = e.clientY; placeCam();
});
canvas.addEventListener("wheel", (e) => {
  e.preventDefault();
  radius *= (e.deltaY > 0 ? 1.1 : 0.9);
  radius = Math.max(1.4, Math.min(8, radius)); placeCam();
}, { passive: false });
window.addEventListener("resize", () => {
  [W, H] = size(); renderer.setSize(W, H, false); camera.aspect = W / H; camera.updateProjectionMatrix();
  if (onResizeHook) onResizeHook();                       // the GPU draw needs its point-size scale refreshed
});

// --- live controls: size + opacity + edge sharpness + DENSITY (§9.5, the O3c unlock). Each dispatches to
// the ACTIVE advection path (GPU → a uniform; CPU → the material / a var), so the sliders work the same
// whichever path the browser picked. Density is the second consumer's first real control knob — the ocean
// producer is why §9.5's control-surface seam opened. Trail length arrives with the O3b trails; colour-ramp
// *menus* and shape menus stay a deferred seam (the ramp DEFAULT is now producer-chosen, but its live UI is
// not built — speculative, no consumer yet). --- #
const sizeR = document.getElementById("sizeRange"), opacR = document.getElementById("opacityRange");
const sharpR = document.getElementById("sharpRange"), densR = document.getElementById("densityRange");
if (sizeR) sizeR.addEventListener("input", () => applySize(parseFloat(sizeR.value)));
if (opacR) opacR.addEventListener("input", () => applyOpacity(parseFloat(opacR.value)));
if (sharpR) sharpR.addEventListener("input", () => applySharp(parseFloat(sharpR.value)));
if (densR) densR.addEventListener("input", () => applyDensity(parseFloat(densR.value)));

let last = performance.now();
function loop(now) {
  let dt = (now - last) / 1000; last = now; if (dt > 0.1) dt = 0.1;
  tick(dt); renderer.render(scene, camera); requestAnimationFrame(loop);
}
requestAnimationFrame(loop);
"""


def flow_globe_html(field: FlowField, *, title: str = "planet-sim — eddy flow-globe (showcase)",
                    subtitle: str = "the emergent eddy life cycle as a particle flow on a real, rotatable "
                                    "planet (§9.5 Rung C — the showcase)",
                    n_particles: int = DEFAULT_N_PARTICLES,
                    particle_size: float = DEFAULT_PARTICLE_SIZE,
                    particle_opacity: float = DEFAULT_PARTICLE_OPACITY,
                    particle_sharpness: float = DEFAULT_PARTICLE_SHARPNESS,
                    crossing_seconds: float = _BAND_CROSSING_SECONDS,
                    colormap: str = "RdBu_r") -> str:
    """Render ``field`` as one deterministic, self-contained three.js HTML page (data + three.js inlined).

    The disclaimer (``field.honesty``) is written into a **visible** ``<div class="disclaimer">`` in the
    static DOM — not merely a JS comment — because under the honest-by-disclosure carve-out it *is* the
    entire license, and is the one thing machine-checked. Opens straight off ``file://`` (no network).

    ``particle_size`` / ``particle_opacity`` / ``particle_sharpness`` set the shipped *defaults*; each
    flows to both the material init and its live slider's initial position (one source, no drift), so a
    viewer can fine-tune appearance in the browser without regenerating, and a notebook can ship different
    defaults. Sharpness ∈ [0,1] is the sprite's opaque-core radius (0 = soft bloom, 1 = near-hard disc).
    ``crossing_seconds`` sets the visual pace (wall-clock seconds for the fastest particle to cross the
    field's longitude span); the default is the eddy band's tuning, so pre-O2 artifacts are unchanged.
    ``colormap`` picks the particle ramp: ``"RdBu_r"`` (default, a **diverging** blue→white→red for a
    signed θ field) or ``"speed"`` (a **sequential** blue→cyan→green→yellow for a 0→max speed field — a
    diverging map bleaches mid-speed, so a speed scalar wants its own monotone ramp). It sets the shipped
    default only; ``FlowField`` is untouched (the §9.3 win) and the ocean demo opts in (§9.6 O3c).
    """
    sequential = colormap == "speed"
    data_json = json.dumps(
        _build_data(field, n_particles, particle_size, particle_opacity, particle_sharpness,
                    crossing_seconds, sequential),
        separators=(",", ":"), ensure_ascii=False)
    disclaimer = html.escape(field.honesty)
    return (
        "<!doctype html>\n<html lang=\"en\">\n<head>\n"
        "<meta charset=\"utf-8\">\n"
        "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
        f"<title>{html.escape(title)}</title>\n"
        f"<style>\n{_CSS}</style>\n</head>\n<body>\n"
        "<div id=\"stage\">\n"
        "  <canvas id=\"globe\"></canvas>\n"
        f"  <div class=\"title\"><h1>{html.escape(title)}</h1><p>{html.escape(subtitle)}</p></div>\n"
        "  <div class=\"hint\">drag to rotate · scroll to zoom</div>\n"
        "  <div class=\"controls\">\n"
        f"    <label>Particle size<input id=\"sizeRange\" type=\"range\" min=\"0.01\" max=\"0.1\" "
        f"step=\"0.005\" value=\"{particle_size}\"></label>\n"
        f"    <label>Opacity<input id=\"opacityRange\" type=\"range\" min=\"0.1\" max=\"1\" "
        f"step=\"0.05\" value=\"{particle_opacity}\"></label>\n"
        f"    <label>Edge sharpness<input id=\"sharpRange\" type=\"range\" min=\"0\" max=\"1\" "
        f"step=\"0.05\" value=\"{particle_sharpness}\"></label>\n"
        f"    <label>Density<input id=\"densityRange\" type=\"range\" min=\"0.1\" max=\"1\" "
        f"step=\"0.05\" value=\"1\"></label>\n"
        "  </div>\n"
        f"  <div class=\"disclaimer\" id=\"disclaimer\"><strong>Illustrative showcase — read this.</strong> "
        f"{disclaimer}</div>\n"
        "</div>\n"
        f"<script>{_three_js()}</script>\n"
        f"<script>\nwindow.FLOW_DATA = {data_json};\n</script>\n"
        f"<script>{_APP_JS}</script>\n"
        "</body>\n</html>\n"
    )


def save_flow_globe_html(field: FlowField, path, **kwargs) -> Path:
    """Render ``field`` and write the standalone HTML showcase (LF newlines, UTF-8). Returns the path."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(flow_globe_html(field, **kwargs), encoding="utf-8", newline="\n")
    return path
