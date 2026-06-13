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
DEFAULT_PARTICLE_SIZE = 0.035    # point size in sphere-radius units (the size slider's initial value)
DEFAULT_PARTICLE_OPACITY = 0.95  # overall particle opacity ceiling (the opacity slider's initial value)


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


def _build_data(field: FlowField, n_particles: int,
                particle_size: float, particle_opacity: float) -> dict:
    """Pack a :class:`FlowField` into the JSON the in-browser renderer consumes (compact, rounded)."""
    lat = np.asarray(field.lat, dtype=float)
    lon = np.asarray(field.lon, dtype=float)
    u = np.asarray(field.u, dtype=float)
    v = np.asarray(field.v, dtype=float)
    scalar = None if field.scalar is None else np.asarray(field.scalar, dtype=float)

    center_lat = float(lat.mean())
    center_lon = float(lon.mean())
    span_lon = float(lon.max() - lon.min()) or 1.0
    umax = float(np.max(np.abs(u))) or 1.0
    # `accel` (model-seconds per wall-clock second) is auto-scaled so the *fastest* particle crosses the
    # band's longitude width in ~`_BAND_CROSSING_SECONDS` — readable streaming regardless of the field's
    # actual speeds. In-browser step: Δλ_deg = deg(u/(a·cosφ)) · accel · dt_real (the honest metric × a
    # purely-visual time-acceleration; the showcase's physics-fidelity is relaxed, ADR 0002 note).
    dlamdt_max = math.degrees(umax / (field.radius_m * math.cos(math.radians(center_lat))))
    accel = (span_lon / _BAND_CROSSING_SECONDS) / dlamdt_max if dlamdt_max > 0 else 1.0

    def flat(arr, nd):
        return [round(float(x), nd) for x in np.asarray(arr).ravel(order="C")]

    return {
        "lat": flat(lat, 4), "lon": flat(lon, 4),
        "u": flat(u, 3), "v": flat(v, 3),
        "scalar": (flat(scalar, 3) if scalar is not None else None),
        "scalar_min": round(float(scalar.min()), 3) if scalar is not None else 0.0,
        "scalar_max": round(float(scalar.max()), 3) if scalar is not None else 1.0,
        "scalar_label": field.scalar_label,
        "radius_m": field.radius_m,
        "accel": accel,
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
# this block. The particle advection (lat/lon metric integration + respawn) and the spherical orbit
# camera are ORIGINAL; three.js supplies only the scene / camera / WebGL renderer (vendored inline).
_APP_JS = r"""
const D = window.FLOW_DATA, T = window.THREE;
const lat = D.lat, lon = D.lon, ny = lat.length, nx = lon.length;
const U = D.u, V = D.v, S = D.scalar, smin = D.scalar_min, smax = D.scalar_max;
const A = D.radius_m, ACCEL = D.accel, cov = D.coverage;
const DEG = 180 / Math.PI, RAD = Math.PI / 180;

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
const base = new T.Mesh(new T.SphereGeometry(1.0, 64, 48),
                        new T.MeshPhongMaterial({ color: 0x1f2c4d, shininess: 6, specular: 0x0c1326 }));
scene.add(base);
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

// RdBu_r colour ramp: 0 = cool blue, 1 = warm red (matches the Rung-A/B theta scale).
function cmap(t) {
  const lo = [0.23, 0.31, 0.75], mid = [0.96, 0.96, 0.96], hi = [0.79, 0.16, 0.18];
  if (t < 0.5) { const s = t / 0.5; return [lo[0] + (mid[0] - lo[0]) * s, lo[1] + (mid[1] - lo[1]) * s, lo[2] + (mid[2] - lo[2]) * s]; }
  const s = (t - 0.5) / 0.5; return [mid[0] + (hi[0] - mid[0]) * s, mid[1] + (hi[1] - mid[1]) * s, mid[2] + (hi[2] - mid[2]) * s];
}

// deterministic particles (a fixed-seed LCG → a reproducible artifact even without a byte-golden test).
const N = D.n_particles;
const pos = new Float32Array(N * 3), col = new Float32Array(N * 4);   // colour is RGBA — alpha carries the fade
const pLat = new Float32Array(N), pLon = new Float32Array(N), pAge = new Float32Array(N), pLife = new Float32Array(N);
let seed = 1234567;
function rnd() { seed = (1103515245 * seed + 12345) % 2147483648; return seed / 2147483648; }
function spawn(i) {
  pLat[i] = cov.lat_min + rnd() * (cov.lat_max - cov.lat_min);
  pLon[i] = cov.lon_min + rnd() * (cov.lon_max - cov.lon_min);
  pAge[i] = 0; pLife[i] = 2.0 + rnd() * 4.0;       // seconds before respawn (staggered so trails don't blink together)
}
for (let i = 0; i < N; i++) spawn(i);

// a soft ROUND particle sprite (a radial alpha falloff): square GL points are the amateur tell, and a
// round sprite is the single biggest "showcase" upgrade. The sprite is white, so the per-vertex
// temperature colour survives (final rgb = vColor·white); its radial alpha rounds off the dot.
function particleSprite() {
  const s = 64, cvs = document.createElement("canvas"); cvs.width = cvs.height = s;
  const g = cvs.getContext("2d");
  const grd = g.createRadialGradient(s / 2, s / 2, 0, s / 2, s / 2, s / 2);
  grd.addColorStop(0.0, "rgba(255,255,255,1)");
  grd.addColorStop(0.45, "rgba(255,255,255,0.85)");
  grd.addColorStop(1.0, "rgba(255,255,255,0)");
  g.fillStyle = grd; g.fillRect(0, 0, s, s);
  return new T.CanvasTexture(cvs);
}
const geo = new T.BufferGeometry();
geo.setAttribute("position", new T.BufferAttribute(pos, 3));
geo.setAttribute("color", new T.BufferAttribute(col, 4));            // RGBA: the 4th channel is the spawn/death fade
const pmat = new T.PointsMaterial({ size: D.particle_size, map: particleSprite(), vertexColors: true,
                                    transparent: true, opacity: D.particle_opacity, depthWrite: false });
const points = new T.Points(geo, pmat);
scene.add(points);

function step(dt) {
  for (let i = 0; i < N; i++) {
    const la = pLat[i], lo = pLon[i];
    const uu = sample(U, la, lo), vv = sample(V, la, lo);
    const cosp = Math.max(0.05, Math.cos(la * RAD));
    pLon[i] += DEG * (uu / (A * cosp)) * ACCEL * dt;   // dλ/dt = u/(a cosφ)
    pLat[i] += DEG * (vv / A) * ACCEL * dt;            // dφ/dt = v/a
    pAge[i] += dt;
    const out = pLat[i] < cov.lat_min || pLat[i] > cov.lat_max || pLon[i] < cov.lon_min || pLon[i] > cov.lon_max;
    if (pAge[i] > pLife[i] || out) spawn(i);
    const xyz = sph(pLat[i], pLon[i], 1.012);
    pos[i * 3] = xyz[0]; pos[i * 3 + 1] = xyz[1]; pos[i * 3 + 2] = xyz[2];
    let t = S ? (sample(S, pLat[i], pLon[i]) - smin) / ((smax - smin) || 1) : 0.5;
    t = Math.max(0, Math.min(1, t));
    const c = cmap(t);                                 // full-brightness colour, ALWAYS (never dimmed toward black)
    // the fade lives in ALPHA, not RGB: a fresh particle fades in from transparent (not from a dark dot),
    // and a dying one fades fully out — both ends clean against whatever sits behind them.
    const alpha = Math.min(1, pAge[i] / 0.3) * Math.min(1, (pLife[i] - pAge[i]) / 0.5);
    col[i * 4] = c[0]; col[i * 4 + 1] = c[1]; col[i * 4 + 2] = c[2]; col[i * 4 + 3] = alpha;
  }
  geo.attributes.position.needsUpdate = true;
  geo.attributes.color.needsUpdate = true;
}

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
});

// --- live appearance controls: size + opacity are live-mutable on the material (no re-render of the data).
// The open-ended rest — colour ramps, particle density, trail length, shape menus — is a deliberately
// deferred seam: speculative, with no second consumer yet, so it stays named-not-built. --- #
const sizeR = document.getElementById("sizeRange"), opacR = document.getElementById("opacityRange");
if (sizeR) sizeR.addEventListener("input", () => { pmat.size = parseFloat(sizeR.value); });
if (opacR) opacR.addEventListener("input", () => { pmat.opacity = parseFloat(opacR.value); });

let last = performance.now();
function loop(now) {
  let dt = (now - last) / 1000; last = now; if (dt > 0.1) dt = 0.1;
  step(dt); renderer.render(scene, camera); requestAnimationFrame(loop);
}
requestAnimationFrame(loop);
"""


def flow_globe_html(field: FlowField, *, title: str = "planet-sim — eddy flow-globe (showcase)",
                    subtitle: str = "the emergent eddy life cycle as a particle flow on a real, rotatable "
                                    "planet (§9.5 Rung C — the showcase)",
                    n_particles: int = DEFAULT_N_PARTICLES,
                    particle_size: float = DEFAULT_PARTICLE_SIZE,
                    particle_opacity: float = DEFAULT_PARTICLE_OPACITY) -> str:
    """Render ``field`` as one deterministic, self-contained three.js HTML page (data + three.js inlined).

    The disclaimer (``field.honesty``) is written into a **visible** ``<div class="disclaimer">`` in the
    static DOM — not merely a JS comment — because under the honest-by-disclosure carve-out it *is* the
    entire license, and is the one thing machine-checked. Opens straight off ``file://`` (no network).

    ``particle_size`` / ``particle_opacity`` set the shipped *defaults*; they flow to both the material
    init and the live size/opacity sliders' initial positions (one source, no drift), so a viewer can
    fine-tune appearance in the browser without regenerating, and a notebook can ship a different default.
    """
    data_json = json.dumps(_build_data(field, n_particles, particle_size, particle_opacity),
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
