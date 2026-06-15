"""Generate ``docs/interactive/index.html`` — the **no-install** browser what-if.

The user-facing front door for *experimentation*: drag three knobs (the Sun's brightness, the
greenhouse strength, and the axial tilt) and watch the planet's climate, its biome bands, and a
plain-language *what changed + why* explanation update **instantly** — no Jupyter, no Python, no
install. It opens straight off disk and serves from GitHub Pages the same way the existing globes do.

How it stays honest (the repo's whole character):

* **Every number is the real model.** A grid of :func:`planet.demo_biomes.compute` runs is
  precomputed here (:func:`compute_grid`); the page only *looks up* the nearest cell. That trades
  off-grid continuity for an instant, deterministic, shareable artifact — the same bargain the
  banked figures make. The one thing a lookup can't show — the Snowball's path-dependence
  (hysteresis) — is *named* in the explanation and demonstrated live in the notebook's §2.
* **The explanation is computed once, in Python** (:mod:`planet.explain`) and baked into the page,
  so the browser carries no climate rules of its own and can never drift from the notebook.
* **Self-contained & deterministic.** The data is *inlined* into the HTML (``file://`` blocks
  ``fetch``, so a separate JSON would not open off disk); the markup, CSS and JS are static; there
  is no build timestamp. So the shell is golden-testable like ``planet/site.py`` and the full
  artifact is reproducible from the model.

Run it (recomputes the grid — tens of seconds — then writes the page)::

    python -m planet interactive
"""
from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from planet.albedo import A_OLR, EBMParams, S0_EARTH
from planet.biomes import BIOME_COLORS, BIOME_NAMES, Biome
from planet.catalog import _REPO_ROOT
from planet.demo_biomes import compute
from planet.explain import Knobs, diagnose, explain
from planet.obliquity import OBLIQUITY_EARTH, OBLIQUITY_FAITHFUL_MAX, obliquity_params

APP_PATH = _REPO_ROOT / "docs" / "interactive" / "index.html"

# --- the grid (the three browser knobs) ------------------------------------------------------- #
# S0 spans dim→bright and lands an exact Earth detent on 1365 (1365 = 1235 + 13·10); the low end
# crosses the Snowball cliff (~1250) so dragging the Sun can freeze the planet. CO2 warming is a
# cut in the OLR offset A (demo_biomes' convention), 0→16 W/m² (present-day → ice-free hothouse).
# Obliquity (axial tilt) is the third axis: an untilted 0° world (s₂ = −5/8, the steepest insolation
# gradient, coldest poles) up to the single-P₂ knob's faithful cap (OBLIQUITY_FAITHFUL_MAX = 45°,
# beyond which the dropped s₄ grows — the obliquity.py "scope edge"). Unlike S0 it is smooth and
# cliff-free, so a coarse 9-value axis captures the whole trend (it triples the grid, not 24×). The
# exact OBLIQUITY_EARTH float is included so Earth's tilt recovers the model bit-for-bit (the
# obliquity factor is exactly 1 there → s₂ unchanged → the (1365, 0, 23.44°) cell is the baseline).
S0_VALUES = [1235.0 + 10.0 * i for i in range(24)]        # 1235 … 1465, 24 steps, includes 1365
CO2_VALUES = [float(i) for i in range(0, 17)]             # 0 … 16 W/m², 17 steps, includes 0
OBLIQUITY_VALUES = sorted(                                # 9 steps, 0 … 45°, includes Earth's 23.44°
    {0.0, 6.0, 12.0, 18.0, OBLIQUITY_EARTH, 30.0, 35.0, 40.0, OBLIQUITY_FAITHFUL_MAX})
_LAT_STRIDE = 6                                           # 180 model latitudes → 30 (equator→pole)


@dataclass(frozen=True)
class Cell:
    """One grid point: the diagnostics, the half-hemisphere profiles, and the baked explanation."""

    Tbar: float
    ice: float
    rainforest: float
    tundra: float
    desert: float
    biome: str          # one digit (0–8) per output latitude, equator→pole
    temp: list[float]   # °C per output latitude, equator→pole
    headline: str
    oneline: str
    paragraph: str


def _axis_default_index(values: list[float], target: float) -> int:
    """The slider's home index — the exact ``target`` if present, else the nearest value."""
    if target in values:
        return values.index(target)
    return min(range(len(values)), key=lambda i: abs(values[i] - target))


def compute_grid(s0_values: list[float] = S0_VALUES,
                 co2_values: list[float] = CO2_VALUES,
                 obliquity_values: list[float] = OBLIQUITY_VALUES) -> dict:
    """Run the validated model over the (S0 × CO2 × tilt) grid → a JSON-ready dict (the slow step).

    The cell list is flattened in this exact nesting order — ``s0`` outermost, ``co2``, then
    ``obliquity`` innermost — so the page decodes it as ``cells[(i·nCo2 + j)·nObl + k]``; the loop
    order and that JS index math must move together.
    """
    base_result = compute(EBMParams())
    base_diag = diagnose(base_result)
    lat_half = [round(float(v), 2) for v in base_result.state.latitude_deg()[::_LAT_STRIDE]]

    cells: list[dict] = []
    for s0 in s0_values:
        for co2 in co2_values:
            for obl in obliquity_values:
                params = obliquity_params(obl, EBMParams(S0=s0, A=A_OLR - co2))
                result = compute(params)
                diag = diagnose(result)
                ex = explain(Knobs(S0=s0, A=A_OLR - co2, obliquity_deg=obl), base_diag, diag)
                codes = result.codes[::_LAT_STRIDE]
                temp = result.state.T[::_LAT_STRIDE]
                cells.append(Cell(
                    Tbar=round(diag.global_mean_T, 2),
                    ice=round(diag.ice_line_lat, 1),
                    rainforest=round(diag.rainforest_pct, 1),
                    tundra=round(diag.tundra_pct, 1),
                    desert=round(diag.desert_pct, 1),
                    biome="".join(str(int(c)) for c in codes),
                    temp=[round(float(t), 1) for t in temp],
                    headline=ex.headline, oneline=ex.oneline, paragraph=ex.paragraph,
                ).__dict__)

    return {
        "axes": {
            "s0": {"label": "Sun — stellar flux S₀", "unit": "W/m²",
                   "values": [round(v, 0) for v in s0_values],
                   "default_index": _axis_default_index(s0_values, S0_EARTH)},
            "co2": {"label": "Greenhouse — added CO₂", "unit": "W/m²",
                    "values": [round(v, 0) for v in co2_values], "default_index": 0},
            "obl": {"label": "Tilt — axial obliquity", "unit": "°",
                    "values": [round(v, 0) for v in obliquity_values],
                    "default_index": _axis_default_index(obliquity_values, OBLIQUITY_EARTH)},
        },
        "lat": lat_half,
        "palette": {str(int(b)): BIOME_COLORS[b] for b in Biome},
        "names": {str(int(b)): BIOME_NAMES[b] for b in Biome},
        "baseline": {"Tbar": round(base_diag.global_mean_T, 2),
                     "ice": round(base_diag.ice_line_lat, 1)},
        "cells": cells,
    }


# --- the page (static shell + inlined data) --------------------------------------------------- #
_CSS = """\
:root { color-scheme: dark; }
* { box-sizing: border-box; }
body { margin: 0; font: 16px/1.55 system-ui, -apple-system, Segoe UI, Roboto, sans-serif;
       color: #e8ecf4; background: #0b1020; }
a { color: #8ab4ff; }
header { padding: 2.2rem 1.5rem 1rem; text-align: center;
         background: radial-gradient(1200px 400px at 50% -10%, #1b2647 0%, #0b1020 70%); }
header h1 { margin: 0; font-size: 2rem; letter-spacing: -.02em; }
header p { max-width: 44rem; margin: .5rem auto 0; color: #aeb7cc; }
main { max-width: 60rem; margin: 0 auto; padding: 1rem 1.2rem 4rem; }
.controls { display: grid; grid-template-columns: repeat(auto-fit, minmax(13rem, 1fr));
            align-items: start; gap: 1rem 1.6rem; background: #131a30; border: 1px solid #232c49;
            border-radius: 12px; padding: 1.1rem 1.3rem; margin-bottom: 1.2rem; }
/* Pin the label box height so every slider starts at the same y — the three knobs sit on one
   aligned row — and so a value string that re-wraps the label as you drag can't jog the slider. */
.knob label { display: block; min-height: 3.4rem; font-weight: 600; margin: 0 0 .35rem; }
.knob .val { color: #8ab4ff; font-variant-numeric: tabular-nums; }
.knob .hint { color: #8b95ad; font-size: .82rem; margin-top: .15rem; }
input[type=range] { -webkit-appearance: none; appearance: none; width: 100%; height: 8px;
   border-radius: 6px; background: #283355; outline: none; }
input[type=range]::-webkit-slider-thumb { -webkit-appearance: none; appearance: none;
   width: 22px; height: 22px; border-radius: 50%; background: #ffd166; border: 2px solid #0b1020;
   box-shadow: 0 0 0 3px #ffd16655; cursor: pointer; }
input[type=range]::-moz-range-thumb { width: 22px; height: 22px; border-radius: 50%;
   background: #ffd166; border: 2px solid #0b1020; box-shadow: 0 0 0 3px #ffd16655; cursor: pointer; }
.knob.greenhouse input[type=range]::-webkit-slider-thumb { background: #7ee0a0; box-shadow: 0 0 0 3px #7ee0a055; }
.knob.greenhouse input[type=range]::-moz-range-thumb { background: #7ee0a0; box-shadow: 0 0 0 3px #7ee0a055; }
.knob.tilt input[type=range]::-webkit-slider-thumb { background: #b69cff; box-shadow: 0 0 0 3px #b69cff55; }
.knob.tilt input[type=range]::-moz-range-thumb { background: #b69cff; box-shadow: 0 0 0 3px #b69cff55; }
.stage { display: flex; flex-wrap: wrap; gap: 1.2rem; }
/* The left column is the live readout (planet disk, temperature curve, stats, biome legend); the
   right panel is pure prose. Keeping the stats and legend OUT of the prose flow is what makes the
   layout steady — a longer one-liner, or expanding "Why", can only push the footer, never the
   numbers or the legend the eye is parked on. So the no-reflow goal is met structurally, not by
   reserving a fixed prose height (see .oneline). */
.vizcol { flex: 1 1 22rem; display: flex; flex-direction: column; gap: 1rem; }
.viz { background: #0c1226; border: 1px solid #232c49; border-radius: 12px; padding: 1rem;
       display: flex; gap: 1rem; justify-content: center; }
.panel { flex: 1 1 18rem; }
/* No reserved height. The one-liner grows a clause per moved knob (a Snowball adds a long one), so
   any fixed min-height is wrong both ways — too tall for the common single-knob case, too short for
   the worst three-knob/Snowball one. Instead the elements that must stay put (stats, legend) live in
   .vizcol, leaving the prose free to grow downward into the footer alone. */
.headline { font-size: 1.2rem; font-weight: 700; margin: .1rem 0 .5rem; }
.oneline { color: #d6def0; }
.more { margin-top: .6rem; }
.more summary { cursor: pointer; color: #8ab4ff; font-size: .9rem; }
.more p { color: #c2cbe0; margin: .5rem 0 0; }
.stats { display: grid; grid-template-columns: 1fr 1fr; gap: .3rem .8rem; margin: 0;
         font-size: .92rem; }
.stats div span { color: #8b95ad; }
.stats div b { color: #e8ecf4; font-variant-numeric: tabular-nums; }
.legend { margin-top: 0; display: grid; grid-template-columns: 1fr 1fr; gap: .25rem .8rem;
          font-size: .82rem; color: #b9c1d6; }
.legend i { display: inline-block; width: .8rem; height: .8rem; border-radius: 2px;
            margin-right: .4rem; vertical-align: -1px; }
#disk { cursor: crosshair; }
/* Hover read-out for the disk: which biome band is under the cursor. position:fixed + pointer
   events off, so following the mouse never reflows the layout (keeps the page steady, like #2). */
.tip { position: fixed; z-index: 10; pointer-events: none; display: none; white-space: nowrap;
       background: #0b1020; color: #e8ecf4; border: 1px solid #2d3a63; border-radius: 6px;
       padding: .25rem .5rem; font-size: .82rem; box-shadow: 0 2px 8px rgba(0,0,0,.45); }
.tip i { display: inline-block; width: .7rem; height: .7rem; border-radius: 2px;
         margin-right: .4rem; vertical-align: -1px; }
footer { color: #8b95ad; font-size: .88rem; margin-top: 2rem; border-top: 1px solid #232c49;
         padding-top: 1rem; }
"""

_BODY = """\
<header>
  <h1>Build a climate — turn a knob, watch a world</h1>
  <p>Brighten or dim the Sun, add greenhouse gas, or tilt the axis. The planet's temperature, its
  polar ice, and its bands of life respond instantly — and the panel tells you <em>what changed and
  why</em>. Every number is the real energy-balance model; this page just looks it up.</p>
</header>
<main>
  <div class="controls">
    <div class="knob sun">
      <label>☀ <span id="s0-label"></span>: <span class="val" id="s0-val"></span></label>
      <input type="range" id="s0" />
      <div class="hint">Dim it far enough and the planet snaps into a Snowball.</div>
    </div>
    <div class="knob greenhouse">
      <label>🏭 <span id="co2-label"></span>: <span class="val" id="co2-val"></span></label>
      <input type="range" id="co2" />
      <div class="hint">More greenhouse gas traps outgoing heat — a warmer world.</div>
    </div>
    <div class="knob tilt">
      <label>🌍 <span id="obl-label"></span>: <span class="val" id="obl-val"></span></label>
      <input type="range" id="obl" />
      <div class="hint">Earth tilts 23°. More tilt spreads the year's sunlight toward the poles.</div>
    </div>
  </div>

  <div class="stage">
    <div class="vizcol">
      <div class="viz">
        <canvas id="disk" width="300" height="300" aria-label="planet disk coloured by biome"></canvas>
        <canvas id="curve" width="320" height="300" aria-label="temperature by latitude"></canvas>
      </div>
      <div class="stats">
        <div><span>Global mean</span> <b id="st-tbar"></b></div>
        <div><span>Ice line</span> <b id="st-ice"></b></div>
        <div><span>Rain forest</span> <b id="st-rf"></b></div>
        <div><span>Tundra</span> <b id="st-tu"></b></div>
      </div>
      <div class="legend" id="legend"></div>
    </div>
    <!-- Pure prose, so it can grow freely: a longer one-liner or an expanded "Why" pushes only the
         footer below — the stats and legend live in .vizcol and never move (the no-reflow goal). -->
    <div class="panel">
      <div class="headline" id="headline"></div>
      <div class="oneline" id="oneline"></div>
      <details class="more"><summary>Why — the fuller mechanism</summary>
        <p id="paragraph"></p></details>
    </div>
  </div>

  <div class="tip" id="tip"></div>

  <footer>
    <p>A lookup over a grid of real <code>planet.demo_biomes.compute</code> runs — the same
    validated model behind the figures. Want continuous knobs, live re-runs, and the Snowball's
    two stable states (it won't re-melt at today's Sun — the climate is path-dependent)? That lives
    in the teaching notebook: <code>python&nbsp;-m&nbsp;planet&nbsp;notebook</code>.
    Back to the <a href="../index.html">gallery</a>.</p>
  </footer>
</main>
"""

# Plain JS (no f-string: braces are JS). DATA is concatenated in ahead of this block.
_APP_JS = r"""
const D = window.PLANET_DATA;
const nS0 = D.axes.s0.values.length, nCo2 = D.axes.co2.values.length, nObl = D.axes.obl.values.length;
// cells are flattened s0-outermost, co2, then obliquity innermost (compute_grid's loop order)
const cell = (i, j, k) => D.cells[(i * nCo2 + j) * nObl + k];
let current = null;   // the cell on screen — the disk-hover read-out names a band from it

const $ = id => document.getElementById(id);
const s0El = $("s0"), co2El = $("co2"), oblEl = $("obl");
s0El.min = 0; s0El.max = nS0 - 1; s0El.step = 1; s0El.value = D.axes.s0.default_index;
co2El.min = 0; co2El.max = nCo2 - 1; co2El.step = 1; co2El.value = D.axes.co2.default_index;
oblEl.min = 0; oblEl.max = nObl - 1; oblEl.step = 1; oblEl.value = D.axes.obl.default_index;
$("s0-label").textContent = D.axes.s0.label;
$("co2-label").textContent = D.axes.co2.label;
$("obl-label").textContent = D.axes.obl.label;

// legend (all biomes, in the equator→pole order they appear)
$("legend").innerHTML = Object.keys(D.names).map(k =>
  `<div><i style="background:${D.palette[k]}"></i>${D.names[k]}</div>`).join("");

// |latitude| (deg) → index of the band whose stored latitude is truly nearest. The model grid is
// equal-area (uniform in sin φ), so the stored latitudes are NOT evenly spaced — ~2° steps near the
// equator widening to ~6° near the pole, and the last one is only ~76°, not 90°. A linear index
// (a/latMax·(N−1)) reads a too-equatorward, too-warm band and drags the warm bands poleward —
// painting forest outside the ice-line ring (which is drawn at its true latitude). Search instead.
// Shared by the disk fill and the hover read-out so the two can never disagree.
function bandForLat(phi) {
  const lat = D.lat, N = lat.length, a = Math.abs(phi);
  let k = 0, best = Infinity;
  for (let m = 0; m < N; m++) { const d = Math.abs(lat[m] - a); if (d < best) { best = d; k = m; } }
  return k;
}

// --- the planet disk: horizontal biome bands by |latitude|, clipped to a circle, lit for 3D --- #
function drawDisk(c) {
  const cv = $("disk"), ctx = cv.getContext("2d"), W = cv.width, H = cv.height;
  const R = Math.min(W, H) / 2 - 6, cx = W / 2, cy = H / 2;
  ctx.clearRect(0, 0, W, H);
  ctx.save();
  ctx.beginPath(); ctx.arc(cx, cy, R, 0, 2 * Math.PI); ctx.clip();
  for (let y = -R; y <= R; y++) {
    const phi = (y / R) * 90;                 // |lat| band, symmetric top/bottom
    const code = c.biome[bandForLat(phi)];
    ctx.fillStyle = D.palette[code] || "#888";
    const half = Math.sqrt(Math.max(0, R * R - y * y));
    ctx.fillRect(cx - half, cy + y, 2 * half, 1.2);
  }
  ctx.restore();
  // ice-line rings (both hemispheres) when there is a finite cap
  if (c.ice < 89 && c.ice > 1) {
    ctx.strokeStyle = "#dff1ff"; ctx.lineWidth = 1; ctx.globalAlpha = .8;
    for (const s of [1, -1]) {
      const y = s * (c.ice / 90) * R, half = Math.sqrt(Math.max(0, R * R - y * y));
      ctx.beginPath(); ctx.moveTo(cx - half, cy + y); ctx.lineTo(cx + half, cy + y); ctx.stroke();
    }
    ctx.globalAlpha = 1;
  }
  // a soft top-left highlight so the disk reads as a sphere
  const g = ctx.createRadialGradient(cx - R * .35, cy - R * .35, R * .1, cx, cy, R * 1.05);
  g.addColorStop(0, "rgba(255,255,255,.18)"); g.addColorStop(.5, "rgba(255,255,255,0)");
  g.addColorStop(1, "rgba(0,0,0,.35)");
  ctx.fillStyle = g; ctx.beginPath(); ctx.arc(cx, cy, R, 0, 2 * Math.PI); ctx.fill();
}

// --- the temperature curve: T(latitude), mirrored to both hemispheres --- #
function drawCurve(c) {
  const cv = $("curve"), ctx = cv.getContext("2d"), W = cv.width, H = cv.height;
  const padL = 38, padB = 24, padT = 12, padR = 8;
  const x0 = padL, x1 = W - padR, y0 = H - padB, y1 = padT;
  ctx.clearRect(0, 0, W, H);
  // temperature range (fixed-ish, clamps to data)
  let tmin = 99, tmax = -99;
  for (const t of c.temp) { tmin = Math.min(tmin, t); tmax = Math.max(tmax, t); }
  tmin = Math.min(tmin, -10); tmax = Math.max(tmax, 30);
  tmin = Math.floor(tmin / 10) * 10; tmax = Math.ceil(tmax / 10) * 10;
  const X = phi => x0 + (phi + 90) / 180 * (x1 - x0);
  const Y = t => y0 + (t - tmin) / (tmax - tmin) * (y1 - y0);
  // axes + gridlines
  ctx.strokeStyle = "#26304f"; ctx.fillStyle = "#7b87a6"; ctx.font = "11px system-ui"; ctx.lineWidth = 1;
  for (let t = tmin; t <= tmax; t += 10) {
    ctx.beginPath(); ctx.moveTo(x0, Y(t)); ctx.lineTo(x1, Y(t)); ctx.stroke();
    ctx.fillText(t + "°", 6, Y(t) + 3);
  }
  ctx.fillText("S", x0 - 2, y0 + 16); ctx.fillText("eq", (x0 + x1) / 2 - 6, y0 + 16);
  ctx.fillText("N", x1 - 6, y0 + 16);
  // 0°C line
  if (tmin < 0 && tmax > 0) { ctx.strokeStyle = "#3a4a78"; ctx.beginPath();
    ctx.moveTo(x0, Y(0)); ctx.lineTo(x1, Y(0)); ctx.stroke(); }
  // the curve (mirror: stored equator→pole over [0,latMax])
  const lat = D.lat, N = lat.length, latMax = lat[N - 1];
  const pts = [];
  for (let k = N - 1; k >= 0; k--) pts.push([-lat[k], c.temp[k]]);   // south pole→equator
  for (let k = 0; k < N; k++) pts.push([lat[k], c.temp[k]]);          // equator→north pole
  ctx.strokeStyle = "#ffb454"; ctx.lineWidth = 2; ctx.beginPath();
  pts.forEach((p, i) => { const x = X(p[0]), y = Y(p[1]); i ? ctx.lineTo(x, y) : ctx.moveTo(x, y); });
  ctx.stroke();
  ctx.fillStyle = "#7b87a6"; ctx.fillText("temperature vs latitude", x0, y1 + 2);
}

function render() {
  const i = +s0El.value, j = +co2El.value, k = +oblEl.value, c = cell(i, j, k);
  current = c;
  $("s0-val").textContent = D.axes.s0.values[i] + " " + D.axes.s0.unit;
  $("co2-val").textContent = "+" + D.axes.co2.values[j] + " " + D.axes.co2.unit;
  $("obl-val").textContent = D.axes.obl.values[k] + D.axes.obl.unit;
  $("headline").textContent = c.headline;
  $("oneline").textContent = c.oneline;
  $("paragraph").textContent = c.paragraph;
  $("st-tbar").textContent = c.Tbar.toFixed(1) + " °C";
  $("st-ice").textContent = c.ice >= 89.5 ? "ice-free" : (c.ice <= 1 ? "frozen over" : c.ice.toFixed(0) + "°");
  $("st-rf").textContent = c.rainforest.toFixed(0) + "%";
  $("st-tu").textContent = c.tundra.toFixed(0) + "%";
  drawDisk(c); drawCurve(c);
}
s0El.addEventListener("input", render);
co2El.addEventListener("input", render);
oblEl.addEventListener("input", render);

// --- hover the disk → name the biome band under the cursor (same |lat|→band map as drawDisk) --- #
const diskEl = $("disk"), tipEl = $("tip");
diskEl.addEventListener("mousemove", e => {
  if (!current) return;
  const rect = diskEl.getBoundingClientRect();                  // CSS px → canvas px (in case scaled)
  const mx = (e.clientX - rect.left) * (diskEl.width / rect.width);
  const my = (e.clientY - rect.top) * (diskEl.height / rect.height);
  const W = diskEl.width, H = diskEl.height;
  const R = Math.min(W, H) / 2 - 6, cx = W / 2, cy = H / 2;
  const dx = mx - cx, dy = my - cy;
  if (dx * dx + dy * dy > R * R) { tipEl.style.display = "none"; return; }   // outside the globe
  const phi = (dy / R) * 90;                                     // drawDisk's mapping, exactly
  const a = Math.abs(phi);
  const code = current.biome[bandForLat(phi)];                   // same nearest-latitude band as the fill
  tipEl.innerHTML = `<i style="background:${D.palette[code]}"></i>${D.names[code]} · ${a.toFixed(0)}°`;
  tipEl.style.display = "block";
  tipEl.style.left = (e.clientX + 14) + "px";
  tipEl.style.top = (e.clientY + 14) + "px";
});
diskEl.addEventListener("mouseleave", () => { tipEl.style.display = "none"; });

render();
"""


def build_app_html(grid: dict) -> str:
    """Render the whole page as one deterministic, self-contained HTML string (data inlined)."""
    data = json.dumps(grid, separators=(",", ":"), ensure_ascii=False)
    return (
        "<!doctype html>\n<html lang=\"en\">\n<head>\n"
        "<meta charset=\"utf-8\">\n"
        "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
        "<title>planet-sim — build a climate</title>\n"
        f"<style>\n{_CSS}</style>\n</head>\n<body>\n"
        f"{_BODY}"
        f"<script>\nwindow.PLANET_DATA = {data};\n</script>\n"
        f"<script>{_APP_JS}</script>\n"
        "</body>\n</html>\n"
    )


def write_app(path: Path = APP_PATH, grid: dict | None = None) -> Path:
    """Recompute the grid (unless supplied) and write the page (LF, UTF-8)."""
    if grid is None:
        grid = compute_grid()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(build_app_html(grid), encoding="utf-8", newline="\n")
    return path


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    print("Running the model over the knob grid (a moment)…")
    saved = write_app()
    print(f"Interactive what-if written → {saved.relative_to(_REPO_ROOT)}")
    print("  open it in a browser (works straight off disk), or publish docs/ via GitHub Pages.")


if __name__ == "__main__":
    main()
