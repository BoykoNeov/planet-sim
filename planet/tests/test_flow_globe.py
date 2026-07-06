"""Planet Rung-C: the emergent eddy life cycle as a **particle flow-globe** (§9.5, the showcase).

The Rung-A peer (``test_demo_eddy_life.py``) and the Rung-B peer (``test_eddy_globe.py``) cover the
matplotlib GIF and the Plotly scalar globe. This covers the three.js particle showcase. Under the
**honest-by-disclosure** carve-out (ADR 0002 status note, 2026-06-12) the discipline is *asymmetric*:

* **physics-fidelity relaxes** — no byte-golden, no numerical-transport proof (the figure was never in
  the correctness path; the science layer keeps the full validation triad);
* **documentation verification tightens** — the on-screen disclaimer **is the entire license**, so it is
  the one thing machine-checked: it is a *visible* DOM element carrying *both* honesty clauses.

Plus the always-green guards every viz module carries: the headless-import discipline, and the
band-vs-globe honesty edge carried *in the data* (the coverage extent). The browser play-through itself
is the one thing not headlessly self-verifiable (no WebGL here) — handed to the user to eyeball.
"""
import re
import shutil
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from planet import flow_globe as fg

REPO_ROOT = Path(__file__).resolve().parents[2]


def _synthetic_eddy(ny=12, nx=20, nf=5):
    """A minimal stand-in for an EddyFlux+EddyFrames — just the geometry/fields the builder consumes.

    ``y`` is the linear β-plane metric ``y = a·rad(φ − φ_ref)`` (so the radius is recovered from it); the
    fields are deterministic noise (no sim is run — this keeps the test in the fast lane)."""
    a = 6.371e6
    phi = np.linspace(30.0, 50.0, ny)
    y = a * np.radians(phi - phi.mean())
    x = np.linspace(0.0, 3.0e6, nx)
    times = np.linspace(0.0, 60.0, nf)
    rng = np.random.default_rng(0)
    frames = SimpleNamespace(times=times, u=rng.standard_normal((nf, ny, nx)),
                             v=rng.standard_normal((nf, ny, nx)),
                             theta=rng.standard_normal((nf, ny, nx)), phi=phi, x=x, y=y)
    return SimpleNamespace(frames=frames, saturation_period=30.0, irreversible_fraction=0.08)


def test_importing_flow_globe_stays_headless():
    # The discipline: flow_globe is NumPy-only at import — it builds an HTML *string*, pulling no plotly /
    # matplotlib / ipywidgets (the eddy_globe import is local to the builder). Checked in a CLEAN
    # subprocess (an in-process sys.modules check is fragile once other tests import plotly).
    code = (
        "import sys, planet.flow_globe\n"
        "print(','.join(m for m in ('plotly', 'ipywidgets', 'matplotlib') if m in sys.modules))\n"
    )
    out = subprocess.run([sys.executable, "-c", code], cwd=str(REPO_ROOT),
                         capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip() == "", f"a heavy dep was pulled at import time: {out.stdout.strip()}"


def test_flow_field_from_eddy_requires_frames():
    """``flow_field_from_eddy`` fails fast with a clear error when the frames side-channel is absent."""
    with pytest.raises(ValueError):
        fg.flow_field_from_eddy(SimpleNamespace(frames=None))


def test_coverage_is_a_bounded_midlat_band_not_a_global_field():
    # The honesty edge, carried IN THE DATA: the coverage is one midlatitude band at its true longitude
    # width — never global, never mirrored to the other hemisphere. Particles seed only within it.
    field = fg.flow_field_from_eddy(_synthetic_eddy())
    cov = field.coverage
    assert cov.is_global is False
    assert np.isclose(cov.lat_min, 30.0) and np.isclose(cov.lat_max, 50.0)
    assert cov.lat_min > 0.0                                      # NH only — not mirrored to the south
    # a bounded sector at the spherical metric width Δlon = Lx/(a·cosφ_c) ≈ 3000 km at 40° ≈ 35°.
    lon_span = cov.lon_max - cov.lon_min
    assert 0.0 < lon_span < 90.0                                 # a band patch, NOT a 360° planet wrap
    # the contract arrays line up: u/v/scalar are (ny, nx) on the (lat, lon) grid.
    ny, nx = field.lat.size, field.lon.size
    assert field.u.shape == (ny, nx) and field.v.shape == (ny, nx)
    assert field.scalar.shape == (ny, nx)


def test_honesty_string_carries_both_clauses():
    # The disclaimer must carry BOTH edges (the carve-out's minimum content): the one-band/illustrative
    # edge AND the ~90%-reversible / κ-residual edge. A band-confined showcase still implies currents, so
    # the reversibility clause is non-negotiable (it is NOT downgraded to honest-by-construction).
    h = fg.flow_field_from_eddy(_synthetic_eddy()).honesty
    assert "band" in h and "planet-wide flow" in h               # one-band / illustrative
    assert "reversible" in h and "sloshes" in h and "κ" in h     # ~90%-reversible / net is the κ residual


def test_artifact_is_self_contained_with_three_js_inlined():
    # The §6 deliverable + the file:// property: three.js is vendored INLINE (its @license banner travels
    # with the artifact), the data is inlined (no fetch), and there is NO external reference — the page
    # opens straight off disk with no network.
    html = fg.flow_globe_html(fg.flow_field_from_eddy(_synthetic_eddy()))
    assert "Three.js Authors" in html and "SPDX-License-Identifier: MIT" in html   # banner intact
    assert "window.FLOW_DATA" in html                            # data inlined, not fetched
    assert '<canvas id="globe">' in html
    assert "<script src" not in html.lower()                     # no external <script src=…>
    assert 'src="http' not in html and 'href="http' not in html  # nothing loaded over the network


def test_disclaimer_is_a_visible_dom_element_carrying_both_clauses():
    # THE machine-checked guarantee (the disclaimer IS the entire license under the carve-out): it lives
    # in a real, VISIBLE DOM element — not a code comment, not buried in the JS data — and carries both
    # honesty clauses, legible to a casual viewer.
    html = fg.flow_globe_html(fg.flow_field_from_eddy(_synthetic_eddy()))
    m = re.search(r'<div class="disclaimer"[^>]*>(.*?)</div>', html, re.S)
    assert m, 'no visible <div class="disclaimer"> in the artifact'
    body = m.group(1)
    assert "illustrative" in body.lower() and "band" in body     # one-band / illustrative edge, on screen
    assert "reversible" in body and "sloshes" in body            # ~90%-reversible edge, on screen
    # the one machine-checked guarantee: the disclaimer is never hidden from the viewer — guard the
    # three ways CSS can silently hide an element (a future edit could otherwise gut the entire license).
    disclaimer_rule = re.search(r"\.disclaimer\s*\{([^}]*)\}", fg._CSS)
    assert disclaimer_rule, "no .disclaimer CSS rule"
    style = disclaimer_rule.group(1).replace(" ", "")
    assert "display:none" not in style and "visibility:hidden" not in style and "opacity:0" not in style


def test_html_carries_both_gpu_pipeline_and_cpu_fallback():
    # The advection runs on the GPU (a float-texture ping-pong pass) by default, but we cannot run WebGL
    # in CI — so the original CPU step() loop ships as a RUNTIME fallback (a GPU failure must degrade to a
    # working globe, never a blank one, and say why in the console). The fallback IS the verification
    # substitute, so pin that BOTH paths survive in the artifact — a future edit can't silently gut it.
    html = fg.flow_globe_html(fg.flow_field_from_eddy(_synthetic_eddy()))
    # the GPU pipeline: the ping-pong state advection (shaders + float render targets + the feature gate).
    assert "buildGPU" in html and "WebGLRenderTarget" in html
    assert "gl_FragColor = vec4(lon, lat, age, life)" in html         # the update-pass state shader (UPDATE_FS)
    assert "EXT_color_buffer_float" in html                           # the float-render feature gate
    # the CPU fallback: the JS step loop is built and wired, with a console reason when the GPU path is out.
    assert "buildCPU" in html and "CPU advection fallback active" in html


def test_emitted_app_js_parses(tmp_path):
    # The in-browser app is one big inlined <script>. A single JS syntax slip kills the ENTIRE block →
    # a blank canvas: no planet, no particles (the WebGL setup never runs). This bit once — a stray
    # backtick inside a GLSL template literal closed the literal early. We cannot run WebGL here, but we
    # CAN parse the JS: `node --check` catches exactly this class of error headlessly, closing the gap
    # between "renders an HTML string" and "that string's script is valid". Skips where node is absent
    # (the browser play-through stays the backstop), so it never blocks a node-less CI.
    node = shutil.which("node")
    if node is None:
        pytest.skip("node not available — JS syntax check skipped")
    html = fg.flow_globe_html(fg.flow_field_from_eddy(_synthetic_eddy()))
    app = re.findall(r"<script>(.*?)</script>", html, re.S)[-1]   # the app block is emitted last
    js = tmp_path / "app.js"
    js.write_text("var window = {};\n" + app, encoding="utf-8")
    r = subprocess.run([node, "--check", str(js)], capture_output=True, text=True)
    assert r.returncode == 0, f"emitted app JS failed to parse:\n{r.stderr}"


def _masked_field() -> fg.FlowField:
    """A tiny FlowField carrying an O1 validity mask (True = data) — a land strip inside the box."""
    lat = np.linspace(-40.0, 40.0, 9)
    lon = np.linspace(-60.0, 60.0, 13)
    ny, nx = lat.size, lon.size
    mask = np.ones((ny, nx), dtype=bool)
    mask[:, (lon >= 0.0) & (lon <= 30.0)] = False
    cov = fg.Coverage(-40.0, 40.0, -60.0, 60.0, is_global=False)
    return fg.FlowField(lat=lat, lon=lon, u=np.ones((ny, nx)), v=np.zeros((ny, nx)), coverage=cov,
                        honesty="a masked probe field — illustrative test data, one band, "
                                "not a planet-wide flow; reversible-flux clause n/a", mask=mask)


def test_mask_packs_into_the_data_and_none_stays_null():
    # The O1 mask rides the payload as flat 0/1 (row-major, like u/v) — and a mask-less field ships
    # `mask: null`, keeping every pre-O1 producer's payload semantics unchanged (default-off).
    d = fg._build_data(_masked_field(), 100, 0.03, 0.9, 0.5)
    assert d["mask"] is not None and len(d["mask"]) == 9 * 13
    assert set(d["mask"]) == {0, 1}
    eddy = fg._build_data(fg.flow_field_from_eddy(_synthetic_eddy()), 100, 0.03, 0.9, 0.5)
    assert eddy["mask"] is None


def test_particles_reject_masked_cells_in_both_advection_paths():
    # The renderer hook (§9.6 O1): the mask rides the velocity texture's formerly-free 4th channel, so
    # BOTH paths reject land — the GPU shader recycles a particle that drifts onto (or respawns on) a
    # masked texel, and the CPU path rejection-samples spawns + recycles in step(). Structural pins, the
    # same style as the GPU/CPU-fallback test above (we cannot run WebGL here).
    html = fg.flow_globe_html(_masked_field())
    assert "toHalf(MASK ? MASK[k] : 1)" in html                       # the mask IS the 4th texture channel
    assert "texture2D(uVel, velUV(lon, lat)).w < 0.5" in html          # GPU: advection kill on a masked texel
    assert "if (rv.w < 0.5) age = life + 1.0;" in html                # GPU: respawn rejects a masked texel (O3c refactor)
    assert "function validAt(" in html                                # CPU: the shared validity sampler…
    assert "if (!validAt(pLat[i], pLon[i])) continue;" in html        # …rejection-sampled at spawn (O3c)…
    assert "!validAt(pLat[i], pLon[i]);     // drifted onto a masked (land) cell" in html  # …and in step()
    assert '"mask":[' in html                                         # the mask data actually shipped


def test_land_ocean_base_layer_is_mask_driven_and_alignment_is_by_construction():
    # §9.6 O3a: the base sphere gets a two-tone land/ocean skin from the O1 mask. The one thing we CAN pin
    # headlessly (no WebGL) is the honesty-by-construction invariant + the graceful degrade: the base
    # fragment shader inverts each surface point to (lat, lon) with the SAME mapping the particles use
    # (`atan(n.z, n.x)` = the inverse of `sph()`), samples the mask on the SAME 0.5 coastline rule, and a
    # no-mask / compile-fail field falls back to the plain solid sphere.
    html = fg.flow_globe_html(_masked_field())
    assert "function buildBase(" in html                        # the mask-driven base path is wired
    assert "m >= 0.5 ? uOcean : uLand" in html                  # two-tone, 0.5 = the coast (as validAt())
    assert "atan(n.z, n.x)" in html                             # lat/lon inverse of sph() → alignment by construction
    assert "function solidBase(" in html and "MeshPhongMaterial" in html   # the no-mask / compile-fail fallback
    # the base only lights up where there IS a mask: a mask-less field (the eddy band) ships mask:null, so
    # buildBase() takes the solidBase() branch at runtime (same static app, the DATA is what differs).
    assert fg._build_data(fg.flow_field_from_eddy(_synthetic_eddy()), 100, 0.03, 0.9, 0.5)["mask"] is None


def test_speed_colormap_is_an_opt_in_default_leaving_flowfield_untouched():
    # §9.6 O3c: a diverging RdBu_r bleaches a 0→max speed field, so the speed scalar gets a SEQUENTIAL ramp
    # — but only as an opt-in default via the flow_globe_html kwarg (FlowField stays untouched, the §9.3
    # win). RdBu_r stays the default (θ is a signed field). The choice rides the payload as a `sequential`
    # flag both cmaps branch on.
    field = _masked_field()
    assert '"sequential":true' in fg.flow_globe_html(field, colormap="speed")
    assert '"sequential":false' in fg.flow_globe_html(field)                 # default = diverging RdBu_r
    # both cmap implementations (JS + GLSL) carry the two ramps and branch on the flag.
    html = fg.flow_globe_html(field, colormap="speed")
    assert "if (SEQ)" in html                                                # JS cmap branch
    assert "vec3 cmap(float t, float seq)" in html                          # GLSL cmap takes the flag


def test_speed_weighted_seeding_lives_in_the_respawn_path_in_both_advection_paths():
    # §9.6 O3c: particles respawn with acceptance ∝ |u,v|/speed_max so fast western-boundary currents
    # dominate. It MUST live in the RESPAWN path (weighting only the initial seed relaxes back to uniform as
    # particles age out), floored so calm water keeps an ambient fill, and it composes with the O1 mask
    # reject via the same invisible-retry idiom. Structural pins (no WebGL here), both paths.
    html = fg.flow_globe_html(_masked_field())
    assert '"speed_max":' in html                                           # the normaliser rides the payload
    # GPU respawn: land rejected outright, THEN valid cells accepted ∝ speed (both inside the age>life branch).
    assert "length(rv.xy) / max(1e-6, uSpeedMax)" in html
    assert "uSeedFloor" in html                                             # the calm-water floor
    # CPU respawn: the same two criteria in spawn() — reject land, then accept ∝ speed.
    assert "if (rnd() < Math.max(SEED_FLOOR, spd / SPEEDMAX)) break;" in html


def test_density_knob_is_wired_in_both_paths():
    # §9.5 unlock (the O3c deliverable): the ocean producer is the "second consumer," so particle density
    # becomes a live control knob. GPU cuts by a per-particle rank uniform; CPU hides the tail in tick.
    html = fg.flow_globe_html(_masked_field())
    assert 'id="densityRange"' in html                                      # the slider is in the DOM
    assert "aSeq > uDensity" in html                                        # GPU: rank-cut discard
    assert "if (i >= nActive) alpha = 0;" in html                           # CPU: tail hidden
    assert "applyDensity" in html                                           # dispatched to the active path


def test_trails_are_an_opt_in_default_off_kwarg():
    # §9.6 O3b: motion trails are default-OFF — no WebGL CI + a blind hand-off means the ocean globe the
    # user eyeballs is the first thing to exercise them, and the already-shipped eddy artifact can't
    # silently regress. The choice rides the payload; trail_decay seeds the trail-length knob.
    field = _masked_field()
    assert '"trails":false' in fg.flow_globe_html(field)                    # default off
    assert '"trails":true' in fg.flow_globe_html(field, trails=True)        # opt-in
    assert '"trail_decay":0.96' in fg.flow_globe_html(field, trails=True)   # the knob's start rides along


def test_trail_pipeline_is_gpu_only_with_the_rotation_smear_fix_and_a_plain_fallback():
    # The trail architecture, pinned structurally (no WebGL here): an accumulate-and-fade feedback buffer
    # gated on the GPU advection path, the depth-only occluder prepass that keeps back-side particles out
    # of the buffer, the additive composite, and the load-bearing rotation-smear fix (decay=0 while
    # dragging). Any miss degrades to the plain single-pass render — the CPU fallback is never touched.
    html = fg.flow_globe_html(_masked_field(), trails=True)
    assert "function buildTrails(" in html                                  # the trail build
    assert "if (useGPU && TRAILS) {" in html                                # gated on the GPU path only
    assert "renderer.render(occScene, camera);          // globe depth only" in html   # the occluder prepass
    assert "const decay = drag ? 0.0 : trailDecay;" in html                 # THE rotation-smear fix
    assert "T.CustomBlending" in html and "blendSrc: T.OneFactor" in html   # additive One+One composite
    # the degrade: renderFrame defaults to the plain single-pass, and trails ride ON TOP of GPU advection
    # (never the CPU fallback) — so a trail miss can only drop back to the O3a/O3c globe, never blank it.
    assert "let renderFrame = () => renderer.render(scene, camera);" in html
    assert "CPU advection fallback active" in html                          # the fade-only fallback survives intact


def test_trail_length_knob_and_resize_realloc_are_wired():
    # The §9.5 trail-length knob (the second control the ocean producer unlocks) + the omission the advisor
    # flagged: the screen-sized trail targets must be reallocated on resize or they misalign with the globe.
    html = fg.flow_globe_html(_masked_field(), trails=True)
    assert 'id="trailRange"' in html                                        # the slider (present only when trails on)
    assert "applyTrail = (v) => { trailDecay = v; };" in html               # dispatches to the decay uniform
    assert "trailResize = () =>" in html and "if (trailResize) trailResize();" in html   # realloc on resize
    # the slider is omitted when trails are off (no dead control on the eddy globe).
    assert 'id="trailRange"' not in fg.flow_globe_html(_masked_field())


@pytest.mark.slow
def test_demo_eddy_particles_banks_the_artifact(tmp_path):
    # ADR 0002: an execution smoke-test, not a physics check (test_eddy_flux validates the numbers). Runs
    # the (short) eddy sim, so it is slow-marked.
    from planet import demo_eddy_particles as demo

    r = demo.compute(nx=40, ny=40, n_frames=8)
    assert r.eddy.frames is not None and r.eddy.frames.times.size == 8

    field = fg.flow_field_from_eddy(r.eddy)
    out = tmp_path / "eddy-particles.html"
    fg.save_flow_globe_html(field, out)
    assert out.exists() and out.stat().st_size > 0
    text = out.read_text(encoding="utf-8")
    assert "Three.js Authors" in text and 'class="disclaimer"' in text
