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
