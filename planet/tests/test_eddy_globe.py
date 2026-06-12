"""Planet Rung-B: the emergent eddy life cycle animated **on the globe** (§9.5).

The Rung-A peer (``test_demo_eddy_life.py``) covers the matplotlib two-panel; this covers the Plotly
globe. Three checks, by the ADR-0002 "reach not correctness" discipline:

* **the layering guard** (always-green): importing :mod:`planet.eddy_globe` must NOT pull plotly /
  matplotlib / ipywidgets — the headless-import discipline the whole webviz stack rides on;
* **the geometry pin** (always-green, the advisor's honesty gate): the band is a **bounded midlatitude
  longitude sector** — its true ~55° width in a *single* hemisphere — NOT a pole-to-pole / 360° global
  field (the overclaim a naive globe would commit);
* **the no-frames guard + the render smoke-test**: a clear ``ValueError`` without the frames
  side-channel, and the figure builds + the HTML banks (``importorskip``-gated on ``[webviz]``).

The frame-fidelity + diagnostic-purity physics lives in ``test_eddy_flux.py``; this is an *execution*
smoke-test, not a physics one (ADR 0002).
"""
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from planet import eddy_globe as eg

REPO_ROOT = Path(__file__).resolve().parents[2]


def test_importing_eddy_globe_stays_headless():
    # The discipline: eddy_globe is NumPy-only at import; plotly is imported lazily inside the figure
    # builder (matplotlib never), so it loads on a bare core install — exactly like planetmap. Checked
    # in a CLEAN subprocess (an in-process sys.modules check is fragile once other tests import plotly).
    code = (
        "import sys, planet.eddy_globe\n"
        "print(','.join(m for m in ('plotly', 'ipywidgets', 'matplotlib') if m in sys.modules))\n"
    )
    out = subprocess.run([sys.executable, "-c", code], cwd=str(REPO_ROOT),
                         capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip() == "", f"a heavy dep was pulled at import time: {out.stdout.strip()}"


def _synthetic_frames(phi_lo=30.0, phi_hi=50.0, ny=24, nx=48, lx_km=3000.0):
    """A minimal stand-in for EddyFrames — just the (x, y, phi) geometry the band-builder consumes.

    ``y`` is the linear β-plane metric ``y = a·rad(φ − φ_ref)`` (so the builder recovers Earth's radius
    from it); ``x`` is a Cartesian periodic span of physical width ``lx_km``. No sim is run."""
    a = 6.371e6
    phi = np.linspace(phi_lo, phi_hi, ny)
    y = a * np.radians(phi - phi.mean())
    x = np.linspace(0.0, lx_km * 1e3, nx)
    return SimpleNamespace(x=x, y=y, phi=phi)


def test_band_is_a_bounded_midlat_sector_not_a_global_field():
    # The advisor's honesty gate: the channel maps to ONE band at its TRUE longitude width — never
    # pole-to-pole, never a 360° wrap, never mirrored into the other hemisphere.
    fr = _synthetic_frames(phi_lo=30.0, phi_hi=50.0, lx_km=3000.0)
    (x, y, z), lat2d, lon2d, lon1d = eg._band_geometry(fr)

    # latitude is the channel's OWN band (single hemisphere), not the whole planet.
    assert np.isclose(lat2d.min(), 30.0) and np.isclose(lat2d.max(), 50.0)
    assert lat2d.min() > 0.0                                   # NH only — not mirrored to the south

    # longitude is a bounded sector centred on 0, NOT a 360° wrap. Its width is the spherical metric
    # Δlon = Lx/(a·cosφ_c): ~3000 km at 40° ≈ 35°, an order of magnitude short of a full wrap.
    a = 6.371e6
    expected = np.degrees(3.0e6 / (a * np.cos(np.radians(40.0))))
    lon_span = lon1d.max() - lon1d.min()
    assert lon_span < 90.0                                     # a band patch, NOT a planet wrap
    assert np.isclose(lon1d.mean(), 0.0, atol=1e-9)            # centred on the camera meridian
    assert np.isclose(lon_span, expected, rtol=0.02)           # the TRUE physical width, not stretched

    # the mesh sits just outside the unit sphere (a patch raised onto the planet, not buried in it).
    assert np.allclose(np.sqrt(x**2 + y**2 + z**2), eg.BAND_RAISE)


def test_eddy_globe_figure_requires_frames():
    """``eddy_globe_figure`` fails fast with a clear error (not an AttributeError, and *before* it even
    reaches for plotly) when the frames side-channel is absent."""
    with pytest.raises(ValueError):
        eg.eddy_globe_figure(SimpleNamespace(frames=None))


@pytest.mark.slow
def test_demo_eddy_globe_renders(tmp_path):
    # ADR 0002: an execution smoke-test, not a physics check (test_eddy_flux validates the numbers).
    pytest.importorskip("plotly")
    from planet import demo_eddy_globe as demo

    r = demo.compute(nx=40, ny=40, n_frames=8)
    assert r.eddy.frames is not None and r.eddy.frames.times.size == 8
    fr = r.eddy.frames

    fig = eg.eddy_globe_figure(r.eddy)
    kinds = [type(t).__name__ for t in fig.data]
    assert kinds.count("Surface") == 2                         # the bare planet + the eddy band
    assert kinds.count("Scatter") >= 4                         # throughput, net, cursor, two markers
    assert len(fig.frames) == 8                                # one frame per banked snapshot

    # the frame wiring is the load-bearing part of an animation — pin it structurally (a figure that
    # builds but whose frames are empty / mis-indexed renders a static globe). Each frame must (a) carry
    # the matching θ snapshot as the band's surfacecolor, distinct frame-to-frame, with the fixed colour
    # range re-stated, and (b) target the band + cursor + the two markers by index (not the base sphere).
    band_i = 1
    for k, frame in enumerate(fig.frames):
        assert np.allclose(np.asarray(frame.data[0].surfacecolor), fr.theta[k])
        assert frame.traces[0] == band_i and len(frame.traces) == 4
        assert frame.data[0].cmin == fr.theta.min() and frame.data[0].cmax == fr.theta.max()
    # the band actually changes across frames (a constant surfacecolor would be a dead animation).
    assert not np.allclose(np.asarray(fig.frames[0].data[0].surfacecolor),
                           np.asarray(fig.frames[-1].data[0].surfacecolor))

    # the two honesty edges are ON the figure, not just in the demo's prose.
    assert "reversible" in fig.layout.title.text              # the ~90%-reversible edge
    assert "band" in fig.layout.title.text.lower()            # the band-not-globe edge
    # the band lat extent tracks the REAL channel (a midlat band), not pole-to-pole.
    (_x, _y, _z), lat2d, _lon2d, lon1d = eg._band_geometry(fr)
    assert np.isclose(lat2d.min(), fr.phi.min()) and np.isclose(lat2d.max(), fr.phi.max())
    assert (lon1d.max() - lon1d.min()) < 90.0                  # the bounded sector, not a 360° wrap

    out = tmp_path / "eddy-globe.html"
    eg.save_eddy_globe_html(r.eddy, out)
    assert out.exists() and out.stat().st_size > 0
