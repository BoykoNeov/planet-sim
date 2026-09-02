"""Planet Rung-5B.3 integration: the banked seasonal-ice-map artifact (still + GIF) end to end.

Runs the demo the figure banks (the 2-D ice-albedo limit cycle on the idealized-Earth mask) and asserts the
headline map story: the continental interior freezes for a good part of the year where the ocean at the same
latitude stays open, land ice is seasonal, and the annual mean now sees the mask (the continent colder than the
ocean at its latitude). ``slow``-marked; the structural physics is covered by ``test_seasonal_ice_map``. The
render smoke-tests are ``importorskip``-gated on ``[viz]`` and write to ``tmp_path`` (ADR 0002).
"""
import numpy as np
import pytest

from planet import demo_seasonal_ice_map as demo
from planet.ebm import T_FREEZE


@pytest.fixture(scope="module")
def result():
    return demo.compute()


@pytest.mark.slow
def test_demo_seasonal_ice_map_story(result):
    c = result.climate
    i = result.band_index()
    interior, ocean = result.sample_columns()
    assert c.land_mask[i, interior] and not c.land_mask[i, ocean]
    frac = c.ice_fraction()
    assert frac[i, interior] > 0.2 and frac[i, ocean] < 0.05
    assert not (c.T[i, interior] < T_FREEZE).all()             # land ice melts every summer
    anom = c.zonal_anomaly()
    assert anom[i, interior] < anom[i, ocean] - 0.5             # the mask is visible in the annual mean
    assert c.converged


@pytest.mark.slow
def test_demo_renders_still_and_animation(result, tmp_path, monkeypatch):
    pytest.importorskip("matplotlib")
    monkeypatch.setattr(demo, "DOCS_FIGURE", tmp_path / "planet-seasonal-ice-map.png")
    monkeypatch.setattr(demo, "OUTPUT_FIGURE", tmp_path / "out" / "planet-seasonal-ice-map.png")
    monkeypatch.setattr(demo, "DOCS_ANIMATION", tmp_path / "planet-seasonal-ice-map.gif")
    monkeypatch.setattr(demo, "OUTPUT_ANIMATION", tmp_path / "out" / "planet-seasonal-ice-map.gif")
    assert demo.save_figure(result).stat().st_size > 0
    assert demo.save_animation(result, fps=8, dpi=50).stat().st_size > 0


@pytest.mark.slow
def test_demo_renders_globe(result, tmp_path, monkeypatch):
    pytest.importorskip("plotly")
    monkeypatch.setattr(demo, "DOCS_GLOBE", tmp_path / "planet-seasonal-ice-globe.html")
    saved = demo.save_globe(result)
    assert saved.exists() and saved.stat().st_size > 0
