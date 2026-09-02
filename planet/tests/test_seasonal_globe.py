"""Render smoke-tests for the seasonal-ice month-slider globe (:mod:`planet.seasonal_globe`, rung 5B.3).

Structural checks on a tiny, fast 2-D ice march (ADR 0002 — execution, not physics): twelve frames, a
slider with twelve steps, pointer spikes off (the standing preference), frozen cells pinned to the ice
sentinel, and a standalone HTML written. ``importorskip``-gated on Plotly (the ``[webviz]`` extra).
"""
import numpy as np
import pytest

from planet import seasonal as sea
from planet import seasonal_map as sm
from planet.seasonal import SeasonalEBM


@pytest.fixture(scope="module")
def tiny_climate():
    x = SeasonalEBM(n_cells=16).x
    lon = (np.arange(12) + 0.5) * 2.0 * np.pi / 12
    m = sm.SeasonalMapEBM(land_mask=sm.box_mask(x, lon, (10.0, 70.0), (30.0, 150.0)),
                          n_cells=16, n_lon=12, n_steps=48)
    return m.march(coalbedo_fn=sea.ice_coalbedo, T_init=10.0, tol=1e-3, max_years=30)


def test_globe_structure(tiny_climate):
    pytest.importorskip("plotly")
    from planet import seasonal_globe as sg
    fig = sg.seasonal_ice_globe(tiny_climate)
    assert len(fig.frames) == 12
    assert len(fig.layout.sliders[0].steps) == 12
    assert fig.layout.scene.xaxis.showspikes is False
    surf = fig.data[0]
    frozen = tiny_climate.frozen(sg.month_steps(tiny_climate.days.size)[0])
    assert frozen.any()
    colour = np.asarray(surf.surfacecolor)[1:-1, :-1]       # strip the pole caps + the seam wrap column
    assert np.asarray(surf.surfacecolor).shape == (tiny_climate.T.shape[0] + 2, tiny_climate.T.shape[1] + 1)
    assert np.all(colour[frozen] == surf.cmax)               # frozen cells sit at the ice sentinel
    assert np.all(colour[~frozen] < surf.cmax)


def test_globe_writes_html(tiny_climate, tmp_path):
    pytest.importorskip("plotly")
    from planet import seasonal_globe as sg
    out = sg.save_seasonal_ice_globe(tiny_climate, tmp_path / "globe.html")
    assert out.exists() and out.stat().st_size > 1000
