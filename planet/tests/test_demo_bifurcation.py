"""Planet Rung-0+ integration: the banked equilibrium-diagram artifact end to end.

Runs the demo the figure banks (the exact S-curve + the Phase-1 sweep + the D-sweep + the step-bias sweep)
and asserts the headline story: two folds, present-day inside and near the top of the finite-cap window,
the sweep jumping at the folds, θ_c growing with D, and the relaxation bias falling with the step.
``slow``-marked (it runs the Phase-1 sweep and a dozen relaxations); the physics is covered fast by
``test_bifurcation``. The figure smoke-test is ``importorskip``-gated on ``[viz]`` (ADR 0002).
"""
import numpy as np
import pytest

from planet import demo_bifurcation as demo
from planet.ebm import S0_EARTH


@pytest.fixture(scope="module")
def result():
    return demo.compute()


@pytest.mark.slow
def test_demo_equilibrium_diagram_story(result):
    c = result.curve
    lo, hi = c.snowball_fold, c.small_ice_cap_fold
    assert 25.0 < lo.latitude_deg < 40.0 and 75.0 < hi.latitude_deg < 85.0
    assert lo.S0 < S0_EARTH < hi.S0 and (hi.S0 - S0_EARTH) / S0_EARTH < 0.01
    step = float(np.abs(np.diff(result.loop.S0_up)).max())
    assert abs(result.loop.freeze_S0 - lo.S0) < step
    assert abs(result.loop.melt_S0 - c.snowball_threshold_S0) < step
    ok = ~np.isnan(result.theta_c)
    assert ok.sum() >= 5 and (~ok).any()                        # a branch that exists… then vanishes
    assert np.all(np.diff(result.theta_c[ok]) > 0.0)
    gaps = np.abs(result.relaxed_ice_line - result.present.latitude_deg)
    assert gaps[0] > gaps[-1] and gaps[-1] < 2.5


@pytest.mark.slow
def test_demo_figure_renders(result, tmp_path, monkeypatch):
    pytest.importorskip("matplotlib")
    monkeypatch.setattr(demo, "DOCS_FIGURE", tmp_path / "planet-bifurcation.png")
    monkeypatch.setattr(demo, "OUTPUT_FIGURE", tmp_path / "out" / "planet-bifurcation.png")
    saved = demo.save_figure(result)
    assert saved.exists() and saved.stat().st_size > 0
