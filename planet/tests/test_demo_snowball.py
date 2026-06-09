"""Planet Phase-1 integration: the banked Snowball artifact end to end.

Runs the full continuation sweep the demo banks (present-day climate + the hysteresis loop) and
asserts the headline numbers against the cited climlab reference bands
(:mod:`projects.planet.climate_reference`). ``slow``-marked — the full-resolution sweep relaxes the
EBM at ~120 solar constants — so the fast lane deselects it (the structural physics is covered fast
by ``test_albedo``); it runs in the full gate. The figure smoke-test is ``importorskip``-gated on the
optional ``[viz]`` extra and is an *execution* check (ADR 0002), not a physics one.
"""
import numpy as np
import pytest

from projects.planet import demo_snowball as demo
from projects.planet.albedo import EBMParams, HysteresisLoop
from projects.planet.climate_reference import REFERENCE


@pytest.mark.slow
def test_demo_banked_numbers():
    r = demo.compute()
    # present-day: the temperate finite-cap branch
    lo, hi = REFERENCE.present_ice_line_band
    assert lo < r.present_ice_line < hi
    glo, ghi = REFERENCE.present_global_mean_band
    assert glo < r.present.global_mean_T < ghi
    assert r.present.net_toa == pytest.approx(0.0, abs=0.5)             # feedback state: albedo-step limited
    # the hysteresis loop: positive width, freeze within the cited dimming band, a frozen Snowball
    assert r.melt_S0 > r.freeze_S0 and r.hysteresis_width > 50.0
    dlo, dhi = REFERENCE.snowball_dimming_pct_band
    assert dlo < r.freeze_dimming_pct < dhi
    assert r.snowball_Tbar < REFERENCE.snowball_global_mean_max_C


@pytest.mark.slow
def test_demo_figure_renders():
    # ADR 0002: the figure is an execution smoke-test, not a physics check (the triads validate the
    # numbers). Gated on the optional [viz] extra; a headless/clean checkout skips rather than errors.
    pytest.importorskip("matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    from projects.planet.albedo import present_day_climate, snowball_hysteresis
    from projects.planet.plots import snowball_figure

    params = EBMParams(n_cells=60)
    present = present_day_climate(params, n_tau=0.05)
    loop = snowball_hysteresis(params=params, n_steps=8, n_tau=0.15, S0_min=1100.0, S0_max=1900.0)
    assert isinstance(loop, HysteresisLoop)
    fig = snowball_figure(loop, present)
    assert len(fig.axes) >= 3
    import matplotlib.pyplot as plt
    plt.close(fig)
