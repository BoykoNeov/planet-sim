"""Planet Phase-3 integration: the banked shallow-water artifact end to end.

Runs the full-resolution demo the figure banks (geostrophic adjustment + Rossby wave) and asserts
the headline numbers — the adjusted state matches the analytic Helmholtz state, most of the bump
radiates away, and the Rossby wave goes westward. ``slow``-marked (the full 96² runs integrate many
inertial periods); the structural physics is covered fast by ``test_circulation`` and the engine
seal (``engines/fluid/tests/``). The figure smoke-test is ``importorskip``-gated on ``[viz]`` and is
an *execution* check (ADR 0002), not a physics one.
"""
import numpy as np
import pytest

from planet import demo_shallowwater as demo


@pytest.mark.slow
def test_demo_banked_numbers():
    r = demo.compute()
    # geostrophic adjustment: most of the bump radiates, the remnant is the analytic Helmholtz state
    assert r.drawdown_fraction > 0.7
    assert r.helmholtz_rel_error < 0.05
    assert np.abs(r.adjustment.mass).max() < 1e-10                  # mass machine-exact through the run
    assert r.adjustment.L_R == pytest.approx(960e3, rel=0.05)      # the cited ~1000 km scale
    # Rossby wave: westward, near the analytic phase speed
    assert r.rossby.c_measured < 0.0
    assert r.rossby.c_measured / r.rossby.c_analytic == pytest.approx(1.0, abs=0.1)


@pytest.mark.slow
def test_demo_figure_renders():
    # ADR 0002: an execution smoke-test, not a physics check. Gated on the optional [viz] extra.
    pytest.importorskip("matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    from planet.circulation import geostrophic_adjustment, rossby_wave
    from planet.plots import shallowwater_figure

    adj = geostrophic_adjustment(nx=48, ny=48, n_periods=6.0)
    ros = rossby_wave(nx=48, ny=48, frac_period=0.3)
    fig = shallowwater_figure(adj, ros)
    assert len(fig.axes) >= 4
    import matplotlib.pyplot as plt
    plt.close(fig)
