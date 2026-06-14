"""Planet rung-4-completion integration: the banked tropical-amplification artifact end to end.

Runs the demo the figure banks (gray OLR wired per latitude at the climlab-matched loading, present +
warmed) and asserts the headline story: the local OLR slope is smallest at the warm equator, and the
warming is tropically amplified — the mirror of rung-2.5's polar amplification. Fast (a handful of Newton
solves, no live external solver). The figure smoke-test is ``importorskip``-gated on ``[viz]`` and is an
*execution* check (ADR 0002), not a physics one (the triad in ``test_radiative_ebm`` validates the numbers).
"""
import numpy as np

from planet import demo_radiative_ebm as demo


def test_demo_banked_tropical_amplification_story():
    r = demo.compute()
    ta = r.ta
    # the discriminator: the local OLR slope is smallest at the warm equator (rung-0 assumes a flat 2).
    assert int(np.argmin(ta.B_loc_present)) == 0
    assert 0.30 < r.wv_fraction < 0.40                              # the climlab-matched loading
    # the headline: tropical amplification, against a uniform-slope null that warms uniformly.
    assert ta.amp_gray < 0.9 and ta.amp_gray_band < 0.9
    assert abs(ta.amp_null - 1.0) < 1e-6
    # present-day: a Jensen warm shift over rung-0, contrast essentially unchanged.
    assert ta.gray_present.global_mean_T > r.rung0_present_T.mean() + 1.0
    # the mean warming is NOT pinned at ΔA/B (Jensen + WV feedback amplify it).
    assert ta.mean_delta_T_gray > ta.dA_over_B_tan


def test_demo_figure_renders():
    # ADR 0002: an execution smoke-test, not a physics check.
    import pytest
    pytest.importorskip("matplotlib")
    import matplotlib
    matplotlib.use("Agg")

    r = demo.compute()
    saved = demo.save_figure(r)
    assert saved.exists()
