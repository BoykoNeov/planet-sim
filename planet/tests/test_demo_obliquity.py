"""Planet §9.1 integration: the banked obliquity-knob artifact end to end.

Runs the demo the figure banks (the geometric s₂(ε) curve + relaxed T(φ) at a range of tilts) and
asserts the two headline stories: the mechanism curve hits its exact anchors and reverses sign, and the
relaxed climate flattens (ice cap retreats) as tilt rises. ``slow``-marked (it relaxes several climates);
the structural physics is covered fast by ``test_obliquity``. The figure smoke-test is
``importorskip``-gated on ``[viz]`` and is an *execution* check (ADR 0002), not a physics one.
"""
import numpy as np
import pytest

from projects.planet import demo_obliquity as demo
from projects.planet.ebm import S2_INSOLATION


@pytest.mark.slow
def test_demo_obliquity_headlines():
    r = demo.compute(n_cells=72, curve_n=28, climate_n_tau=0.04)
    # Mechanism — the s₂(ε) curve: exact −5/8 at no tilt, climlab ≈−0.48 at Earth, monotone rise, sign flip.
    assert r.s2_grid[0] == pytest.approx(-0.625, abs=1e-3)              # ε=0 → exactly −5/8
    assert r.s2_earth == pytest.approx(S2_INSOLATION, abs=0.005)        # Earth → the climlab cross-check
    assert np.all(np.diff(r.s2_grid) > 0.0)                            # rises monotonically toward zero
    assert r.s2_grid[-1] > 0.0                                          # reverses sign by 90° (poles warmer)
    # Consequence — the relaxed climate flattens with tilt: the ice cap retreats poleward (CLIMATE_TILTS
    # is (0, 10, 23.44, 40), increasing tilt), Earth lands near the ~70° benchmark.
    icelines = [st.ice_line_lat for st in r.climate_states]
    assert icelines[0] < icelines[1] < icelines[-1]                    # more tilt → ice line poleward
    earth_idx = r.tilts.index(min(r.tilts, key=lambda t: abs(t - 23.44)))
    assert 60.0 < r.climate_states[earth_idx].ice_line_lat < 80.0       # Earth's ~70° ice line


@pytest.mark.slow
def test_demo_obliquity_figure_renders():
    # ADR 0002: an execution smoke-test, not a physics check (the triad validates the numbers).
    pytest.importorskip("matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    from projects.planet.plots import obliquity_figure

    r = demo.compute(n_cells=48, curve_n=20, climate_n_tau=0.06)
    fig = obliquity_figure(r.eps_grid, r.s2_grid, r.s2_earth, r.climate_states, r.tilts)
    assert len(fig.axes) >= 2                                           # the s₂ curve + the climate panel
    import matplotlib.pyplot as plt
    plt.close(fig)
