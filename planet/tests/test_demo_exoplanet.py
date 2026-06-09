"""Planet §9.1 integration: the banked exoplanet-knob artifact end to end.

Runs the demo the figure banks (Sun-vs-M-dwarf Snowball loops + size-scaled T(φ) profiles) and asserts
the two headline stories: a redder star narrows the hysteresis loop and lowers the freeze threshold, and
a bigger planet sharpens the equator-to-pole gradient. ``slow``-marked (it traces two continuation
sweeps); the structural physics is covered fast by ``test_exoplanet``. The figure smoke-test is
``importorskip``-gated on ``[viz]`` and is an *execution* check (ADR 0002), not a physics one.
"""
import numpy as np
import pytest

from projects.planet import demo_exoplanet as demo


@pytest.mark.slow
def test_demo_exoplanet_headlines():
    r = demo.compute(n_cells=72, n_steps=12, sweep_n_tau=0.12)
    # Knob 1 — a redder star is harder to snowball: the M-dwarf loop is much narrower and freezes lower.
    assert r.mdwarf_loop.hysteresis_width < 0.5 * r.sun_loop.hysteresis_width
    assert r.mdwarf_loop.freeze_S0 < r.sun_loop.freeze_S0
    assert r.stellar_ai["M5V"] < r.stellar_ai["G2V (Sun)"]
    # Knob 2 — a bigger planet sharpens the gradient: the ice line marches equatorward with size, while
    # the global mean barely moves (the 0-D mean is size-invariant; only the feedback shifts it).
    icelines = [st.ice_line_lat for st in r.size_states]
    assert icelines[0] > icelines[1] > icelines[2]                  # sizes (0.5, 1.0, 2.0) → ice equatorward
    means = np.array([st.global_mean_T for st in r.size_states])
    assert means.max() - means.min() < 13.0                          # mean moves only via the feedback (loose)


@pytest.mark.slow
def test_demo_exoplanet_figure_renders():
    # ADR 0002: an execution smoke-test, not a physics check (the triad validates the numbers).
    pytest.importorskip("matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    from projects.planet.plots import exoplanet_figure

    r = demo.compute(n_cells=48, n_steps=8, sweep_n_tau=0.15)
    fig = exoplanet_figure(r.sun_loop, r.mdwarf_loop, r.mdwarf_label, r.stellar_ai,
                           r.size_states, r.sizes)
    assert len(fig.axes) >= 3                                        # stellar + albedo + size panels
    import matplotlib.pyplot as plt
    plt.close(fig)
