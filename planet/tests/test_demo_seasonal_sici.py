"""Planet Rung-5B.4 integration: the banked seasonal-SICI artifact end to end.

Runs :mod:`planet.demo_seasonal_sici` and asserts the headline story — the annual-mean parent has a fold,
the seasonal sweep walks through it one cell at a time with no hysteresis, and the bistability returns as
the mixed layer deepens. ``slow``-marked: it marches many limit cycles. The demo's *shipped* constants
sweep 720 and 1440 cells and take several minutes, so the fixture shrinks every axis to the cheapest
configuration that still carries the story — the resolution and fineness of the banked figure are the
demo's business, the *logic* is what is under test here (the physics is covered fast by
``test_seasonal_sici``). The figure smoke-test is ``importorskip``-gated on ``[viz]`` (ADR 0002).
"""
import numpy as np
import pytest

from planet import demo_seasonal_sici as demo


@pytest.fixture(scope="module")
def result():
    # Shrink the demo to a test-sized run: a coarse grid, a coarse S₀ step, two depths (one seasonal, one
    # near-seasonless), two resolutions. Module-scoped, so the constants are swapped and restored by hand
    # rather than with the (function-scoped) monkeypatch fixture.
    saved = {k: getattr(demo, k) for k in
             ("SWEEP_N_CELLS", "SWEEP_N_STEPS", "SWEEP_S0", "DEPTHS", "DEPTH_N_CELLS",
              "RES_N_CELLS", "RES_S0")}
    demo.SWEEP_N_CELLS = 180
    demo.SWEEP_N_STEPS = 90
    demo.SWEEP_S0 = np.arange(1368.0, 1382.01, 2.0)
    demo.DEPTHS = (50.0, 800.0)
    # 360 for the depth leg: the physics test carries the 720-cell version of this control, and here the
    # only claim is that the demo's depth axis wires up and changes sign — but keep it above 180 so the
    # surviving cap is still a resolved object rather than a lone cell.
    demo.DEPTH_N_CELLS = 360
    demo.RES_N_CELLS = (360, 720)
    demo.RES_S0 = np.arange(1374.0, 1380.01, 2.0)
    try:
        yield demo.compute()
    finally:
        for k, v in saved.items():
            setattr(demo, k, v)


@pytest.mark.slow
def test_demo_seasonal_sici_story(result):
    # The annual-mean parent must actually have the fold this rung is about — otherwise the comparison is
    # vacuous and every "no loop" below would be trivially true.
    assert result.parent.small_ice_cap_fold is not None
    assert result.theta_c > 3.0 and result.parent_loop_width > 3.0

    # The seasonal sweep: no hysteresis, and the cap grows one cell at a time (the un-interpolated read).
    assert result.loop.all_converged
    assert not result.loop.detected and result.loop.width == 0.0
    assert result.loop.down.max_cell_jump == 1 and result.loop.up.max_cell_jump == 1
    assert result.forbidden_caps.size >= 1                       # caps the annual mean cannot hold

    # The depth axis: seasonal Earth is one climate, the near-seasonless deep ocean is two — and the
    # seasonal swing that separates them falls with depth (the mechanism, not just the outcome).
    shallow, deep = result.seed_points
    assert not shallow.bistable and deep.bistable
    assert deep.warm_cap_deg == 0.0 and deep.survived_cap_deg > 0.0
    assert deep.warm_polar_amplitude_K < shallow.warm_polar_amplitude_K

    # The resolution check: the verdict does not move as the polar cell shrinks. Only the LOOP WIDTH is
    # asserted here — it is the resolution-robust signature. The cell-count claim is deliberately NOT made
    # at this scale: on a coarse grid a θ_c cap spans about one cell, so "grew by one cell" and "a fold
    # switched the whole cap on" are the same observation; `test_seasonal_sici` makes that claim at 720,
    # where the cap spans several. What is checked here is that the fold's own scale is tracked and grows
    # with the grid — the guard that keeps the comparison honest as resolution changes.
    assert all(w == 0.0 for w in result.res_loop_width)
    assert result.res_polar_cell[1] < result.res_polar_cell[0]
    assert result.res_fold_cells[1] > result.res_fold_cells[0]
    assert result.res_max_jump[-1] <= result.res_fold_cells[-1]


@pytest.mark.slow
def test_demo_figure_renders(result, tmp_path):
    pytest.importorskip("matplotlib")
    saved_docs, saved_out = demo.DOCS_FIGURE, demo.OUTPUT_FIGURE
    demo.DOCS_FIGURE = tmp_path / "planet-seasonal-sici.png"
    demo.OUTPUT_FIGURE = tmp_path / "out" / "planet-seasonal-sici.png"
    try:
        saved = demo.save_figure(result)
        assert saved.exists() and saved.stat().st_size > 0
    finally:
        demo.DOCS_FIGURE, demo.OUTPUT_FIGURE = saved_docs, saved_out
