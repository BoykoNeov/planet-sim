"""Planet Rung-5B.1 integration: the banked seasonal-cycle / continentality artifact end to end.

Runs the demo the figure banks (the seasonal limit cycle, both solvers) and asserts the headline story:
continentality — the land tile swings far more than the ocean tile at the same latitude and the ocean
lags more — plus the two solvers agreeing (the anti-damping cross-check). ``slow``-marked (it marches the
EBM to a limit cycle), so the fast lane deselects it; the structural physics is covered fast by
``test_seasonal``. The figure smoke-test is ``importorskip``-gated on ``[viz]`` and is an *execution*
check (ADR 0002), not a physics one.
"""
import numpy as np
import pytest

from planet import demo_seasonal as demo


@pytest.mark.slow
def test_demo_continentality_story():
    r = demo.compute()
    m, c = r.model, r.climate
    i = r.band_index()
    ampL, ampO = c.amplitude("land")[i], c.amplitude("ocean")[i]
    lagL, lagO = m.phase_lag_days(c.T_land)[i], m.phase_lag_days(c.T_ocean)[i]
    assert ampL > 4.0 * ampO                                    # land swings several× the ocean
    assert lagO > lagL                                          # the ocean lags the sun more
    # the core insight: the annual means are identical (continentality is all amplitude)
    assert np.max(np.abs(c.annual_mean("land") - c.annual_mean("ocean"))) < 1e-9
    # the two solvers agree (backward-Euler transport is not damping the swing)
    gap = max(np.max(np.abs(r.marched.T_land - c.T_land)), np.max(np.abs(r.marched.T_ocean - c.T_ocean)))
    assert r.marched.converged and gap < 0.15


@pytest.mark.slow
def test_demo_figure_renders():
    pytest.importorskip("matplotlib")
    r = demo.compute()
    saved = demo.save_figure(r)
    assert saved.exists() and saved.stat().st_size > 0
