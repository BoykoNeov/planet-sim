"""Planet Phase-4 integration: the banked emergent-jet artifact end to end.

Runs the demo the figure banks (present-day EBM → forced shallow-water → emergent jet) and asserts
the headline story: a geostrophically-balanced westerly jet at midlatitudes, mass machine-exact under
forcing, and the release test re-confirming the engine's invariants while the jet persists.
``slow``-marked (it spins the engine up many inertial periods); the structural physics is covered fast
by ``test_coupler`` / the engine seal. The figure smoke-test is ``importorskip``-gated on ``[viz]`` and
is an *execution* check (ADR 0002), not a physics one.
"""
import numpy as np
import pytest

from planet import demo_coupler as demo


@pytest.mark.slow
def test_demo_banked_jet_story():
    r = demo.compute(nx=48, ny=48)
    j = r.jet
    assert 30.0 <= j.jet_lat <= 50.0                          # midlatitude jet (benchmark)
    assert j.jet_speed > 5.0                                  # tens of m/s
    assert j.core_balance_residual < 0.05                     # geostrophically balanced (anchor)
    assert np.abs(j.mass).max() < 1e-10                       # mass machine-exact under forcing
    # release re-confirms the engine's invariants AND the jet persists (the reframed conservation leg)
    assert np.abs(j.energy_release).max() < 1e-6
    assert j.u_profile_release.max() == pytest.approx(j.jet_speed, rel=0.05)


@pytest.mark.slow
def test_demo_figure_renders():
    # ADR 0002: an execution smoke-test, not a physics check (the triad validates the numbers).
    pytest.importorskip("matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    from planet.plots import coupler_figure

    r = demo.compute(nx=48, ny=48)
    fig = coupler_figure(r.jet, r.state)
    assert len(fig.axes) >= 4                                 # jet + chain(+twin) + field(+cbar) + conserve
    import matplotlib.pyplot as plt
    plt.close(fig)
