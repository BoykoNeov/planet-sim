"""Planet rung-A integration: the banked animated-eddy-life artifact end to end (§9.5).

Runs the demo the GIF banks (the released eddy life cycle with the ``n_frames`` side-channel) and
renders the two-panel animation. ``slow``-marked (it spins the engine up and releases it); the
frame-fidelity + diagnostic-purity physics is covered in ``test_eddy_flux.py``. The render is an
*execution* smoke-test (ADR 0002), ``importorskip``-gated on ``[viz]`` — not a physics one. The
no-frames guard is a fast check (it never runs the sim).
"""
import pytest


def test_eddy_life_animation_requires_frames():
    """``eddy_life_animation`` fails fast with a clear error if the frames side-channel is absent
    (computed with ``n_frames=0``) — not an ``AttributeError`` deep in the update loop."""
    pytest.importorskip("matplotlib")
    from types import SimpleNamespace

    from planet.plots import eddy_life_animation

    with pytest.raises(ValueError):
        eddy_life_animation(SimpleNamespace(frames=None))


@pytest.mark.slow
def test_demo_eddy_life_renders(tmp_path):
    # ADR 0002: an execution smoke-test, not a physics check (test_eddy_flux validates the numbers).
    pytest.importorskip("matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.animation import PillowWriter
    import matplotlib.pyplot as plt

    from planet import demo_eddy_life as demo
    from planet.plots import eddy_life_animation

    r = demo.compute(nx=40, ny=40, n_frames=8)
    assert r.eddy.frames is not None
    assert r.eddy.frames.times.size == 8

    anim = eddy_life_animation(r.eddy)
    assert len(anim._fig.axes) >= 2                       # the two panels (+ the θ colorbar)

    out = tmp_path / "eddy-life.gif"
    anim.save(out, writer=PillowWriter(fps=4))
    assert out.exists() and out.stat().st_size > 0
    plt.close("all")
