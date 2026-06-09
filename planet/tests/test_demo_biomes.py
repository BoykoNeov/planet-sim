"""Planet Phase-2 integration: the banked biome-map artifact end to end.

Runs the demo the figure banks (present-day climate → precip → biomes, plus the warmed-planet
comparison) and asserts the headline story: the present-day band ordering, the partition tiling the
planet, and the **poleward migration under warming** (the tropics expand, the ice/tundra retreat).
``slow``-marked — it relaxes the EBM at the demo's fine ``n_tau`` — so the fast lane deselects it (the
structural physics is covered fast by ``test_biomes`` / ``test_precip``); it runs in the full gate.
The figure smoke-test is ``importorskip``-gated on the optional ``[viz]`` extra and is an *execution*
check (ADR 0002), not a physics one.
"""
import pytest

from projects.planet import demo_biomes as demo
from projects.planet.biomes import Biome


@pytest.mark.slow
def test_demo_banked_biome_story():
    r = demo.compute()
    targets = (0, 15, 30, 45, 60, 75, 90)
    bands = {t: b for t, (_lat, _T, _P, b) in zip(targets, r.bands(targets))}
    assert bands[0] == Biome.TROPICAL_RAIN_FOREST
    assert bands[60] == Biome.BOREAL_FOREST
    assert bands[90] == Biome.TUNDRA
    assert sum(r.area_fractions.values()) == pytest.approx(1.0)             # the map tiles the planet
    assert r.fraction(Biome.TROPICAL_RAIN_FOREST) > 0.0


@pytest.mark.slow
def test_demo_warming_migrates_biomes_poleward():
    r = demo.compute()
    warm = demo.warmed()
    assert warm.state.global_mean_T > r.state.global_mean_T                 # the knob warms the planet
    # tropics expand, ice/tundra retreat (the poleward migration)
    assert warm.fraction(Biome.TROPICAL_RAIN_FOREST) > r.fraction(Biome.TROPICAL_RAIN_FOREST)
    assert warm.fraction(Biome.TUNDRA) < r.fraction(Biome.TUNDRA)
    assert warm.state.ice_line_lat >= r.state.ice_line_lat                  # ice retreats poleward


@pytest.mark.slow
def test_demo_figure_renders():
    # ADR 0002: an execution smoke-test, not a physics check (the triads validate the numbers).
    pytest.importorskip("matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    from projects.planet.plots import biomes_figure

    r = demo.compute()
    fig = biomes_figure(r.state, r.precip_cm, r.codes)
    assert len(fig.axes) >= 4                                               # map + whittaker + profile(+twin) + legend
    import matplotlib.pyplot as plt
    plt.close(fig)
