"""Planet Rung-5B.2 integration: the banked 2-D continentality-map artifact end to end.

Runs the demo the figure banks (the 2-D seasonal limit cycle on a coarse idealized-Earth mask) and asserts
the headline map story: at one latitude the continental interior swings far more than the coast, which
swings more than the open ocean (continentality resolved *within* a latitude), while the annual-mean map
stays zonally flat (the mask is invisible in the mean). ``slow``-marked (it marches a 2-D field to a limit
cycle); the structural physics is covered fast by ``test_seasonal_map``. The figure smoke-test is
``importorskip``-gated on ``[viz]`` and is an *execution* check (ADR 0002), not a physics one.
"""
import numpy as np
import pytest

from planet import demo_seasonal_map as demo


@pytest.mark.slow
def test_demo_continentality_map_story():
    r = demo.compute()
    m, c = r.model, r.climate
    i = r.band_index()
    interior, coast, ocean = r.sample_columns()
    rng = c.seasonal_range()
    # the sample points are what the story claims (two land, one ocean)
    assert c.land_mask[i, interior] and c.land_mask[i, coast] and not c.land_mask[i, ocean]
    # continentality resolved within a latitude: interior ≫ coast ≫ ocean
    assert rng[i, interior] > rng[i, coast] > rng[i, ocean]
    assert rng[i, interior] > 3.0 * rng[i, ocean]
    # the NMS headline: the annual-mean map is zonally flat (the mask is invisible in the mean)
    amean = c.annual_mean()
    assert float(np.max(amean.max(axis=1) - amean.min(axis=1))) < 0.5
    assert c.converged


@pytest.mark.slow
def test_demo_figure_renders():
    pytest.importorskip("matplotlib")
    r = demo.compute()
    saved = demo.save_figure(r)
    assert saved.exists() and saved.stat().st_size > 0
