"""Planet Rung-5A integration: the banked rain-shadow artifact end to end.

Runs the demo the figure banks (a meridional ridge on a regional patch, the Smith & Barstad windward
rain + the lee moisture depletion) and asserts the headline story: the windward slope wettens far above
the zonal-mean baseline, the lee falls *below* it (a real rain-shadow desert), and the biome map is
re-classified over a large fraction of the patch. Uses ``use_jet=False`` (the S&B reference westerly) so
the check is fast and needs no shallow-water spin-up — the wind is essentially the same ~15 m/s the jet
gives at this latitude. The structural physics is covered by ``test_orographic`` / ``test_orographic_scene``
/ ``test_orographic_depletion``; the figure smoke-test is ``importorskip``-gated on ``[viz]`` and is an
*execution* check (ADR 0002), not a physics one.
"""
import numpy as np

from planet import demo_orographic as demo


def test_demo_rain_shadow_story():
    scene = demo.compute(use_jet=False)
    base_mean = float(scene.baseline_precip_cm.mean())
    crest_row = int(np.argmax(scene.elevation_m.max(axis=1)))
    tot_row = scene.precip_cm[crest_row]
    # windward slope wettens well above the zonal-mean baseline
    assert tot_row.max() > 2.0 * base_mean
    # the lee falls BELOW the baseline: a real rain-shadow desert (the 5A.3 depletion, not enhancement-only)
    assert tot_row.min() < base_mean
    assert scene.lee_desert_fraction > 0.0
    # the payoff: the mountain re-classifies a large fraction of the biome map
    assert scene.biome_changed_fraction > 0.2


def test_demo_figure_renders(tmp_path, monkeypatch):
    import pytest
    pytest.importorskip("matplotlib")
    # Render into a temp dir so the fast (use_jet=False) render never clobbers the committed figure,
    # which is the honest emergent-jet version banked by `python -m planet.demo_orographic`.
    monkeypatch.setattr(demo, "DOCS_FIGURE", tmp_path / "planet-orographic.png")
    monkeypatch.setattr(demo, "OUTPUT_FIGURE", tmp_path / "planet-orographic-out.png")
    scene = demo.compute(use_jet=False)
    saved = demo.save_figure(scene)
    assert saved.exists() and saved.stat().st_size > 0
