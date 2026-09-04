"""Planet Rung-5A.4 integration: the banked alpine-biomes artifact end to end.

Runs the demo the figure banks — the same patch, the same Rung-5A.3 rainfall, classified with and
without the terrain temperature correction — and asserts the headline story *and* the negative:

* **the payoff** — a 2500 m crest at a temperate-forest latitude cools ~16 K and re-classifies to a
  cold-limited biome, and the cooling **alone** (same rain, uncooled control) moves a large fraction of
  the patch;
* **the negative that chose the default** — the emergent moist adiabat *confirms* the pinned 6.5 K/km
  constant at mid-latitudes rather than retiring it, and in the deep tropics it puts the freezing level
  ~45 % **above** the observed ~4.5–5 km band while the constant lands just below it. The assertion is
  the **ordering** (one close, one far), not "the constant is in the band" — see
  :mod:`planet.tests.test_elevation_temperature` for why that distinction is deliberate.

The structural physics is covered by ``test_elevation_temperature`` / ``test_orographic_scene``; the
figure smoke-test is ``importorskip``-gated on ``[viz]`` and is an *execution* check (ADR 0002), not a
physics one. ``use_jet=False`` keeps it fast — the wind is irrelevant to the cooling either way.
"""
import numpy as np

from planet import biomes, demo_alpine_biomes as demo, elevation_temperature as elev


def test_demo_alpine_story():
    scene, control, diagnostic = demo.compute(use_jet=False)
    mid = scene.temperature_C.shape[0] // 2
    crest = int(np.argmax(scene.elevation_m[mid, :]))

    # the crest cooled by ~Γ·z and crossed into a cold-limited biome the uncooled control does not have
    assert scene.elevation_cooling_K.max() > 10.0
    assert scene.biome_codes[mid, crest] in (biomes.Biome.TUNDRA, biomes.Biome.BOREAL_FOREST)
    assert control.biome_codes[mid, crest] not in (biomes.Biome.TUNDRA, biomes.Biome.BOREAL_FOREST)
    # the cooling ALONE re-classifies a large fraction (the control isolates it from the rain shadow)
    assert scene.alpine_fraction > 0.2
    assert scene.biome_changed_fraction > control.biome_changed_fraction


def test_demo_lapse_rate_verdict_is_the_negative():
    _, _, d = demo.compute(use_jet=False)
    lo, hi = elev.OBSERVED_TROPICAL_FREEZING_LEVEL_M
    mid = int(np.argmin(np.abs(d.latitude_deg - 47.0)))
    # at mid-latitudes the emergent rate reproduces the pinned constant instead of retiring it
    assert abs(d.gamma_moist[mid] - d.gamma_constant[mid]) * 1e3 < 0.5
    # in the deep tropics the constant lands in the observed band and the moist adiabat overshoots it
    assert lo * 0.9 < d.freezing_constant_m[0] < hi
    assert d.freezing_moist_m[0] > hi * 1.2


def test_demo_figure_renders(tmp_path, monkeypatch):
    import pytest
    pytest.importorskip("matplotlib")
    monkeypatch.setattr(demo, "DOCS_FIGURE", tmp_path / "planet-alpine-biomes.png")
    monkeypatch.setattr(demo, "OUTPUT_FIGURE", tmp_path / "planet-alpine-biomes-out.png")
    scene, control, diagnostic = demo.compute(use_jet=False)
    saved = demo.save_figure(scene, control, diagnostic)
    assert saved.exists() and saved.stat().st_size > 0
