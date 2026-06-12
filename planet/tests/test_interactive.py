"""Guards for the no-install browser what-if (:mod:`planet.interactive`).

Two tiers, like the rest of the repo:

* **fast** — the page *shell* is deterministic and self-contained (data inlined, no CDN), the grid
  has the right shape, and the committed page is well-formed. These build only a tiny 2×2 grid (a
  few model runs), so they run in the always-on lane.
* **slow** — the committed ``docs/interactive/index.html`` equals a fresh full-grid generation, so a
  change to the model or the page that isn't re-banked fails the gate. Marked ``slow`` because it
  reruns the whole knob grid (tens of seconds).
"""
from __future__ import annotations

import json
import re

import pytest

from planet import interactive

# A tiny grid whose (0,0) cell is the exact Earth detent (S0=1365, CO2=0) → the baseline message.
_SMALL = dict(s0_values=[1365.0, 1375.0], co2_values=[0.0, 2.0])


def test_compute_grid_shape():
    grid = interactive.compute_grid(**_SMALL)
    assert len(grid["cells"]) == 4                         # 2 × 2
    assert grid["axes"]["s0"]["values"] == [1365.0, 1375.0]
    assert grid["axes"]["co2"]["default_index"] == 0
    assert len(grid["lat"]) == len(grid["cells"][0]["temp"]) == len(grid["cells"][0]["biome"])
    assert set(grid["palette"]) == set(grid["names"])      # a colour for every named biome
    cell = grid["cells"][0]
    assert all(ch in "012345678" for ch in cell["biome"])  # biome string is single-digit codes


def test_earth_detent_is_the_baseline():
    grid = interactive.compute_grid(**_SMALL)
    earth = grid["cells"][0]                                # (s0=1365, co2=0)
    assert "baseline" in earth["headline"].lower()
    assert 13.0 < earth["Tbar"] < 16.0


def test_page_is_self_contained_and_deterministic():
    grid = interactive.compute_grid(**_SMALL)
    html = interactive.build_app_html(grid)
    assert html == interactive.build_app_html(grid)        # deterministic (no timestamp / set order)
    # data is inlined (file:// blocks fetch), and the page pulls in nothing external
    assert "window.PLANET_DATA = " in html
    assert "http://" not in html and "https://" not in html
    assert 'src="http' not in html
    assert '<canvas id="disk"' in html and '<canvas id="curve"' in html


def test_inlined_data_round_trips():
    grid = interactive.compute_grid(**_SMALL)
    html = interactive.build_app_html(grid)
    blob = re.search(r"window\.PLANET_DATA = (\{.*?\});\n", html, re.S).group(1)
    assert json.loads(blob)["axes"]["s0"]["values"] == grid["axes"]["s0"]["values"]


def test_committed_page_is_well_formed():
    """The shipped page must exist and be self-contained — cheap guard, no regeneration."""
    html = interactive.APP_PATH.read_text(encoding="utf-8")
    assert "window.PLANET_DATA = " in html
    assert "http://" not in html and "https://" not in html
    assert '<canvas id="disk"' in html


@pytest.mark.slow
def test_committed_page_is_up_to_date():
    """docs/interactive/index.html must equal a fresh full-grid build (re-run `python -m planet.interactive`)."""
    expected = interactive.build_app_html(interactive.compute_grid())
    actual = interactive.APP_PATH.read_text(encoding="utf-8")
    assert actual == expected, "docs/interactive/index.html is stale — regenerate it and commit"
