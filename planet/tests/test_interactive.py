"""Guards for the no-install browser what-if (:mod:`planet.interactive`).

Two tiers, like the rest of the repo:

* **fast** — the page *shell* is deterministic and self-contained (data inlined, no CDN), the grid
  has the right shape, and the committed page is well-formed. These build only a tiny 2×2 grid (a
  few model runs), so they run in the always-on lane.
* **slow** — the committed ``docs/interactive/index.html`` equals a fresh full-grid generation, so a
  change to the model or the page that isn't re-banked fails the gate. Marked ``slow`` because it
  reruns the whole (S0 × CO2 × tilt) knob grid (a few thousand model solves, minutes).
"""
from __future__ import annotations

import json
import os
import re

import pytest

from planet import interactive
from planet.biomes import Biome
from planet.obliquity import OBLIQUITY_EARTH

# The slow golden regenerates the whole (S0 × CO2 × tilt) grid of EBM solves and compares the page
# byte-for-byte. That is safe as a *local* drift guard (it runs on the machine that banked the page),
# but fragile cross-platform: a last-bit LAPACK difference on the Linux CI runner — especially one
# that flips a digit in a biome string near a Whittaker threshold — would fail the comparison for a
# non-bug. So, exactly like the notebook smoke-test, it is skipped in CI; the fast structural tests
# below cover correctness there. The local full gate still runs it. (REMOVE if the page is ever
# generated deterministically cross-platform, e.g. from a committed grid artifact rather than a live
# re-solve.)
_SKIP_IN_CI = os.environ.get("CI", "").lower() in {"true", "1"}

# A tiny grid whose (0,0,0) cell is the exact Earth detent (S0=1365, CO2=0, tilt=23.44°) → the
# baseline message. The obliquity axis must carry the exact OBLIQUITY_EARTH float as its first value
# (not a typed 23.44) so the obliquity factor is exactly 1 there and s₂ is bit-identical to the model.
_SMALL = dict(s0_values=[1365.0, 1375.0], co2_values=[0.0, 2.0],
              obliquity_values=[OBLIQUITY_EARTH, 45.0])


def test_compute_grid_shape():
    grid = interactive.compute_grid(**_SMALL)
    assert len(grid["cells"]) == 8                         # 2 × 2 × 2
    assert grid["axes"]["s0"]["values"] == [1365.0, 1375.0]
    assert grid["axes"]["co2"]["default_index"] == 0
    assert grid["axes"]["obl"]["default_index"] == 0       # OBLIQUITY_EARTH is the first tilt value
    assert len(grid["lat"]) == len(grid["cells"][0]["temp"]) == len(grid["cells"][0]["biome"])
    assert set(grid["palette"]) == set(grid["names"])      # a colour for every named biome
    cell = grid["cells"][0]
    assert all(ch in "012345678" for ch in cell["biome"])  # biome string is single-digit codes


def test_earth_detent_is_the_baseline():
    grid = interactive.compute_grid(**_SMALL)
    earth = grid["cells"][0]                                # (s0=1365, co2=0, tilt=Earth)
    assert "baseline" in earth["headline"].lower()
    assert 13.0 < earth["Tbar"] < 16.0


def test_obliquity_axis_moves_the_climate():
    """A flatter sky (less tilt) starves the poles → more ice and a cooler planet than Earth's tilt.

    Holds the Sun and greenhouse at present-day and walks only the tilt, so the (S0, CO2) detent is
    fixed and the response is the obliquity knob alone — the s₂(ε) gradient steepening at low tilt.
    """
    grid = interactive.compute_grid(s0_values=[1365.0], co2_values=[0.0],
                                    obliquity_values=[0.0, OBLIQUITY_EARTH])
    flat, earth = grid["cells"][0], grid["cells"][1]       # tilt = 0°, then Earth's 23.44°
    assert flat["ice"] < earth["ice"]                      # an untilted world has more polar ice
    assert flat["Tbar"] < earth["Tbar"]                    # and runs colder (stronger ice-albedo)
    assert "tilt" in flat["paragraph"].lower()             # the prose names the obliquity knob
    assert "baseline" in earth["headline"].lower()         # Earth's tilt at the detent stays baseline


def test_page_is_self_contained_and_deterministic():
    grid = interactive.compute_grid(**_SMALL)
    html = interactive.build_app_html(grid)
    assert html == interactive.build_app_html(grid)        # deterministic (no timestamp / set order)
    # data is inlined (file:// blocks fetch), and the page pulls in nothing external
    assert "window.PLANET_DATA = " in html
    assert "http://" not in html and "https://" not in html
    assert 'src="http' not in html
    assert '<canvas id="disk"' in html and '<canvas id="curve"' in html
    assert all(f'id="{knob}"' in html for knob in ("s0", "co2", "obl"))   # all three knob sliders
    assert 'id="tip"' in html and "mousemove" in html                     # the disk biome-hover read-out


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
    assert all(f'id="{knob}"' in html for knob in ("s0", "co2", "obl"))   # all three knob sliders
    assert 'id="tip"' in html and "mousemove" in html                     # the disk biome-hover read-out


def test_disk_hover_reads_the_biome_band():
    """The disk has a hover read-out: a tooltip element, a mousemove handler, and a name lookup.

    Pins the wiring (element + listener + ``D.names`` lookup) without a browser — the actual
    cursor→band mapping reuses ``drawDisk``'s math and is eyeballed in the play-through.
    """
    html = interactive.build_app_html(interactive.compute_grid(**_SMALL))
    assert '<div class="tip" id="tip">' in html        # the tooltip element
    assert 'addEventListener("mousemove"' in html      # follows the cursor over the disk
    assert "D.names[code]" in html                     # names the band it lands on
    assert "getBoundingClientRect" in html             # maps CSS px → canvas px before the |lat| math


def test_disk_band_lookup_is_nearest_latitude_not_linear():
    """The disk must place biome bands by *true nearest* latitude, not a linear index.

    The model grid is equal-area (uniform in sin φ), so the stored latitudes are non-uniform (~2° at
    the equator → ~6° at the pole, last band ~76°). A linear ``a/latMax·(N−1)`` index drags the warm
    bands poleward of the true-latitude ice-line ring (forest on the ice cap). Both the fill and the
    hover must go through one shared nearest-latitude helper so they cannot drift.
    """
    html = interactive.build_app_html(interactive.compute_grid(**_SMALL))
    assert "function bandForLat(" in html              # the shared nearest-latitude helper
    assert html.count("bandForLat(phi)") >= 2          # used by BOTH the disk fill and the hover
    assert "latMax) * (N - 1)" not in html             # the broken linear index is gone


def _band_for_lat(lat, phi):
    """Python replica of the page's ``bandForLat``: index of the nearest stored latitude."""
    a = abs(phi)
    return min(range(len(lat)), key=lambda m: abs(lat[m] - a))


@pytest.mark.parametrize("s0, co2, obl", [
    (1265.0, 2.0, 35.0),            # the user-flagged dim/tilted cell that painted boreal past the ring
    (1265.0, 0.0, OBLIQUITY_EARTH), # dim Sun, Earth tilt — ice line drops to ~40°, boreal just inside it
    (1265.0, 9.0, OBLIQUITY_EARTH), # warmed back up — ice line ~59°, boreal hugging it
])
def test_disk_paints_no_forest_on_the_ice_cap(s0, co2, obl):
    """Render-path guard: replay the pixel→band map and assert no warm biome sits poleward of the ring.

    The bug hid from a *model* scan (the field is correct) because it lived in the *render*. Poleward
    of the ice line the surface is below −10 °C, so the only valid biome is tundra (boreal ends at the
    warmer −5 °C isotherm, which is equatorward of the ring). This replays ``drawDisk``'s fixed
    nearest-latitude lookup over the disk's pixel rows and pins that physics on the screen.
    """
    grid = interactive.compute_grid(s0_values=[s0], co2_values=[co2], obliquity_values=[obl])
    cell, lat, ice = grid["cells"][0], grid["lat"], grid["cells"][0]["ice"]
    assert 1 < ice < 89, "this cell must have a finite ice cap for the ring to be drawn"
    R = 144                                            # the shipped disk radius (canvas 300 → R≈144)
    boreal_phis = []
    for y in range(R + 1):                             # disk pixel rows, equator (0) → pole (90°)
        phi = (y / R) * 90.0
        code = int(cell["biome"][_band_for_lat(lat, phi)])
        if phi > ice:
            assert code == int(Biome.TUNDRA), (
                f"biome {Biome(code).name} painted at |lat|≈{phi:.0f}° — poleward of the {ice:.0f}° "
                f"ice ring (only tundra is valid below −10 °C)")
        if code == int(Biome.BOREAL_FOREST):
            boreal_phis.append(phi)
    assert max(boreal_phis) < ice                      # the boreal band lies entirely equatorward of the ring


@pytest.mark.slow
@pytest.mark.skipif(
    _SKIP_IN_CI,
    reason="byte-exact over the full (S0 × CO2 × tilt) grid of live EBM solves — fragile cross-platform "
    "(LAPACK last-bit near a Whittaker biome threshold); a local-only drift guard. Fast structural "
    "tests cover CI.",
)
def test_committed_page_is_up_to_date():
    """docs/interactive/index.html must equal a fresh full-grid build (re-run `python -m planet.interactive`)."""
    expected = interactive.build_app_html(interactive.compute_grid())
    actual = interactive.APP_PATH.read_text(encoding="utf-8")
    assert actual == expected, "docs/interactive/index.html is stale — regenerate it and commit"
