"""Guards for the no-install browser what-if (:mod:`planet.interactive`).

Two tiers, like the rest of the repo:

* **fast** — the page *shell* is deterministic and self-contained (data inlined, no CDN), the grid
  has the right shape, and the committed page is well-formed. These build only a tiny 2×2 grid (a
  few model runs), so they run in the always-on lane.
* **slow** — the committed ``docs/interactive/index.html`` equals a fresh full-grid generation, so a
  change to the model or the page that isn't re-banked fails the gate. Marked ``slow`` because it
  reruns the whole (S0 × CO2 × tilt × ocean) knob grid (~ten thousand model solves, minutes).
"""
from __future__ import annotations

import json
import os
import re

import numpy as np
import pytest

from planet import interactive
from planet.biomes import Biome
from planet.obliquity import OBLIQUITY_EARTH
from planet.ocean import OCEAN_FRACTION_EARTH

# The slow golden regenerates the whole (S0 × CO2 × tilt × ocean) grid of EBM solves and compares the page
# byte-for-byte. That is safe as a *local* drift guard (it runs on the machine that banked the page),
# but fragile cross-platform: a last-bit LAPACK difference on the Linux CI runner — especially one
# that flips a digit in a biome string near a Whittaker threshold — would fail the comparison for a
# non-bug. So, exactly like the notebook smoke-test, it is skipped in CI; the fast structural tests
# below cover correctness there. The local full gate still runs it. (REMOVE if the page is ever
# generated deterministically cross-platform, e.g. from a committed grid artifact rather than a live
# re-solve.)
_SKIP_IN_CI = os.environ.get("CI", "").lower() in {"true", "1"}

# A tiny grid whose (0,0,0,0) cell is the exact Earth detent (S0=1365, CO2=0, tilt=23.44°, ocean=0.71)
# → the baseline message. The obliquity and ocean axes must carry the exact OBLIQUITY_EARTH and
# OCEAN_FRACTION_EARTH floats as their first value (not a typed 23.44 / 0.71) so both knobs are exactly
# the identity there and s₂/a0/D are bit-identical to the model.
_SMALL = dict(s0_values=[1365.0, 1375.0], co2_values=[0.0, 2.0],
              obliquity_values=[OBLIQUITY_EARTH, 45.0], ocean_values=[OCEAN_FRACTION_EARTH, 1.0])


def test_compute_grid_shape():
    grid = interactive.compute_grid(**_SMALL)
    assert len(grid["cells"]) == 16                        # 2 × 2 × 2 × 2
    assert grid["axes"]["s0"]["values"] == [1365.0, 1375.0]
    assert grid["axes"]["co2"]["default_index"] == 0
    assert grid["axes"]["obl"]["default_index"] == 0       # OBLIQUITY_EARTH is the first tilt value
    assert grid["axes"]["ocean"]["values"] == [71, 100]    # stored as whole-percent of surface
    assert grid["axes"]["ocean"]["default_index"] == 0     # OCEAN_FRACTION_EARTH is the first value
    assert len(grid["lat"]) == len(grid["cells"][0]["temp"]) == len(grid["cells"][0]["biome"])
    assert set(grid["palette"]) == set(grid["names"])      # a colour for every named biome
    cell = grid["cells"][0]
    assert all(ch in "012345678" for ch in cell["biome"])  # biome string is single-digit codes
    assert len(grid["cold_cells"]) == 4                    # 2 S0 × 2 CO2 — the cold branch ignores tilt/ocean


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
                                    obliquity_values=[0.0, OBLIQUITY_EARTH],
                                    ocean_values=[OCEAN_FRACTION_EARTH])   # pin ocean so only tilt moves
    flat, earth = grid["cells"][0], grid["cells"][1]       # tilt = 0°, then Earth's 23.44°
    assert flat["ice"] < earth["ice"]                      # an untilted world has more polar ice
    assert flat["Tbar"] < earth["Tbar"]                    # and runs colder (stronger ice-albedo)
    assert "tilt" in flat["paragraph"].lower()             # the prose names the obliquity knob
    assert "baseline" in earth["headline"].lower()         # Earth's tilt at the detent stays baseline


def test_ocean_axis_moves_the_climate():
    """More ocean darkens the surface → a warmer planet than a drier world at the same Sun/CO₂/tilt.

    Holds the other three knobs at present-day and walks only the ocean fraction from a land world to
    Earth's sea fraction, so the response is the ocean knob alone (lower a0 + a touch more transport).
    """
    grid = interactive.compute_grid(s0_values=[1365.0], co2_values=[0.0],
                                    obliquity_values=[OBLIQUITY_EARTH],
                                    ocean_values=[0.0, OCEAN_FRACTION_EARTH])
    land, earth = grid["cells"][0], grid["cells"][1]       # 0% ocean, then Earth's 71%
    assert land["Tbar"] < earth["Tbar"]                    # a drier, brighter world runs colder
    assert "ocean" in land["paragraph"].lower()            # the prose names the ocean knob
    assert "baseline" in earth["headline"].lower()         # Earth's sea fraction stays the baseline


def test_cold_branch_is_the_snowball_toggle():
    """The cold (Snowball) sub-grid is the §2 hysteresis: the same knobs, a frozen start.

    It is keyed by (S0, CO2) only — on a Snowball the mean T and ice cover are exactly tilt/ocean-
    independent (uniform ice albedo) — and every cell is frozen over. At Earth's Sun the warm branch
    is temperate while the cold branch is a Snowball: two stable climates for one knob (the bistability
    the toggle exposes).
    """
    grid = interactive.compute_grid(**_SMALL)
    assert len(grid["cold_cells"]) == 4                    # 2 S0 × 2 CO2
    for c in grid["cold_cells"]:
        assert c["ice"] <= 1.0                             # every cold cell is frozen over (Snowball)
        assert "snowball" in c["headline"].lower()
    warm, cold = grid["cells"][0], grid["cold_cells"][0]   # both at (S0=1365, CO2=0)
    assert warm["ice"] > 60 and cold["ice"] <= 1.0         # temperate vs frozen at the SAME Sun
    assert "started warm sits temperate" in cold["oneline"]  # the path-dependence, named


def test_cold_branch_mean_is_tilt_independent_but_curve_is_not():
    """The invariant that justifies the lean 2-D cold sub-grid — verified, not just asserted in prose.

    On a Snowball the global-mean T and ice cover are *exactly* tilt/ocean-independent (uniform ice
    albedo → mean = [(S₀/4)(1−α_ice)−A]/B, no obliquity/transport term; frozen to the equator either
    way), which is why the cold data is stored at Earth's tilt/ocean only. The per-latitude *profile*
    is NOT invariant — obliquity reshapes the insolation — so the page's hint scopes its claim to the
    mean/ice and labels the shown curve as Earth's tilt (this test pins both halves of that honesty).
    """
    from planet.demo_biomes import compute
    from planet.obliquity import obliquity_params
    from planet.albedo import EBMParams, A_OLR, S0_EARTH

    earth = compute(EBMParams(S0=S0_EARTH, A=A_OLR), ic_equator=-40.0, ic_pole=-40.0)
    tilted = compute(obliquity_params(45.0, EBMParams(S0=S0_EARTH, A=A_OLR)),
                     ic_equator=-40.0, ic_pole=-40.0)
    # The mean is invariant to the continuum identity B·⟨T⟩ = (S₀/4)(1−α_ice)−A; the discrete grid
    # leaks ~1e-4 °C, far below the cell's 2-dp rounding, so the *stored/displayed* value is identical.
    assert round(earth.state.global_mean_T, 2) == round(tilted.state.global_mean_T, 2)   # mean: invariant
    assert round(earth.state.ice_line_lat, 1) == round(tilted.state.ice_line_lat, 1)     # ice: invariant
    curve_gap = float(np.max(np.abs(earth.state.T - tilted.state.T)))
    assert curve_gap > 1.0     # but the profile shifts several °C — why the curve is labeled Earth's


def test_starting_climate_toggle_is_wired():
    """The page ships the warm/frozen segmented control, the cold sub-grid, and the render branch."""
    html = interactive.build_app_html(interactive.compute_grid(**_SMALL))
    assert 'id="warm-btn"' in html and 'id="cold-btn"' in html   # the segmented control
    assert '"cold_cells":' in html                               # the Snowball sub-grid is inlined
    assert "coldCell(i, j)" in html                              # render() reads the cold branch when toggled
    # the disclosure hint: only the frozen mean/ice are tilt/ocean-independent; the curve is Earth's
    assert "ignore tilt & ocean" in html and "shown at Earth's tilt" in html


def test_page_is_self_contained_and_deterministic():
    grid = interactive.compute_grid(**_SMALL)
    html = interactive.build_app_html(grid)
    assert html == interactive.build_app_html(grid)        # deterministic (no timestamp / set order)
    # data is inlined (file:// blocks fetch), and the page pulls in nothing external
    assert "window.PLANET_DATA = " in html
    assert "http://" not in html and "https://" not in html
    assert 'src="http' not in html
    assert '<canvas id="disk"' in html and '<canvas id="curve"' in html
    assert all(f'id="{knob}"' in html for knob in ("s0", "co2", "obl", "ocean"))  # all four knob sliders
    assert 'id="warm-btn"' in html and 'id="cold-btn"' in html             # the starting-climate toggle
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
    assert all(f'id="{knob}"' in html for knob in ("s0", "co2", "obl", "ocean"))  # all four knob sliders
    assert 'id="warm-btn"' in html and 'id="cold-btn"' in html             # the starting-climate toggle
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
    grid = interactive.compute_grid(s0_values=[s0], co2_values=[co2], obliquity_values=[obl],
                                    ocean_values=[OCEAN_FRACTION_EARTH])   # Earth sea fraction
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
    reason="byte-exact over the full (S0 × CO2 × tilt × ocean) grid of live EBM solves — fragile cross-platform "
    "(LAPACK last-bit near a Whittaker biome threshold); a local-only drift guard. Fast structural "
    "tests cover CI.",
)
def test_committed_page_is_up_to_date():
    """docs/interactive/index.html must equal a fresh full-grid build (re-run `python -m planet interactive`)."""
    expected = interactive.build_app_html(interactive.compute_grid())
    actual = interactive.APP_PATH.read_text(encoding="utf-8")
    assert actual == expected, "docs/interactive/index.html is stale — regenerate it and commit"
