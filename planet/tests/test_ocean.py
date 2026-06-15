"""Guards for the ocean-fraction knob (:mod:`planet.ocean`).

Like the obliquity knob, this is a *parameter derivation* with one hard contract — **Earth's sea
fraction recovers the model exactly** (a clean perturbation) — plus the cited directions and
magnitudes of its two channels (the firm albedo leg, the loose transport leg). The end-to-end check
that a wetter world actually runs warmer pins the sign through a real EBM solve.
"""
from __future__ import annotations

import pytest

from planet.albedo import EBMParams, present_day_climate
from planet.ebm import ALBEDO_A0, D_TRANSPORT
from planet.obliquity import obliquity_params
from planet.ocean import (
    A0_LAND_MINUS_OCEAN, OCEAN_FRACTION_EARTH,
    ocean_albedo_a0, ocean_params, ocean_transport_D,
)


def test_earth_fraction_is_a_clean_perturbation():
    """At Earth's sea fraction both channels are the identity → params equal the model defaults."""
    assert ocean_albedo_a0(OCEAN_FRACTION_EARTH) == ALBEDO_A0          # exactly, no offset
    assert ocean_transport_D(OCEAN_FRACTION_EARTH) == D_TRANSPORT      # exactly, unit factor
    assert ocean_params(OCEAN_FRACTION_EARTH) == EBMParams()           # the whole bundle unchanged


def test_more_ocean_darkens_and_carries_more_heat():
    """Monotone: a wetter world has a lower a0 (darker, warmer) and a larger D (more transport)."""
    assert ocean_albedo_a0(1.0) < ocean_albedo_a0(OCEAN_FRACTION_EARTH) < ocean_albedo_a0(0.0)
    assert ocean_transport_D(1.0) > ocean_transport_D(OCEAN_FRACTION_EARTH) > ocean_transport_D(0.0)


def test_albedo_swing_matches_the_cited_magnitude():
    """The all-land minus all-ocean planetary-albedo swing is the pinned A0_LAND_MINUS_OCEAN."""
    assert ocean_albedo_a0(0.0) - ocean_albedo_a0(1.0) == pytest.approx(A0_LAND_MINUS_OCEAN)


def test_fraction_is_clamped_to_the_unit_interval():
    assert ocean_albedo_a0(-0.5) == ocean_albedo_a0(0.0)
    assert ocean_albedo_a0(1.5) == ocean_albedo_a0(1.0)
    assert ocean_transport_D(2.0) == ocean_transport_D(1.0)


def test_params_touch_only_a0_and_D():
    """The knob replaces a0 and D and nothing else (the disjoint-channel contract)."""
    p = ocean_params(0.30)
    base = EBMParams()
    assert p.a0 != base.a0 and p.D != base.D
    assert (p.S0, p.s2, p.A, p.B, p.T_freeze, p.a2, p.ai) == (
        base.S0, base.s2, base.A, base.B, base.T_freeze, base.a2, base.ai)


def test_commutes_with_obliquity():
    """Ocean (a0/D) and obliquity (s2) are disjoint, so the two knobs compose in either order."""
    eps, w = 35.0, 0.30
    assert ocean_params(w, obliquity_params(eps)) == obliquity_params(eps, ocean_params(w))


def test_a_wetter_world_runs_warmer_end_to_end():
    """The sign through a real EBM solve: more ocean → a warmer equilibrium global mean."""
    land = present_day_climate(ocean_params(0.0)).global_mean_T
    sea = present_day_climate(ocean_params(1.0)).global_mean_T
    assert sea > land
