"""Planet Phase-1 validation: the frozen climlab reference table + the live climlab cross-check.

The benchmark leg's two halves: (1) the **frozen** reference table keeps the triad green without the
``[climate]`` extra — here we check it is internally sane and that its pinned climlab defaults are the
*same numbers* :mod:`planet.ebm` computes with (single source of truth); (2) the **live**
climlab cross-check is ``slow`` / ``importorskip`` — it runs only where the opt-in ``[climate]`` extra
is installed, comparing climlab's own present-day climate to the reference bands.
"""
import pytest

from planet import ebm
from planet.climate_reference import REFERENCE, climlab_present_day


def test_reference_climlab_defaults_match_the_model_constants():
    # The frozen table's climlab defaults ARE the ebm.py constants — so the benchmark and the model
    # cannot silently disagree on what "climlab defaults" means.
    assert REFERENCE.climlab_S0 == ebm.S0_EARTH
    assert REFERENCE.climlab_A == ebm.A_OLR
    assert REFERENCE.climlab_B == ebm.B_OLR
    assert REFERENCE.climlab_D == ebm.D_TRANSPORT
    assert REFERENCE.climlab_Tf == ebm.T_FREEZE
    assert REFERENCE.climlab_a0 == ebm.ALBEDO_A0
    assert REFERENCE.climlab_a2 == ebm.ALBEDO_A2
    assert REFERENCE.climlab_ai == ebm.ALBEDO_ICE
    assert REFERENCE.climlab_s2 == ebm.S2_INSOLATION


def test_reference_bands_are_sane():
    lo, hi = REFERENCE.present_ice_line_band
    assert lo < REFERENCE.present_ice_line_deg < hi
    glo, ghi = REFERENCE.present_global_mean_band
    assert glo < REFERENCE.present_global_mean_C < ghi
    assert REFERENCE.snowball_dimming_pct_band[0] < REFERENCE.snowball_dimming_pct_band[1]
    assert REFERENCE.snowball_global_mean_max_C < 0.0
    assert REFERENCE.hysteresis_positive


@pytest.mark.slow
def test_live_climlab_present_day_cross_check():
    # The pycalphad pattern: climlab is the opt-in reference tool — skip cleanly when it is absent
    # (the committed reality). Where present, its present-day climate must land in the same bands.
    pytest.importorskip("climlab")
    global_mean, ice_line = climlab_present_day()
    glo, ghi = REFERENCE.present_global_mean_band
    assert glo < global_mean < ghi
    lo, hi = REFERENCE.present_ice_line_band
    assert lo < ice_line < hi
