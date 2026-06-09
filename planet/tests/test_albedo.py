"""Planet Phase-1 validation: ice-albedo feedback & the Snowball hysteresis — the banked-payoff triad.

The nonlinear half of the plan §3 triad (planet-earth-system.md). What is asserted **tight** is
*structural/qualitative* — the existence and shape of the bistability (two stable branches, a
catastrophic freezing jump, a re-melt at a brighter sun than the freeze); what is **loose** is the
*calibrated thresholds* (the exact freeze/melt S₀, the ice-line latitude), bounded against the cited
climlab reference bands ([[ebm-radiation-source]] / :mod:`projects.planet.climate_reference`). Nothing
but the feedback produces the hysteresis — it is emergent, so it is asserted firmly.
"""
import numpy as np
import pytest

from projects.planet import ebm, albedo
from projects.planet.climate_reference import REFERENCE


# --------------------------------------------------------------------------- #
# The step-function albedo and the absorbed shortwave
# --------------------------------------------------------------------------- #
def test_planetary_albedo_is_a_step_at_the_freeze_isotherm():
    x = np.linspace(0.0, 1.0, 5)
    warm = albedo.planetary_albedo(x, np.full_like(x, 20.0))             # ice-free
    cold = albedo.planetary_albedo(x, np.full_like(x, -20.0))            # iced
    assert np.allclose(cold, ebm.ALBEDO_ICE)
    assert np.allclose(warm, ebm.ALBEDO_A0 + ebm.ALBEDO_A2 * ebm.legendre_P2(x))
    # a mixed field: warm cells keep the ice-free value, frozen cells jump to the ice value
    a = albedo.planetary_albedo(x, np.array([20.0, 20.0, -20.0, -20.0, -20.0]))
    assert a[0] < ebm.ALBEDO_ICE and a[-1] == ebm.ALBEDO_ICE


def test_absorbed_shortwave_is_insolation_times_coalbedo():
    x = np.linspace(0.0, 1.0, 10)
    got = albedo.absorbed_shortwave(x, np.full_like(x, 20.0))
    want = ebm.insolation(x) * (1.0 - (ebm.ALBEDO_A0 + ebm.ALBEDO_A2 * ebm.legendre_P2(x)))
    assert np.allclose(got, want)


# --------------------------------------------------------------------------- #
# Benchmark (loose): the present-day finite-cap branch + the bistability
# --------------------------------------------------------------------------- #
def test_present_day_is_the_finite_cap_branch():
    pd = albedo.present_day_climate(n_tau=0.02)
    lo, hi = REFERENCE.present_ice_line_band
    assert lo < pd.ice_line_lat < hi                                    # ice line ~70° (loose band)
    glo, ghi = REFERENCE.present_global_mean_band
    assert glo < pd.global_mean_T < ghi                                 # global mean ~14–15 °C
    # Conservation in a feedback state: net-TOA is small but the albedo *discontinuity* at the ice
    # line limits exact closure (the machine-exact conservation is the no-feedback test_ebm leg).
    assert pd.net_toa == pytest.approx(0.0, abs=0.5)


def test_present_day_warm_start_finds_the_ice_free_branch():
    # At the SAME present S₀, a warm-uniform start settles on the (equally valid) ice-free branch —
    # the very bistability the Snowball loop traces, visible already at present insolation.
    params = albedo.EBMParams()
    ice_free = params.model().equilibrate(params.absorbed_fn(), 40.0, n_tau=0.02)
    assert ice_free.ice_line_lat == 90.0                               # no cap (warm branch)
    assert ice_free.global_mean_T > albedo.present_day_climate(n_tau=0.02).global_mean_T  # warmer


# --------------------------------------------------------------------------- #
# Analytical/structural (tight): the hysteresis exists and has the right shape
# --------------------------------------------------------------------------- #
def test_snowball_hysteresis_structure():
    # A coarse continuation sweep is enough to expose the STRUCTURE (the exact thresholds shift with
    # resolution/dt and are only asserted loosely below). Nothing but the feedback produces this.
    loop = albedo.snowball_hysteresis(
        params=albedo.EBMParams(n_cells=90), S0_min=1050.0, S0_max=1900.0, n_steps=14, n_tau=0.1,
    )
    # (1) a positive-width loop: the white planet re-melts only at a BRIGHTER sun than it froze
    assert loop.melt_S0 > loop.freeze_S0
    assert loop.hysteresis_width > 50.0
    assert REFERENCE.hysteresis_positive
    # (2) a catastrophic jump on the dimming branch (a large discontinuous drop in T̄)
    assert np.min(np.diff(loop.Tbar_down)) < -20.0
    # (3) the Snowball state is deeply frozen
    assert loop.Tbar_down.min() < REFERENCE.snowball_global_mean_max_C
    # (4) bistability: at an S₀ inside the loop the dimming (warm) and brightening (cold) branches differ
    mid = 0.5 * (loop.freeze_S0 + loop.melt_S0)
    di = int(np.argmin(np.abs(loop.S0_down - mid)))
    ui = int(np.argmin(np.abs(loop.S0_up - mid)))
    assert loop.Tbar_down[di] - loop.Tbar_up[ui] > 20.0


def test_snowball_branch_is_fully_iced_and_warm_branch_is_not():
    loop = albedo.snowball_hysteresis(
        params=albedo.EBMParams(n_cells=90), S0_min=1050.0, S0_max=1900.0, n_steps=14, n_tau=0.1,
    )
    assert loop.iceline_down.min() == 0.0          # the down sweep reaches a full Snowball (ice to equator)
    assert loop.iceline_up.max() == 90.0           # the up sweep reaches an ice-free planet again
