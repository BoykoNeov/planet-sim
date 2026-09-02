"""Triad for the complete equilibrium diagram + the small-ice-cap instability (:mod:`planet.bifurcation`).

*Tight (analytic).* (a) The finite-volume inverse curve converges to North's **Legendre-mode** solution
(even modes, exact piecewise Gauss–Legendre projection of the albedo step) at ~2nd order in Δx, in both
face modes. (b) **Reduction to Phase 1**: at a given sun the curve's stable finite-cap crossing is what the
nonlinear relaxation (``present_day_climate``) lands on, to within the marcher's *cell-quantized* ice edge
(the gap halves as the grid doubles — a Δx effect, not a Δt one).
*Tight (structural).* (c) The **slope-stability theorem** (Cahalan & North 1979) checked by marching: a
relaxation seeded on a stable segment returns, one seeded on an unstable segment departs to another
branch. (d) The Phase-1 **continuation sweep** — an independent method — freezes and re-melts within one
sweep step of the fold / branch-end the curve reads exactly. (e) Every curve equilibrium balances net TOA
exactly (conservation), and the theorem's verdict is consistent with the fold structure.
*Loose (the payoff).* Two folds at Earth parameters — the Snowball fold near 33° and the small-ice-cap
fold near 80° (``θ_c ≈ 10°``); present-day sits *inside and near the top of* the finite-cap window; ``θ_c``
grows with ``D`` and the window closes (the branch vanishes) at strong transport — North 1984.
"""
import numpy as np
import pytest
from dataclasses import replace

from planet import bifurcation as bf
from planet.albedo import EBMParams, present_day_climate, snowball_hysteresis
from planet.ebm import S0_EARTH, insolation


@pytest.fixture(scope="module")
def curve():
    return bf.equilibrium_curve()


# --------------------------------------------------------------------------- #
# TIGHT (a) — convergence to the Legendre-mode analytic solution.
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("face", ["harmonic", "exact"])
def test_fv_curve_converges_to_legendre_solution_at_second_order(face):
    # The inverse curve is a linear solve on the engine-pinned operator; North's even-Legendre expansion
    # is the continuous model's exact solution. The FV curve must approach it at ~2nd order in Δx — the
    # error quarters per grid doubling. Read on 0.05 < x_s < 0.95: the end faces are one-sided, and the
    # polar-most cells of the coarse grids are ~4° wide in latitude where the curve is steepest — the
    # asymptotic rate lives inside that window.
    errs = []
    for n in (45, 90, 180):
        c = bf.equilibrium_curve(EBMParams(n_cells=n, face=face))
        sel = (c.x_ice > 0.05) & (c.x_ice < 0.95)
        ref = bf.legendre_equilibrium_curve(c.x_ice[sel], c.params, n_modes=200)
        errs.append(float(np.max(np.abs(c.S0[sel] - ref) / ref)))
    assert errs[-1] < 1e-4                                   # tight at n=180
    assert errs[0] / errs[1] > 2.8 and errs[1] / errs[2] > 2.8   # ~2nd order (measured ratios 3.3–4.2)


def test_norths_two_mode_truncation_is_the_n2_special_case():
    # North 1975's original closed form is the P₀+P₂ truncation. It carries a ~1 % truncation error against
    # the many-mode solution (the step albedo needs the higher modes) — visible, and converging away with n.
    p = EBMParams()
    xs = np.array([0.3, 0.6, 0.9])
    two = bf.legendre_equilibrium_curve(xs, p, n_modes=2)
    many = bf.legendre_equilibrium_curve(xs, p, n_modes=150)
    more = bf.legendre_equilibrium_curve(xs, p, n_modes=300)
    assert np.max(np.abs(two - many) / many) > 2e-3          # the truncation is a real (small) error…
    assert np.max(np.abs(more - many) / many) < 1e-5         # …and the expansion has converged by 150 modes


# --------------------------------------------------------------------------- #
# TIGHT (b) — reduction to Phase 1's relaxation: the exact curve is what the marcher lands on.
# --------------------------------------------------------------------------- #
def test_relaxation_lands_on_the_curve_to_within_one_cell_and_the_gap_shrinks_with_resolution():
    # At S₀ = 1330 the finite-cap branch is steep enough that the ice line is a robust read. The nonlinear
    # relaxation quantizes the ice edge to a CELL (a whole cell flips albedo), the curve places it at a face
    # by interpolation — so they agree to within ~one cell's latitude width, and the gap halves as the grid
    # doubles: a Δx effect, not a Δt one (the O(Δt) fixed-point bias is below it at n_tau = 0.02).
    gaps = []
    for n in (60, 240):
        c = bf.equilibrium_curve(EBMParams(n_cells=n))
        e = c.stable_finite_cap_at(1330.0)
        st = present_day_climate(replace(c.params, S0=1330.0), n_tau=0.02, tol=1e-11, max_iter=400000)
        cell_deg = np.degrees(1.0 / n / np.cos(np.radians(e.latitude_deg)))
        gap = abs(st.ice_line_lat - e.latitude_deg)
        assert gap < 1.5 * cell_deg
        gaps.append(gap)
    assert gaps[1] < 0.5 * gaps[0]


# --------------------------------------------------------------------------- #
# TIGHT (c) — the slope-stability theorem, checked by marching.
# --------------------------------------------------------------------------- #
def _index_at(curve, lat_deg):
    return int(np.argmin(np.abs(curve.latitude_deg() - lat_deg)))


@pytest.mark.parametrize("lat_deg", [50.0, 65.0])
def test_stable_segment_returns_under_perturbation(curve, lat_deg):
    k = _index_at(curve, lat_deg)
    assert curve.stable[k]                                   # dS₀/dx_s > 0 here
    for dK in (-1.0, +1.0):
        st = bf.relax_from_curve(curve, k, dK, n_tau=0.02, tol=1e-12, max_iter=400000)
        assert abs(st.ice_line_lat - curve.latitude_deg()[k]) < 1.5


def test_unstable_middle_branch_departs_to_snowball_or_a_warm_branch(curve):
    # Between the Snowball fold and the equator the curve bends back (dS₀/dx_s < 0): a cooling nudge runs
    # away to the Snowball, a warming nudge to the finite-cap / ice-free branch — the same sun, the
    # equilibrium cannot hold. That IS the theorem's unstable verdict, seen by the nonlinear model.
    for lat_deg in (15.0, 25.0):
        k = _index_at(curve, lat_deg)
        assert not curve.stable[k]
        cold = bf.relax_from_curve(curve, k, -1.0, n_tau=0.02, tol=1e-12, max_iter=400000)
        warm = bf.relax_from_curve(curve, k, +1.0, n_tau=0.02, tol=1e-12, max_iter=400000)
        assert cold.ice_line_lat < 1.0                       # → Snowball
        assert warm.ice_line_lat > lat_deg + 10.0            # → a much smaller cap (a different branch)


def test_unstable_small_cap_departs_to_the_finite_cap_or_ice_free(curve):
    # Poleward of the small-ice-cap fold the tiny cap is unstable: warmed, it vanishes (ice-free); cooled,
    # it GROWS to the stable finite cap the same sun supports — the SICI jump, run forwards and backwards.
    f = curve.small_ice_cap_fold
    k = _index_at(curve, 85.0)
    assert curve.latitude_deg()[k] > f.latitude_deg and not curve.stable[k]
    warm = bf.relax_from_curve(curve, k, +1.0, n_tau=0.02, tol=1e-12, max_iter=400000)
    cold = bf.relax_from_curve(curve, k, -1.0, n_tau=0.02, tol=1e-12, max_iter=400000)
    assert warm.ice_line_lat > 89.0
    assert cold.ice_line_lat < f.latitude_deg              # crossed the fold onto the stable branch


def test_verdict_is_consistent_with_the_fold_structure(curve):
    # Stable exactly between the Snowball fold (min) and the small-ice-cap fold (max); unstable equatorward
    # of the min and poleward of the max — the two folds bracket the whole stable finite-cap branch.
    lo, hi = curve.snowball_fold, curve.small_ice_cap_fold
    lat = curve.latitude_deg()
    inside = (lat > lo.latitude_deg + 1.0) & (lat < hi.latitude_deg - 1.0)
    outside = (lat < lo.latitude_deg - 1.0) | ((lat > hi.latitude_deg + 1.0) & (lat < 89.0))
    assert curve.stable[inside].all()
    assert not curve.stable[outside].any()


# --------------------------------------------------------------------------- #
# TIGHT (d) — the independent Phase-1 continuation sweep jumps at the folds the curve reads.
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_continuation_sweep_jumps_within_one_step_of_the_exact_folds(curve):
    loop = snowball_hysteresis(n_steps=60, n_tau=0.05)
    step = float(np.abs(np.diff(loop.S0_up)).max())
    assert abs(loop.freeze_S0 - curve.snowball_fold.S0) < step          # the freeze = the Snowball fold
    assert abs(loop.melt_S0 - curve.snowball_threshold_S0) < step       # the re-melt = the Snowball branch end
    # the sweep rides the STABLE branch only: every finite-cap point it visits has a stable curve twin.
    # Where the branch is steep (below ~70°) the twin is within a few degrees; on the flat polar stretch
    # near the small-ice-cap fold a few W/m² of relaxation bias slide the ice line by many degrees, so only
    # existence is asserted there.
    for S0, lat in zip(loop.S0_down, loop.iceline_down):
        if 5.0 < lat < 85.0:
            e = curve.stable_finite_cap_at(S0)
            assert e is not None
            if lat < 70.0:
                assert abs(e.latitude_deg - lat) < 4.0


# --------------------------------------------------------------------------- #
# CONSERVATION — every equilibrium on the curve balances the top of the atmosphere exactly.
# --------------------------------------------------------------------------- #
def test_every_curve_point_balances_net_toa(curve):
    p = curve.params
    x = curve.x_cells
    ice_free = p.a0 + p.a2 * (0.5 * (3 * x * x - 1))
    for k in range(0, curve.S0.size, 15):
        alb = ice_free.copy(); alb[k:] = p.ai
        absorbed = insolation(x, curve.S0[k], p.s2) * (1.0 - alb)
        net = float(np.mean(absorbed) - p.A - p.B * curve.global_mean_T[k])
        assert abs(net) < 1e-9


# --------------------------------------------------------------------------- #
# LOOSE (the payoff) — two folds, a narrow window, present-day near its top, θ_c grows with D.
# --------------------------------------------------------------------------- #
def test_two_folds_and_the_finite_cap_window(curve):
    lo, hi = curve.snowball_fold, curve.small_ice_cap_fold
    assert lo is not None and hi is not None
    assert 25.0 < lo.latitude_deg < 40.0                     # the Snowball fold (Phase 1's ~8 % dimming)
    assert 0.90 < lo.S0 / S0_EARTH < 0.95
    assert 75.0 < hi.latitude_deg < 85.0                     # the small-ice-cap fold: θ_c ≈ 10°
    assert 7.0 < hi.cap_radius_deg < 14.0
    win = curve.finite_cap_window
    assert win is not None and win[0] == lo.S0 and win[1] == hi.S0
    # present day sits inside the window, close to its upper (SICI) edge — a few W/m² from losing its cap
    assert win[0] < S0_EARTH < win[1]
    assert (win[1] - S0_EARTH) / S0_EARTH < 0.01
    # the ice-free branch begins BELOW the SICI fold: between them, ice-free and finite-cap coexist
    assert curve.ice_free_threshold_S0 < hi.S0
    e = curve.stable_finite_cap_at(S0_EARTH)
    assert e is not None and 68.0 < e.latitude_deg < 80.0   # Earth's finite-cap branch (~70° benchmark)


def test_five_equilibria_at_present_day(curve):
    # At today's sun the model holds FIVE equilibria: ice-free (stable), Snowball (stable), the finite
    # cap (stable, Earth's), and two unstable separators (the mid-latitude one and the sub-critical small
    # cap) — the full bistable picture Phase 1's sweep could only sample two-thirds of.
    eqs = curve.equilibria_at(S0_EARTH)
    kinds = sorted(e.kind for e in eqs)
    assert kinds == ["finite-cap", "finite-cap", "finite-cap", "ice-free", "snowball"]
    assert sum(e.stable for e in eqs) == 3


def test_critical_cap_grows_with_transport_and_the_branch_vanishes():
    theta, lo, hi = bf.critical_cap_sweep([0.3, 0.555, 0.8, 1.0, 1.6])
    assert np.all(np.diff(theta[:4]) > 0.0)                  # θ_c grows with D over the Earth-relevant range
    assert np.all(np.diff(hi[:4] - lo[:4]) < 0.0)            # …and the finite-cap window narrows
    assert np.isnan(theta[4]) and np.isnan(hi[4])            # strong transport: no stable finite cap at all


def test_relaxation_bias_is_quantified_and_falls_with_the_step():
    # The Phase-1 O(Δt) fixed-point bias: at the default n_tau = 0.5 the relaxed ice line sits ~10° off the
    # exact curve; shrinking the step brings it onto the curve (to the cell quantization above).
    c = bf.equilibrium_curve()
    exact = c.stable_finite_cap_at(S0_EARTH).latitude_deg
    lats = bf.relaxation_bias_sweep([0.5, 0.05, 0.02])
    gaps = np.abs(lats - exact)
    assert gaps[0] > 5.0 and gaps[1] < 3.0 and gaps[2] < 2.5
    assert gaps[0] > gaps[1] > gaps[2]


def test_legendre_rejects_callable_D():
    with pytest.raises(ValueError, match="uniform scalar D"):
        bf.legendre_equilibrium_curve([0.5], EBMParams(D=lambda x: 0.5 + 0 * x))
