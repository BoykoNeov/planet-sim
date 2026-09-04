"""Triad for the **seasonal** small-ice-cap instability sweep (:mod:`planet.seasonal_sici`, rung 5B.4).

Rung 0+ found a second fold in the *annual-mean* ice-albedo EBM — a polar cap below a critical radius
``θ_c`` that no sun can hold, so the cap vanishes in a **jump** and leaves a hysteresis loop in ``S₀``.
Rung 5B.1+ put the same ice-albedo feedback on the *seasonal* marcher and named the bifurcation study as
deferred. This is it, and the question it asks is whether a cap that melts every summer is even the object
that fold is about.

*Tight (structural, machine).* The **annual-mean parent** :func:`~planet.seasonal_sici.annual_mean_curve` —
the seasonal model's own operator and its own annual-mean insolation, inverse-solved — satisfies the steady
EBM it claims to solve, hits ``T_f`` exactly at the prescribed ice line, and is hemispherically symmetric;
its forcing equals the pinned :func:`planet.obliquity.annual_mean_insolation` kernel to machine precision.
*Tight (cross-model).* That parent finds the **same fold** as rung 0+'s independent hemisphere solve
(:func:`planet.bifurcation.equilibrium_curve`), differing only by the ``P₂`` truncation of the insolation —
so ~1° in ``θ_c``, **not** machine (stated in advance, not discovered).
*Convergence (the reduction to the parent).* With the seasons off (``ε = 0``) the marched ice line converges
onto the parent curve's exact ice line at **first order in dt** — the Strang-splitting rate, the same shape
as 5B.1+'s ε=0 self-consistency check. This is what ties child to parent.
*Conservation.* A swept point still balances global-and-annual net TOA with the realized co-albedo.
*Unit (the new numerical objects).* The loop-width detector and the cap-planting seed are pinned on
synthetic input — including the below-threshold case, so *"no loop"* and *"a loop below detection"* cannot
silently merge.
*Positive control (the one that makes the null result readable).* Every other marched call to the loop
detector here returns zero, which on its own is also what a **blind** detector returns. So the same detector
and the same sweep are pointed at a deep (800 m) mixed layer, where the fold genuinely survives, and must
find it at the parent's own scale.
*Loose (the payoff, slow-marked).* At Earth's tilt with a 50 m mixed layer the perennial cap grows **one
cell at a time** through the parent's critical radius with **no detectable hysteresis**, where the parent's
own loop is a clear ~7 W m⁻²; deepening the mixed layer brings the bistability back. Direction banked, the
magnitudes ride the calibrated heat capacities.
"""
import math

import numpy as np
import pytest

from dataclasses import replace

from planet import seasonal_sici as sici
from planet.albedo import EBMParams
from planet.bifurcation import equilibrium_curve
from planet.ebm import T_FREEZE
from planet.obliquity import annual_mean_insolation, insolation_s2
from planet.seasonal import ice_coalbedo, ice_edge_latitude


# --------------------------------------------------------------------------- #
# TIGHT (structural) — the parent curve solves the equations it says it solves.
# --------------------------------------------------------------------------- #
def test_parent_curve_satisfies_the_steady_ebm_and_pins_the_ice_line():
    # Every point on the inverse-solved curve is an exact equilibrium of the annual-mean model: the residual
    # L_T·T + S₀⟨s⟩(1−α) − A − B·T vanishes to machine precision (it is a linear solve, so nothing weaker
    # would be acceptable), and the ice-line condition T(x_s) = T_f holds at the prescribed face. Together
    # these pin the whole inverse solve — the operator, the two-hemisphere albedo mask, and the face read.
    cfg = sici.SICIConfig(n_cells=180, n_steps=180)
    c = sici.annual_mean_curve(cfg)
    m = cfg.model()
    s_unit = m.insolation_series().mean(axis=1) / m.S0
    ice_free = 1.0 - m.coalbedo()
    a_ice = EBMParams().ai
    n = m.n_cells
    for j in (5, 30, 60, 89):                                     # equatorward → poleward ice lines
        k = n // 2 + j
        alb = ice_free.copy()
        alb[k:] = a_ice
        alb[:n - k] = a_ice
        T, S0 = c.T[j], c.S0[j]
        residual = m._apply_LT(T) + S0 * s_unit * (1.0 - alb) - m.A - m.B * T
        assert np.max(np.abs(residual)) < 1e-8                     # ~1e-10 relative to the ~340 W/m² forcing
        T_face = 0.5 * (T[k - 1] + T[k])
        assert abs(T_face - T_FREEZE) < 1e-9                       # the ice-line condition, exactly


def test_parent_curve_is_hemispherically_symmetric():
    # The annual mean is symmetric about the equator and the prescribed cap is mirrored, so every profile on
    # the curve must be even in x. A one-sided albedo mask (the easy slip in going full-sphere) reddens here.
    cfg = sici.SICIConfig(n_cells=180, n_steps=180)
    c = sici.annual_mean_curve(cfg)
    for j in (5, 40, 88):
        assert np.max(np.abs(c.T[j] - c.T[j][::-1])) < 1e-9


def test_parent_forcing_is_the_pinned_annual_mean_insolation_kernel():
    # The parent integrates the MARCHER's own time-mean insolation — not the P₂-truncated insolation() — so
    # the reduction is exact rather than truncation-limited (the 5B.1 lesson). That time mean is the pinned
    # obliquity kernel: at n_steps = 720 both are the same uniform-in-orbit sampling, so they agree to
    # machine precision, and a coarser year converges onto it.
    cfg = sici.SICIConfig(n_cells=90, n_steps=720)
    m = cfg.model()
    own = m.insolation_series().mean(axis=1) / m.S0
    pinned = annual_mean_insolation(m.phi, cfg.obliquity_deg) / math.pi
    assert np.max(np.abs(own - pinned)) < 1e-14
    coarse = sici.SICIConfig(n_cells=90, n_steps=180).model()
    assert np.max(np.abs(coarse.insolation_series().mean(axis=1) / coarse.S0 - pinned)) < 1e-5


# --------------------------------------------------------------------------- #
# TIGHT (cross-model) — the same fold rung 0+ found, by an independent implementation.
# --------------------------------------------------------------------------- #
def test_parent_curve_finds_the_rung0_small_ice_cap_fold():
    # Two independent solves of the same continuous problem: rung 0+'s hemisphere curve (P₂-truncated
    # insolation, EnergyBalanceModel's operator) and this full-sphere one (the seasonal model's operator and
    # its true annual-mean insolation). Both must find a polar maximum — the SICI — at a comparable critical
    # radius and sun. They are NOT expected to agree to machine precision, and this is stated in advance:
    # the true annual-mean insolation carries moments beyond P₂ that the truncation drops, which moves θ_c
    # by ~1° and S₀ by a few W/m² (~0.3 %). A disagreement of that size is the truncation, not a bug.
    cfg = sici.SICIConfig(n_cells=1440, n_steps=720)
    ours = sici.annual_mean_curve(cfg)
    theirs = equilibrium_curve(replace(EBMParams(), s2=insolation_s2(cfg.obliquity_deg), n_cells=720))
    f_ours, f_theirs = ours.small_ice_cap_fold, theirs.small_ice_cap_fold
    assert f_ours is not None and f_theirs is not None            # both find the second fold
    assert abs(f_ours.cap_radius_deg - f_theirs.cap_radius_deg) < 1.5      # the P₄+ moments, ~1°
    assert abs(f_ours.S0 - f_theirs.S0) / f_theirs.S0 < 0.01               # ~0.3 % in the sun
    assert ours.finite_cap_window is not None and theirs.finite_cap_window is not None


def test_parent_fold_converges_with_resolution():
    # θ_c and the loop width must CONVERGE, not drift — otherwise the reference the seasonal sweep is read
    # against is itself resolution noise (which is exactly what happens at lower obliquity, and why tilt was
    # rejected as the sweep axis). Refining 720 → 1440 → 2880 must move both by successively less.
    folds, loops = [], []
    for n in (720, 1440, 2880):
        c = sici.annual_mean_curve(sici.SICIConfig(n_cells=n, n_steps=720))
        f = c.small_ice_cap_fold
        assert f is not None
        folds.append(f.cap_radius_deg)
        loops.append(f.S0 - c.ice_free_threshold_S0)
    assert abs(folds[2] - folds[1]) < abs(folds[1] - folds[0])
    assert abs(loops[2] - loops[1]) < abs(loops[1] - loops[0])
    assert abs(folds[2] - folds[1]) < 0.1                          # θ_c settled to ~0.1°
    assert loops[2] > 5.0                                          # and the parent's loop is a real ~7 W/m²


# --------------------------------------------------------------------------- #
# CONVERGENCE — the reduction to the parent: seasons off ⟹ the exact curve, at first order in dt.
# --------------------------------------------------------------------------- #
def _epsilon0_ice_line_gap(n_steps: int) -> float:
    """|marched − exact| ice-line latitude at ε=0, seeded from the parent's own equilibrium (degrees)."""
    cfg = sici.SICIConfig(n_cells=180, n_steps=n_steps, obliquity_deg=0.0)
    c = sici.annual_mean_curve(cfg)
    idx = int(np.argmin(np.abs(c.latitude_deg() - 70.0)))          # well inside the stable finite-cap branch
    m = cfg.model(float(c.S0[idx]))
    cl = m.march(coalbedo_fn=ice_coalbedo, T_init=np.asarray(c.T[idx]), tol=1e-11, max_years=2000)
    assert cl.converged
    return abs(ice_edge_latitude(m.x, cl.T_ocean, kind="perennial") - float(c.latitude_deg()[idx]))


def test_epsilon0_marched_ice_line_converges_to_the_parent_at_first_order():
    # With the tilt at zero there is no seasonal cycle to average over, so the marcher's limit cycle IS an
    # equilibrium and must sit on the parent curve — offset only by the O(Δt) Strang-splitting bias that
    # rung 0+ had to quantify. Halving dt must halve the gap. This is the reduction that ties this module's
    # new parent to the marcher it is used to judge; a fixed tolerance would pass on a wrong operator, a
    # first-order rate will not.
    g1, g2, g3 = (_epsilon0_ice_line_gap(n) for n in (45, 90, 180))
    assert 1.7 < g1 / g2 < 2.3
    assert 1.7 < g2 / g3 < 2.3


def test_continuation_warm_start_agrees_with_a_cold_start_on_a_stable_branch():
    # The continuation's whole method is warm-starting, so it has to be shown that the warm start selects a
    # branch WITHOUT biasing the state on it: well inside the finite-cap branch, a swept point and an
    # independently cold-started march at the same sun must reach the same limit cycle.
    cfg = sici.SICIConfig(n_cells=180, n_steps=90)
    swept = sici.continuation_sweep(cfg, [1400.0, 1380.0, 1360.0, 1340.0], seed=25.0, max_years=600)
    assert swept.all_converged
    m = cfg.model(1340.0)
    alone = m.march(coalbedo_fn=ice_coalbedo, T_init=25.0, tol=1e-7, max_years=600)
    assert alone.converged
    cap_alone = 90.0 - ice_edge_latitude(m.x, alone.T_ocean, kind="perennial")
    assert abs(swept.points[-1].perennial_cap_deg - cap_alone) < cfg.cap_resolution_deg


# --------------------------------------------------------------------------- #
# CONSERVATION — a swept point still closes the energy budget.
# --------------------------------------------------------------------------- #
def test_swept_point_conserves_global_annual_energy():
    # The sweep changes nothing about the physics — it only chooses seeds — so a converged point must still
    # balance net TOA in the global-and-annual mean with the co-albedo it ACTUALLY realized. A seed that
    # quietly corrupted the state (e.g. a broadcast/aliasing slip in the new tuple T_init branch) would
    # break this before it broke anything else.
    cfg = sici.SICIConfig(n_cells=120, n_steps=180)
    swept = sici.continuation_sweep(cfg, [1380.0, 1370.0], seed=25.0, tol=1e-8, max_years=600)
    assert swept.all_converged
    m = cfg.model(1370.0)
    c = m.march(coalbedo_fn=ice_coalbedo, T_init=swept.final_state, tol=1e-8, max_years=600)
    S, xc = m.insolation_series(), m.x[:, None]
    aL = S * ice_coalbedo(xc, c.T_land)
    aO = S * ice_coalbedo(xc, c.T_ocean)
    net = m.f_land * (aL - m.A - m.B * c.T_land) + m.f_ocean * (aO - m.A - m.B * c.T_ocean)
    assert abs(float(net.mean())) < 5e-3


# --------------------------------------------------------------------------- #
# UNIT — the two new numerical objects, pinned on synthetic input.
# --------------------------------------------------------------------------- #
def _fake_continuation(cfg, S0, caps, direction):
    pts = tuple(sici.SweepPoint(S0=float(s), perennial_cap_deg=float(c), seasonal_cap_deg=float(c) + 5.0,
                                n_perennial_cells=int(c // 4), global_mean_T=10.0, polar_amplitude_K=3.0,
                                converged=True, years=10)
                for s, c in zip(S0, caps))
    return sici.Continuation(cfg, pts, direction, (np.zeros(1), np.zeros(1)))


def test_loop_width_detector_on_synthetic_branches():
    # The loop width is the headline metric, so its definition is pinned rather than trusted: the span of
    # samples whose two legs disagree by more than half a polar cell, PLUS one sampling interval (a single
    # disagreeing sample already means the branches part somewhere inside the neighbouring gaps).
    cfg = sici.SICIConfig(n_cells=720)
    S0 = np.arange(1360.0, 1372.1, 2.0)                            # 7 samples, dS0 = 2
    thr = cfg.cap_resolution_deg
    same = [10.0, 8.0, 6.0, 4.0, 2.0, 0.0, 0.0]
    loop = sici.HysteresisLoop(_fake_continuation(cfg, S0[::-1], same[::-1], "down"),
                               _fake_continuation(cfg, S0, same, "up"), thr)
    assert not loop.detected and loop.width == 0.0 and loop.interval is None
    assert loop.dS0 == pytest.approx(2.0)

    # Two adjacent samples disagree by well over the threshold → a loop spanning them plus one interval.
    differ = list(same)
    differ[4] = 2.0 + 3.0 * thr
    differ[5] = 0.0 + 3.0 * thr
    loop2 = sici.HysteresisLoop(_fake_continuation(cfg, S0[::-1], same[::-1], "down"),
                                _fake_continuation(cfg, S0, differ, "up"), thr)
    assert loop2.detected
    assert loop2.interval == (pytest.approx(1368.0), pytest.approx(1370.0))
    assert loop2.width == pytest.approx(4.0)                       # 2 W/m² span + one 2 W/m² sample

    # A disagreement UNDER the threshold is noise, not a loop — the guard that keeps "no loop" and "a loop
    # below detection" from silently merging.
    small = list(same)
    small[4] = 2.0 + 0.4 * thr
    loop3 = sici.HysteresisLoop(_fake_continuation(cfg, S0[::-1], same[::-1], "down"),
                                _fake_continuation(cfg, S0, small, "up"), thr)
    assert not loop3.detected and loop3.width == 0.0


def test_max_cell_jump_reads_the_uninterpolated_cap():
    cfg = sici.SICIConfig(n_cells=720)
    S0 = np.arange(1360.0, 1366.1, 2.0)
    one_at_a_time = _fake_continuation(cfg, S0, [0.0, 4.0, 8.0, 12.0], "up")
    assert one_at_a_time.max_cell_jump == 1                        # cells = cap//4 → 0,1,2,3
    jumpy = _fake_continuation(cfg, S0, [0.0, 0.0, 16.0, 20.0], "up")
    assert jumpy.max_cell_jump == 4                                # a fold the radius would have smoothed


def test_plant_cap_cools_exactly_the_polar_cells_in_both_hemispheres():
    # The planted-cap seed is the second observable's whole mechanism: it must ice the cells poleward of the
    # requested radius in BOTH hemispheres and leave every other cell bit-untouched (an off-by-one here
    # would plant a cap of the wrong size and quietly change the answer).
    cfg = sici.SICIConfig(n_cells=180, n_steps=90)
    m = cfg.model()
    base = (np.full(m.n_cells, 12.0), np.full(m.n_cells, 8.0))
    cap_deg = 15.0
    TL, TO = sici.plant_cap(m, base, cap_deg, margin_K=5.0)
    cold = np.abs(m.x) >= math.sin(math.radians(90.0 - cap_deg))
    assert cold.sum() > 0 and cold.sum() < m.n_cells
    assert np.array_equal(cold, np.abs(np.degrees(np.arcsin(np.clip(m.x, -1, 1)))) >= 90.0 - cap_deg - 1e-9)
    assert np.all(TL[cold] == T_FREEZE - 5.0) and np.all(TO[cold] == T_FREEZE - 5.0)
    assert np.all(TL[~cold] == 12.0) and np.all(TO[~cold] == 8.0)  # untouched elsewhere
    assert np.all(base[0] == 12.0)                                 # and the input is not mutated


def test_config_rejects_an_odd_grid_and_an_unknown_tile():
    # An odd cell count puts no face at the equator, so the mirrored cap mask would be asymmetric.
    with pytest.raises(ValueError, match="even"):
        sici.SICIConfig(n_cells=181)
    with pytest.raises(ValueError, match="tile"):
        sici.SICIConfig(tile="atmosphere")


def test_polar_cell_quantum_matches_the_grid():
    # Every cap radius is quoted against this, so the arithmetic is pinned: the polar cell spans
    # 90° − asin(1 − Δx), and it shrinks only as √Δx (quadrupling the cells halves it).
    for n in (360, 720, 1440):
        cfg = sici.SICIConfig(n_cells=n)
        assert cfg.polar_cell_deg == pytest.approx(math.degrees(math.acos(1.0 - 2.0 / n)))
        assert cfg.cap_resolution_deg == pytest.approx(0.5 * cfg.polar_cell_deg)
    assert sici.SICIConfig(n_cells=1440).polar_cell_deg == pytest.approx(
        0.5 * sici.SICIConfig(n_cells=360).polar_cell_deg, rel=0.02)


# --------------------------------------------------------------------------- #
# LOOSE (the payoff) — what the seasonal cycle does to the fold. Slow: these march.
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_earth_mixed_layer_dissolves_the_fold_cap_grows_cell_by_cell():
    # THE HEADLINE. The annual-mean parent of this very model says a polar cap smaller than θ_c cannot be
    # held by any sun: brightening, the cap must vanish in a jump, leaving a hysteresis loop of ~9 W/m².
    # March the same model with the seasons on, at Earth's tilt and a 50 m mixed layer, and none of that
    # survives:
    #   (a) the up- and down-sweeps retrace each other — no loop clears half a polar cell;
    #   (b) the perennial cap grows ONE CELL AT A TIME, which is the un-interpolated read (the cap radius
    #       alone could look smooth over a stepping state — this cannot); and
    #   (c) stable caps exist at radii the annual-mean model forbids, i.e. strictly below θ_c.
    # Direction banked; the numbers ride the calibrated heat capacities.
    #
    # RESOLUTION IS PART OF THE CLAIM, not a cost knob. Claim (b) only means something on a grid where a
    # θ_c-sized cap spans SEVERAL cells: at 360 cells it spans exactly ONE, so "the count grew by one cell"
    # and "a fold switched the whole cap on at once" are literally the same observation and the assertion
    # would be vacuous. Hence 720 (θ_c ≈ 5 cells) — and the guard below refuses the claim on any grid too
    # coarse to separate the two, so this cannot silently rot back into meaninglessness.
    cfg = sici.SICIConfig(n_cells=720, n_steps=180)
    parent = sici.annual_mean_curve(cfg)
    fold = parent.small_ice_cap_fold
    assert fold is not None
    parent_loop = fold.S0 - parent.ice_free_threshold_S0
    assert parent_loop > 5.0                                       # the parent's fold is unmistakable

    loop = sici.hysteresis_loop(cfg, np.arange(1364.0, 1382.1, 1.0))
    assert loop.all_converged                                      # never read a fold off a drifting march
    assert not loop.detected                                       # (a) the legs retrace each other
    assert loop.width == 0.0
    assert np.max(loop.cap_gap_deg) < cfg.cap_resolution_deg
    fold_cells = cfg.cells_in_cap(fold.cap_radius_deg)
    assert fold_cells >= 4, "grid too coarse for claim (b) to mean anything — see the comment above"
    assert loop.down.max_cell_jump == 1 and loop.up.max_cell_jump == 1     # (b) one cell at a time …
    assert loop.down.max_cell_jump * 4 <= fold_cells                       # … vs the whole cap a fold flips

    caps = loop.down.perennial_cap_deg
    forbidden = caps[(caps > 0.0) & (caps < fold.cap_radius_deg)]
    assert forbidden.size >= 2                                     # (c) caps the annual mean cannot hold
    assert caps.max() > fold.cap_radius_deg                        # and the sweep spans θ_c, not stops short


@pytest.mark.slow
def test_deepening_the_mixed_layer_brings_the_bistability_back():
    # The reduction on the payoff side, and the control that (a) above is about SEASONS and not about this
    # module's machinery: deepen the mixed layer and the seasonal swing dies as 1/(ωC), so the marcher walks
    # back toward its own annual-mean parent — and the parent's bistability reappears. Planting a cap of the
    # parent's critical size at the parent's fold sun: at 50 m it just relaxes to whatever a warm start
    # found (ONE climate), at depth it survives where the warm start grows no cap at all (TWO climates).
    # The approach is asymptotic, never attained — the swing falls as 1/(ωC) while the spin-up cost grows
    # linearly in C — so what is asserted is the ORDERING and the sign change, not the parent's numbers.
    #
    # 720 cells for the same reason the test above uses it: on a grid where a θ_c cap spans one cell, "a
    # planted cap survived" would mean one frozen cell versus none, and the surviving cap would be narrower
    # than the polar cell it sits in. The guard below keeps that from creeping back.
    cfg = sici.SICIConfig(n_cells=720, n_steps=180)
    parent = sici.annual_mean_curve(cfg)
    fold = parent.small_ice_cap_fold
    assert cfg.cells_in_cap(fold.cap_radius_deg) >= 4     # the grid must resolve a critical-sized cap
    gaps = {}
    for depth in (50.0, 200.0, 800.0):
        sd = sici.seed_dependence(cfg, fold.S0, curve=parent, ocean_mixed_depth=depth, max_years=6000)
        assert sd.converged
        gaps[depth] = sd
    assert not gaps[50.0].bistable                                 # seasonal Earth: one climate at this sun
    assert gaps[50.0].gap_deg < cfg.cap_resolution_deg
    assert gaps[800.0].bistable                                    # near-seasonless: two climates, the fold
    assert gaps[800.0].warm_cap_deg == 0.0                         # the warm start finds NO cap here …
    assert gaps[800.0].survived_cap_deg > cfg.polar_cell_deg       # … yet a planted one survives, and it is
    #                                                                wider than a single polar cell, so the
    #                                                                surviving climate is resolved, not a
    #                                                                one-cell artifact of the quantization
    assert gaps[50.0].gap_deg < gaps[200.0].gap_deg < gaps[800.0].gap_deg   # monotone in depth


@pytest.mark.slow
def test_the_loop_detector_reads_a_real_loop_when_one_exists():
    # THE INSTRUMENT'S POSITIVE CONTROL, and the reason the "loop width = 0" headline can be believed.
    # Every other marched call to `hysteresis_loop` in this module returns zero, so on its own that result
    # rests on a detector only ever shown to read non-zero on hand-built arrays (the unit test above). Here
    # the SAME detector, the SAME sweep machinery, is pointed at a configuration that genuinely has a fold:
    # an 800 m mixed layer, where the seasonal swing is damped to ~0.3 K and the marcher sits close to its
    # own annual-mean parent. It must find the loop — a wide band of suns where the dimming leg has no
    # perennial cap at all while the brightening leg still carries one — and its width must land near the
    # parent's own 8.9 W/m², not at some unrelated scale.
    #
    # This is the asymmetry worth stating: a detector that always reads zero is indistinguishable from a
    # broken one. Only this test separates "the seasons dissolved the fold" from "the instrument is blind".
    cfg = sici.SICIConfig(n_cells=360, n_steps=90)
    parent = sici.annual_mean_curve(cfg)
    fold = parent.small_ice_cap_fold
    parent_loop = fold.S0 - parent.ice_free_threshold_S0

    deep = sici.hysteresis_loop(cfg, np.arange(1368.0, 1378.1, 2.5),
                                ocean_mixed_depth=800.0, max_years=8000)
    assert deep.all_converged
    assert deep.detected                                          # the detector fires …
    assert deep.width > 0.5 * parent_loop                         # … at the parent's scale, not some other
    assert deep.width < 2.0 * parent_loop
    lo, hi = deep.interval
    assert parent.ice_free_threshold_S0 - 5.0 < lo and hi < fold.S0 + 5.0   # and in the right band of suns

    # The shape of a fold, not just its width: somewhere inside the loop the dimming leg carries NO
    # year-round ice while the brightening leg still does — two branches at one sun, traced by sweeping.
    inside = (deep.S0 >= lo) & (deep.S0 <= hi)
    caps_down = dict(zip(np.round(deep.down.S0, 9), deep.down.perennial_cap_deg))
    down_inside = np.array([caps_down[k] for k in np.round(deep.S0[inside], 9)])
    up_inside = np.array([p.perennial_cap_deg for p in deep.up.points])[inside]
    assert np.any((down_inside == 0.0) & (up_inside > cfg.polar_cell_deg))
