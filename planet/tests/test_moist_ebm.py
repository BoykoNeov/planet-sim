"""Triad for the rung-2.5 MSE-diffusing moist EBM (:mod:`planet.moist_ebm`).

Re-classed for honesty (plan §10 / the module docstring):

* **tight** — ``β(T)`` is the *exact* analytic Clausius–Clapeyron derivative; the transport coefficient
  ``D_eff`` sits **inside** the conservative divergence (a varying-``D_eff`` transport step conserves
  ``∫T dx``); the global-mean warming is **pinned** ``⟨δT⟩ = ΔA/B`` to machine precision; and the
  **frozen-D_eff null** warms *exactly uniformly* (PA = 1) — proving the emergent PA is the
  ``dD_eff/dT`` feedback, not the recalibrated ``D``-shape.
* **real-but-loose (the unlock)** — the **polar amplification**: the dt-free moist EBM warms the poles
  ~1.8–2.05× the tropics while the dry EBM warms uniformly; direction banked, magnitude loose
  (RH-dependent). The headline is the **dt-free** ``moist_steady_direct``; the Strang relaxation carries an
  O(Δt) shape bias that suppresses it to ~1.5 at the default step (one test pins that).
* **plumbing (by construction)** — RH = 0 **and** ``D_s = D_TRANSPORT`` reduce **both** paths to the
  genuine ``EnergyBalanceModel`` rung-0 solve **bit-for-bit**: the relaxation to ``equilibrate``, the
  dt-free direct solve to ``steady_linear``.
* **named modeling choices** — the recalibration matches the present contrast (``D_s < 0.555``), and the
  PA factor is **invariant** to that target (contrast vs the P₂ amplitude).
"""
import numpy as np
import pytest

from engines.diffusion import Diffusion1D, Neumann, uniform_grid
from planet import moist, moist_ebm as me
from planet.albedo import EBMParams
from planet.ebm import D_TRANSPORT, EnergyBalanceModel, legendre_P2


@pytest.fixture(scope="module")
def params():
    return EBMParams()


@pytest.fixture(scope="module")
def pa(params):
    # The headline experiment once (recalibrates D_s, warms dry + moist by ΔA=10) — reused across tests.
    return me.polar_amplification(params, RH=moist.RH_DEFAULT, dA=10.0)


# --------------------------------------------------------------------------- #
# TIGHT — β is the exact C–C derivative; D_eff is inside the conservative operator.
# --------------------------------------------------------------------------- #
def test_dqsat_dT_is_the_exact_analytic_cc_derivative():
    # The analytic derivative must match a centered finite difference of the rung-2 q_sat function …
    T = np.array([-20.0, 0.0, 15.0, 30.0])
    h = 1e-4
    fd = (moist.saturation_specific_humidity(T + h) - moist.saturation_specific_humidity(T - h)) / (2 * h)
    assert np.allclose(me.dqsat_dT(T), fd, rtol=1e-5)
    # … and its log-slope is the Clausius–Clapeyron ~7 %/K (distinct from the energy-constrained rate).
    rate = me.dqsat_dT(15.0) / moist.saturation_specific_humidity(15.0)
    assert 0.06 < rate < 0.075


def test_moisture_amplification_is_large_warm_small_cold(params):
    # β(T) = (L/c_p)·RH·dq_sat/dT: monotone increasing, large in the warm tropics, ≈ 0 at the cold pole.
    beta = me.moisture_amplification(np.array([-30.0, 0.0, 30.0]))
    assert beta[0] < beta[1] < beta[2]
    assert beta[2] > 1.0                       # warm equator strongly amplified
    assert beta[0] < 0.2                        # cold pole barely amplified
    # RH = 0 ⟹ no amplification (the dry limit)
    assert me.moisture_amplification(15.0, RH=0.0) == 0.0


def test_Deff_sits_inside_the_conservative_divergence(params):
    # The inside-form proof: a transport-only step with a strongly VARYING D_eff still conserves ∫T dx
    # (only ∂/∂x[(1−x²)·D_eff·∂T/∂x] — coefficient inside — conserves under insulated ends; the outside
    # form D_eff·∂/∂x[…] would not). This is exactly how moist_equilibrium builds its operator.
    grid = uniform_grid(1.0, params.n_cells)
    x = grid.centers
    C = me.RHO_WATER * me.CW_WATER * params.water_depth
    T = legendre_P2(x) + 0.3 * x ** 3                          # arbitrary non-uniform field
    D_eff = me.effective_diffusivity(20.0 - 50.0 * x, D_TRANSPORT, RH=0.8)   # varies ~3× across the grid
    assert D_eff[0] / D_eff[-1] > 2.0                          # genuinely varying
    Dcells = (D_eff / C) * (1.0 - x ** 2)
    stepped = Diffusion1D(grid, Dcells, Neumann(0.0), Neumann(0.0)).step(T, 0.3 * C / params.B)
    assert abs(float(np.mean(stepped)) - float(np.mean(T))) < 1e-12   # ∫T dx conserved (inside form)


def test_global_mean_warming_is_pinned_to_dA_over_B(params, pa):
    # TIGHT conservation: the transport conserves ∫T dx, so ⟨δT⟩ = ΔA/B for BOTH models to machine
    # precision — moisture REDISTRIBUTES a pinned global mean poleward, it does not add net warming.
    expected = pa.dA / params.B
    assert pa.mean_delta_T == pytest.approx(expected, abs=1e-9)
    assert float(np.mean(pa.delta_T_dry)) == pytest.approx(expected, abs=1e-9)
    assert float(np.mean(pa.delta_T_moist)) == pytest.approx(expected, abs=1e-9)


def test_frozen_Deff_warming_is_exactly_uniform(params, pa):
    # THE ATTRIBUTION NULL: freeze D_eff at the present moist profile (a T-INDEPENDENT array-D EBM — the
    # genuine rung-1 callable-D path) and warm by the same ΔA. Uniform δT solves the perturbation for any
    # frozen D(x), so PA must be exactly 1 — proving the moist PA is 100 % the dD_eff/dT feedback and
    # 0 % the recalibrated D-shape. (If this is not ~1, D_eff is leaking an update somewhere.)
    absorbed = me.constant_albedo_absorbed(params)
    D_frozen = me.effective_diffusivity(pa.moist_present.T, pa.D_s, pa.RH)
    Dcall = lambda xx: np.interp(np.asarray(xx, dtype=float), pa.moist_present.x, D_frozen)
    common = dict(B=params.B, D=Dcall, T_freeze=params.T_freeze,
                  water_depth=params.water_depth, n_cells=params.n_cells, face="harmonic")
    fz_present = EnergyBalanceModel(A=params.A, **common).equilibrate(absorbed, pa.moist_present.T)
    fz_warm = EnergyBalanceModel(A=params.A - pa.dA, **common).equilibrate(absorbed, fz_present.T)
    dT = fz_warm.T - fz_present.T
    assert dT[-1] / dT[0] == pytest.approx(1.0, abs=1e-3)
    assert float(np.ptp(dT)) < 1e-3


# --------------------------------------------------------------------------- #
# PLUMBING (by construction) — RH=0 AND D_s=D_TRANSPORT reduce to the genuine rung-0 engine solve.
# --------------------------------------------------------------------------- #
def test_RH0_default_D_reduces_to_rung0_engine_bit_for_bit(params):
    # The self-contained moist relaxation, with no moisture (RH=0) and the rung-0 D, executes the
    # IDENTICAL operations as the validated EnergyBalanceModel relaxation (β≡0 ⟹ the per-step operator
    # is the dry one every step; face="harmonic" makes the cells match). Bit-for-bit, not just close.
    absorbed = me.constant_albedo_absorbed(params)
    eng = EnergyBalanceModel(A=params.A, B=params.B, D=D_TRANSPORT, T_freeze=params.T_freeze,
                             water_depth=params.water_depth, n_cells=params.n_cells, face="harmonic")
    T_init = 30.0 + (-30.0 - 30.0) * eng.x
    ref = eng.equilibrate(absorbed, T_init)
    mine = me.moist_equilibrium(params, D_TRANSPORT, absorbed, RH=0.0, T_init=T_init)
    assert np.array_equal(mine.T, ref.T)


def test_direct_solve_RH0_reduces_to_steady_linear_to_machine_precision(params):
    # The DT-FREE companion reduction: at RH=0 (constant D_eff = D_s) the Picard direct solve IS the dry
    # LINEAR EBM, so it must reproduce EnergyBalanceModel.steady_linear — the headline path's plumbing
    # anchor, pinning the lean tridiagonal assembly to the engine's. (To machine precision, not bit-for-bit:
    # the engine forms the cell coefficient as (D/C)(1−x²)·C, a different float order than D·(1−x²).)
    absorbed = me.constant_albedo_absorbed(params)
    eng = EnergyBalanceModel(A=params.A, B=params.B, D=D_TRANSPORT, T_freeze=params.T_freeze,
                             water_depth=params.water_depth, n_cells=params.n_cells, face="harmonic")
    ref = eng.steady_linear(absorbed)
    mine = me.moist_steady_direct(params, D_TRANSPORT, absorbed, RH=0.0)
    assert np.allclose(mine.T, ref.T, rtol=0.0, atol=1e-9)
    assert mine.converged


# --------------------------------------------------------------------------- #
# REAL-BUT-LOOSE (the unlock) — emergent polar amplification; the dry model is the uniform null.
# --------------------------------------------------------------------------- #
def test_moist_polar_amplifies_while_dry_warms_uniformly(params, pa):
    # The headline: the moist EBM warms the pole MORE than the equator (PA > 1); the dry EBM warms
    # EXACTLY uniformly (the clean null — transport of a uniform field is zero).
    assert pa.moist_present.converged and pa.moist_warm.converged   # the relaxation reached a fixed point
    assert pa.dry_present.converged and pa.dry_warm.converged
    assert pa.delta_T_moist[-1] > pa.delta_T_moist[0]          # pole warms more than equator
    assert pa.pa_moist > 1.2                                   # genuine amplification
    assert pa.pa_dry == pytest.approx(1.0, abs=1e-3)           # dry null
    assert float(np.ptp(pa.delta_T_dry)) < 1e-3               # dry warming is uniform
    assert pa.polar_excess > 0.0                               # pole warms above the pinned mean


def test_polar_amplification_factor_in_loose_band_both_metrics(pa):
    # Loose benchmark (Earth defaults, RH 0.8), read off the DT-FREE direct-solve climates. NAME THE
    # METRIC: the single-endpoint ratio (the headline ~2.05, the most generous, polar cell on the
    # harmonic-face bias) and the area-band ratio mean(>=60°)/mean(<=30°) (~1.80, less generous) are BOTH
    # honest polar amplification — assert both in the loose band, with endpoint ≥ band (the endpoint reads
    # the warmest/coldest extremes). (Earlier these were ~1.5 / ~1.4 — the n_tau=0.5 Strang relaxation's
    # O(Δt) shape bias; test_relaxation_underestimates_dt_free_PA pins that artifact.)
    assert 1.9 < pa.pa_moist < 2.2                            # endpoint ratio ~2.05 (dt-free)
    assert 1.7 < pa.pa_moist_band < 1.9                       # area-band ratio ~1.80 (dt-free)
    assert pa.pa_moist >= pa.pa_moist_band                     # endpoint is the more generous metric
    assert pa.pa_dry_band == pytest.approx(1.0, abs=1e-3)      # dry null both ways (band metric too)


def test_relaxation_underestimates_dt_free_PA_and_converges_to_it(params):
    # THE SPLITTING-ARTIFACT PIN. The headline PA is the dt-FREE direct solve (~2.05 endpoint). The Strang
    # relaxation carries an O(Δt) SHAPE bias (backward-Euler transport split against the exact radiation
    # half-step) that SUPPRESSES the amplification: at the default n_tau=0.5 it reads ~1.5, and it climbs
    # toward the direct value as n_tau shrinks. The global mean stays exact at every step — only the shape,
    # hence the RATIO, is biased. Pins the finding (and the direction of the bias) so it cannot silently
    # regress to the old under-converged headline.
    absorbed = me.constant_albedo_absorbed(params)
    RH = moist.RH_DEFAULT
    D_s = me.recalibrate_sensible_D(params, RH)                # the dt-free recalibrated D_s (~0.28)
    A_warm = params.A - 10.0

    def endpoint_pa(present, warm):
        dT = warm.T - present.T
        return dT[-1] / dT[0], float(np.mean(dT))

    dp = me.moist_steady_direct(params, D_s, absorbed, RH=RH)
    dw = me.moist_steady_direct(params, D_s, absorbed, RH=RH, A=A_warm, T_init=dp.T)
    pa_direct, mean_direct = endpoint_pa(dp, dw)

    cp = me.moist_equilibrium(params, D_s, absorbed, RH=RH, n_tau=0.5)
    cw = me.moist_equilibrium(params, D_s, absorbed, RH=RH, A=A_warm, T_init=cp.T, n_tau=0.5)
    pa_coarse, mean_coarse = endpoint_pa(cp, cw)

    fp = me.moist_equilibrium(params, D_s, absorbed, RH=RH, n_tau=0.1, tol=1e-11, max_iter=200000)
    fw = me.moist_equilibrium(params, D_s, absorbed, RH=RH, A=A_warm, T_init=fp.T, n_tau=0.1,
                              tol=1e-11, max_iter=200000)
    pa_fine, _ = endpoint_pa(fp, fw)

    assert pa_direct > 1.9                                     # dt-free headline ~2.05
    assert pa_coarse < 1.65                                    # the n_tau=0.5 artifact ~1.5 (suppressed)
    assert pa_coarse < pa_fine < pa_direct                    # shrinking dt climbs toward the dt-free value
    # the global-mean warming is the pinned ΔA/B at BOTH steps — only the SHAPE drifts with dt
    assert mean_coarse == pytest.approx(mean_direct, abs=1e-6)
    assert mean_direct == pytest.approx(10.0 / params.B, abs=1e-9)


def test_pa_direction_is_robust_grows_with_RH(params):
    # Direction banked, magnitude loose: more moisture (higher RH) ⟹ stronger amplification, but PA > 1
    # at every RH. Each RH recalibrated to its own present contrast so the comparison is fair.
    pa_lo = me.polar_amplification(params, RH=0.6, dA=10.0)
    pa_hi = me.polar_amplification(params, RH=0.8, dA=10.0)
    assert pa_lo.moist_warm.converged and pa_hi.moist_warm.converged   # both relaxations converged
    assert pa_lo.pa_moist > 1.0 and pa_hi.pa_moist > 1.0
    assert pa_hi.pa_moist > pa_lo.pa_moist


# --------------------------------------------------------------------------- #
# NAMED MODELING CHOICES — the recalibration (the wall) and its target-invariance.
# --------------------------------------------------------------------------- #
def test_recalibration_matches_the_present_contrast(params, pa):
    # The moist present-day climate reproduces the dry present-day equator-pole CONTRAST (the
    # recalibration target) — that is what re-deriving D_s buys, so latent transport is not double-counted.
    assert me.equator_pole_contrast(pa.moist_present) == pytest.approx(
        me.equator_pole_contrast(pa.dry_present), abs=0.05)


def test_recalibrated_Ds_is_below_the_dry_default(pa):
    # D_s < rung-0's 0.555: because the moisture-amplified D_eff transports MORE heat, less dry D is
    # needed to match the same present contrast (the explicit-latent-transport correction = the wall).
    assert pa.D_s < D_TRANSPORT
    assert 0.2 < pa.D_s < 0.45                                 # ≈ 0.28 for Earth + RH 0.8 (dt-free)


def test_pa_is_invariant_to_the_recalibration_target(params):
    # NAMED choice: matching the endpoint contrast vs the P₂ amplitude T₂ moves the PA factor < 5 %,
    # because PA is set by the SHAPE of β(T), not by the overall D_s scaling — so the target is a
    # modeling choice, not a tuning of the result.
    D_contrast = me.recalibrate_sensible_D(params, RH=0.8, target="contrast")
    D_T2 = me.recalibrate_sensible_D(params, RH=0.8, target="T2")
    pa_contrast = me.polar_amplification(params, RH=0.8, D_s=D_contrast)
    pa_T2 = me.polar_amplification(params, RH=0.8, D_s=D_T2)
    assert pa_T2.pa_moist == pytest.approx(pa_contrast.pa_moist, rel=0.05)
