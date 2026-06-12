"""Triad for the rung-2.5 MSE-diffusing moist EBM (:mod:`planet.moist_ebm`).

Re-classed for honesty (plan §10 / the module docstring):

* **tight** — ``β(T)`` is the *exact* analytic Clausius–Clapeyron derivative; the transport coefficient
  ``D_eff`` sits **inside** the conservative divergence (a varying-``D_eff`` transport step conserves
  ``∫T dx``); the global-mean warming is **pinned** ``⟨δT⟩ = ΔA/B`` to machine precision; and the
  **frozen-D_eff null** warms *exactly uniformly* (PA = 1) — proving the emergent PA is the
  ``dD_eff/dT`` feedback, not the recalibrated ``D``-shape.
* **real-but-loose (the unlock)** — the **polar amplification**: the moist EBM warms the poles ~1.4–1.5×
  the tropics while the dry EBM warms uniformly; direction banked, magnitude loose (RH-dependent).
* **plumbing (by construction)** — RH = 0 **and** ``D_s = D_TRANSPORT`` reduce the moist relaxation to
  the genuine ``EnergyBalanceModel`` rung-0 solve **bit-for-bit**.
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


# --------------------------------------------------------------------------- #
# REAL-BUT-LOOSE (the unlock) — emergent polar amplification; the dry model is the uniform null.
# --------------------------------------------------------------------------- #
def test_moist_polar_amplifies_while_dry_warms_uniformly(params, pa):
    # The headline: the moist EBM warms the pole MORE than the equator (PA > 1); the dry EBM warms
    # EXACTLY uniformly (the clean null — transport of a uniform field is zero).
    assert pa.delta_T_moist[-1] > pa.delta_T_moist[0]          # pole warms more than equator
    assert pa.pa_moist > 1.2                                   # genuine amplification
    assert pa.pa_dry == pytest.approx(1.0, abs=1e-3)           # dry null
    assert float(np.ptp(pa.delta_T_dry)) < 1e-3               # dry warming is uniform
    assert pa.polar_excess > 0.0                               # pole warms above the pinned mean


def test_polar_amplification_factor_in_loose_band(pa):
    # Loose benchmark (Earth defaults, RH 0.8): poles warm ~1.5× the tropics from moisture alone. The
    # magnitude is loose — the observed ~2–3× also needs ice-albedo + lapse-rate feedbacks (out of scope).
    assert 1.3 < pa.pa_moist < 1.7


def test_pa_direction_is_robust_grows_with_RH(params):
    # Direction banked, magnitude loose: more moisture (higher RH) ⟹ stronger amplification, but PA > 1
    # at every RH. Each RH recalibrated to its own present contrast so the comparison is fair.
    pa_lo = me.polar_amplification(params, RH=0.6, dA=10.0)
    pa_hi = me.polar_amplification(params, RH=0.8, dA=10.0)
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
    assert 0.2 < pa.D_s < 0.45                                 # ≈ 0.30 for Earth + RH 0.8


def test_pa_is_invariant_to_the_recalibration_target(params):
    # NAMED choice: matching the endpoint contrast vs the P₂ amplitude T₂ moves the PA factor < 5 %,
    # because PA is set by the SHAPE of β(T), not by the overall D_s scaling — so the target is a
    # modeling choice, not a tuning of the result.
    D_contrast = me.recalibrate_sensible_D(params, RH=0.8, target="contrast")
    D_T2 = me.recalibrate_sensible_D(params, RH=0.8, target="T2")
    pa_contrast = me.polar_amplification(params, RH=0.8, D_s=D_contrast)
    pa_T2 = me.polar_amplification(params, RH=0.8, D_s=D_T2)
    assert pa_T2.pa_moist == pytest.approx(pa_contrast.pa_moist, rel=0.05)
