"""Rung 5A validation: linear orographic precipitation, the first step off the zonal mean (plan §12.5).

The rain-shadow half of the "north star" (regional climate from geography). What is asserted **tight**
is *exact* and *structural*: the model converges to the **closed-form triangle-ridge solution** (which
pins the *reduced* transfer function — the ``C_w`` scaling, the upslope gradient ``iσ``, the fallout
``τ_f``; in that limit ``H_w`` and the vertical-wavenumber **branch** drop out), recovers the
**upslope limit** to the analytic terrain gradient, and vanishes over flat ground. The ``sgn(σ)``
branch is guarded *only* by the **rain-shadow direction** test (windward wetter than lee, peak upwind);
the wind-reversal test is a reflection self-consistency check, not a branch test. What is **loose** is
the absolute mm/hr magnitude — set by the cited Smith & Barstad constants
([[smith-barstad-orographic-source]]).

Scope (the honesty flag): a *diagnostic* on a *prescribed* uniform wind over a *regional Cartesian
patch* — it makes the *precipitation* 2-D, not the engine. Sphere placement, jet wiring and
serialization are 5A.2 (deferred).
"""
import numpy as np
import pytest

from planet import orographic as og


def _ridge_field(x: np.ndarray, h1: np.ndarray) -> np.ndarray:
    """Tile a 1-D cross-ridge profile ``h1(x)`` across latitude → a y-invariant 2-D terrain ``[y, x]``."""
    return np.tile(h1, (len(x), 1))


# --------------------------------------------------------------------------- #
# Pinned constants ([[smith-barstad-orographic-source]]) — cited Smith & Barstad (2004) values
# --------------------------------------------------------------------------- #
def test_pinned_constants_are_cited_values():
    assert og.TAU_C_S == pytest.approx(1000.0)                 # conversion time
    assert og.TAU_F_S == pytest.approx(1000.0)                 # fallout time
    assert og.NM_PER_S == pytest.approx(0.005)                 # moist stability frequency
    assert og.HW_M == pytest.approx(2500.0)                    # water-vapour scale height
    assert og.U_REF_M_S == pytest.approx(15.0)                 # reference wind speed
    # C_w = ρ_Sref·Θ_m/γ, both lapse rates negative → C_w > 0
    assert og.CW_KG_M3 == pytest.approx(7.4e-3 * -6.5 / -5.8)
    assert og.CW_KG_M3 > 0.0


def test_westerly_wind_blows_eastward():
    # meteorological convention: direction is where the wind comes FROM; 270° (west) → +x (eastward)
    u, v = og.wind_components(15.0, 270.0)
    assert u == pytest.approx(15.0) and v == pytest.approx(0.0, abs=1e-12)


# --------------------------------------------------------------------------- #
# Analytical (tight, *exact*): convergence to the closed-form triangle-ridge solution
# --------------------------------------------------------------------------- #
def _triangle_run(dx: float):
    """FFT model vs the exact solution on a triangle ridge, in the reduced limit H_w = τ_c = 0, f = 0."""
    x, _ = og.make_grid(half_width_m=100e3, dx=dx)
    oro = _ridge_field(x, og.triangle_ridge(x))
    P = og.orographic_precip(oro, dx, dx, speed=15.0, direction_deg=270.0,
                             tau_c=0.0, Hw=0.0, latitude_deg=0.0, truncate=True)
    p_model = P[len(x) // 2, :]
    p_exact = og.triangle_ridge_exact(x, u=15.0, Cw=og.CW_KG_M3, tau=og.TAU_F_S)
    return np.max(np.abs(p_model - p_exact)), float(p_exact.max())


def test_triangle_ridge_matches_exact_solution():
    # the tight anchor: at 1 km spacing the FFT model matches the analytic rain shadow to < 0.5 % of peak
    err, peak = _triangle_run(dx=1e3)
    assert err / peak < 0.005
    # the crest value is the upslope rate C = C_w·U·A/d reduced by the finite fallout advection Uτ:
    # P(crest) = C·(1 − e^{−d/Uτ}), the windward-branch limit at x → 0⁻.
    C = og.CW_KG_M3 * 15.0 * 500.0 / 50e3 * og.SECONDS_PER_HOUR
    Ut = 15.0 * og.TAU_F_S
    assert peak == pytest.approx(C * (1.0 - np.exp(-50e3 / Ut)), rel=0.01)


def test_triangle_ridge_error_converges_with_resolution():
    # refining the grid drives the error toward zero (validates it's a discretization error, not a bug)
    err_coarse, peak = _triangle_run(dx=5e3)
    err_fine, _ = _triangle_run(dx=1e3)
    assert err_fine < 0.3 * err_coarse                         # ~O(dx²): 5× refinement cuts error ≫ 3×
    assert err_fine / peak < 0.005


# --------------------------------------------------------------------------- #
# Structural (tight): the upslope limit
# --------------------------------------------------------------------------- #
def test_upslope_limit_recovers_Cw_U_grad_h():
    # H_w = τ_c = τ_f = 0 collapses the transfer function to the classic upslope model
    # P = C_w · max(0, U·∂h/∂x). Compare to the *analytic* Gaussian gradient ∂h/∂x = −x/σ² · h.
    dx = 1e3
    x, _ = og.make_grid(half_width_m=200e3, dx=dx)
    A, sig = 1000.0, 15e3
    h1 = og.gaussian_ridge(x, amplitude_m=A, sigma_m=sig)
    P = og.orographic_precip(_ridge_field(x, h1), dx, dx, speed=15.0, direction_deg=270.0,
                             tau_c=0.0, tau_f=0.0, Hw=0.0, latitude_deg=0.0, truncate=True)
    p_model = P[len(x) // 2, :]
    dhdx = -x / sig**2 * h1                                    # analytic gradient of the Gaussian
    upslope = np.maximum(og.CW_KG_M3 * 15.0 * dhdx, 0.0) * og.SECONDS_PER_HOUR
    assert np.max(np.abs(p_model - upslope)) / upslope.max() < 0.01


# --------------------------------------------------------------------------- #
# Structural (tight): the rain shadow — the SOLE guard on the sgn(σ) vertical-wavenumber branch.
# In the triangle-ridge anchor above (H_w = 0) the branch drops out of the transfer function, so it is
# *not* validated there; nor by the wind-reversal mirror (invariant under a branch flip). This test is
# the only thing that reddens if the branch sign is wrong (the shadow then flips to the windward side).
# --------------------------------------------------------------------------- #
def test_rain_shadow_windward_wet_lee_dry_is_the_branch_anchor():
    dx = 2e3
    x, _ = og.make_grid(half_width_m=200e3, dx=dx)
    h1 = og.gaussian_ridge(x)                                  # crest at x = 0, westerly wind → +x
    P = og.orographic_precip(_ridge_field(x, h1), dx, dx, speed=15.0, direction_deg=270.0)
    p = P[len(x) // 2, :]
    crest = len(x) // 2
    windward = p[:crest].sum()                                 # upwind of the crest (x < 0)
    lee = p[crest:].sum()                                      # downwind (x > 0)
    # load-bearing branch assertions (a flipped branch drives windward ≪ lee and the peak into the lee):
    assert lee < 0.5 * windward                                # the lee is materially drier — not just <
    assert x[int(np.argmax(p))] <= -4e3                        # the peak sits several km *upwind* of the crest


def test_wind_reversal_is_reflection_symmetry_not_a_branch_test():
    # reversing the wind (westerly 270° → easterly 90°) must mirror the precip in x. This is a reflection
    # self-consistency check (catches k-grid / fftfreq asymmetry); it holds for *either* branch sign, so
    # it does NOT discriminate the sgn(σ) branch — that is the rain-shadow test's job.
    dx = 2e3
    x, _ = og.make_grid(half_width_m=200e3, dx=dx)
    oro = _ridge_field(x, og.gaussian_ridge(x))
    p_west = og.orographic_precip(oro, dx, dx, speed=15.0, direction_deg=270.0, latitude_deg=0.0)[len(x) // 2, :]
    p_east = og.orographic_precip(oro, dx, dx, speed=15.0, direction_deg=90.0, latitude_deg=0.0)[len(x) // 2, :]
    assert np.max(np.abs(p_west - p_east[::-1])) < 1e-9        # exact mirror (machine precision)


# --------------------------------------------------------------------------- #
# Structural (tight): degenerate cases — flat ground, non-negativity, no NaNs
# --------------------------------------------------------------------------- #
def test_flat_terrain_gives_no_orographic_anomaly():
    P = og.orographic_precip(np.zeros((40, 40)), 2e3, 2e3)
    assert np.max(np.abs(P)) == 0.0


def test_background_precip_is_added_uniformly():
    P = og.orographic_precip(np.zeros((20, 20)), 2e3, 2e3, background_mm_hr=1.5)
    assert np.allclose(P, 1.5)


def test_precip_is_nonnegative_and_finite():
    dx = 2e3
    x, _ = og.make_grid(half_width_m=200e3, dx=dx)
    oro = _ridge_field(x, og.gaussian_ridge(x, amplitude_m=2000.0, sigma_m=8e3))  # steep/narrow → evanescent modes
    P = og.orographic_precip(oro, dx, dx, speed=15.0, direction_deg=270.0)
    assert np.all(np.isfinite(P))                              # σ=0 locus and evanescent modes don't NaN
    assert np.all(P >= 0.0)                                    # truncation holds
