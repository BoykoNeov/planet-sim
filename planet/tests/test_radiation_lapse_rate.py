"""Triad for the rung-4 lapse-rate-feedback slice (moist adiabat → emergent λ_LR, :mod:`planet.radiation`).

The within-rung upgrade that makes the lapse-rate feedback **emergent**: swap the fixed convective ``Γ``
for a temperature-dependent moist adiabat that flattens as it warms, so surface warming amplifies in the
upper troposphere and ``OLR(Ts)`` steepens. The §12 scoping guess ("supplies ``λ_LR ≈ 0.84``, closing the
gap to climlab's 2") was **OVERTURNED**: the emergent value is ``≈ +1.5`` and the column *overshoots*.

* **tight** — ``moist_adiabat=False`` is bit-for-bit the rung-4-core column (the reduction); the moist
  adiabat is **derived by its limits** (dry ``→ g/c_p ≈ 9.8 K/km``, flattens when warm); the Soden & Held
  kernel split **closes** (Planck + LR + WV = B_total to ~1e-3) and is resolution-converged.
* **discriminator (sign & kind, banked)** — the moist adiabat warms the *upper troposphere more than the
  surface* (``ΔT_aloft/ΔTs > 1``), which fixed ``Γ`` cannot do; the kernel lapse-rate term is ``> 0``.
* **unlock (loose)** — the magnitude is ``≈ 1.5``: the *right order*, but it **overshoots** the global-mean
  Soden & Held ``λ_LR = 0.84`` (single global column = the tropical branch only) and pushes the with-WV
  ``B`` *above* climlab's 2. Banked as order + sign, not as a matched magnitude.
"""
import numpy as np
import pytest

from planet import radiation as rad
from planet.ebm import B_OLR


@pytest.fixture(scope="module")
def fixed_column():
    return rad.calibrate_column()                       # the rung-4-core column (fixed Γ, default)


@pytest.fixture(scope="module")
def moist_column():
    return rad.calibrate_column(moist_adiabat=True)     # the opt-in moist-adiabat column


# --------------------------------------------------------------------------- #
# TIGHT — the moist adiabat is derived (by its limits), and the off-switch reduces bit-for-bit.
# --------------------------------------------------------------------------- #
def test_moist_adiabat_recovers_the_dry_adiabat_in_the_dry_limit():
    """As r_s → 0 (very cold / dry) the moist adiabat → the dry adiabat g/c_p ≈ 9.8 K/km."""
    dry = rad.GRAVITY / rad.C_P_AIR * 1000.0           # K/km
    cold = rad.moist_adiabatic_lapse_rate(180.0, 5.0e4) * 1000.0
    assert cold == pytest.approx(dry, rel=0.01)        # within 1% of the dry adiabat
    assert dry == pytest.approx(9.77, abs=0.1)


def test_moist_adiabat_flattens_when_warm():
    """The lapse rate decreases with temperature (latent-heat release) — the whole mechanism.

    A warm, moist surface (~300 K) sits near the ~4 K/km tropical value, well below the dry adiabat.
    """
    g = [rad.moist_adiabatic_lapse_rate(T, 1.0e5) * 1000.0 for T in (180.0, 250.0, 300.0)]
    assert g[0] > g[1] > g[2]                           # monotone: warmer ⇒ flatter
    assert 3.0 < g[2] < 5.0                             # tropical surface ≈ 4 K/km


def test_moist_adiabat_off_is_bit_for_bit_the_fixed_gamma_column(fixed_column):
    """``moist_adiabat=False`` (default) is identical to a hand-built fixed-Γ profile — the reduction.

    Guards that adding the flag and factoring out ``_olr_from`` did not perturb the rung-4-core column.
    """
    col = fixed_column
    for Ts in (270.0, 288.0, 300.0):
        p = np.linspace(0.02 * rad.P_SURFACE, rad.P_SURFACE, col.n_levels)
        sh = rad.R_DRY_AIR * rad.emission_temperature() / rad.GRAVITY
        z = sh * np.log(rad.P_SURFACE / p)
        T = np.maximum(Ts - col.lapse_rate * z, col.strat_T)
        tau = col._optical_depth(Ts, 1.0, True)
        assert col.outgoing_longwave(Ts) == col._olr_from(T, tau)        # exact, same computation


def test_kernel_closure_and_convergence(moist_column):
    """Planck + LR + WV = B_total to ~1e-3 (first order), and λ_LR is resolution-converged."""
    k = moist_column.feedback_kernel()
    assert abs(k.closure_residual) < 3e-3                               # the kernels sum to the slope
    lr = [rad.calibrate_column(moist_adiabat=True, n_levels=n).feedback_kernel().lapse_rate
          for n in (100, 400)]
    assert abs(lr[0] - lr[1]) < 1e-2                                    # converged (not a discretization artifact)


# --------------------------------------------------------------------------- #
# DISCRIMINATOR — the sign and kind (upper-troposphere amplification) are banked.
# --------------------------------------------------------------------------- #
def test_moist_adiabat_amplifies_upper_tropospheric_warming(fixed_column, moist_column):
    """Warming Ts warms the upper troposphere MORE than the surface for a moist adiabat; fixed Γ shifts uniformly.

    This is the sign of the lapse-rate feedback, measured (not assumed): the emission level sits aloft, so
    amplified warming aloft steepens OLR(Ts).
    """
    dT = 4.0
    p, _ = moist_column._profile(rad.PRESENT_SURFACE_T)
    level = int(np.argmin(np.abs(p - 0.4 * rad.P_SURFACE)))            # a mid-upper-tropospheric level

    def amp(col):
        _, hi = col._profile(rad.PRESENT_SURFACE_T + dT)
        _, lo = col._profile(rad.PRESENT_SURFACE_T - dT)
        return (hi[level] - lo[level]) / (2.0 * dT)

    assert amp(fixed_column) == pytest.approx(1.0, abs=0.02)           # uniform shift (above the strat floor)
    assert amp(moist_column) > 1.10                                    # upper-trop amplification


def test_lapse_rate_term_is_a_positive_feedback_on_B(fixed_column, moist_column):
    """The kernel lapse-rate term is clearly positive for the moist adiabat, ~null for fixed Γ.

    Fixed Γ has only a small tropopause-migration residual (the column's null is not perfectly clean);
    the moist adiabat adds a large positive lapse-rate term to B.
    """
    assert abs(fixed_column.feedback_kernel().lapse_rate) < 0.4        # the null (small residual only)
    assert moist_column.feedback_kernel().lapse_rate > 1.0             # a real, large feedback


# --------------------------------------------------------------------------- #
# UNLOCK (loose) — the magnitude is the right order but OVERSHOOTS (the overturn).
# --------------------------------------------------------------------------- #
def test_emergent_lapse_rate_feedback_overshoots_the_global_mean(moist_column):
    """λ_LR ≈ 1.5: the right order, but it OVERSHOOTS the global-mean Soden & Held 0.84 (the overturn).

    A single global moist-adiabat column captures the *tropical* branch only; the magnitude is loose
    (rides the τ shape + WV loading, the wall). Banked as order + sign, NOT as a match to 0.84.
    """
    lr = moist_column.feedback_kernel().lapse_rate
    assert 1.0 < lr < 2.2                                              # order-1, ~1.5
    assert lr > rad.SH_LAPSE_RATE + 0.3                               # clearly above the global-mean 0.84 → overshoot


def test_moist_adiabat_pushes_B_above_climlabs_two(fixed_column, moist_column):
    """With the emergent lapse-rate feedback the with-WV B sits ABOVE climlab's 2 (overshoot), not at it.

    The fixed-Γ default sits below 2 (omits the feedback); making it emergent overshoots — the headline
    reconciliation with the existing decomposition (0.84 is the global mean, this is the tropical column).
    """
    assert fixed_column.feedback_kernel().total < B_OLR               # default omits LR → below 2
    assert moist_column.feedback_kernel().total > B_OLR               # emergent (tropical) LR → above 2


def test_present_operating_point_is_preserved_after_recalibration(moist_column):
    """The moist-adiabat column is recalibrated so OLR(288)=239 — same operating point, only the slope changes."""
    assert moist_column.outgoing_longwave(rad.PRESENT_SURFACE_T) == pytest.approx(rad.PRESENT_OLR, abs=0.5)


# --------------------------------------------------------------------------- #
# The banked-figure guard (slow).
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_demo_reproduces_the_lapse_rate_headline():
    # Guards docs/figures/planet-lapse-rate.png: the emergent λ_LR is positive, order ~1.5, overshoots the
    # global-mean 0.84, and the moist-adiabat column's B sits above climlab's 2.
    from planet import demo_lapse_rate as demo
    r = demo.compute()
    assert r.lapse_rate_feedback > 1.0
    assert r.lapse_rate_feedback > rad.SH_LAPSE_RATE
    assert r.b_total_moist > B_OLR
    assert abs(r.kernel_residual) < 3e-3
