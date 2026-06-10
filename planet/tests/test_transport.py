"""Planet rung-1 step-2 validation: the two-way feedback *machinery* (Phase A, synthetic flux).

:mod:`planet.transport` closes the EBM ⇄ circulation loop — the resolved meridional eddy heat flux
``⟨v'θ'⟩`` is diagnosed as an effective diffusivity ``D_eff`` and fed back into the EBM's transport.
Phase A validates the **feedback machinery** driven by a *synthetic* exactly-down-gradient flux
(the Phase-B eddy simulation is not built yet — it plugs into the ``flux_fn`` seam). What is
asserted, and its honesty class:

* **THE anchor — reduction-to-EBM (tight).** The down-gradient closure ``⟨v'θ'⟩ = −D_eff·∂θ̄/∂y``
  has the **same form** as the EBM's transport term, so the two-way model with a constant
  flow-diagnosed ``D_eff`` *is* a rung-0 diffusive EBM with that ``D``: re-equilibrating at the
  diagnosed ``D_eff`` reproduces a rung-0 EBM run directly at the bridge-implied ``D``
  (``test_reduction_to_ebm``), and rung-0 is a **fixed point** of the two-way map when the flux is
  the EBM's own diffusive flux (``test_rung0_is_a_fixed_point``). This threads the *entire* pipeline
  (flux → gradient → κ → D → climate) — a sign, an ``a²``, the ``C_atm`` value, or the band mask
  being wrong all break it. (Distinct from the *refactor-hygiene* array≡scalar ``D`` test in
  ``test_ebm`` — that has no eddy content.)
* **The κ→D bridge (physical, tight).** ``D = C_atm·κ/a²`` round-trips exactly, and rung-0's
  ``D = 0.555`` maps to ``κ ≈ 2.2×10⁶ m²/s`` — the observed midlatitude eddy-diffusivity order.
* **Right-signed climate response (the headline-supporting leg).** A *stronger* eddy flux (larger
  ``D_eff``) flattens the equator-to-pole contrast; a weaker one steepens it — monotone. (One
  feedback pass; the converged self-consistent fixed point and the *emergent* flux are Phase B.)
* **Diagnosis (tight).** ``bulk_diffusivity`` recovers the prescribed ``κ`` from an exactly
  down-gradient flux and reports the correct **sign** (down-gradient → positive).
"""
from dataclasses import replace

import numpy as np
import pytest

from planet import transport as tr
from planet.albedo import EBMParams, present_day_climate
from planet.ebm import D_TRANSPORT

# A coarse EBM (fewer cells) so the repeated relaxations stay fast; the physics is grid-converged
# enough for these checks (the reduction is a self-consistency, independent of resolution).
COARSE = EBMParams(n_cells=60)


# --------------------------------------------------------------------------- #
# The κ→D bridge — physical, with the magnitude sanity check.
# --------------------------------------------------------------------------- #
def test_bridge_roundtrips_and_recovers_observed_eddy_diffusivity():
    # Round-trip is exact (the bridge is a single multiplicative constant C_atm/a²).
    for D in (0.2, 0.555, 1.0):
        assert float(tr.kappa_to_ebm_D(tr.ebm_D_to_kappa(D))) == pytest.approx(D, rel=1e-12)
    # Magnitude sanity: rung-0 D = 0.555 ⟺ κ in the observed midlatitude range (~1–5×10⁶ m²/s).
    kappa0 = float(tr.ebm_D_to_kappa(D_TRANSPORT))
    assert 1.0e6 < kappa0 < 5.0e6
    # C_atm is the textbook tropospheric column heat capacity (~10⁷ J m⁻² K⁻¹).
    assert 0.8e7 < tr.C_ATM < 1.2e7


# --------------------------------------------------------------------------- #
# Diagnosing the diffusivity — recovery + the down-gradient sign.
# --------------------------------------------------------------------------- #
def test_bulk_diffusivity_recovers_kappa_and_signs_down_gradient_positive():
    y = np.linspace(0.0, 3.0e6, 64)
    theta = 30.0 - 1.0e-5 * y                       # a smooth poleward-decreasing tracer (g < 0)
    g = np.gradient(theta, y)
    kappa = 1.7e6
    F_down = -kappa * g                             # exactly down-gradient
    assert tr.bulk_diffusivity(F_down, g) == pytest.approx(kappa, rel=1e-10)
    assert tr.bulk_diffusivity(F_down, g) > 0.0     # down-gradient ⇒ positive
    assert tr.bulk_diffusivity(+kappa * g, g) < 0.0  # an up-gradient flux ⇒ negative
    # pointwise companion agrees where the gradient is non-zero
    assert tr.pointwise_diffusivity(F_down, g) == pytest.approx(kappa, rel=1e-10)


def test_bulk_diffusivity_raises_on_zero_gradient():
    g = np.zeros(10)
    with pytest.raises(ValueError):
        tr.bulk_diffusivity(np.ones(10), g)


# --------------------------------------------------------------------------- #
# THE anchor — reduction-to-EBM.
# --------------------------------------------------------------------------- #
def test_reduction_to_ebm():
    """With a constant flow-diagnosed D_eff, the two-way model *is* a rung-0 diffusive EBM with that
    D: re-equilibrating at the diagnosed D_eff reproduces a rung-0 EBM run directly at the bridge-
    implied D. Driven by an INDEPENDENT κ (not the bridge's own κ₀), so it is not circular."""
    kappa_test = 1.5e6                               # an independent physical diffusivity
    res = tr.two_way_pass(params=COARSE,
                          flux_fn=lambda theta, y: tr.diffusive_flux(theta, y, kappa_test))
    # the machinery recovers the prescribed κ from the synthetic flux …
    assert res.kappa_eff == pytest.approx(kappa_test, rel=1e-9)
    # … and its re-equilibrated climate is a rung-0 EBM at the bridge-implied D (the reduction)
    ref = present_day_climate(replace(COARSE, D=float(tr.kappa_to_ebm_D(kappa_test))))
    assert np.allclose(res.climate_after.T, ref.T, atol=1e-6)


def test_rung0_is_a_fixed_point():
    """The down-gradient limit recovers rung 0: when the eddy flux equals the EBM's own diffusive
    flux (the default flux_fn), the two-way map returns rung-0 — D_eff = D and the climate unchanged."""
    res = tr.two_way_pass(params=COARSE)             # default flux = rung-0's diffusive flux
    assert res.D_eff == pytest.approx(D_TRANSPORT, rel=1e-9)
    assert res.kappa_eff == pytest.approx(float(tr.ebm_D_to_kappa(D_TRANSPORT)), rel=1e-9)
    assert np.allclose(res.climate_after.T, res.climate_before.T, atol=1e-6)
    assert res.contrast_after == pytest.approx(res.contrast_before, abs=1e-4)


# --------------------------------------------------------------------------- #
# Right-signed climate response (the headline-supporting leg, one pass).
# --------------------------------------------------------------------------- #
def test_stronger_flux_flattens_the_gradient():
    """More eddy transport (larger D_eff) ⇒ a flatter equator-to-pole contrast; less ⇒ steeper.
    The right-signed feedback that, in Phase B, makes the two-way loop a real negative feedback."""
    k0 = float(tr.ebm_D_to_kappa(D_TRANSPORT))
    strong = tr.two_way_pass(params=COARSE, flux_fn=lambda th, y: tr.diffusive_flux(th, y, 2.0 * k0))
    weak = tr.two_way_pass(params=COARSE, flux_fn=lambda th, y: tr.diffusive_flux(th, y, 0.5 * k0))
    base_contrast = strong.contrast_before          # identical present-day start for both
    assert strong.D_eff > D_TRANSPORT > weak.D_eff
    assert strong.contrast_after < base_contrast    # stronger transport flattens
    assert weak.contrast_after > base_contrast      # weaker transport steepens
    assert strong.contrast_after < weak.contrast_after   # monotone in D_eff


# --------------------------------------------------------------------------- #
# The flux_fn seam + the diagnosed flux is net down-gradient (the Phase-B handoff point).
# --------------------------------------------------------------------------- #
def test_default_flux_is_net_down_gradient_over_the_interior():
    res = tr.two_way_pass(params=COARSE)
    assert res.kappa_eff > 0.0                       # net down-gradient over the window-flat interior
    assert res.interior.sum() > 0                    # the flat top is a non-empty band
    assert res.interior.sum() < res.phi.size         # but a strict subset (the taper is excluded)
