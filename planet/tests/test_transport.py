"""Planet rung-1 step-2 validation: the two-way feedback *machinery* (Phase A, synthetic flux).

:mod:`planet.transport` closes the EBM ⇄ circulation loop — the resolved meridional eddy heat flux
``⟨v'θ'⟩`` is diagnosed as an effective diffusivity ``D_eff`` and fed back into the EBM's transport.
Phase A validates the **feedback machinery** driven by a *synthetic* exactly-down-gradient flux
(the Phase-B eddy simulation is not built yet — it plugs into the ``flux_fn`` seam). What is
asserted, and its honesty class:

* **The κ→D bridge — physical, pinned (tight).** ``D = C_atm·κ/a²`` (``C_atm = c_p·p_s/g``, the
  atmospheric column heat capacity) is pinned to its **computed value** (``test_bridge…``) — not
  just round-tripped, which would let a wrong ``a²``/``C_atm`` cancel — so the "physical/citable
  bridge" claim is actually backed: rung-0's ``D = 0.555`` maps to ``κ ≈ 2.17×10⁶ m²/s``, the
  observed midlatitude eddy-diffusivity order.
* **Right-signed climate response (tight, non-tautological — the headline-supporting leg).** A
  *stronger* eddy flux (larger ``D_eff``) flattens the equator-to-pole contrast; a weaker one
  steepens it (``test_stronger_flux…``) — the EBM's genuine *physical* response to transport (it
  varies ``D`` and checks the climate, so it cannot pass by construction). One feedback pass; the
  converged self-consistent fixed point and the *emergent* flux are Phase B.
* **Feedback plumbing (tight).** ``bulk_diffusivity`` recovers a prescribed ``κ`` with the correct
  **sign** (down-gradient → positive; up-gradient → negative → the EBM rejects it), and
  ``two_way_pass`` *routes* the diagnosed ``D_eff`` into the EBM re-equilibration
  (``test_reduction_to_ebm``, ``test_rung0_is_a_fixed_point``).
* **What Phase A does NOT yet anchor (honest scope).** The design anchor "reduction-to-EBM"
  *reduces to rung-0 by construction* here — ``two_way_pass`` literally re-runs the scalar-``D``
  rung-0 EBM at ``D_eff``, so it cannot fail and is plumbing, not an independent test. The genuinely
  tight reduction (an *independent* two-way budget whose flux-divergence must match the EBM operator
  ``D·∂/∂x[(1−x²)∂T/∂x]``, plus the Cartesian-channel ↔ spherical-EBM geometry correspondence) needs
  the emergent flux and arrives in **Phase B**. (Distinct from the *refactor-hygiene* array≡scalar
  ``D`` test in ``test_ebm`` — that has no eddy content either.)
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
    # PIN THE BRIDGE ABSOLUTELY (not just the round-trip, which would let a wrong a²/C_atm cancel):
    # C_atm = c_p·p_s/g ≈ 1.037×10⁷ J m⁻² K⁻¹, and rung-0 D = 0.555 ⟺ κ ≈ 2.17×10⁶ m²/s — the
    # computed value, which also lands in the observed midlatitude eddy-diffusivity range (~1–5×10⁶).
    assert tr.C_ATM == pytest.approx(1.037e7, rel=1e-2)
    assert float(tr.ebm_D_to_kappa(D_TRANSPORT)) == pytest.approx(2.17e6, rel=2e-2)


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
    """PLUMBING (not the tight anchor — see the module docstring): the machinery recovers a
    prescribed κ from the synthetic flux and *routes* the diagnosed D_eff into re-equilibration —
    the re-equilibrated climate is a rung-0 EBM at the diagnosed D. Note this is correct *by
    construction* (two_way_pass re-runs the scalar-D rung-0 EBM) and the bridge cancels on both
    sides, so it does NOT test the bridge value (that is pinned absolutely above); the genuinely
    tight reduction (independent flux-divergence = EBM operator) is Phase B."""
    kappa_test = 1.5e6
    res = tr.two_way_pass(params=COARSE,
                          flux_fn=lambda theta, y: tr.diffusive_flux(theta, y, kappa_test))
    assert res.kappa_eff == pytest.approx(kappa_test, rel=1e-9)   # κ recovered from the flux
    ref = present_day_climate(replace(COARSE, D=res.D_eff))       # re-equilibration routes D_eff …
    assert np.allclose(res.climate_after.T, ref.T, atol=1e-6)     # … into a rung-0 EBM at that D


def test_rung0_is_a_fixed_point():
    """The down-gradient limit recovers rung 0: when the eddy flux equals the EBM's own diffusive
    flux (the default flux_fn), the two-way map returns rung-0 — D_eff = D (a κ→D round-trip) and the
    climate unchanged (by construction, the re-equilibration re-runs rung-0). A consistency check on
    the pipeline, not an independent bridge/structure test (those are pinned/deferred respectively)."""
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
