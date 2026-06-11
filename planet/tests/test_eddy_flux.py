"""Planet rung-1 step-2 **Phase-B** validation: the EMERGENT eddy heat flux (``slow`` — runs the sim).

:mod:`planet.eddy_flux` diagnoses ``⟨v'θ'⟩(φ)`` from a passive tracer advected on the *released*
barotropically-unstable jet — the real flux that fills Phase A's ``flux_fn`` seam. Every test here
runs the shallow-water life cycle, so the whole module is ``slow``; the **fast** geometry-correspondence
legs (the spherical operator + the order-unity ``cos φ`` metric, pure math) live in ``test_transport.py``.
What is asserted, and its honesty class (see :mod:`planet.eddy_flux`):

* **(headline — DIRECTION banked) state-dependent diffusivity.** A *flatter* EBM gradient → a weaker
  jet → weaker eddies → a **smaller** ``κ_eff`` (``α`` held fixed). This non-circularity is what makes
  the two-way loop a real, right-signed feedback rather than a re-labelling of a fixed ``D``.
* **(sign) the emergent flux is net down-gradient** (``κ_bulk > 0`` — up-gradient would be negative)
  on a genuinely **unstable** jet (Rayleigh–Kuo), and is **mostly reversible** (the irreversible
  fraction is ``~10⁻³``: the instantaneous flux oscillates in sign).
* **(FINDING, not a manufactured match) the flux does NOT tightly reduce to the EBM operator at
  rung 1.** The resolved flux-divergence is jet-localised / structurally non-diffusive (low shape
  correlation with smooth down-gradient diffusion) — the geometry correspondence is *in place* for the
  strong baroclinic flux of rung 3, where the reduction becomes non-vacuous.
* **(seam) ``close_loop`` routes the emergent ``D_eff``** through the Phase-A bridge + re-equilibration
  with the **right sign** (weaker transport ⇒ steeper contrast); the absolute climate is degenerate
  (the named magnitude edge) and not banked.

**Magnitude is named, not asserted:** no test pins ``κ_eff``'s value — only its sign, ordering, and
reversibility, the legs that carry validation.
"""
import numpy as np
import pytest

from planet import eddy_flux as ef
from planet.albedo import EBMParams, present_day_climate
from planet.ebm import D_TRANSPORT

pytestmark = pytest.mark.slow

# A reduced grid so the two life cycles stay within a reasonable slow-test budget. The diagnosed
# magnitude is named-not-banked anyway; the legs asserted here (sign, climate ordering, reversibility,
# the non-diffusive reduction) are robust to resolution (verified at nx=80 vs 96 at build).
NX = 64


@pytest.fixture(scope="module")
def steep():
    """Present-day (steep gradient, s₂=−0.48): the strong-jet reference life cycle."""
    return ef.eddy_life_cycle(present_day_climate(EBMParams(s2=-0.48)), nx=NX, ny=NX)


@pytest.fixture(scope="module")
def flat():
    """A flatter, high-obliquity-like world (s₂=−0.32): weaker jet, weaker eddies."""
    return ef.eddy_life_cycle(present_day_climate(EBMParams(s2=-0.32)), nx=NX, ny=NX)


def test_emergent_flux_is_unstable_down_gradient_and_mostly_reversible(steep):
    assert steep.rayleigh_kuo                       # the jet IS barotropically unstable (emergent eddies)
    assert steep.jet_speed > 5.0                    # a real jet spun up
    assert steep.kappa_bulk > 0.0                   # net DOWN-gradient (up-gradient ⇒ negative)
    assert steep.D_eff > 0.0
    assert 0.0 < steep.irreversible_fraction < 0.3  # mostly reversible (measured ~0.1: ~90% sloshing)


def test_state_dependent_diffusivity_is_the_non_circularity(steep, flat):
    # THE headline: a flatter climate makes a weaker jet → weaker eddies → a smaller κ_eff. Without
    # this state dependence the two-way loop would be cosmetic (a fixed D in disguise).
    assert flat.jet_speed < steep.jet_speed
    assert flat.rayleigh_kuo                        # the flat jet is still unstable (a fair comparison)
    assert flat.kappa_bulk < steep.kappa_bulk


def test_reduction_to_ebm_operator_is_not_diffusive_at_rung1(steep):
    # The honest Phase-B finding: the resolved barotropic flux-divergence is jet-localised and
    # structurally non-diffusive, so it does NOT match the smooth down-gradient operator's shape. The
    # geometry correspondence (test_transport.py) is in place for rung 3, where the flux is strong.
    red = ef.reduction_to_ebm_operator(steep)
    assert red["shape_correlation"] < 0.8           # ~0.6: partial, NOT tight (tight would be ≳0.9)


def test_close_loop_routes_emergent_Deff_with_the_right_sign(steep):
    out = ef.close_loop(steep)
    assert out["D_eff"] < D_TRANSPORT               # ~1000× below rung-0 — the named magnitude edge
    # weaker transport ⇒ steeper equator-to-pole contrast (the Phase-A right-signed leg, emergent κ)
    assert out["steeper"]
