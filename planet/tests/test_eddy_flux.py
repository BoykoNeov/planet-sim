"""Planet rung-1 step-2 **Phase-B** validation: the EMERGENT eddy heat flux (``slow`` — runs the sim).

:mod:`planet.eddy_flux` diagnoses ``⟨v'θ'⟩(φ)`` from a passive tracer advected on the *released*
barotropically-unstable jet — the real flux that fills Phase A's ``flux_fn`` seam. Running the
shallow-water life cycle is expensive, so this is one ``slow`` test (the **fast**
geometry-correspondence legs — the spherical operator + the order-unity ``cos φ`` metric, pure math —
live in ``test_transport.py``). What is asserted, and its honesty class (see :mod:`planet.eddy_flux`):

* **(headline — DIRECTION banked) state-dependent diffusivity.** A *flatter* EBM gradient → a weaker
  jet → weaker eddies → a **smaller** ``κ_eff`` (``α`` held fixed). This non-circularity is what makes
  the two-way loop a real, right-signed feedback rather than a re-labelling of a fixed ``D``.
* **(sign) the emergent flux is net down-gradient** (``κ_bulk > 0`` — up-gradient would be negative)
  on a genuinely **unstable** jet (Rayleigh–Kuo), and is **mostly reversible** (the irreversible
  fraction is ``~0.1``: the instantaneous flux oscillates in sign).
* **(FINDING, not a manufactured match) the flux does NOT tightly reduce to the EBM operator at
  rung 1.** The resolved flux-divergence is only *partially* down-gradient-shaped (correlation ``~0.6``
  with smooth diffusion, and near-vacuous on the near-linear gradient) — the tight reduction becomes
  non-vacuous only at rung 3, where the geometry correspondence (``test_transport.py``) is in place.
* **(seam) ``close_loop`` routes the emergent ``D_eff``** through the Phase-A bridge + re-equilibration
  with the **right sign** (weaker transport ⇒ steeper contrast); the absolute climate is degenerate
  (the named magnitude edge) and not banked.

**Magnitude is named, not asserted:** no assertion pins ``κ_eff``'s value — only its sign, ordering,
and reversibility, the legs that carry validation.

**One test, by design:** both life cycles run inside a single test function (not split across
fixtures) so the gate's default ``-n auto`` (``--dist load``, which does *not* share module-scoped
fixtures across xdist workers) computes exactly **two** life cycles, never re-running them per worker.
"""
import pytest

from planet import eddy_flux as ef
from planet.albedo import EBMParams, present_day_climate
from planet.ebm import D_TRANSPORT

pytestmark = pytest.mark.slow

# A reduced grid so the two life cycles stay within a reasonable slow-test budget. The diagnosed
# magnitude is named-not-banked anyway; the legs asserted here (sign, climate ordering, reversibility,
# the non-diffusive reduction) are robust to resolution (verified at nx=80 vs 96 at build).
NX = 64


def test_emergent_eddy_flux_phase_b():
    """The whole Phase-B eddy-sim seal in one run-pair: sign + state-dependence + non-tight reduction
    + close_loop direction. (Two life cycles total — see the module docstring on why it is one test.)"""
    steep = ef.eddy_life_cycle(present_day_climate(EBMParams(s2=-0.48)), nx=NX, ny=NX)
    flat = ef.eddy_life_cycle(present_day_climate(EBMParams(s2=-0.32)), nx=NX, ny=NX)

    # -- the emergent flux is unstable, down-gradient, and mostly reversible -- #
    assert steep.rayleigh_kuo                        # the jet IS barotropically unstable (emergent eddies)
    assert steep.jet_speed > 5.0                     # a real jet spun up
    assert steep.kappa_bulk > 0.0                    # net DOWN-gradient (up-gradient ⇒ negative)
    assert steep.D_eff > 0.0
    assert 0.0 < steep.irreversible_fraction < 0.3   # mostly reversible (measured ~0.1: ~90% sloshing)

    # -- THE headline: state-dependent diffusivity (the non-circularity). A flatter climate makes a
    #    weaker jet → weaker eddies → a smaller κ_eff. Without it the loop would be cosmetic. -- #
    assert flat.jet_speed < steep.jet_speed
    assert flat.rayleigh_kuo                         # the flat jet is still unstable (a fair comparison)
    assert flat.kappa_bulk < steep.kappa_bulk

    # -- the tight reduction is a FINDING, not a match: the resolved divergence is jet-localised and
    #    does NOT match the smooth operator's shape (the geometry is delivered for rung 3, not here). -- #
    red = ef.reduction_to_ebm_operator(steep)
    assert red["shape_correlation"] < 0.8            # ~0.6: partial, NOT tight (tight would be ≳0.9)

    # -- the seam: the emergent D_eff routes through the Phase-A bridge with the right sign -- #
    out = ef.close_loop(steep)
    assert out["D_eff"] < D_TRANSPORT                # ~1000× below rung-0 — the named magnitude edge
    assert out["steeper"]                            # weaker transport ⇒ steeper contrast (right sign)
