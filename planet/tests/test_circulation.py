"""Planet Phase-3 validation: the planetary instantiation of the frozen shallow-water engine.

The *engine* seal (``engines/fluid/tests/``) validates the generic solver (wave speeds,
geostrophic balance, PV at finite amplitude). *Here* the **planetary numbers** are pinned and
checked — the consumer-side validation, exactly as Steel/Chip validate their diffusion-spine
constants against published data:

* the β-plane constants are Earth's (``f₀ = 2Ω sin φ``, ``β = 2Ω cos φ / a``);
* the equivalent depth yields the cited **extratropical deformation radius ``L_R ≈ 1000 km``**
  (calibrated, loose — [[shallow-water-source]]);
* **geostrophic adjustment** of a height anomaly settles to the analytic Helmholtz state over
  ``L_R`` (the published Rossby benchmark), conserving mass exactly; and
* a balanced mode propagates **westward** (Rossby) at ~the analytic phase speed.
"""
import numpy as np
import pytest

from projects.planet import circulation as circ


# --------------------------------------------------------------------------- #
# Pinned planetary constants & scales ([[shallow-water-source]])
# --------------------------------------------------------------------------- #
def test_coriolis_constants():
    assert circ.OMEGA_EARTH == pytest.approx(7.292e-5)
    assert circ.R_EARTH == pytest.approx(6.371e6)
    # f0 = 2Ω sin45, β = 2Ω cos45 / a
    assert circ.coriolis_f0(45.0) == pytest.approx(2 * 7.292e-5 * np.sin(np.radians(45)))
    assert circ.coriolis_beta(45.0) == pytest.approx(2 * 7.292e-5 * np.cos(np.radians(45)) / 6.371e6)
    # f0 grows toward the pole, β shrinks
    assert circ.coriolis_f0(60) > circ.coriolis_f0(30)
    assert circ.coriolis_beta(60) < circ.coriolis_beta(30)


def test_deformation_radius_is_extratropical_scale():
    grid, sw = circ.midlatitude_beta_plane(32, 32)
    # the cited extratropical deformation radius, ~1000 km (loose band)
    assert 0.7e6 < sw.rossby_radius < 1.3e6
    # the domain spans several deformation radii and L_R is resolved by several cells
    assert grid.Lx / sw.rossby_radius == pytest.approx(6.0, rel=0.01)
    assert grid.dx < 0.25 * sw.rossby_radius


def test_beta_plane_is_an_fplane_when_beta_zero():
    _, sw = circ.midlatitude_beta_plane(16, 16, beta=0.0)
    f = sw.f_at_corners()
    assert np.allclose(f, circ.coriolis_f0())                 # constant f, no β


def test_deformation_radius_value():
    # L_R = √(gH)/f0 from the pinned constants (no integration — the planetary-number check)
    _, sw = circ.midlatitude_beta_plane(16, 16)
    assert sw.rossby_radius == pytest.approx(np.sqrt(circ.G_EARTH * circ.H_EQUIV) / circ.coriolis_f0())


# --------------------------------------------------------------------------- #
# The integration demos (slow — the engine's own fast suite already seals the physics;
# these validate the *planetary* instantiation end to end and bank the figure's numbers).
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_adjustment_matches_helmholtz_over_LR():
    adj = circ.geostrophic_adjustment(nx=96, ny=96, n_periods=25.0)
    rel = np.max(np.abs(adj.eta_balanced - adj.eta_helmholtz)) / np.max(np.abs(adj.eta_helmholtz))
    assert rel < 0.05                                         # the balanced remnant IS the analytic state
    assert adj.eta_balanced.max() == pytest.approx(adj.eta_helmholtz.max(), rel=0.05)
    assert adj.eta_balanced.max() < 0.3 * adj.eta_init.max()  # most of the bump radiated away
    assert np.abs(adj.mass).max() < 1e-10                     # mass conserved throughout


@pytest.mark.slow
def test_rossby_wave_westward_near_analytic():
    ros = circ.rossby_wave(nx=96, ny=96, frac_period=0.5)
    assert ros.c_measured < 0.0 and ros.c_analytic < 0.0      # westward
    assert ros.c_measured / ros.c_analytic == pytest.approx(1.0, abs=0.1)
