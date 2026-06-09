"""Conservation seals — mass (machine-exact) + the discriminating finite-amplitude PV/enstrophy leg.

The contract's conservation invariants (plan §3, Phase-3 triad):

* **Mass ``∫h`` — machine precision.** Flux-form continuity telescopes on the periodic
  domain, so ``Σ h ΔxΔy`` changes only by floating-point round-off — true for *any*
  state and step (the shallow-water analogue of the diffusion engine's no-flux seal).
* **Energy (KE+PE) — bounded, dt→0 convergent drift.** The symmetric vector-invariant
  scheme conserves energy *semi-discretely*, so the drift is purely the RK3 time-truncation
  and shrinks as ``dt³``. Asserted as a convergent bound, not machine-exact (the honest
  claim for a centered explicit solver — see the module docstring).
* **Potential vorticity / enstrophy — bounded at finite amplitude (THE Coriolis seal).**
  This is the *discriminating* leg, and — per the design review — it is run at **finite
  amplitude** (a balanced vortex at Rossby number ~0.5 that genuinely advects PV around),
  where a wrong Coriolis averaging would show up as spurious PV-extrema growth or gross
  potential-enstrophy drift. At small amplitude this test would be near-vacuous (advection
  ≈ 0); the finite-amplitude requirement is what gives it teeth. Enstrophy is conserved to
  a small *spatial*-discretization-limited bound (it shrinks with Δx, not dt — this scheme
  conserves energy, not enstrophy, semi-discretely), and PV extrema barely move.
* **``∫ζ dA = 0``** — the relative vorticity is a discrete curl, so its area integral
  vanishes to machine precision (the structural fact behind circulation conservation).
"""
import numpy as np
import pytest

from engines.fluid import ShallowWater, SWState, uniform_grid


# --------------------------------------------------------------------------- #
# Helpers: build a geostrophically-balanced Gaussian vortex (finite amplitude).
# --------------------------------------------------------------------------- #
def balanced_vortex(grid, sw, amp_eta, sigma):
    """A geostrophically-balanced Gaussian height anomaly: u = −(g/f)∂η/∂y, v = (g/f)∂η/∂x,
    sampled analytically at the C-grid face locations (exact geostrophy for the smooth η)."""
    f0 = sw.f0
    cx, cy = grid.Lx / 2, grid.Ly / 2
    eta = lambda x, y: amp_eta * np.exp(-(((x - cx) ** 2 + (y - cy) ** 2)) / (2 * sigma ** 2))
    Xc, Yc = grid.center_mesh()
    xU = grid.x_corners()[None, :] * np.ones((grid.ny, 1)); yU = grid.y_centers()[:, None] * np.ones((1, grid.nx))
    xV = grid.x_centers()[None, :] * np.ones((grid.ny, 1)); yV = grid.y_corners()[:, None] * np.ones((1, grid.nx))
    u = -(sw.g / f0) * eta(xU, yU) * (-(yU - cy) / sigma ** 2)
    v = (sw.g / f0) * eta(xV, yV) * (-(xV - cx) / sigma ** 2)
    return SWState(h=sw.H + eta(Xc, Yc), u=u, v=v)


# --------------------------------------------------------------------------- #
# Mass — machine precision, any state, any step
# --------------------------------------------------------------------------- #
def test_mass_conserved_to_machine_precision():
    g = uniform_grid(4e6, 4e6, 48, 48)
    sw = ShallowWater(g, 9.81, 1000.0, f0=1e-4, beta=1.6e-11)
    s = balanced_vortex(g, sw, amp_eta=40.0, sigma=4e5)
    m0 = sw.mass(s)
    dt = sw.max_dt(s)
    for _ in range(300):
        s = sw.step(s, dt)
    assert sw.mass(s) == pytest.approx(m0, rel=1e-12)


def test_relative_vorticity_integrates_to_zero():
    g = uniform_grid(3e6, 3e6, 40, 40)
    sw = ShallowWater(g, 9.81, 800.0, f0=1e-4)
    s = balanced_vortex(g, sw, amp_eta=30.0, sigma=3e5)
    zeta = sw.relative_vorticity(s)
    assert np.sum(zeta) * g.cell_area == pytest.approx(0.0, abs=1e-6 * abs(sw.f0) * g.Lx * g.Ly)


# --------------------------------------------------------------------------- #
# The discriminating leg — finite-amplitude PV / potential enstrophy
# --------------------------------------------------------------------------- #
def test_finite_amplitude_pv_and_enstrophy_bounded():
    """The Coriolis seal: a balanced vortex at Rossby number ~0.5, integrated several inertial
    periods, must conserve mass exactly, hold potential enstrophy to a small bound, and not grow
    spurious PV extrema (PV is materially conserved → extrema cannot be created without forcing)."""
    g = uniform_grid(4e6, 4e6, 48, 48)
    f0 = 1e-4
    sw = ShallowWater(g, 9.81, 1000.0, f0=f0)
    s0 = balanced_vortex(g, sw, amp_eta=40.0, sigma=4e5)

    # genuinely finite amplitude — the whole point of this leg
    Ro = np.max(np.abs(sw.relative_vorticity(s0))) / f0
    assert Ro > 0.3

    Z0 = sw.potential_enstrophy(s0)
    m0 = sw.mass(s0)
    q0 = sw.potential_vorticity(s0)
    dt = sw.max_dt(s0)
    n = int(3 * 2 * np.pi / f0 / dt)            # ~3 inertial periods
    s = s0
    for _ in range(n):
        s = sw.step(s, dt)
    qf = sw.potential_vorticity(s)

    assert sw.mass(s) == pytest.approx(m0, rel=1e-11)
    assert abs(sw.potential_enstrophy(s) - Z0) / Z0 < 1e-4        # bounded enstrophy drift
    # PV materially conserved: no spurious extrema creation (no forcing/dissipation). The bound is
    # deliberately loose (~1%) — the seal catches a GROSS Coriolis error (which blows PV extrema out
    # by orders of magnitude), not this margin, so it must not be BLAS/numpy-fragile in the fast lane.
    assert qf.max() <= q0.max() * (1 + 1e-2)
    assert qf.min() >= q0.min() * (1 - 1e-2)


@pytest.mark.slow
def test_energy_drift_converges_third_order_in_dt():
    """Energy is conserved semi-discretely, so its drift is RK3 time-truncation: halving dt cuts
    the drift by ~2³=8 (third order). Confirms energy is the scheme's semi-discrete invariant."""
    g = uniform_grid(4e6, 4e6, 64, 64)
    f0 = 1e-4
    sw = ShallowWater(g, 9.81, 1000.0, f0=f0)
    s0 = balanced_vortex(g, sw, amp_eta=40.0, sigma=4e5)
    E0 = sw.energy(s0)
    T_end = 2 * 2 * np.pi / f0
    drifts = []
    base = sw.max_dt(s0)
    for dt in (base, base / 2, base / 4):
        s = s0
        for _ in range(int(T_end / dt)):
            s = sw.step(s, dt)
        drifts.append(abs(sw.energy(s) - E0) / E0)
    # each halving of dt reduces the drift by a factor ≳ 5 (≈ 2³, allowing margin)
    assert drifts[0] / drifts[1] > 5
    assert drifts[1] / drifts[2] > 5
