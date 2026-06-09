"""Geostrophy anchors — balanced-state persistence + Rossby's geostrophic-adjustment benchmark.

* **Geostrophic balance is a steady state.** A zonal jet ``u(y)``, ``v=0``, in geostrophic
  balance ``f·u = −g ∂η/∂y`` is an *exact* steady solution of the full nonlinear system (with
  ``v=0`` and no x-dependence, the gradient-wind / vorticity-flux corrections cancel identically),
  so it must sit still: ``v`` stays ~0 and ``u`` barely changes. The analytical-limit leg's
  "balance" half, valid at any amplitude.

* **Geostrophic adjustment (the published benchmark — Rossby's problem).** An unbalanced
  localized height bump (``u=v=0``) radiates inertia-gravity waves and **settles to a
  geostrophically-balanced remnant of scale ``L_R``**. Linear potential vorticity is conserved
  pointwise, so the adjusted height field solves the Helmholtz relation
  ``(1 − L_R²∇²) η_adj = η_init`` — an exact analytic target the time-averaged late state must
  reproduce (it does, to ~few %), with the central height drawn down as energy radiates away.
  The discriminating content is the **``L_R`` scale**: a wrong ``f₀`` (wrong ``L_R``) misses the
  adjusted state's width and amplitude.
"""
import numpy as np
import pytest

from engines.fluid import ShallowWater, SWState, uniform_grid


# --------------------------------------------------------------------------- #
# Geostrophic balance — a balanced zonal jet is steady
# --------------------------------------------------------------------------- #
def test_balanced_zonal_jet_is_steady():
    g = uniform_grid(4e6, 4e6, 64, 64)
    f0 = 1e-4
    sw = ShallowWater(g, 9.81, 1000.0, f0=f0)
    ly = 2 * np.pi / g.Ly
    eta0 = 20.0                                   # finite amplitude
    Xc, Yc = g.center_mesh()
    eta = eta0 * np.cos(ly * Yc)
    yU = g.y_centers()[:, None] * np.ones((1, g.nx))
    u = (sw.g * eta0 * ly / f0) * np.sin(ly * yU)   # f u = -g dη/dy
    v = np.zeros((g.ny, g.nx))
    s0 = SWState(h=sw.H + eta, u=u, v=v)

    jet_speed = np.max(np.abs(u))
    s = s0
    dt = sw.max_dt(s0)
    for _ in range(400):
        s = sw.step(s, dt)
    # the jet holds: v stays a tiny fraction of the jet speed; u barely moves
    assert np.max(np.abs(s.v)) < 1e-3 * jet_speed
    assert np.max(np.abs(s.u - u)) < 1e-3 * jet_speed
    assert sw.mass(s) == pytest.approx(sw.mass(s0), rel=1e-12)


# --------------------------------------------------------------------------- #
# Geostrophic adjustment — bump → balanced state over L_R (Helmholtz target)
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_geostrophic_adjustment_matches_helmholtz_over_LR():
    g = uniform_grid(6e6, 6e6, 96, 96)
    f0 = 1e-4
    sw = ShallowWater(g, 9.81, 1000.0, f0=f0)
    LR = sw.rossby_radius
    Xa, Ya = g.center_mesh()
    sigma = LR / 3.0                              # narrower than L_R → adjustment is visible
    eta0 = 1.0 * np.exp(-(((Xa - g.Lx / 2) ** 2 + (Ya - g.Ly / 2) ** 2)) / (2 * sigma ** 2))
    s0 = SWState(h=sw.H + eta0, u=np.zeros((g.ny, g.nx)), v=np.zeros((g.ny, g.nx)))

    dt = sw.max_dt(s0) * 0.6
    T_total = 25 * 2 * np.pi / f0
    n = int(T_total / dt)
    nwin = n // 3
    eta_acc = np.zeros((g.ny, g.nx))
    s = s0
    for i in range(n):
        s = sw.step(s, dt)
        if i >= n - nwin:
            eta_acc += (s.h - sw.H)
    eta_bal = eta_acc / nwin

    # analytic adjusted state: (1 − L_R² ∇²) η_adj = η_init, solved spectrally on the periodic grid
    kx = 2 * np.pi * np.fft.fftfreq(g.nx, d=g.dx)
    ky = 2 * np.pi * np.fft.fftfreq(g.ny, d=g.dy)
    KX, KY = np.meshgrid(kx, ky)
    eta_helm = np.real(np.fft.ifft2(np.fft.fft2(eta0) / (1.0 + LR ** 2 * (KX ** 2 + KY ** 2))))

    # the balanced remnant matches Helmholtz to a few %, and the centre has drawn down a lot
    assert np.max(np.abs(eta_bal - eta_helm)) / np.max(np.abs(eta_helm)) < 0.05
    assert eta_bal.max() < 0.25 * eta0.max()                 # most of the bump radiated away
    assert eta_bal.max() == pytest.approx(eta_helm.max(), rel=0.05)
    assert sw.mass(s) == pytest.approx(sw.mass(s0), rel=1e-11)
