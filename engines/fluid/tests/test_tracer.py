"""Passive tracer-advection seal — the rung-1 engine extension (ADR 0005; GCM climb §5).

:mod:`engines.fluid` advects an optional scalar tracer ``θ`` in **flux form**
(``∂(hθ)/∂t = −∇·(hθ u)``), carried through the *same* SSP-RK3 as the dry dynamics and
reusing the same mass fluxes. The tracer is **strictly passive** — it does not feed back on
``(h, u, v)`` — so the extension is additive and the dry trajectory is unchanged. What the
suite guarantees (the honesty class mirrors the engine's mass / energy / enstrophy split):

* **Passivity / re-seal.** A state carrying a tracer evolves ``(h, u, v)`` **byte-identically**
  to the same state with ``tracer=None`` (``np.array_equal`` over many steps) — the concrete
  proof the dry dynamics did not move under the extension (the regression guard that, post
  ADR 0005, replaces the old "re-run the frozen seal").
* **Tracer mass ``∫hθ`` — machine precision.** Flux-form advection telescopes on the periodic
  domain exactly like the fluid mass ``∫h`` — the clean anchor of this extension.
* **Consistency.** A *uniform* tracer stays uniform under any flow (the scheme is consistent:
  for ``θ=const`` the tracer tendency reduces to ``θ·∂h/∂t``), so no spurious source/sink.
* **Analytic limit.** A uniform flow (``f=0``, flat ``h``) translates a smooth tracer blob at
  the flow speed: its centre of mass advects at ``U₀`` to <1 %, the field to a few %
  (the centered scheme's dispersive phase error — loose, named).
* **Tracer variance ``∫½hθ²`` — bounded** (NOT machine-exact like mass, NOT dt-convergent
  like energy): like potential enstrophy it is a near-invariant the centered scheme holds to a
  small, *spatially*-limited drift. **NOT monotone** — there is no flux limiter, so the scheme
  **over/undershoots** on sharp tracer gradients (Gibbs ripples). That is the named scope edge
  (a TVD/WENO limiter or hyperdiffusion is the unbuilt upgrade); we *assert the overshoot
  exists* rather than claim boundedness of ``θ``.
"""
import numpy as np
import pytest

from engines.fluid import ShallowWater, SWState, uniform_grid


# --------------------------------------------------------------------------- #
# Helper: a geostrophically-balanced Gaussian vortex (as in test_conservation),
# optionally carrying a tracer — a finite-amplitude flow that genuinely advects.
# --------------------------------------------------------------------------- #
def balanced_vortex(grid, sw, amp_eta, sigma, tracer=None):
    f0 = sw.f0
    cx, cy = grid.Lx / 2, grid.Ly / 2
    eta = lambda x, y: amp_eta * np.exp(-(((x - cx) ** 2 + (y - cy) ** 2)) / (2 * sigma ** 2))
    Xc, Yc = grid.center_mesh()
    xU = grid.x_corners()[None, :] * np.ones((grid.ny, 1)); yU = grid.y_centers()[:, None] * np.ones((1, grid.nx))
    xV = grid.x_centers()[None, :] * np.ones((grid.ny, 1)); yV = grid.y_corners()[:, None] * np.ones((1, grid.nx))
    u = -(sw.g / f0) * eta(xU, yU) * (-(yU - cy) / sigma ** 2)
    v = (sw.g / f0) * eta(xV, yV) * (-(xV - cx) / sigma ** 2)
    return SWState(h=sw.H + eta(Xc, Yc), u=u, v=v, tracer=tracer)


# --------------------------------------------------------------------------- #
# tracer=None: the dry path is unchanged
# --------------------------------------------------------------------------- #
def test_tracer_none_runs_the_dry_dynamics():
    g = uniform_grid(4e6, 4e6, 32, 32)
    sw = ShallowWater(g, 9.81, 1000.0, f0=1e-4, beta=1.6e-11)
    s = balanced_vortex(g, sw, amp_eta=30.0, sigma=4e5)        # tracer=None
    s2 = sw.step(s, sw.max_dt(s))
    assert s2.tracer is None


def test_passivity_dry_dynamics_bit_for_bit():
    """Carrying a tracer must not perturb (h, u, v) by a single bit (the re-seal)."""
    g = uniform_grid(4e6, 4e6, 48, 48)
    sw = ShallowWater(g, 9.81, 1000.0, f0=1e-4, beta=1.6e-11)
    X, Y = g.center_mesh()
    theta = 5.0 + np.sin(2 * np.pi * X / g.Lx) * np.cos(2 * np.pi * Y / g.Ly)
    s_dry = balanced_vortex(g, sw, amp_eta=40.0, sigma=4e5)
    s_trc = balanced_vortex(g, sw, amp_eta=40.0, sigma=4e5, tracer=theta)
    dt = sw.max_dt(s_dry)
    for _ in range(200):
        s_dry = sw.step(s_dry, dt)
        s_trc = sw.step(s_trc, dt)
    assert np.array_equal(s_dry.h, s_trc.h)
    assert np.array_equal(s_dry.u, s_trc.u)
    assert np.array_equal(s_dry.v, s_trc.v)


# --------------------------------------------------------------------------- #
# Tracer mass — machine precision (the clean anchor)
# --------------------------------------------------------------------------- #
def test_tracer_mass_conserved_to_machine_precision():
    g = uniform_grid(4e6, 4e6, 48, 48)
    sw = ShallowWater(g, 9.81, 1000.0, f0=1e-4, beta=1.6e-11)
    X, Y = g.center_mesh()
    cx, cy = g.Lx / 2, g.Ly / 2
    theta = 2.0 + np.exp(-(((X - cx) ** 2 + (Y - cy) ** 2) / (2 * 5e5 ** 2)))   # strictly positive ⇒ m0 = O(1)
    s = balanced_vortex(g, sw, amp_eta=40.0, sigma=4e5, tracer=theta)
    m0 = sw.tracer_mass(s)
    dt = sw.max_dt(s)
    for _ in range(300):
        s = sw.step(s, dt)
    assert sw.tracer_mass(s) == pytest.approx(m0, rel=1e-12)


def test_tracer_mass_raises_without_a_tracer():
    g = uniform_grid(2e6, 2e6, 16, 16)
    sw = ShallowWater(g, 9.81, 1000.0, f0=1e-4)
    s = balanced_vortex(g, sw, amp_eta=10.0, sigma=3e5)       # tracer=None
    with pytest.raises(ValueError):
        sw.tracer_mass(s)
    with pytest.raises(ValueError):
        sw.tracer_variance(s)


# --------------------------------------------------------------------------- #
# Consistency — a uniform tracer is preserved exactly (no spurious source)
# --------------------------------------------------------------------------- #
def test_uniform_tracer_stays_uniform():
    g = uniform_grid(4e6, 4e6, 40, 40)
    sw = ShallowWater(g, 9.81, 1000.0, f0=1e-4, beta=1.6e-11)
    s = balanced_vortex(g, sw, amp_eta=50.0, sigma=4e5, tracer=np.full((40, 40), 7.0))
    dt = sw.max_dt(s)
    for _ in range(200):
        s = sw.step(s, dt)
    assert s.tracer.max() - s.tracer.min() == pytest.approx(0.0, abs=1e-9)
    assert s.tracer.mean() == pytest.approx(7.0, rel=1e-12)


# --------------------------------------------------------------------------- #
# Analytic limit — uniform flow translates a smooth blob at U0
# --------------------------------------------------------------------------- #
def test_uniform_flow_translates_a_smooth_blob():
    """f=0, flat h, uniform u=U₀ ⇒ pure advection. The flow is an exact steady solution
    (so (h,u,v) are untouched) and the smooth blob translates at U₀ — its centre of mass to
    <1 %, the field to a few % (the centered scheme's dispersive error, loose)."""
    ny = nx = 64
    L = 2.0e6
    g = uniform_grid(L, L, nx, ny)
    sw = ShallowWater(g, 9.81, 1000.0, f0=0.0, beta=0.0)
    X, Y = g.center_mesh()
    U0 = 20.0
    sigma = L / 12
    cx, cy = L / 2, L / 2
    blob = lambda c: np.exp(-(((X - c) ** 2 + (Y - cy) ** 2) / (2 * sigma ** 2)))
    theta0 = blob(cx)
    s = SWState(h=sw.H * np.ones((ny, nx)), u=U0 * np.ones((ny, nx)), v=np.zeros((ny, nx)), tracer=theta0.copy())
    dt = sw.max_dt(s) * 0.5
    n = int(0.2 * L / U0 / dt)
    t = n * dt
    for _ in range(n):
        s = sw.step(s, dt)
    shift = (U0 * t) % L
    # the uniform translating flow is an exact steady solution → (h, u, v) untouched
    assert np.allclose(s.u, U0) and np.allclose(s.v, 0.0)
    # analytic shifted blob (sum the periodic images), and the centre of mass (circular mean in x)
    theta_exact = sum(blob((cx + shift) % L + k * L) for k in (-1, 0, 1))

    def com_x(th):
        w = th.sum(axis=0)
        ang = 2 * np.pi * g.x_centers() / L
        return np.arctan2((w * np.sin(ang)).sum(), (w * np.cos(ang)).sum()) % (2 * np.pi) * L / (2 * np.pi)

    disp = (com_x(s.tracer) - com_x(theta0)) % L
    assert disp == pytest.approx(shift, rel=0.01)
    assert np.max(np.abs(s.tracer - theta_exact)) < 0.04


# --------------------------------------------------------------------------- #
# Tracer variance — bounded; and the named not-monotone scope edge
# --------------------------------------------------------------------------- #
def test_tracer_variance_bounded():
    """∫½hθ² is a near-invariant the centered scheme holds to a small, spatially-limited drift —
    bounded (NOT machine-exact like mass, NOT dt-convergent like energy). Asserted loose."""
    g = uniform_grid(4e6, 4e6, 48, 48)
    sw = ShallowWater(g, 9.81, 1000.0, f0=1e-4, beta=1.6e-11)
    X, Y = g.center_mesh()
    cx = g.Lx / 2
    theta = 3.0 + np.tanh((X - cx) / 3e5)            # a front the vortex filaments
    s = balanced_vortex(g, sw, amp_eta=40.0, sigma=4e5, tracer=theta)
    V0 = sw.tracer_variance(s)
    dt = sw.max_dt(s)
    for _ in range(int(3 * 2 * np.pi / sw.f0 / dt)):
        s = sw.step(s, dt)
    assert abs(sw.tracer_variance(s) - V0) / V0 < 1e-3


def test_sharp_gradient_overshoots_not_monotone():
    """The named scope edge: no flux limiter ⇒ the centered scheme over/undershoots a sharp
    tracer front (Gibbs ripples). We ASSERT the overshoot exists — monotonicity is NOT claimed —
    while the tracer mass stays machine-exact through it."""
    ny = nx = 64
    L = 2.0e6
    g = uniform_grid(L, L, nx, ny)
    sw = ShallowWater(g, 9.81, 1000.0, f0=0.0, beta=0.0)
    X, Y = g.center_mesh()
    U0 = 20.0
    theta0 = np.where(X < L / 2, 1.0, 0.0)           # a sharp step
    s = SWState(h=sw.H * np.ones((ny, nx)), u=U0 * np.ones((ny, nx)), v=np.zeros((ny, nx)), tracer=theta0.copy())
    m0 = sw.tracer_mass(s)
    dt = sw.max_dt(s) * 0.5
    for _ in range(120):
        s = sw.step(s, dt)
    assert s.tracer.max() > 1.0 + 0.05               # overshoot above the initial max
    assert s.tracer.min() < 0.0 - 0.05               # undershoot below the initial min
    assert sw.tracer_mass(s) == pytest.approx(m0, rel=1e-11)   # mass still exact through the ripples
