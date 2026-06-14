"""TVD flux-limiter seal — the opt-in rung-1 tracer upgrade (GCM climb §12.2).

The default tracer advection in :mod:`engines.fluid` is the *unlimited centered* flux
(second-order, but it over/undershoots sharp fronts — its Gibbs ripples are *asserted* in
:mod:`engines.fluid.tests.test_tracer`). Constructing :class:`ShallowWater` with a
``tracer_limiter`` swaps the centered face value for a **TVD-limited** one
``θ_face = θ_up + ½·ψ(r)·(θ_down − θ_up)``. What this suite pins (the triad):

* **Reduction (byte-for-bit).** ``tracer_limiter=None`` is the original centered scheme to the
  last bit (``np.array_equal``) — the limiter is a pure opt-in; the dry dynamics never move.
* **Monotone (the tight anchor).** A uniform x-flow advecting a sharp step — the exact mirror of
  ``test_sharp_gradient_overshoots_not_monotone`` — produces **no new extrema** (θ stays in the
  IC range, and positive) under every limiter, where the centered scheme overshoots. This is
  rigorous 1-D TVD (grid-aligned flow, flat h frozen by an exact steady state).
* **Conservation kept.** The scheme stays in conservative flux form, so ``∫hθ`` is still
  machine-exact with the limiter on; and a uniform tracer stays uniform (consistency).

**Honest scope edge (advisor):** strict TVD is a *1-D* property. In genuinely 2-D flow the
dimension-split limiting (Goodman–LeVeque) no longer guarantees a maximum principle — so the 2-D
test asserts only that the limiter *reduces* the overshoot vs centered, not that it eliminates it.
The other standard trade is first-order clipping at smooth extrema (a gentle limiter like van Leer
clips little); the limited scheme is also *dissipative*, so the tracer variance only decreases.
"""
import numpy as np
import pytest

from engines.fluid import ShallowWater, SWState, uniform_grid

LIMITERS = ["minmod", "vanleer", "mc", "superbee"]


def balanced_vortex(grid, sw, amp_eta, sigma, tracer=None):
    """A geostrophically-balanced Gaussian vortex (as in test_tracer / test_conservation)."""
    f0 = sw.f0
    cx, cy = grid.Lx / 2, grid.Ly / 2
    eta = lambda x, y: amp_eta * np.exp(-(((x - cx) ** 2 + (y - cy) ** 2)) / (2 * sigma ** 2))
    Xc, Yc = grid.center_mesh()
    xU = grid.x_corners()[None, :] * np.ones((grid.ny, 1)); yU = grid.y_centers()[:, None] * np.ones((1, grid.nx))
    xV = grid.x_centers()[None, :] * np.ones((grid.ny, 1)); yV = grid.y_corners()[:, None] * np.ones((1, grid.nx))
    u = -(sw.g / f0) * eta(xU, yU) * (-(yU - cy) / sigma ** 2)
    v = (sw.g / f0) * eta(xV, yV) * (-(xV - cx) / sigma ** 2)
    return SWState(h=sw.H + eta(Xc, Yc), u=u, v=v, tracer=tracer)


def uniform_step_state(g, sw, speed, axis="x"):
    """A sharp tracer step in an exact uniform-translation steady state (f=0, flat h).

    ``axis`` ∈ {"x", "y"} sets the flow/step direction; ``speed`` is **signed** — a negative speed
    drives the ``U<0`` / ``V<0`` upwind branch of the limiter (the quadrant a +x-only anchor never
    probes), so parametrizing over ±x, ±y is what makes the monotonicity anchor sign-complete."""
    X, Y = g.center_mesh()
    zeros = np.zeros((g.ny, g.nx))
    if axis == "x":
        theta0 = np.where(X < g.Lx / 2, 1.0, 0.0)
        u, v = speed * np.ones((g.ny, g.nx)), zeros
    else:
        theta0 = np.where(Y < g.Ly / 2, 1.0, 0.0)
        u, v = zeros, speed * np.ones((g.ny, g.nx))
    return SWState(h=sw.H * np.ones((g.ny, g.nx)), u=u, v=v, tracer=theta0.copy())


# --------------------------------------------------------------------------- #
# Validation — an unknown limiter name is rejected at construction
# --------------------------------------------------------------------------- #
def test_unknown_limiter_raises():
    g = uniform_grid(2e6, 2e6, 16, 16)
    with pytest.raises(ValueError):
        ShallowWater(g, 9.81, 1000.0, f0=1e-4, tracer_limiter="bogus")
    # the documented names and None all construct fine
    for name in LIMITERS + [None]:
        ShallowWater(g, 9.81, 1000.0, f0=1e-4, tracer_limiter=name)


# --------------------------------------------------------------------------- #
# Reduction — tracer_limiter=None is the centered scheme, byte-for-byte
# --------------------------------------------------------------------------- #
def test_limiter_none_reproduces_centered_bit_for_bit():
    g = uniform_grid(4e6, 4e6, 48, 48)
    X, Y = g.center_mesh()
    theta = 2.0 + np.exp(-(((X - g.Lx / 2) ** 2 + (Y - g.Ly / 2) ** 2) / (2 * 5e5 ** 2)))
    sw_default = ShallowWater(g, 9.81, 1000.0, f0=1e-4, beta=1.6e-11)               # no arg
    sw_none = ShallowWater(g, 9.81, 1000.0, f0=1e-4, beta=1.6e-11, tracer_limiter=None)
    sw_vl = ShallowWater(g, 9.81, 1000.0, f0=1e-4, beta=1.6e-11, tracer_limiter="vanleer")
    s0 = balanced_vortex(g, sw_default, amp_eta=40.0, sigma=4e5, tracer=theta)
    sd, sn, sv = s0, s0, s0
    for _ in range(30):
        sd = sw_default.step(sd, sw_default.max_dt(sd))
        sn = sw_none.step(sn, sw_none.max_dt(sn))
        sv = sw_vl.step(sv, sw_vl.max_dt(sv))
    # None is the default is the original centered scheme — to the last bit
    assert np.array_equal(sd.tracer, sn.tracer)
    # the limiter actually changes the tracer (it is not a silent no-op)
    assert not np.array_equal(sd.tracer, sv.tracer)
    # ...yet the dry (h, u, v) trajectory is byte-identical even WITH the limiter on (passivity):
    # _limited_face only touches the tracer fluxes Fx/Fy, never dh/du/dv.
    assert np.array_equal(sd.h, sv.h) and np.array_equal(sd.u, sv.u) and np.array_equal(sd.v, sv.v)


# --------------------------------------------------------------------------- #
# The tight anchor — a uniform flow advecting a step is monotone (no new extrema)
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("limiter", LIMITERS)
@pytest.mark.parametrize("axis,speed", [("x", 20.0), ("x", -20.0), ("y", 20.0), ("y", -20.0)])
def test_limited_step_is_monotone(limiter, axis, speed):
    """The mirror of test_sharp_gradient_overshoots_not_monotone: with a limiter the sharp
    step develops NO new extrema (θ ∈ [0, 1], stays positive) while ∫hθ stays machine-exact.

    Parametrized over ±x AND ±y so the monotonicity (no-new-extrema) guarantee is verified in all
    four directions, not just the +x quadrant a single anchor probes. (A pure step stays monotone
    even under a *buggy* limiter — it has no interior extremum to over-amplify — so the sign of the
    smoothness ratio is pinned by the separate direction-symmetry test below, not here.)"""
    ny = nx = 64
    L = 2.0e6
    g = uniform_grid(L, L, nx, ny)
    sw = ShallowWater(g, 9.81, 1000.0, f0=0.0, beta=0.0, tracer_limiter=limiter)
    s = uniform_step_state(g, sw, speed, axis=axis)
    m0 = sw.tracer_mass(s)
    dt = sw.max_dt(s) * 0.5
    for _ in range(120):
        s = sw.step(s, dt)
    assert s.tracer.max() <= 1.0 + 1e-9        # no overshoot above the IC max
    assert s.tracer.min() >= 0.0 - 1e-9        # no undershoot below the IC min (positivity)
    assert sw.tracer_mass(s) == pytest.approx(m0, rel=1e-11)
    # and the front genuinely moved (the scheme is not trivially frozen)
    assert not np.array_equal(s.tracer, uniform_step_state(g, sw, speed, axis=axis).tracer)


def test_centered_default_still_overshoots_the_same_step():
    """Belt-and-suspenders: the SAME step on the default (no-limiter) solver DOES overshoot —
    so the monotonicity above is the limiter's doing, not a too-gentle test setup."""
    ny = nx = 64
    L = 2.0e6
    g = uniform_grid(L, L, nx, ny)
    sw = ShallowWater(g, 9.81, 1000.0, f0=0.0, beta=0.0)        # default centered
    s = uniform_step_state(g, sw, 20.0)
    dt = sw.max_dt(s) * 0.5
    for _ in range(120):
        s = sw.step(s, dt)
    assert s.tracer.max() > 1.0 + 0.05
    assert s.tracer.min() < 0.0 - 0.05


def _vanleer_blob_peak(g, speed, axis):
    """Advect a well-resolved Gaussian a fixed distance under uniform flow; return the peak."""
    sw = ShallowWater(g, 9.81, 1000.0, f0=0.0, beta=0.0, tracer_limiter="vanleer")
    X, Y = g.center_mesh()
    sigma = g.Lx / 12
    blob = np.exp(-(((X - g.Lx / 2) ** 2 + (Y - g.Ly / 2) ** 2) / (2 * sigma ** 2)))
    u = (speed if axis == "x" else 0.0) * np.ones((g.ny, g.nx))
    v = (speed if axis == "y" else 0.0) * np.ones((g.ny, g.nx))
    s = SWState(h=sw.H * np.ones((g.ny, g.nx)), u=u, v=v, tracer=blob.copy())
    dt = sw.max_dt(s) * 0.5
    for _ in range(int(0.2 * g.Lx / abs(speed) / dt)):
        s = sw.step(s, dt)
    return s.tracer.max()


def test_limiter_is_direction_symmetric():
    """The discriminator with teeth for the upwind-ratio SIGN: the limiter must be
    reflection-equivariant — a smooth Gaussian advects equally well leftward as rightward (and
    down as up). A sign error in the ``r`` numerator for ``U<0``/``V<0`` silently degrades the
    negative-direction flux to first order (peak retention ~0.83 vs ~0.96) without ever
    overshooting, so a +x-only blob test cannot see it. Peak retention is **bit-identical** for ±
    on the correct scheme, making this a sharp catch; we also assert the peak stays > 0.9 (gentle
    van Leer barely clips a resolved Gaussian) — the buggy first-order branch fails both asserts."""
    g = uniform_grid(2e6, 2e6, 64, 64)
    for axis in ("x", "y"):
        peak_plus = _vanleer_blob_peak(g, 20.0, axis)
        peak_minus = _vanleer_blob_peak(g, -20.0, axis)
        assert peak_plus == pytest.approx(peak_minus, abs=1e-9)   # reflection-symmetric (bug: ~0.13 gap)
        assert peak_minus > 0.9                                   # well-retained both directions


# --------------------------------------------------------------------------- #
# Conservation & consistency hold with the limiter on
# --------------------------------------------------------------------------- #
def test_limiter_mass_conserved_to_machine_precision():
    g = uniform_grid(4e6, 4e6, 48, 48)
    sw = ShallowWater(g, 9.81, 1000.0, f0=1e-4, beta=1.6e-11, tracer_limiter="vanleer")
    X, Y = g.center_mesh()
    theta = 2.0 + np.tanh((X - g.Lx / 2) / 3e5)                 # sharp front, strictly positive ⇒ m0=O(big)
    s = balanced_vortex(g, sw, amp_eta=40.0, sigma=4e5, tracer=theta)
    m0 = sw.tracer_mass(s)
    dt = sw.max_dt(s)
    for _ in range(300):
        s = sw.step(s, dt)
    assert sw.tracer_mass(s) == pytest.approx(m0, rel=1e-11)


def test_uniform_tracer_stays_uniform_with_limiter():
    g = uniform_grid(4e6, 4e6, 40, 40)
    sw = ShallowWater(g, 9.81, 1000.0, f0=1e-4, beta=1.6e-11, tracer_limiter="superbee")
    s = balanced_vortex(g, sw, amp_eta=50.0, sigma=4e5, tracer=np.full((40, 40), 7.0))
    dt = sw.max_dt(s)
    for _ in range(150):
        s = sw.step(s, dt)
    assert s.tracer.max() - s.tracer.min() == pytest.approx(0.0, abs=1e-9)
    assert s.tracer.mean() == pytest.approx(7.0, rel=1e-12)


# --------------------------------------------------------------------------- #
# 2-D — the honest claim: the limiter REDUCES the overshoot (not a maximum principle)
# --------------------------------------------------------------------------- #
def test_limiter_reduces_overshoot_in_2d():
    """In 2-D the dimension-split limiting no longer guarantees a maximum principle
    (Goodman–LeVeque), so we assert only that van Leer pulls the over/undershoot toward the
    IC range vs the centered scheme on the same filamenting front — and keeps ∫hθ exact."""
    g = uniform_grid(4e6, 4e6, 48, 48)
    X, _ = g.center_mesh()
    theta = 1.0 + np.tanh((X - g.Lx / 2) / 3e5)                 # front over [0, 2], strictly positive
    swc = ShallowWater(g, 9.81, 1000.0, f0=1e-4, beta=1.6e-11)                       # centered
    swl = ShallowWater(g, 9.81, 1000.0, f0=1e-4, beta=1.6e-11, tracer_limiter="vanleer")
    sc = balanced_vortex(g, swc, amp_eta=40.0, sigma=4e5, tracer=theta.copy())
    sl = balanced_vortex(g, swl, amp_eta=40.0, sigma=4e5, tracer=theta.copy())
    m0 = swl.tracer_mass(sl)
    dt = swc.max_dt(sc)
    for _ in range(int(2 * 2 * np.pi / swc.f0 / dt)):
        sc = swc.step(sc, dt)
        sl = swl.step(sl, dt)
    # the centered scheme overshoots the IC range [0, 2]; the limiter is strictly closer to it
    assert sc.tracer.max() > 2.0 + 1e-3
    assert sl.tracer.max() < sc.tracer.max()
    assert sl.tracer.min() > sc.tracer.min()
    assert swl.tracer_mass(sl) == pytest.approx(m0, rel=1e-11)


def test_limiter_variance_is_one_sided_dissipative():
    """A TVD limiter is dissipative — unlike the centered near-invariant, the tracer variance
    only *decreases* (one-sided). Asserted on the filamenting front."""
    g = uniform_grid(4e6, 4e6, 48, 48)
    sw = ShallowWater(g, 9.81, 1000.0, f0=1e-4, beta=1.6e-11, tracer_limiter="vanleer")
    X, _ = g.center_mesh()
    theta = np.tanh((X - g.Lx / 2) / 3e5)
    s = balanced_vortex(g, sw, amp_eta=40.0, sigma=4e5, tracer=theta)
    V0 = sw.tracer_variance(s)
    dt = sw.max_dt(s)
    for _ in range(int(2 * 2 * np.pi / sw.f0 / dt)):
        s = sw.step(s, dt)
    Vf = sw.tracer_variance(s)
    assert Vf <= V0 * (1.0 + 1e-9)             # only decreases (dissipative), within round-off
    assert Vf < V0 * 0.999                     # and it genuinely dissipates (not a no-op)


def test_vanleer_translates_smooth_blob_without_overshoot():
    """The gentle-limiter trade: van Leer carries a well-resolved Gaussian at the flow speed
    with NO overshoot and only mild peak clipping (it does not flatten the smooth extremum)."""
    ny = nx = 64
    L = 2.0e6
    g = uniform_grid(L, L, nx, ny)
    sw = ShallowWater(g, 9.81, 1000.0, f0=0.0, beta=0.0, tracer_limiter="vanleer")
    X, Y = g.center_mesh()
    U0 = 20.0
    sigma = L / 12
    cx, cy = L / 2, L / 2
    theta0 = np.exp(-(((X - cx) ** 2 + (Y - cy) ** 2) / (2 * sigma ** 2)))
    s = SWState(h=sw.H * np.ones((ny, nx)), u=U0 * np.ones((ny, nx)), v=np.zeros((ny, nx)), tracer=theta0.copy())
    dt = sw.max_dt(s) * 0.5
    n = int(0.2 * L / U0 / dt)
    for _ in range(n):
        s = sw.step(s, dt)
    assert s.tracer.max() <= theta0.max() + 1e-9      # no spurious overshoot
    assert s.tracer.min() >= -1e-9                     # stays positive
    assert s.tracer.max() > 0.8                        # peak only mildly clipped (gentle limiter)
