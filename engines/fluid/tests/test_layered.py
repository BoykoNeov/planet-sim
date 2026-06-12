"""N-layer shallow-water engine — the rung-3 baroclinic extension of ``engines.fluid``.

:class:`LayeredShallowWater` stacks N free-surface SW layers coupled *only* through the
Montgomery pressure. What the suite guarantees (the honesty classes mirror the single-layer
engine, with the linear growth rate as the new tight anchor):

* **Single-layer reduction is byte-identical.** ``nl=1``, no background ⇒ the trajectory is
  ``np.array_equal`` to :class:`ShallowWater` over many steps (the by-construction rung-0
  reduction; the single-layer engine is untouched, so this is a *meaningful* cross-engine check).
* **Per-layer mass — machine precision** on the ``background=None`` path (flux-form continuity
  telescopes per layer, exactly like the single-layer mass).
* **Two-layer Poincaré dispersions** (no basic state): the engine reproduces the external
  ``√(gH_tot)`` and internal ``√(g'H_e)`` modes — the tight check that the Montgomery coupling is
  right (mirrors ``test_waves.py`` for two layers). *[slow — external-mode CFL integration.]*
* **Baroclinic growth rate — the tight anchor.** With a supercritical thermal-wind background, a
  small ``l=0`` perturbation at the analytic most-unstable wavenumber grows at the
  ``stability.TwoLayerStability`` rate within a few %, converging with resolution. *[slow.]*

Conservation is asserted **only** with ``background=None``: a background extracts available
potential energy (that is the instability), so the perturbation energy *grows* — the signal, not
a drift to bound.
"""
import numpy as np
import pytest

from engines.fluid import (
    ShallowWater,
    SWState,
    LayeredShallowWater,
    LayeredState,
    TwoLayerStability,
    uniform_grid,
)


# Idealized rung-3 parameters (the spike's).
P = dict(f0=1.0e-4, g=10.0, gp=0.2, H1=500.0, H2=500.0)


# --------------------------------------------------------------------------- #
# Single-layer reduction — byte-identical to ShallowWater (the rung-0 reduction)
# --------------------------------------------------------------------------- #
def test_single_layer_reduction_is_bit_for_bit():
    """nl=1 with no background reproduces ShallowWater to the bit — the structural extension did
    not perturb the single-layer trajectory by a single ULP (M₁ = g(h_b+h) is bit-commutative
    with the single-layer Bernoulli pressure; a length-1 leading axis is the same float ops)."""
    ny, nx = 24, 32
    grid = uniform_grid(4e6, 4e6, nx, ny)
    rng = np.random.default_rng(0)
    H = 1000.0
    h2d = H + 20.0 * rng.standard_normal((ny, nx))
    u2d = 0.5 * rng.standard_normal((ny, nx))
    v2d = 0.5 * rng.standard_normal((ny, nx))
    sw = ShallowWater(grid, 9.81, H, f0=1e-4, beta=1.6e-11)
    lay = LayeredShallowWater(grid, 9.81, [H], [], f0=1e-4, beta=1.6e-11)
    s_sw = SWState(h=h2d.copy(), u=u2d.copy(), v=v2d.copy())
    s_la = LayeredState(h=h2d[None].copy(), u=u2d[None].copy(), v=v2d[None].copy())
    dt = 0.3 * min(grid.dx, grid.dy) / np.sqrt(9.81 * H)
    for _ in range(60):
        s_sw = sw.step(s_sw, dt)
        s_la = lay.step(s_la, dt)
        assert np.array_equal(s_sw.h, s_la.h[0])
        assert np.array_equal(s_sw.u, s_la.u[0])
        assert np.array_equal(s_sw.v, s_la.v[0])


def test_single_layer_reduction_with_topography_bit_for_bit():
    """The Montgomery h_b enters every layer; for nl=1 it must still reduce to g(h+h_b) exactly."""
    ny, nx = 16, 20
    grid = uniform_grid(2e6, 2e6, nx, ny)
    rng = np.random.default_rng(1)
    H = 800.0
    X, Y = grid.center_mesh()
    h_b = 50.0 * np.sin(2 * np.pi * X / grid.Lx)
    h2d = H + 10.0 * rng.standard_normal((ny, nx))
    z = np.zeros((ny, nx))
    sw = ShallowWater(grid, 9.81, H, f0=1e-4, bottom=h_b)
    lay = LayeredShallowWater(grid, 9.81, [H], [], f0=1e-4, bottom=h_b)
    s_sw = SWState(h=h2d.copy(), u=z.copy(), v=z.copy())
    s_la = LayeredState(h=h2d[None].copy(), u=z[None].copy(), v=z[None].copy())
    dt = 0.3 * min(grid.dx, grid.dy) / np.sqrt(9.81 * H)
    for _ in range(40):
        s_sw = sw.step(s_sw, dt)
        s_la = lay.step(s_la, dt)
        assert np.array_equal(s_sw.h, s_la.h[0])
        assert np.array_equal(s_sw.u, s_la.u[0])
        assert np.array_equal(s_sw.v, s_la.v[0])


# --------------------------------------------------------------------------- #
# Per-layer mass — machine precision (background=None, flux-form telescoping)
# --------------------------------------------------------------------------- #
def test_per_layer_mass_conserved_to_machine_precision():
    ny, nx = 32, 40
    grid = uniform_grid(2e6, 2e6, nx, ny)
    lay = LayeredShallowWater(grid, P["g"], [P["H1"], P["H2"]], [P["gp"]], f0=P["f0"])
    rng = np.random.default_rng(2)
    h = np.stack([P["H1"] + 5.0 * rng.standard_normal((ny, nx)),
                  P["H2"] + 5.0 * rng.standard_normal((ny, nx))])
    z = np.zeros((2, ny, nx))
    s = LayeredState(h=h, u=0.2 * rng.standard_normal((2, ny, nx)), v=z.copy())
    m0 = lay.layer_mass(s)
    dt = lay.max_dt(s)
    for _ in range(200):
        s = lay.step(s, dt)
    assert np.allclose(lay.layer_mass(s), m0, rtol=1e-12, atol=0.0)


# --------------------------------------------------------------------------- #
# thermal_wind injects exactly the analytic G_k (the engine ↔ anchor cross-link)
# --------------------------------------------------------------------------- #
def test_thermal_wind_matches_analytic_gradients():
    grid = uniform_grid(1e6, 1e6, 8, 8)
    lay = LayeredShallowWater(grid, P["g"], [P["H1"], P["H2"]], [P["gp"]], f0=P["f0"])
    st = lay.stability()
    Us = 4.0
    bg = lay.thermal_wind([0.5 * Us, -0.5 * Us])
    G1, G2 = st.basic_state_gradients(0.5 * Us, -0.5 * Us)
    assert bg.G[0] == pytest.approx(G1, rel=1e-12)
    assert bg.G[1] == pytest.approx(G2, rel=1e-12)
    assert np.array_equal(bg.U, np.array([0.5 * Us, -0.5 * Us]))


# --------------------------------------------------------------------------- #
# Construction guards
# --------------------------------------------------------------------------- #
def test_construction_validation():
    grid = uniform_grid(1e6, 1e6, 8, 8)
    with pytest.raises(ValueError):                       # nl-1 reduced gravities required
        LayeredShallowWater(grid, 10.0, [500.0, 500.0], [], f0=1e-4)
    with pytest.raises(ValueError):                       # unstable stratification
        LayeredShallowWater(grid, 10.0, [500.0, 500.0], [-0.2], f0=1e-4)
    with pytest.raises(ValueError):                       # positive thickness
        LayeredShallowWater(grid, 10.0, [500.0, -1.0], [0.2], f0=1e-4)
    lay = LayeredShallowWater(grid, 10.0, [500.0, 500.0], [0.2], f0=1e-4)
    s = LayeredState(h=np.full((2, 8, 8), 500.0), u=np.zeros((2, 8, 8)), v=np.zeros((2, 8, 8)))
    with pytest.raises(ValueError):                       # non-positive dt
        lay.step(s, 0.0)
    with pytest.raises(ValueError):                       # dt above the external-mode CFL
        lay.step(s, 1e9)
    # a background on a β-plane is inconsistent (the mean Coriolis would be unbalanced) → rejected
    bg = lay.thermal_wind([2.0, -2.0])
    with pytest.raises(ValueError):
        LayeredShallowWater(grid, 10.0, [500.0, 500.0], [0.2], f0=1e-4, beta=1e-11, background=bg)


# --------------------------------------------------------------------------- #
# Two-layer Poincaré dispersions reproduced by the engine (no basic state) [slow]
# --------------------------------------------------------------------------- #
def _measure_omega(lay, s0, project, dt, n):
    series, s = [], s0
    for _ in range(n):
        s = lay.step(s, dt)
        series.append(project(s))
    series = np.asarray(series) - np.mean(series)
    freqs = np.fft.rfftfreq(len(series), d=dt)
    return 2 * np.pi * freqs[np.argmax(np.abs(np.fft.rfft(series)))]


@pytest.mark.slow
@pytest.mark.parametrize("mode", ["external", "internal"])
def test_two_layer_poincare_dispersion(mode):
    """A small height mode in phase (external, √(gH_tot)) or in anti-phase (internal, √(g'H_e))
    oscillates at the two-layer Poincaré frequency — the tight check on the Montgomery coupling."""
    ny, nx = 4, 64
    Lx = 2e6
    grid = uniform_grid(Lx, Lx * ny / nx, nx, ny)
    lay = LayeredShallowWater(grid, P["g"], [P["H1"], P["H2"]], [P["gp"]], f0=P["f0"])
    st = lay.stability()
    k = 2 * np.pi / Lx
    x = grid.x_centers()
    seed = 1e-3 * np.cos(k * x)[None, :] * np.ones((ny, 1))
    z = np.zeros((ny, nx))
    sign = +1.0 if mode == "external" else -1.0
    om_an = st.poincare_external(k) if mode == "external" else st.poincare_internal(k)
    s0 = LayeredState(h=np.array([P["H1"] + seed, P["H2"] + sign * seed]),
                      u=np.array([z.copy(), z.copy()]), v=np.array([z.copy(), z.copy()]))
    dt = lay.max_dt(s0) * 0.5
    n = int(3 * 2 * np.pi / om_an / dt)
    om = _measure_omega(lay, s0, lambda s: np.mean((s.h[0] - P["H1"]) * np.cos(k * x)), dt, n)
    assert om == pytest.approx(om_an, rel=5e-3)


# --------------------------------------------------------------------------- #
# Baroclinic growth rate vs the analytic two-layer SW rate — the tight anchor [slow]
# --------------------------------------------------------------------------- #
def _measure_growth(nx, ny=4, n_efold=4.0, Us=4.0):
    """Integrate the most-unstable l=0 mode on the nonlinear engine; return (σ_meas, σ_analytic)."""
    st = TwoLayerStability(**P)
    kstar, sig_an = st.most_unstable(Us)
    Lx = 2 * np.pi / kstar
    grid = uniform_grid(Lx, Lx * ny / nx, nx, ny)
    lay = LayeredShallowWater(grid, P["g"], [P["H1"], P["H2"]], [P["gp"]], f0=P["f0"], beta=0.0)
    lay = LayeredShallowWater(grid, P["g"], [P["H1"], P["H2"]], [P["gp"]], f0=P["f0"], beta=0.0,
                              background=lay.thermal_wind([0.5 * Us, -0.5 * Us]))
    x = grid.x_centers()
    seed = 1e-3 * np.cos(kstar * x)[None, :] * np.ones((ny, 1))
    z = np.zeros((ny, nx))
    s = LayeredState(h=np.array([P["H1"] + seed, P["H2"] - seed]),
                     u=np.array([z.copy(), z.copy()]), v=np.array([z.copy(), z.copy()]))
    dt = lay.max_dt(s, safety=0.3)
    nsteps = int(np.ceil((n_efold / sig_an) / dt))
    ts, lnE = [], []
    for n in range(nsteps + 1):
        if n % max(1, nsteps // 200) == 0:
            ts.append(n * dt)
            lnE.append(0.5 * np.log(lay.perturbation_energy(s) + 1e-300))
        s = lay.step(s, dt)
    ts, lnE = np.array(ts), np.array(lnE)
    i0, i1 = len(ts) // 3, 9 * len(ts) // 10          # the clean exponential window
    sig_meas = float(np.polyfit(ts[i0:i1], lnE[i0:i1], 1)[0])
    return sig_meas, sig_an


@pytest.mark.slow
def test_baroclinic_growth_rate_matches_analytic():
    sig_meas, sig_an = _measure_growth(nx=64)
    assert sig_meas == pytest.approx(sig_an, rel=0.06)    # within a few % (the spike's ~4%)


@pytest.mark.slow
def test_growth_rate_converges_with_resolution():
    """The measured rate approaches the analytic σ as the grid refines (the discretization is the
    error, not a modelling gap) — the route's soundness, not just a single lucky resolution."""
    m32, an = _measure_growth(nx=32)
    m64, _ = _measure_growth(nx=64)
    err32 = abs(m32 - an) / an
    err64 = abs(m64 - an) / an
    assert err64 < err32
    assert err64 < 0.06
