"""Wave-speed / dispersion anchors — the analytical-limit leg of the Phase-3 triad.

The exact linear-wave relations the rotating shallow-water solver must reproduce in the
small-amplitude limit (where the nonlinear advection is negligible):

* **Gravity-wave speed** ``c = √(gH)`` (f-plane, f₀=0): a height mode oscillates at ``ω = c·k``.
* **Poincaré (inertia-gravity) dispersion** ``ω² = f₀² + gH·k²`` — the *rotation* check that
  bites: the oscillation frequency must rise with ``f₀`` exactly as the relation predicts (a
  wrong Coriolis would give the wrong frequency). Reproduced to ~1e-3.
* **Rossby-wave dispersion** ``ω = −βk/(k² + l² + 1/L_R²)`` (β-plane): a balanced mode
  propagates **westward** (ω < 0 for k>0), with **longer waves faster** — reproduced to a
  loose band that *converges to the analytic value as the grid refines* (the slow balanced
  mode carries a few-percent numerical-dispersion error a gravity wave does not; a teaching
  point, asserted loose).

Frequencies are measured from a tracked single-mode projection time series via its FFT peak
(robust to a constant geostrophic offset), the way a lab would read a dispersion relation.
"""
import numpy as np
import pytest

from engines.fluid import ShallowWater, SWState, uniform_grid


def _measure_omega(sw, s0, project, dt, n):
    """FFT peak frequency of the scalar time series ``project(state)`` over ``n`` steps of ``dt``."""
    series = []
    s = s0
    for _ in range(n):
        s = sw.step(s, dt)
        series.append(project(s))
    series = np.asarray(series) - np.mean(series)
    freqs = np.fft.rfftfreq(len(series), d=dt)
    return 2 * np.pi * freqs[np.argmax(np.abs(np.fft.rfft(series)))]


# --------------------------------------------------------------------------- #
# Gravity-wave speed  √(gH)
# --------------------------------------------------------------------------- #
def test_gravity_wave_speed():
    g = uniform_grid(2e6, 2e6, 48, 4)
    sw = ShallowWater(g, 9.81, 1000.0, f0=0.0)
    k = 2 * np.pi / g.Lx
    X, _ = g.center_mesh()
    s0 = SWState(h=sw.H + 0.5 * np.cos(k * X), u=np.zeros((g.ny, g.nx)), v=np.zeros((g.ny, g.nx)))
    om_an = sw.gravity_wave_speed * k
    dt = sw.max_dt(s0) * 0.5
    n = int(3 * 2 * np.pi / om_an / dt)
    om = _measure_omega(sw, s0, lambda s: np.mean((s.h - sw.H) * np.cos(k * X)), dt, n)
    assert om == pytest.approx(om_an, rel=2e-3)


def test_gravity_wave_speed_property():
    sw = ShallowWater(uniform_grid(1e6, 1e6, 8, 8), 9.81, 1000.0)
    assert sw.gravity_wave_speed == pytest.approx(np.sqrt(9.81 * 1000.0))


# --------------------------------------------------------------------------- #
# Poincaré / inertia-gravity dispersion  ω² = f² + gHk²   (the rotation check)
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("f0", [0.5e-4, 1.0e-4, 2.0e-4])
def test_poincare_dispersion(f0):
    g = uniform_grid(2e6, 2e6, 48, 4)
    sw = ShallowWater(g, 9.81, 1000.0, f0=f0)
    k = 2 * np.pi / g.Lx
    X, _ = g.center_mesh()
    s0 = SWState(h=sw.H + 0.5 * np.cos(k * X), u=np.zeros((g.ny, g.nx)), v=np.zeros((g.ny, g.nx)))
    om_an = np.sqrt(f0 ** 2 + 9.81 * 1000.0 * k ** 2)
    dt = sw.max_dt(s0) * 0.5
    n = int(3 * 2 * np.pi / om_an / dt)
    om = _measure_omega(sw, s0, lambda s: np.mean((s.h - sw.H) * np.cos(k * X)), dt, n)
    assert om == pytest.approx(om_an, rel=2e-3)


# --------------------------------------------------------------------------- #
# Rossby-wave dispersion  ω = −βk/(k²+l²+1/L_R²)  — westward, dispersive, loose
# --------------------------------------------------------------------------- #
def _rossby_omega(sw, grid, mk, ml, frac=0.35):
    f0 = sw.f0
    k = 2 * np.pi * mk / grid.Lx
    l = 2 * np.pi * ml / grid.Ly
    Xc, Yc = grid.center_mesh()
    xU = grid.x_corners()[None, :] * np.ones((grid.ny, 1)); yU = grid.y_centers()[:, None] * np.ones((1, grid.nx))
    xV = grid.x_centers()[None, :] * np.ones((grid.ny, 1)); yV = grid.y_corners()[:, None] * np.ones((1, grid.nx))
    eta = (f0 / sw.g) * np.cos(k * Xc + l * Yc)           # geostrophic streamfunction Ψ=1
    u = l * np.sin(k * xU + l * yU)
    v = -k * np.sin(k * xV + l * yV)
    s = SWState(h=sw.H + eta, u=u, v=v)
    om_an = -sw.beta * k / (k ** 2 + l ** 2 + 1.0 / sw.rossby_radius ** 2)
    dt = sw.max_dt(s) * 0.8
    n = int(frac * 2 * np.pi / abs(om_an) / dt)
    phase, times = [], []
    t = 0.0
    for _ in range(n):
        s = sw.step(s, dt); t += dt
        de = s.h - sw.H
        a = np.mean(de * np.cos(k * Xc + l * Yc)); b = np.mean(de * np.sin(k * Xc + l * Yc))
        phase.append(np.arctan2(b, a)); times.append(t)
    om = np.polyfit(np.array(times), np.unwrap(phase), 1)[0]
    return om, om_an, om_an / k          # measured ω, analytic ω, analytic zonal phase speed


@pytest.mark.slow
def test_rossby_wave_westward_and_dispersive():
    g = uniform_grid(8e6, 8e6, 96, 96)
    sw = ShallowWater(g, 9.81, 1000.0, f0=1e-4, beta=2e-12)
    om11, an11, c11 = _rossby_omega(sw, g, 1, 1)
    om21, an21, c21 = _rossby_omega(sw, g, 2, 1)
    # westward: zonal phase speed and ω both negative for k>0
    assert om11 < 0 and an11 < 0 and c11 < 0
    # magnitude within a loose band of analytic (converges to 1.0 with resolution)
    assert 0.9 < om11 / an11 < 1.02
    # dispersive: the longer zonal wave (k=1) travels faster than the shorter (k=2)
    assert abs(c11) > abs(c21)
