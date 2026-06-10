"""Solver mechanics — the frozen data boundary, construction guards, and the explicit-CFL contract.

Not physics (the wave / geostrophy / conservation files carry the triad) but the API
promises: construction validation, the no-mutation step contract, the CFL stability guard
(the explicit analogue of the diffusion engine's unconditional-stability promise — here
*conditional*, so enforced), and the **passively-advected tracer slot** (the rung-1 extension,
ADR 0005 — stepping it advects the tracer; the full triad is `test_tracer.py`).
"""
import numpy as np
import pytest

from engines.fluid import Grid2D, ShallowWater, SWState, uniform_grid


# --------------------------------------------------------------------------- #
# Grid
# --------------------------------------------------------------------------- #
def test_uniform_grid_geometry():
    g = uniform_grid(4e6, 2e6, 80, 40)
    assert g.nx == 80 and g.ny == 40
    assert g.dx == pytest.approx(5e4) and g.dy == pytest.approx(5e4)
    assert g.Lx == pytest.approx(4e6) and g.Ly == pytest.approx(2e6)
    assert g.cell_area == pytest.approx(2.5e9)
    assert g.x_centers().shape == (80,) and g.y_centers().shape == (40,)
    assert g.x_centers()[0] == pytest.approx(0.5 * g.dx)        # cell-centred
    assert g.x_corners()[0] == pytest.approx(0.0)               # corner at the origin


def test_uniform_grid_rejects_tiny():
    with pytest.raises(ValueError):
        uniform_grid(1.0, 1.0, 1, 8)


# --------------------------------------------------------------------------- #
# Construction guards
# --------------------------------------------------------------------------- #
def test_construction_rejects_bad_params():
    g = uniform_grid(1e6, 1e6, 16, 16)
    with pytest.raises(ValueError):
        ShallowWater(g, g=-1.0, mean_depth=100.0)
    with pytest.raises(ValueError):
        ShallowWater(g, g=9.81, mean_depth=0.0)
    with pytest.raises(ValueError):
        ShallowWater(g, 9.81, 100.0, bottom=np.zeros((3, 3)))   # wrong topography shape


def test_derived_scales():
    sw = ShallowWater(uniform_grid(1e6, 1e6, 16, 16), 9.81, 1000.0, f0=1e-4)
    assert sw.gravity_wave_speed == pytest.approx(np.sqrt(9.81e3))
    assert sw.rossby_radius == pytest.approx(np.sqrt(9.81e3) / 1e-4)
    assert sw.f_at_corners().shape == (16, 16)
    # f-plane (f0=0) → infinite deformation radius
    assert ShallowWater(uniform_grid(1e6, 1e6, 8, 8), 9.81, 1000.0, f0=0.0).rossby_radius == np.inf


def test_beta_plane_coriolis_field():
    g = uniform_grid(1e6, 2e6, 8, 16)
    sw = ShallowWater(g, 9.81, 1000.0, f0=1e-4, beta=1e-11)   # y_ref defaults to domain centre
    f = sw.f_at_corners()
    # f varies only in y, linearly, equal to f0 at the domain centre
    assert np.allclose(f, f[:, [0]])                          # no x-variation
    y = g.y_corners()
    assert np.allclose(f[:, 0], 1e-4 + 1e-11 * (y - 0.5 * g.Ly))


# --------------------------------------------------------------------------- #
# The step contract — no mutation, CFL guard, tracer slot
# --------------------------------------------------------------------------- #
def _simple_state(g, sw):
    X, Y = g.center_mesh()
    eta = 1.0 * np.exp(-(((X - g.Lx / 2) ** 2 + (Y - g.Ly / 2) ** 2)) / (2 * (g.Lx / 8) ** 2))
    return SWState(h=sw.H + eta, u=np.zeros((g.ny, g.nx)), v=np.zeros((g.ny, g.nx)))


def test_step_does_not_mutate_input():
    g = uniform_grid(2e6, 2e6, 32, 32)
    sw = ShallowWater(g, 9.81, 1000.0, f0=1e-4)
    s = _simple_state(g, sw)
    h_before = s.h.copy()
    sw.step(s, sw.max_dt(s))
    assert np.array_equal(s.h, h_before)                      # input untouched


def test_step_rejects_nonpositive_dt():
    g = uniform_grid(1e6, 1e6, 16, 16)
    sw = ShallowWater(g, 9.81, 1000.0)
    with pytest.raises(ValueError):
        sw.step(_simple_state(g, sw), 0.0)


def test_step_rejects_dt_above_cfl():
    g = uniform_grid(2e6, 2e6, 32, 32)
    sw = ShallowWater(g, 9.81, 1000.0, f0=1e-4)
    s = _simple_state(g, sw)
    too_big = 10.0 * sw.max_dt(s, safety=1.0)
    with pytest.raises(ValueError, match="CFL"):
        sw.step(s, too_big)


def test_recommended_step_is_well_under_the_cfl_limit():
    g = uniform_grid(2e6, 2e6, 32, 32)
    sw = ShallowWater(g, 9.81, 1000.0, f0=1e-4)
    s = _simple_state(g, sw)
    assert sw.max_dt(s, safety=0.3) < sw.max_dt(s, safety=1.0)
    sw.step(s, sw.max_dt(s))                                  # the recommended step never trips the guard


def test_tracer_slot_is_advected():
    """The rung-1 extension (ADR 0005): a set tracer is advected (no longer raises). The full
    passive-tracer triad — ∫hθ conservation, consistency, bit-for-bit passivity — is test_tracer.py."""
    g = uniform_grid(1e6, 1e6, 16, 16)
    sw = ShallowWater(g, 9.81, 1000.0, f0=1e-4)
    s = _simple_state(g, sw)
    with_tracer = SWState(h=s.h, u=s.u, v=s.v, tracer=np.ones((g.ny, g.nx)))
    stepped = sw.step(with_tracer, sw.max_dt(s))                # advects the tracer; does not raise
    assert stepped.tracer is not None and stepped.tracer.shape == (g.ny, g.nx)


def test_solve_marches_to_t_end():
    g = uniform_grid(2e6, 2e6, 32, 32)
    sw = ShallowWater(g, 9.81, 1000.0, f0=1e-4)
    s0 = _simple_state(g, sw)
    s = sw.solve(s0, t_end=3600.0)                           # one hour
    assert s.h.shape == (g.ny, g.nx)
    assert sw.mass(s) == pytest.approx(sw.mass(s0), rel=1e-12)
    assert not np.array_equal(s.h, s0.h)                      # something actually happened
