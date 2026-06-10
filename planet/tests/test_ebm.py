"""Planet Phase-1 validation: the latitudinal EBM machinery — the transport + Strang-split triad.

Carries the engine-coupled half of the plan §3 triad (planet-earth-system.md). The EBM reuses the
**shared diffusion spine** for the latitudinal transport and Strang-splits the radiation around it
(the Jominy idiom); these tests validate that assembly through the **interchangeable A/B/C trio**:

* **Analytical limit (tight).** (a) the **0-D** global mean ``T̄ = [(S₀/4)(1−ᾱ) − A]/B`` reproduced to
  *machine precision* (the transport is mean-preserving); (b) the **North (1975) two-mode** profile
  ``T₀ + T₂·P₂(x)`` reproduced near-exactly by the dt-free direct solve with exact faces, with the
  harmonic-mean polar floor named, and the Strang relaxation shown to **converge** to it as dt→0.
* **Conservation (tight).** net-TOA ``⟨S(1−α)⟩ − A − B⟨T⟩ = 0`` at equilibrium; pure transport
  preserves ``∫T dx`` exactly (the engine's no-flux invariant).
* **Pinning.** the exact-face construction reproduces the true face coefficient; the direct solve's
  reconstructed operator matches the engine's ``step`` (so mode C cannot drift from the engine).

The radiation/albedo constants are the cited climlab/North defaults ([[ebm-radiation-source]]), pinned
in :mod:`planet.ebm` — not carried from memory.
"""
import math

import numpy as np
import pytest
from scipy.linalg import solve_banded

from planet import ebm

S0 = ebm.S0_EARTH
ALPHA = ebm.ALBEDO_A0


def const_absorbed(alpha: float = ALPHA, S0: float = S0):
    """A constant-albedo (no-feedback) absorbed-shortwave callable — the linear North problem."""
    return lambda x, T: ebm.insolation(x, S0) * (1.0 - alpha)


def _cell_centers(n: int) -> np.ndarray:
    return (np.arange(n) + 0.5) / n               # uniform cell centers on [0,1] → mean = ∫₀¹·dx


# --------------------------------------------------------------------------- #
# Pinned constants ([[ebm-radiation-source]] — the climlab `EBM` defaults)
# --------------------------------------------------------------------------- #
def test_constants_are_the_cited_climlab_values():
    assert ebm.S0_EARTH == pytest.approx(1365.2)
    assert ebm.S2_INSOLATION == pytest.approx(-0.48)
    assert ebm.A_OLR == pytest.approx(210.0)
    assert ebm.B_OLR == pytest.approx(2.0)
    assert ebm.D_TRANSPORT == pytest.approx(0.555)
    assert ebm.T_FREEZE == pytest.approx(-10.0)
    assert ebm.ALBEDO_A0 == pytest.approx(0.30)
    assert ebm.ALBEDO_A2 == pytest.approx(0.078)
    assert ebm.ALBEDO_ICE == pytest.approx(0.62)


# --------------------------------------------------------------------------- #
# Insolation / P₂ geometry
# --------------------------------------------------------------------------- #
def test_legendre_p2_values_and_zero_mean():
    assert ebm.legendre_P2(0.0) == pytest.approx(-0.5)        # equator
    assert ebm.legendre_P2(1.0) == pytest.approx(1.0)         # pole
    x = _cell_centers(2000)
    assert float(np.mean(ebm.legendre_P2(x))) == pytest.approx(0.0, abs=1e-6)   # ∫₀¹ P₂ dx = 0


def test_insolation_global_mean_is_S0_over_4():
    x = _cell_centers(2000)
    assert float(np.mean(ebm.insolation(x))) == pytest.approx(S0 / 4.0, rel=1e-5)
    assert ebm.insolation(0.0) > ebm.insolation(1.0)         # more sun at the equator (s₂ < 0)


# --------------------------------------------------------------------------- #
# Analytical limit (a): the 0-D global mean (machine-exact)
# --------------------------------------------------------------------------- #
def test_zero_d_equilibrium_temperature_is_earthlike():
    assert ebm.equilibrium_temperature_0d() == pytest.approx(14.45, abs=0.05)


def test_relaxed_global_mean_matches_0d_to_the_grid_quadrature_limit():
    # The transport is mean-preserving, so the steady mean obeys the *discrete* energy balance
    # exactly (net-TOA ~ machine, asserted below); vs the *continuous* T₀ it differs only by the
    # O(1/n²) quadrature error of point-sampled insolation (∫₀¹ P₂ dx is not exactly 0 on a grid) —
    # ~5e-4 °C at n=120, and tightening with n.
    T0 = ebm.equilibrium_temperature_0d(S0, ALPHA)
    d_coarse = ebm.EnergyBalanceModel(n_cells=120).equilibrium(const_absorbed(), method="direct")
    d_fine = ebm.EnergyBalanceModel(n_cells=240).equilibrium(const_absorbed(), method="direct")
    assert d_coarse.global_mean_T == pytest.approx(T0, abs=2e-3)
    assert abs(d_fine.global_mean_T - T0) < abs(d_coarse.global_mean_T - T0)     # O(1/n²) convergence
    r = ebm.EnergyBalanceModel(n_cells=120).equilibrium(const_absorbed(), 10.0, n_tau=0.01)
    assert r.global_mean_T == pytest.approx(T0, abs=2e-3)


# --------------------------------------------------------------------------- #
# Analytical limit (b): the North two-mode profile across the A/B/C trio
# --------------------------------------------------------------------------- #
def test_north_two_mode_exact_direct_solve_is_tight():
    # mode C + exact faces: the dt-free linear solve reproduces North to near machine precision —
    # THE tight analytic profile anchor (no splitting error, no polar harmonic bias).
    m = ebm.EnergyBalanceModel(n_cells=180, face="exact")
    d = m.equilibrium(const_absorbed(), method="direct")
    ana = ebm.two_mode_solution(m.x, S0, ALPHA)
    assert np.max(np.abs(d.T - ana)) < 1e-3


def test_north_two_mode_harmonic_names_the_polar_floor():
    # mode C + harmonic faces: reproduces North but with the named ~0.1°C harmonic-mean polar bias
    # (the engine's face averaging where (1−x²)→0) — a bounded scope edge, not zero.
    m = ebm.EnergyBalanceModel(n_cells=180, face="harmonic")
    d = m.equilibrium(const_absorbed(), method="direct")
    err = np.max(np.abs(d.T - ebm.two_mode_solution(m.x, S0, ALPHA)))
    assert 0.05 < err < 0.30


def test_relaxation_converges_to_north_as_dt_shrinks():
    # The Strang-splitting steady-state error → 0 as dt → 0: shrinking n_tau tightens the North
    # match by ≳3× — the honest validation of the SPLITTING (the Jominy-idiom reuse).
    m = ebm.EnergyBalanceModel(n_cells=120, face="exact")
    ana = ebm.two_mode_solution(m.x, S0, ALPHA)
    err_coarse = np.max(np.abs(m.equilibrium(const_absorbed(), 10.0, n_tau=0.1).T - ana))
    err_fine = np.max(np.abs(m.equilibrium(const_absorbed(), 10.0, n_tau=0.01).T - ana))
    assert err_fine < 0.3 * err_coarse


def test_relaxation_converges_to_the_direct_steady_state():
    # The relaxation lands on the EXACT discrete steady state the direct solve computes — so the
    # splitting is validated against the engine's own discretization, independent of the analytic.
    m = ebm.EnergyBalanceModel(n_cells=120, face="exact")
    d = m.equilibrium(const_absorbed(), method="direct")
    r = m.equilibrium(const_absorbed(), 10.0, n_tau=0.005)
    assert np.max(np.abs(r.T - d.T)) < 0.1


# --------------------------------------------------------------------------- #
# Pinning: exact-face construction + the direct operator vs the engine
# --------------------------------------------------------------------------- #
def test_exact_face_helper_reproduces_the_true_face_coefficient():
    # The pre-distorted cell array's harmonic-mean faces equal the true (1−x²) face coefficient to
    # machine precision (and stay positive) — the construction that removes the polar floor.
    m = ebm.EnergyBalanceModel(n_cells=180, face="exact")
    edges = m.grid.edges
    coeff = lambda x: (m.D / m.C) * (1.0 - np.asarray(x, dtype=float) ** 2)
    Dc = m._Dcells
    Dface = 2.0 * Dc[:-1] * Dc[1:] / (Dc[:-1] + Dc[1:])
    assert np.allclose(Dface, coeff(edges[1:-1]), rtol=1e-12)
    assert np.all(Dc > 0.0)


def test_direct_transport_operator_matches_the_engine():
    # The reconstructed L_T is pinned to the engine: the engine's pure-transport backward-Euler step
    # equals solving (I − dt·L_T/C) — so mode C cannot silently drift from the engine.
    for face in ("harmonic", "exact"):
        m = ebm.EnergyBalanceModel(n_cells=90, face=face)
        sub, diag, sup = m._transport_tridiag()
        rng = np.random.default_rng(0)
        u = rng.standard_normal(90)
        n = 90
        for dt in (1e4, 1e6, 1e8):
            eng = m.solver.step(u, dt)
            ab = np.zeros((3, n))
            ab[0, 1:] = -dt * sup[:-1] / m.C
            ab[1, :] = 1.0 - dt * diag / m.C
            ab[2, :-1] = -dt * sub[1:] / m.C
            assert np.allclose(eng, solve_banded((1, 1), ab, u), atol=1e-10)


# --------------------------------------------------------------------------- #
# Conservation
# --------------------------------------------------------------------------- #
def test_net_toa_is_zero_at_equilibrium():
    # Global energy balance: absorbed solar = OLR at the steady state (machine-exact for the direct
    # solve; sub-mK for the relaxation). The transport only redistributes — it cannot change the mean.
    m = ebm.EnergyBalanceModel(n_cells=120)
    assert m.equilibrium(const_absorbed(), method="direct").net_toa == pytest.approx(0.0, abs=1e-9)
    assert m.equilibrium(const_absorbed(), 10.0, n_tau=0.01).net_toa == pytest.approx(0.0, abs=1e-3)


def test_pure_transport_conserves_total_under_no_flux():
    # ∫T dx is conserved exactly by the insulated (Neumann 0) transport — the engine's no-flux
    # invariant, re-confirmed for the Neumann/Neumann pole pair the EBM uses.
    m = ebm.EnergyBalanceModel(n_cells=120)
    T = 10.0 + 5.0 * np.cos(np.pi * m.x)
    total0 = m.solver.total(T)
    for _ in range(50):
        T = m.solver.step(T, 0.3 * m.tau_rad)
    assert m.solver.total(T) == pytest.approx(total0, rel=1e-9)


def test_equilibrium_independent_of_heat_capacity():
    # C (water_depth) sets only the relaxation timescale; the steady state is C-independent.
    a = ebm.EnergyBalanceModel(n_cells=120, water_depth=5.0).equilibrium(const_absorbed(), method="direct")
    b = ebm.EnergyBalanceModel(n_cells=120, water_depth=50.0).equilibrium(const_absorbed(), method="direct")
    assert np.allclose(a.T, b.T, atol=1e-9)


# --------------------------------------------------------------------------- #
# Ice-line diagnostic
# --------------------------------------------------------------------------- #
def test_ice_line_latitude_limits_and_crossing():
    x = _cell_centers(400)
    assert ebm.ice_line_latitude(x, np.full_like(x, 20.0)) == 90.0       # ice-free → cap at the pole
    assert ebm.ice_line_latitude(x, np.full_like(x, -50.0)) == 0.0       # Snowball → ice to the equator
    T = 20.0 - 60.0 * x                                                  # crosses −10 °C at x = 0.5
    assert ebm.ice_line_latitude(x, T) == pytest.approx(math.degrees(math.asin(0.5)), abs=1.0)


# --------------------------------------------------------------------------- #
# Input validation
# --------------------------------------------------------------------------- #
def test_invalid_inputs_raise():
    with pytest.raises(ValueError):
        ebm.EnergyBalanceModel(face="bogus")
    with pytest.raises(ValueError):
        ebm.EnergyBalanceModel(B=0.0)                                    # unstable relaxation
    m = ebm.EnergyBalanceModel(n_cells=60)
    with pytest.raises(ValueError):
        m.equilibrium(const_absorbed(), method="bogus")
    with pytest.raises(ValueError):
        m.equilibrium(const_absorbed(), method="relax")                 # relax needs T_init
    with pytest.raises(ValueError):
        # the ice feedback (state-dependent absorbed) cannot go through the direct linear solve
        m.steady_linear(lambda x, T: ebm.insolation(x, S0) * (1.0 - np.where(np.asarray(T) < -10.0, 0.6, 0.3)))
