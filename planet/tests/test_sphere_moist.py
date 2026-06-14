"""Triad for the full-sphere moisture budget — emergent ITCZ rain co-located with the EFE
(:mod:`planet.sphere_moist`, rung 2.x).

*Tight (structural)* — the cross-model reduction to the hemisphere :mod:`planet.moist` at ``φ_EFE = 0``
(machine-exact), ``∫(P − E) dx = 0`` (machine-exact, even for an asymmetric displaced cell), and a
symmetric climate ⟹ a symmetric ``P − E`` peaking at the equator. *Real-but-loose (the unlock)* — the net
``P − E`` peaks **on** the EFE because the prescribed Hadley cell beats the down-gradient eddy export there
(the falsifiable check, ~2–3× margin — not assumed), and the ITCZ convergence intensifies at the ~C–C rate
under warming (transported from the Hadley fix). *Plumbing / named negatives* — co-location *at all* is
by-construction; the displaced-peak intensification is **geometric** (pinned-edge cell narrowing), **not**
emergent ``q`` (the clean negative result, pinned so it is not silently re-read as a win); the wet/dry
dipole is displacement-driven with a by-construction direction; ``strength = 0`` reduces to the eddy-only
budget.
"""
from dataclasses import replace

import numpy as np
import pytest

from planet import ebm, moist
from planet import sphere_moist as sm
from planet.ebm import ClimateState
from planet.sphere_ebm import SphereEBM


def _absorbed_const(x, T):
    """Constant-albedo absorbed shortwave — the splitting-free linear reference forcing."""
    return ebm.insolation(x) * (1.0 - ebm.ALBEDO_A0)


def _sym_climate(n_cells=360):
    m = SphereEBM(n_cells=n_cells)
    return m, m.steady_linear(_absorbed_const)


def _displaced_climate(q_drive, n_cells=360):
    m = SphereEBM(n_cells=n_cells)
    return m, m.steady_linear(_absorbed_const, Q=q_drive * m.x)


# --------------------------------------------------------------------------- #
# TIGHT — structural anchors.
# --------------------------------------------------------------------------- #
def test_reduces_to_hemisphere_moist_budget_at_efe0():
    # The cross-model reduction (moist.py untouched): at φ_EFE = 0 the full-sphere operators reproduce the
    # hemisphere moist.moisture_budget on the NH to machine precision (a symmetric q zeroes the equatorial
    # face-flux, collapsing the full-sphere stencil to the hemisphere's equatorial-symmetry boundary).
    m, c = _sym_climate()
    assert c.phi_efe == pytest.approx(0.0, abs=1e-6)
    nh = m.x > 0
    hemi = ClimateState(x=m.x[nh], T=c.T[nh], global_mean_T=c.global_mean_T,
                        ice_line_lat=float("nan"), net_toa=0.0, converged=True, iterations=0)
    b = sm.sphere_moisture_budget(c, hadley=True)
    # eddy-only and full both reduce
    eddy_only = sm.sphere_moisture_budget(c, hadley=False)
    assert np.max(np.abs(eddy_only.p_minus_e[nh] - moist.moisture_convergence(hemi))) < 1e-9
    assert np.max(np.abs(b.p_minus_e[nh] - moist.moisture_budget(hemi, hadley=True).p_minus_e)) < 1e-9


def test_conservation_machine_exact_even_when_displaced():
    # ∫(P − E) dx = 0 to machine precision — symmetric AND displaced (the conservative face form holds for
    # the asymmetric two-cell circulation; the hemispheric integrals are equal and opposite).
    _, c0 = _sym_climate()
    b0 = sm.sphere_moisture_budget(c0)
    assert abs(b0.net_p_minus_e) < 1e-9
    for q_drive in (2.0, 4.0, 6.0):
        _, c = _displaced_climate(q_drive)
        b = sm.sphere_moisture_budget(c)
        assert abs(b.net_p_minus_e) < 1e-9
        assert b.nh_p_minus_e == pytest.approx(-b.sh_p_minus_e, abs=1e-9)   # equal-and-opposite dipole


def test_symmetric_climate_has_symmetric_budget_peaking_at_equator():
    # A symmetric climate ⟹ an even P − E, the rain belt on the equator, no hemispheric dipole.
    m, c = _sym_climate()
    b = sm.sphere_moisture_budget(c)
    assert np.max(np.abs(b.p_minus_e - b.p_minus_e[::-1])) < 1e-6     # even (to the linear solve's symmetry)
    assert abs(b.rain_max_lat) < 1.0                                  # ITCZ on the equator (innermost cell)
    assert b.nh_p_minus_e == pytest.approx(0.0, abs=1e-9)
    assert b.sh_p_minus_e == pytest.approx(0.0, abs=1e-9)


def test_streamfunction_anchored_at_efe():
    # Ψ vanishes at the EFE and at the fixed ±edge descents, flips sign across the EFE (two cells), and is
    # zero poleward of the edges — the plumbing of the migrating two-cell circulation.
    m = SphereEBM(n_cells=360)
    efe = 6.0
    psi = sm.hadley_streamfunction(m.x, phi_efe=efe, edge_deg=30.0)
    x_efe = np.sin(np.radians(efe))
    i_efe = int(np.argmin(np.abs(m.x - x_efe)))
    assert abs(psi[i_efe]) < 0.05                                     # ~0 at the ascent
    south = (m.x < x_efe) & (m.x > np.sin(np.radians(-30.0)))
    north = (m.x > x_efe) & (m.x < np.sin(np.radians(30.0)))
    assert np.all(psi[south] >= 0.0) and np.max(psi[south]) > 0.5     # northward low-level flow south of EFE
    assert np.all(psi[north] <= 0.0) and np.min(psi[north]) < -0.5    # southward low-level flow north of EFE
    assert np.all(psi[np.abs(m.x) > np.sin(np.radians(30.0))] == 0.0) # vanishes poleward of the descents


# --------------------------------------------------------------------------- #
# REAL-BUT-LOOSE (the unlock) — the falsifiable co-location + the warming response.
# --------------------------------------------------------------------------- #
def test_net_pme_co_locates_with_efe_hadley_beating_the_eddy_export():
    # THE falsifiable check: the down-gradient eddy term EXPORTS moisture from the warm EFE (backwards), so
    # the net P − E peaks on the EFE only because the prescribed cell BEATS that export at the displaced
    # latitude — not guaranteed a priori, true (~2–3× margin) for the calibrated strength.
    for q_drive in (2.0, 4.0, 6.0):
        m, c = _displaced_climate(q_drive)
        b = sm.sphere_moisture_budget(c)
        i_efe = int(np.argmin(np.abs(b.phi - b.phi_efe)))
        assert b.phi_efe > 0.3                                        # displaced into the warmed hemisphere
        assert b.p_minus_e_eddy[i_efe] < 0.0                          # eddy is backwards (export) at the EFE
        assert b.p_minus_e_hadley[i_efe] > -b.p_minus_e_eddy[i_efe]   # Hadley beats the export
        assert b.p_minus_e_hadley[i_efe] > 2.0 * abs(b.p_minus_e_eddy[i_efe])   # by a comfortable margin
        # The rain belt TRACKS the EFE — null-rejecting (not an absolute degree tol, which is vacuous at a
        # small shift): the net rain max is NEARER the EFE than the equator (rejects "stayed at the equator")
        # and tracks it without wild overshoot (a small geometric poleward skew from the narrowing near-cell).
        assert abs(b.rain_max_lat - b.phi_efe) < abs(b.rain_max_lat)
        assert 0.8 < b.rain_max_lat / b.phi_efe < 1.7


def test_warming_intensifies_itcz_at_cc_rate():
    # The clean emergent q(T) signature (transported from the Hadley fix): warm the climate, hold the cell
    # strength FIXED — only q(T) moves — and the ITCZ convergence intensifies at the ~C–C moisture rate,
    # faster than the energy-constrained global mean.
    m, c = _sym_climate()
    i0 = int(np.argmax(sm.sphere_moisture_budget(c).p_minus_e_hadley))
    dT = 4.0
    base = sm.hadley_moisture_convergence(c)[i0]
    warm = sm.hadley_moisture_convergence(replace(c, T=c.T + dT, global_mean_T=c.global_mean_T + dT))[i0]
    itcz_rate = (warm / base - 1.0) / dT
    cc_rate = moist.L_VAPOR / (moist.R_VAPOR * (c.T[i0] + moist.T0_KELVIN) ** 2)
    assert itcz_rate == pytest.approx(cc_rate, rel=0.2)               # intensifies at the ~C–C rate
    assert itcz_rate > moist.energy_constrained_rate()               # faster than the energy-constrained mean


# --------------------------------------------------------------------------- #
# PLUMBING / NAMED NEGATIVES — what is by-construction, and the clean negative result.
# --------------------------------------------------------------------------- #
def test_displaced_peak_intensification_is_geometric_not_q():
    # THE CLEAN NEGATIVE RESULT (pinned so it is not silently re-read as a win): the displaced-ITCZ peak
    # grows because the pinned-edge near-side cell NARROWS (geometry), NOT because the warm hemisphere is
    # moister. Holding the cell position fixed and replacing q with its hemispheric symmetrization leaves
    # the ITCZ peak essentially unchanged — so the meridional q-asymmetry does not set the intensity.
    m, c = _displaced_climate(6.0)
    b_real = sm.sphere_moisture_budget(c)
    c_symq = replace(c, T=0.5 * (c.T + c.T[::-1]))                    # symmetric q, SAME phi_efe (isolation)
    b_symq = sm.sphere_moisture_budget(c_symq)
    assert b_symq.phi_efe == b_real.phi_efe                           # cell anchored at the same EFE
    assert b_symq.p_minus_e.max() == pytest.approx(b_real.p_minus_e.max(), rel=0.01)
    # and the peak IS above the symmetric-climate peak — but that growth is the geometry, not q
    assert b_real.p_minus_e.max() > sm.sphere_moisture_budget(_sym_climate()[1]).p_minus_e.max()


def test_wet_dry_dipole_points_toward_the_warm_hemisphere():
    # The displacement-driven wet/dry dipole: the displaced ITCZ puts net convergence in the warm (N)
    # hemisphere and the compensating descent in the cold (S) one. Direction is by-construction (the EFE is
    # defined to move toward the warm side); magnitude grows with the displacement.
    asyms = []
    for q_drive in (2.0, 4.0, 6.0):
        _, c = _displaced_climate(q_drive)
        b = sm.sphere_moisture_budget(c)
        assert b.nh_p_minus_e > 0.0 and b.sh_p_minus_e < 0.0
        asyms.append(b.nh_p_minus_e)
    assert asyms[0] < asyms[1] < asyms[2]                            # grows with the imposed displacement


def test_zero_strength_reduces_to_eddy_only():
    # The Hadley term is an independent diff against the eddy default — strength = 0 recovers it bit-for-bit.
    _, c = _displaced_climate(4.0)
    b0 = sm.sphere_moisture_budget(c, strength=0.0)
    eddy = sm.sphere_moisture_budget(c, hadley=False)
    assert np.array_equal(b0.p_minus_e, eddy.p_minus_e)


def test_eddy_only_default_is_backwards_at_the_itcz():
    # The rung-2 trade carried onto the full sphere: pure eddy diffusion EXPORTS moisture from the moist
    # ITCZ (the equator for a symmetric climate) — the backwards sign the mean cell repairs.
    _, c = _sym_climate()
    b = sm.sphere_moisture_budget(c, hadley=False)
    i_eq = int(np.argmin(np.abs(c.x)))
    assert b.p_minus_e[i_eq] < 0.0


@pytest.mark.slow
def test_demo_reproduces_the_banked_headline():
    # Guards the committed figure (planet-sphere-moist.png): a fresh clone reproduces the headline, not just
    # reads it — emergent rain co-located with the EFE (a check: Hadley beats the eddy export), the budget
    # conserves, the dipole points to the warm hemisphere, and the peak intensity is geometric (not q).
    from planet import demo_sphere_moist as demo
    r = demo.compute()
    assert r.phi_efe > 0.5                                       # ITCZ displaced into the warmed hemisphere
    assert abs(r.rain_max_lat - r.phi_efe) < abs(r.rain_max_lat) # net rain belt TRACKS the EFE (not the equator)
    assert 0.8 < r.rain_max_lat / r.phi_efe < 1.7                # tracks it without wild overshoot
    assert r.eddy_at_efe < 0.0 < r.hadley_at_efe                 # eddy exports, Hadley converges
    assert r.hadley_at_efe > 2.0 * abs(r.eddy_at_efe)            # the falsifiable margin
    assert abs(r.net_p_minus_e) < 1e-9                           # conserved
    assert r.nh_grid[-1] > 0.0 > r.sh_grid[-1]                   # warm hemisphere wets, cold dries
    assert r.peak_real == pytest.approx(r.peak_symq, rel=0.01)   # intensity is GEOMETRIC, not emergent q
