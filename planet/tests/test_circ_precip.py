"""Planet rung-1 step-3 validation: circulation-informed precipitation (:mod:`planet.circ_precip`).

Step 3 wires the precip storm-track band's centre to the **emergent jet latitude** instead of the
prescribed constant. What is asserted, and its honesty class (see :mod:`planet.circ_precip`):

* **(banked — the reduction, *by construction*) the default centre recovers the rung-0 field exactly.**
  ``midlat_center_deg`` defaults to the cited :data:`planet.precip.MIDLAT_CENTER_DEG`, so the
  parameterized pattern is the rung-0 pattern bit-for-bit — plumbing, not an independent test (the
  :func:`planet.transport.two_way_pass` honesty), asserted so the seam can never silently drift.
* **(banked — the mechanism) the band tracks the centre / the emergent jet.** Moving the centre moves
  the midlatitude precip peak with it; the ``slow`` composition test drives the centre from a real
  coupled jet on a *synthetic off-centre gradient* and shows the band follows the flow-selected jet —
  *migration*. (The coupler's own "jet tracks the gradient" proof is not re-tested here.)
* **(structure preserved across the realistic range) the band triad survives the moved centre.** ITCZ
  wettest / subtropical trough a local min / poles dry, for the realistic-to-modest centre range — and
  the named **scope edge**: a large *equatorward* displacement shallows the trough toward merging with
  the ITCZ, so the band-tracking is asserted only where the structure holds.

The fast tests use a tiny duck-typed fake jet (only ``jet_lat`` is read); only the migration test runs
the coupler (``slow``).
"""
import types

import numpy as np
import pytest

from planet import circ_precip as cp
from planet import coupler, precip
from planet.albedo import present_day_climate
from planet.ebm import ClimateState


def _fake_jet(jet_lat: float):
    """A duck-typed stand-in for :class:`planet.coupler.CoupledJet` (only ``jet_lat`` is read)."""
    return types.SimpleNamespace(jet_lat=jet_lat)


def _midlat_peak_lat(center: float, lo: float = 15.0, hi: float = 85.0, n: int = 1401) -> float:
    """Latitude (deg) of the midlatitude storm-track precip maximum (poleward of the ITCZ)."""
    phi = np.linspace(lo, hi, n)
    return float(phi[int(np.argmax(precip.precip_pattern(phi, midlat_center_deg=center)))])


# --------------------------------------------------------------------------- #
# Reduction — by construction (the default centre IS the rung-0 field)
# --------------------------------------------------------------------------- #
def test_default_center_reproduces_rung0_pattern_bit_for_bit():
    phi = np.linspace(-90.0, 90.0, 361)
    # explicit default == implicit default == the rung-0 pattern (plumbing, not an independent test)
    assert np.array_equal(precip.precip_pattern(phi, precip.MIDLAT_CENTER_DEG),
                          precip.precip_pattern(phi))
    assert np.array_equal(precip.precipitation(phi, 20.0, precip.MIDLAT_CENTER_DEG),
                          precip.precipitation(phi, 20.0))


def test_circulation_informed_reduces_to_rung0_when_jet_at_prescribed_center():
    st = present_day_climate(n_tau=0.02)
    got = cp.circulation_informed_precip(st, _fake_jet(precip.MIDLAT_CENTER_DEG))
    assert np.array_equal(got, precip.precip_field(st))            # the reduction (jet at 50° → rung 0)


def test_circulation_informed_uses_the_jet_latitude():
    st = present_day_climate(n_tau=0.02)
    got = cp.circulation_informed_precip(st, _fake_jet(43.0))
    want = precip.precipitation(st.latitude_deg(), st.global_mean_T, midlat_center_deg=43.0)
    assert np.array_equal(got, want)
    assert cp.storm_track_center(_fake_jet(-43.0)) == pytest.approx(43.0)   # symmetric (|jet_lat|)


# --------------------------------------------------------------------------- #
# Mechanism — the band tracks the centre (migration)
# --------------------------------------------------------------------------- #
def test_midlat_peak_tracks_the_centre():
    # moving the storm-track centre equatorward moves the midlatitude precip peak equatorward with it
    assert _midlat_peak_lat(50.0) > _midlat_peak_lat(44.0) > _midlat_peak_lat(38.0)
    for center in (38.0, 44.0, 50.0):
        assert _midlat_peak_lat(center) == pytest.approx(center, abs=2.0)   # peak ≈ the centre


def test_relocation_reports_the_equatorward_trade():
    st = present_day_climate(n_tau=0.02)
    rel = cp.relocate(st, _fake_jet(44.0))
    assert rel.center_rung0 == pytest.approx(precip.MIDLAT_CENTER_DEG)
    assert rel.center_circ == pytest.approx(44.0)
    assert rel.displacement == pytest.approx(44.0 - precip.MIDLAT_CENTER_DEG)   # negative = equatorward
    assert rel.displacement < 0.0
    assert np.array_equal(rel.precip_rung0, precip.precip_field(st))
    assert np.array_equal(rel.precip_circ, cp.circulation_informed_precip(st, _fake_jet(44.0)))
    assert not np.array_equal(rel.precip_circ, rel.precip_rung0)                # the band actually moved


# --------------------------------------------------------------------------- #
# Structure preserved across the realistic centre range (+ the named merge edge)
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("center", [38.0, 44.0, 50.0])
def test_band_triad_holds_across_realistic_centers(center):
    phi = np.linspace(0.0, 90.0, 361)
    p = precip.precip_pattern(phi, midlat_center_deg=center)
    assert p[0] == p.max()                                         # ITCZ (equator) wettest
    assert float(precip.precip_pattern(90.0, center)) < float(precip.precip_pattern(center, center))  # poles dry
    belt = (phi > 8.0) & (phi < center)                            # equator → midlat peak
    phi_min = float(phi[belt][int(np.argmin(p[belt]))])
    assert 12.0 < phi_min < center                                 # a subtropical trough equatorward of the band
    pmin = float(precip.precip_pattern(phi_min, center))
    assert pmin < float(precip.precip_pattern(phi_min - 8.0, center))   # local min: drier equatorward (toward ITCZ)
    assert pmin < float(precip.precip_pattern(phi_min + 8.0, center))   # ... and poleward (toward the storm track)


def test_large_equatorward_displacement_shallows_the_trough_named_edge():
    # Named scope edge: pushing the band well equatorward shallows the subtropical trough toward merging
    # with the ITCZ. Prominence = (midlat-peak − trough) shrinks monotonically as the centre moves in.
    def prominence(center):
        phi = np.linspace(0.0, 90.0, 361)
        p = precip.precip_pattern(phi, midlat_center_deg=center)
        belt = (phi > 8.0) & (phi < center)
        trough = float(p[belt].min())
        return float(precip.precip_pattern(center, center)) - trough
    assert prominence(50.0) > prominence(44.0) > prominence(36.0) > prominence(30.0)


# --------------------------------------------------------------------------- #
# The mechanism end-to-end (slow): the band follows a dynamically-selected jet
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_band_follows_the_emergent_jet_on_a_synthetic_gradient():
    """Migration mechanism: an off-centre *synthetic* EBM gradient → an equatorward jet → the storm-track
    rain band follows it (anchored to the flow, displaced from the prescribed 50°).

    The coupler's "jet tracks the gradient" proof is *not* re-tested; this asserts only the new content —
    the precip band tracks the emergent ``jet_lat``. A modest off-centre gradient (35°) keeps the band in
    the range where the subtropical trough survives.
    """
    x = np.linspace(0.01, 0.999, 180)
    phi = np.degrees(np.arcsin(x))
    T = 20.0 - 30.0 / (1.0 + np.exp(-(phi - 35.0) / 3.0))          # smooth step down centred at 35°
    synth = ClimateState(x=x, T=T, global_mean_T=float(T.mean()), ice_line_lat=90.0,
                         net_toa=0.0, converged=True, iterations=0)

    jet = coupler.couple_jet(state=synth, nx=48, ny=48)
    assert jet.jet_lat < 46.0                                      # the flow-selected jet is equatorward of 50°
    rel = cp.relocate(synth, jet)
    assert rel.center_circ == pytest.approx(jet.jet_lat)           # the band centre IS the emergent jet
    assert rel.displacement < -3.0                                 # migration: moved equatorward of the rung-0 band
    peak_circ = float(synth.latitude_deg()[np.argmax(rel.precip_circ * (synth.latitude_deg() > 15.0))])
    assert peak_circ == pytest.approx(jet.jet_lat, abs=3.0)        # the midlat rain peak tracks the jet
