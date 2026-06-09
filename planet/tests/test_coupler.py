"""Planet Phase-4 validation: the one-way EBM → shallow-water coupler (the emergent jet).

The validation triad (plan §3), with the conservation leg **reframed** for a forced–dissipative
system (advisor-blessed; see :mod:`projects.planet.coupler`):

* **Analytical limit (tight, amplitude-independent).** The emergent jet is in **geostrophic
  balance** — the steady zonal-mean ``f·u ≈ −g·∂h/∂y`` to a few percent in the core; and the jet
  **latitude tracks the EBM gradient maximum** (an off-centre *synthetic* gradient makes the jet
  follow it — the decisive emergence proof: nothing imposes the jet location).
* **Conservation (reframed).** Mass is machine-exact under forcing; the **release test** (forcing &
  drag off → bare frozen engine) re-confirms the engine's mass/energy/enstrophy invariants **and the
  jet persists** — a genuine balanced state, not a forcing-propped artifact. (Energy/PV are *not*
  conserved *under* the forcing — that balance is what selects the steady jet; claiming otherwise
  would be false. The leg asserts what is true: mass forced + invariants on release.)
* **Benchmark (loose).** Westerly jet at midlatitudes (~30–45°), tens of m/s.

Fast structural tests cover the builders (the periodic, zero-mean height target; the gradient-peak
diagnostic) always-green; the integration legs are ``slow``-marked (each spins the engine up many
inertial periods) and run at coarse resolution — the engine's own physics is sealed in
``engines/fluid/tests/`` and ``test_circulation``.
"""
import numpy as np
import pytest

from projects.planet import circulation as circ
from projects.planet import coupler
from projects.planet.albedo import present_day_climate
from projects.planet.ebm import ClimateState


# --------------------------------------------------------------------------- #
# Fast structural tests — the builders only (no integration), always-green.
# --------------------------------------------------------------------------- #
def _channel():
    f0 = circ.coriolis_f0(coupler.PHI_REF_DEG)
    L_R = np.sqrt(circ.G_EARTH * circ.H_EQUIV) / f0
    from engines.fluid import ShallowWater, uniform_grid
    grid = uniform_grid(coupler.CHANNEL_N_LR * L_R, coupler.CHANNEL_N_LR * L_R, 48, 48)
    sw = ShallowWater(grid, circ.G_EARTH, circ.H_EQUIV, f0=f0, beta=circ.coriolis_beta(coupler.PHI_REF_DEG))
    return grid, sw


def test_channel_latitudes_bracket_the_reference():
    grid, sw = _channel()
    phi = coupler.channel_latitudes(grid, sw)
    assert phi[0] < coupler.PHI_REF_DEG < phi[-1]            # the channel brackets the reference latitude
    assert np.all(np.diff(phi) > 0)                          # y increases poleward
    assert phi[-1] < 70.0                                    # the poleward edge stays equatorward of the ice cliff


def test_height_target_zero_mean_and_periodic():
    grid, sw = _channel()
    state = present_day_climate(n_tau=0.1)
    _, eta_profile, _ = coupler.height_target(state, grid, sw)
    # discretely zero-mean → the Newtonian relaxation conserves ∫h exactly (the mass leg's basis)
    assert abs(eta_profile.mean()) < 1e-9
    # the Tukey window forces matched VALUE and (near-)zero SLOPE at the y-seam → admissible on the
    # periodic grid (the edge value is a constant −mean offset, which is itself periodic).
    assert eta_profile[0] == pytest.approx(eta_profile[-1], abs=1e-9)
    edge_step = max(abs(eta_profile[1] - eta_profile[0]), abs(eta_profile[-1] - eta_profile[-2]))
    interior_step = np.max(np.abs(np.diff(eta_profile)))
    assert edge_step < 0.5 * interior_step                   # the seam slope is a fraction of the jet's, not a discontinuity


def test_height_target_warm_is_high():
    """Thermal/hydrostatic sign: the height anomaly is high where the EBM is warm (the equatorward side)."""
    grid, sw = _channel()
    state = present_day_climate(n_tau=0.1)
    _, eta_profile, phi = coupler.height_target(state, grid, sw)
    # in the windowed interior, η rises toward the warm (equatorward, smaller φ) side
    interior = slice(eta_profile.size // 4, 3 * eta_profile.size // 4)
    equatorward_half = eta_profile[interior][: (interior.stop - interior.start) // 2]
    poleward_half = eta_profile[interior][(interior.stop - interior.start) // 2:]
    assert equatorward_half.mean() > poleward_half.mean()


def test_gradient_peak_latitude_finds_a_synthetic_step():
    x = np.linspace(0.01, 0.999, 180)
    phi = np.degrees(np.arcsin(x))
    T = 20.0 - 30.0 / (1.0 + np.exp(-(phi - 48.0) / 3.0))   # a smooth step down centred at 48°
    state = ClimateState(x=x, T=T, global_mean_T=float(T.mean()), ice_line_lat=90.0,
                         net_toa=0.0, converged=True, iterations=0)
    assert coupler.gradient_peak_latitude(state, 20.0, 70.0) == pytest.approx(48.0, abs=4.0)


# --------------------------------------------------------------------------- #
# The integration triad — slow (each spins the engine up many inertial periods).
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_emergent_jet_is_geostrophically_balanced_at_midlatitudes():
    """Analytic anchor + loose benchmark: a balanced westerly jet emerges at midlatitudes."""
    r = coupler.couple_jet(nx=48, ny=48)
    assert r.converged
    assert r.jet_speed > 5.0                                 # a real jet (tens of m/s, loose)
    assert 30.0 <= r.jet_lat <= 50.0                         # midlatitude jet (benchmark, loose)
    assert r.core_balance_residual < 0.05                    # f·u ≈ −g ∂h/∂y in the core (the tight anchor)
    # the doubly-periodic channel requires a flanking easterly return (named scope edge, not a bug)
    assert r.u_profile.min() < 0.0


@pytest.mark.slow
def test_jet_latitude_tracks_the_climate_gradient_not_the_channel():
    """Emergence proof: with the gradient placed off-centre (synthetic), the jet follows it.

    The jet is not pinned to the channel centre — moving the EBM gradient maximum moves the jet the
    same way. This is the amplitude-independent, *validated* leg of the non-circularity split.
    """
    def synth(step_lat):
        x = np.linspace(0.01, 0.999, 180)
        phi = np.degrees(np.arcsin(x))
        T = 20.0 - 30.0 / (1.0 + np.exp(-(phi - step_lat) / 3.0))
        return ClimateState(x=x, T=T, global_mean_T=float(T.mean()), ice_line_lat=90.0,
                            net_toa=0.0, converged=True, iterations=0)

    eq = coupler.couple_jet(state=synth(30.0), nx=48, ny=48)    # gradient equatorward
    pole = coupler.couple_jet(state=synth(50.0), nx=48, ny=48)  # gradient poleward
    assert pole.jet_lat > eq.jet_lat + 5.0                      # the jet moved poleward WITH the gradient
    # each jet lands in the neighbourhood of its (independently-computed) gradient maximum
    assert abs(eq.jet_lat - eq.gradient_peak_lat) < 8.0
    assert abs(pole.jet_lat - pole.gradient_peak_lat) < 8.0


@pytest.mark.slow
def test_mass_machine_exact_and_release_reconfirms_engine_invariants():
    """Conservation leg (reframed): mass forced-exact; release re-confirms invariants + the jet persists."""
    r = coupler.couple_jet(nx=48, ny=48)
    assert np.abs(r.mass).max() < 1e-10                        # mass machine-exact under forcing (zero-mean target)
    # release: forcing & drag OFF → the bare frozen engine conserves its invariants ...
    assert np.abs(r.mass_release).max() < 1e-10
    assert np.abs(r.energy_release).max() < 1e-6
    assert np.abs(r.enstrophy_release).max() < 1e-3
    # ... and the jet persists (a genuine balanced state, not propped up by the forcing)
    assert r.u_profile_release.max() == pytest.approx(r.jet_speed, rel=0.05)


@pytest.mark.slow
def test_jet_speed_scales_with_amplitude_but_latitude_does_not():
    """Non-circularity split: speed scales with the calibrated amplitude α (tuning); latitude does not."""
    base = coupler.couple_jet(nx=48, ny=48, alpha=coupler.HEIGHT_PER_KELVIN)
    strong = coupler.couple_jet(nx=48, ny=48, alpha=2.0 * coupler.HEIGHT_PER_KELVIN)
    assert strong.jet_speed / base.jet_speed == pytest.approx(2.0, abs=0.3)   # speed ∝ α (loose, tuning)
    assert strong.jet_lat == pytest.approx(base.jet_lat, abs=2.0)             # latitude amplitude-independent
