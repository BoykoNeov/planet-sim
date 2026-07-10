"""Triad for the seasonal **ice-albedo** feedback on the marcher (:mod:`planet.seasonal`, rung 5B.1+).

The fixed-albedo scope edge of rung 5B.1, lifted: :func:`planet.seasonal.ice_coalbedo` (the rung-0
step-function albedo, reused) passed as ``march``'s ``coalbedo_fn`` makes the absorbed shortwave
**state-dependent per tile** — a *nonlinear* path the spectral solve cannot carry. What it buys:

*Tight (structural, bit-identical).* (a) **Warm-limit reduction** — an ice feedback with an *unreachable*
freezing threshold never trips, so the ice march is **bit-identical** to the fixed ice-free march (hence,
via 5B.1, to the spectral solve). (b) **Frozen-limit reduction** — a threshold above every temperature
ices everywhere → a *constant* ice co-albedo → **bit-identical** to the fixed-uniform-``a_ice`` march.
*Convergence (the ε=0 self-consistency check).* With the seasons off (``ε=0``) the ice march relaxes to a
**steady** field that is a fixed point of the continuous ice EBM ``L_T·T̄ + S(1−α(T)) − A − B·T = 0`` — the
residual falls at **first order in dt** (the Strang-splitting rate, exactly as the fixed-albedo
marcher→spectral check), reading the model against *itself* (no external insolation/branch to match).
*Conservation.* Global-and-annual net TOA with the **realized** co-albedo ≈ 0.
*Loose (the payoff).* The small-``C`` **land** tile freezes to a far lower latitude and for a larger share
of the year than the sluggish **ocean** tile — continentality as an ice asymmetry — and the ice feedback
keeps Phase 1's **bistability** (a warm vs. a cold seed → a finite-ice vs. a snowball climate at one sun).
"""
import numpy as np
import pytest

from planet import seasonal as sea
from planet.ebm import ALBEDO_ICE, T_FREEZE


# --------------------------------------------------------------------------- #
# TIGHT (a) — warm-limit reduction: a never-tripping ice feedback ≡ the fixed ice-free march.
# --------------------------------------------------------------------------- #
def test_warm_limit_bit_identical_to_fixed_albedo_march():
    # An ice-albedo march whose freezing threshold is unreachably low never flips a single cell to ice, so
    # its co-albedo is the ice-free α everywhere, every step — the march must be BIT-identical to the plain
    # fixed-albedo march (same seed → same trajectory). This proves the coalbedo_fn plumbing reduces to the
    # spectral-consistent path exactly (the default coalbedo_fn=None path is literally unchanged code).
    m = sea.SeasonalEBM(n_cells=90, n_steps=180)
    never_ice = lambda x, T: sea.ice_coalbedo(x, T, T_freeze=-1.0e3)   # threshold never reached
    fixed = m.march(T_init=15.0, tol=1e-9, max_years=120)
    icy = m.march(coalbedo_fn=never_ice, T_init=15.0, tol=1e-9, max_years=120)
    assert np.max(np.abs(icy.T_land - fixed.T_land)) < 1e-12
    assert np.max(np.abs(icy.T_ocean - fixed.T_ocean)) < 1e-12


# --------------------------------------------------------------------------- #
# TIGHT (b) — frozen-limit reduction: an all-ice feedback ≡ the fixed constant-a_ice march.
# --------------------------------------------------------------------------- #
def test_frozen_limit_bit_identical_to_fixed_ice_albedo_march():
    # A freezing threshold above every temperature ices EVERY cell every step → a constant co-albedo
    # (1 − a_ice) — so the nonlinear march collapses onto the fixed-albedo march at α = a_ice. Bit-identical
    # (same seed). Confirms the state-dependent branch evaluates the step function correctly at the ice end.
    m = sea.SeasonalEBM(n_cells=90, n_steps=180)
    always_ice = lambda x, T: sea.ice_coalbedo(x, T, T_freeze=1.0e3)   # threshold always exceeded
    fixed_ice = m.march(albedo=ALBEDO_ICE, T_init=-30.0, tol=1e-9, max_years=120)
    icy = m.march(coalbedo_fn=always_ice, T_init=-30.0, tol=1e-9, max_years=120)
    assert np.max(np.abs(icy.T_land - fixed_ice.T_land)) < 1e-12
    assert np.max(np.abs(icy.T_ocean - fixed_ice.T_ocean)) < 1e-12


# --------------------------------------------------------------------------- #
# CONVERGENCE — the ε=0 in-model self-consistency (a fixed point of the continuous ice EBM).
# --------------------------------------------------------------------------- #
def _epsilon0_residual(n_steps: int) -> float:
    """Max |L_T·T̄ + S(1−α(T)) − A − B·T| over both tiles for the converged ε=0 ice march (W m⁻²)."""
    m = sea.SeasonalEBM(n_cells=120, n_steps=n_steps, obliquity_deg=0.0)
    c = m.march(coalbedo_fn=sea.ice_coalbedo, T_init=10.0, tol=1e-10, max_years=400)
    assert c.converged
    # ε=0 ⟹ the limit cycle is a steady state: every time column is identical, take t=0.
    assert np.max(np.abs(c.T_land - c.T_land[:, :1])) < 1e-9      # genuinely steady
    TL, TO = c.T_land[:, 0], c.T_ocean[:, 0]
    Tbar = m.f_land * TL + m.f_ocean * TO
    S = m.insolation_series()[:, 0]                                # time-constant at ε=0
    LT = m._apply_LT(Tbar)                                         # L_T·T̄  (W m⁻²)
    rL = S * sea.ice_coalbedo(m.x, TL) - m.A - m.B * TL + LT
    rO = S * sea.ice_coalbedo(m.x, TO) - m.A - m.B * TO + LT
    return max(np.max(np.abs(rL)), np.max(np.abs(rO)))


def test_epsilon0_steady_state_is_first_order_fixed_point():
    # With the seasons off, the ice march relaxes to a genuine steady state; that state is a fixed point of
    # the continuous ice EBM only up to the Strang-splitting error, which is O(dt) — so halving dt halves
    # the balance residual. This reads the nonlinear scheme against ITSELF (no external annual-mean
    # reference, sidestepping both the insolation-truncation ghost and the ice-albedo branch-matching
    # problem), and pins the transport-operator assembly under the state-dependent forcing.
    r1, r2, r3 = _epsilon0_residual(90), _epsilon0_residual(180), _epsilon0_residual(360)
    assert 1.7 < r1 / r2 < 2.3
    assert 1.7 < r2 / r3 < 2.3                                     # clean first order in dt → converges


# --------------------------------------------------------------------------- #
# CONSERVATION — global + annual energy balance with the realized (state-dependent) albedo.
# --------------------------------------------------------------------------- #
def test_global_annual_net_toa_with_ice_albedo():
    # Over one converged year the global-and-annual-mean net TOA — using the co-albedo ACTUALLY realized on
    # each tile each step, not the ice-free one — is ≈ 0: the transport is untouched by the ice feedback
    # (it still conserves ∫T̄dx via C_a), so absorbed solar balances OLR in the global-annual mean. The
    # small residual is the O(dt) diagnostic-sampling error (dt-limited, not tol-limited).
    m = sea.SeasonalEBM(n_cells=120, n_steps=360)
    c = m.march(coalbedo_fn=sea.ice_coalbedo, T_init=20.0, tol=1e-8, max_years=300)
    assert c.converged
    S, xc = m.insolation_series(), m.x[:, None]
    aL = S * sea.ice_coalbedo(xc, c.T_land)
    aO = S * sea.ice_coalbedo(xc, c.T_ocean)
    net = m.f_land * (aL - m.A - m.B * c.T_land) + m.f_ocean * (aO - m.A - m.B * c.T_ocean)
    assert abs(float(net.mean())) < 5e-3


# --------------------------------------------------------------------------- #
# LOOSE (the payoff) — continentality as an ice asymmetry.
# --------------------------------------------------------------------------- #
def test_land_freezes_lower_and_longer_than_ocean():
    # The headline: the small-C land tile plunges below freezing over a wide winter band while the sluggish
    # ocean tile barely does — so the land grows a wide seasonal ice zone reaching far lower latitude, and
    # spends much more of the year iced, at a midlatitude the ocean keeps open all year. Direction banked;
    # the exact latitudes/fractions ride the calibrated heat capacities (loose bands).
    m = sea.SeasonalEBM(n_cells=180, n_steps=360)
    c = m.march(coalbedo_fn=sea.ice_coalbedo, T_init=20.0, tol=1e-6, max_years=300)
    assert c.converged
    edge_land = sea.ice_edge_latitude(m.x, c.T_land)
    edge_ocean = sea.ice_edge_latitude(m.x, c.T_ocean)
    assert edge_land < edge_ocean - 10.0                          # land ice reaches well equatorward of ocean
    i45 = m.nearest_index(45)
    frac_land = c.ice_fraction("land")[i45]
    frac_ocean = c.ice_fraction("ocean")[i45]
    assert frac_land > 0.25                                       # land at 45° freezes a good part of the year
    assert frac_ocean < 0.05                                      # the ocean at 45° stays essentially open
    assert 20.0 < edge_land < 55.0                                # strong-continental seasonal-ice reach (loose)


def test_land_ice_is_seasonal_ocean_ice_is_perennial():
    # A subtler nugget the heat-capacity contrast forces: land ice is PURELY seasonal — the tiny-C land tile
    # climbs above freezing every summer, so no land cell holds ice year-round (perennial edge at the pole);
    # ocean ice, once formed, is ~perennial — the huge-C ocean barely warms in summer, so its winter ice
    # survives (perennial edge ≈ its seasonal edge). Continentality decides not just where ice forms but
    # whether it melts.
    m = sea.SeasonalEBM(n_cells=180, n_steps=360)
    c = m.march(coalbedo_fn=sea.ice_coalbedo, T_init=20.0, tol=1e-6, max_years=300)
    land_perennial = sea.ice_edge_latitude(m.x, c.T_land, kind="perennial")
    ocean_seasonal = sea.ice_edge_latitude(m.x, c.T_ocean, kind="seasonal")
    ocean_perennial = sea.ice_edge_latitude(m.x, c.T_ocean, kind="perennial")
    assert land_perennial > 85.0                                  # ~no year-round land ice (all melts)
    assert abs(ocean_perennial - ocean_seasonal) < 5.0           # ocean ice, once formed, persists


# --------------------------------------------------------------------------- #
# LOOSE — the feedback keeps Phase 1's bistability (warm vs. cold seed → two climates).
# --------------------------------------------------------------------------- #
def test_warm_vs_cold_start_bistability():
    # The ice-albedo nonlinearity carries the Snowball bistability INTO the seasonal cycle: at one and the
    # same sun a warm seed settles on a finite-ice climate, a cold seed on a frozen snowball (ice to the
    # equator). The path dependence IS the bistability the rung-0 hysteresis traces — here inside the
    # marching seasonal EBM, which the fixed-albedo (linear, unique-limit-cycle) march cannot show.
    m = sea.SeasonalEBM(n_cells=120, n_steps=180)
    warm = m.march(coalbedo_fn=sea.ice_coalbedo, T_init=20.0, tol=1e-6, max_years=300)
    cold = m.march(coalbedo_fn=sea.ice_coalbedo, T_init=-40.0, tol=1e-6, max_years=300)
    assert warm.converged and cold.converged
    Tw = float(warm.annual_mean("mean").mean())
    Tc = float(cold.annual_mean("mean").mean())
    assert Tw - Tc > 30.0                                         # two genuinely distinct climates
    assert sea.ice_edge_latitude(m.x, cold.T_mean) < 5.0         # the cold branch is a snowball (ice to eq)
    assert sea.ice_edge_latitude(m.x, warm.T_mean) > 30.0        # the warm branch has open midlatitudes


# --------------------------------------------------------------------------- #
# API guard — the ice feedback is exclusive with a fixed albedo/absorbed field.
# --------------------------------------------------------------------------- #
def test_coalbedo_fn_exclusive_with_fixed_forcing():
    m = sea.SeasonalEBM(n_cells=30, n_steps=90)
    with pytest.raises(ValueError, match="exclusive"):
        m.march(coalbedo_fn=sea.ice_coalbedo, albedo=0.3)
    with pytest.raises(ValueError, match="exclusive"):
        m.march(coalbedo_fn=sea.ice_coalbedo, absorbed=np.zeros((30, 90)))
