"""One-way EBM → shallow-water coupler: the climate gradient forces an emergent jet (Planet Phase 4).

The capstone's final phase, and the one that **couples the two shared engines**: the frozen
diffusion spine (via the EBM, :mod:`projects.planet.ebm`) hands its equilibrium
**meridional temperature gradient** to the frozen rotating shallow-water engine (via
:mod:`projects.planet.circulation`), and a **geostrophically-balanced midlatitude jet emerges**.
This is the planet analogue of the whole-program payoff — *forcing → climate → circulation* —
the third link in the chain after Steel's *cooling → microstructure* and Chip's *process → device*.

**One-way (v1).** Climate forces circulation; the flow does **not** feed back on the climate.
Two-way coupling (an advected temperature tracer closing the heat budget back onto the EBM) is
**rung 1** of the §5 GCM climb — seamed at the engine's ``tracer`` slot, not built (plan §3–4).

How the two engines are coupled — forcing split *around* the frozen fluid engine
---------------------------------------------------------------------------------
The frozen :mod:`engines.fluid` solves the **bare** shallow-water equations: ``step`` carries no
forcing term. So the thermal forcing is composed *around* it by **operator splitting** — the
**identical idiom** the EBM uses to graft radiation around the diffusion engine
(:meth:`~projects.planet.ebm.EnergyBalanceModel._radiation_half`) and Steel's Jominy used for its
lateral sink. Each coupler step is *half-forcing / full bare-engine-step / half-forcing*, where the
forcing half-step is two **exact-exponential** relaxations (so it is unconditionally stable and
contributes no splitting instability):

* **Newtonian thermal relaxation** of the height field toward the EBM-derived target,
  ``h ← h_target + (h − h_target)·exp(−½Δt/τ_relax)`` — the "EBM hands the flow a height field"
  of the §8 cadence mandate, here a *steady* (one-way, quasi-equilibrium EBM) target;
* **weak Rayleigh drag** on the momentum, ``(u, v) ← (u, v)·exp(−½Δt/τ_drag)``, which lets the
  spun-up flow settle to a steady jet (``τ_drag ≫ 1/f`` so the steady jet is **near-geostrophic**).

The rotation does the rest: a thermally-raised height field at rest is *not* balanced, so it
accelerates a flow that Coriolis turns into a **zonal jet** in geostrophic balance with the forced
height gradient. The jet is **emergent** — nothing imposes ``u``; only ``h_target`` is imposed, and
``u`` is the dynamics' response. (A non-rotating planet would produce *no* jet — just a stagnant
high-pressure field; the jet's very existence is the rotation signature.)

The geometry bridge — embedding the 1-D zonal-mean EBM in the 2-D β-plane channel
---------------------------------------------------------------------------------
The EBM is 1-D zonal-mean on ``x = sin φ ∈ [0, 1]`` (one hemisphere); the engine is a 2-D
doubly-periodic β-plane channel (:mod:`projects.planet.circulation`). The bridge
(:func:`height_target`): sample the EBM temperature at the channel's latitudes, form the anomaly
``θ(y) = T(φ(y)) − ⟨T⟩``, scale it to a height anomaly ``η_target = α·θ`` (warm → high, the
thermal/hydrostatic sign), and **make it admissible for the doubly-periodic engine** — this is the
one non-obvious step:

* **Periodicity (the engine has no walls).** The engine is doubly periodic; rigid channel walls in
  ``y`` are its *named, unbuilt* BC extension. A monotonic ``θ`` (warm equator → cold pole) is **not**
  periodic — its y-seam would jump by the full equator-to-pole contrast and force a spurious giant
  boundary jet. So ``θ`` is **windowed** (a Tukey taper to matched, zero-slope edges) and made
  **discretely zero-mean**, giving a C¹-periodic, mass-neutral target. The price the periodic channel
  exacts is real: near-zero net zonal momentum forces a **compensating easterly return** flanking the
  westerly jet. The *existence* of that return is geometric (the periodic constraint) and physical, and
  its east–west–east **sign** banding qualitatively resembles the general circulation; but this
  single-layer periodic channel does **not** reproduce the observed *westerly-dominant magnitudes* —
  here the poleward easterly is actually the strongest band, and its concentration at the poleward edge
  is **window-construction-dependent**, not observed. So the emergent westerly jet (the validated,
  benchmarked feature) is the claim; the return is a **named scope edge**, not a faithful
  trade-wind/polar-easterly reconstruction (plan §3).
* **The baroclinic zone, not the ice cliff.** The channel is centred at :data:`PHI_REF_DEG` and
  spans the smooth midlatitude gradient, deliberately **excluding the ice-line albedo cliff**
  (~73°, :mod:`projects.planet.albedo`), whose sharp gradient would otherwise dominate the forcing.
* **Emergence, not placement.** Because the window is broad and flat across the interior, the jet
  lands where the **EBM gradient is steepest** (≈ the gradient maximum, :data:`PHI_REF_DEG`-poleward
  for present-day Earth) — *not* hand-placed at the channel centre. ``test_coupler`` proves this with
  a synthetic off-centre gradient whose jet follows it.

Validation triad (plan §3) — what is asserted tight vs loose
------------------------------------------------------------
* **Analytical limit (tight, amplitude-independent).** The emergent jet is in **geostrophic
  balance**: the steady zonal-mean ``f·u ≈ −g·∂h/∂y`` to a few percent in the jet core
  (:func:`geostrophic_balance`). This is *distinct* from Phase 3's seal (geostrophic *adjustment in
  isolation*): here the **coupled** system produces a balanced jet, and the balance is checked on the
  model's *actual* steady field. It does not depend on the forcing amplitude. *Caveat:* because the flow
  is forced toward a (near-)balanced target, this balance leg is **partly self-fulfilling** — so the
  **emergence** and **release** legs below, not the small balance residual, carry the decisive validation.
* **Conservation (reframed — see below).** **Mass to machine precision** throughout the forced run
  (the zero-mean target + the engine's exact mass invariant), and a **release test**: turn the
  forcing & drag *off* and run the bare frozen engine — mass / energy / potential enstrophy are
  conserved (the engine's Phase-3 guarantees re-confirmed in the coupled configuration) **and the jet
  persists**, proving it is a genuine balanced state, not a forcing-propped artifact.
* **Benchmark (loose).** The westerly **jet latitude** (~30–45°) and **strength** (tens of m/s) vs the
  observed midlatitude jet. Latitude is amplitude-independent (validated-ish, emergent); *strength*
  scales with the calibrated forcing amplitude ``α`` (tuning, flagged).

The conservation reframe (build-time honesty call, advisor-blessed)
-------------------------------------------------------------------
The plan's Phase-4 prose lists "(mass, PV, energy) preserved under the steady forcing", but that is
internally inconsistent with its own "*forced* steady flow" — a **forced–dissipative** system does
**not** conserve energy or PV: the forcing–drag balance is *precisely what selects* the steady jet
(inject available potential energy via relaxation; remove kinetic energy via drag). Claiming
otherwise would be false (the same honesty class as Phase 3's "energy *or* enstrophy, as measured").
So the conservation leg is reframed: **mass machine-exact under forcing** + the **release test**
re-confirms the engine's mass/energy/enstrophy invariants once the forcing is removed. This is the
honest reading of "the frozen engine's guarantees re-confirmed in the coupled configuration".

Non-circularity, named scope edge (plan §3)
-------------------------------------------
*Validated tight:* geostrophic balance of the emergent jet (amplitude-independent) and the jet
**latitude** tracking the EBM gradient maximum (emergent — the synthetic off-centre test).
*Calibrated/flagged (loose):* the height-per-Kelvin amplitude ``α`` (→ jet *speed*, tuning), and the
exact wind magnitudes. *Scope edge:* **one-way** (no feedback); **dry single layer** — there is no
advected thermodynamic variable, so this coupler **cannot** claim *poleward heat transport*, the
*reduction-to-diffusive-EBM* limit (both rung 1, needing the tracer), or *thermal-wind* balance
(rung 3, needing vertical shear) — those anchors live on the §5 climb, not here. The geometry bridge
(1-D zonal-mean EBM ↔ 2-D β-plane) is a **reduced coupling** (the channel's meridional structure is
forced from the zonal-mean gradient; **doubly-periodic**, so a flanking easterly return is required).

Units — SI ([[shallow-water-source]], [[ebm-radiation-source]])
---------------------------------------------------------------
Inherits the fluid engine's SI (``h``, ``H`` m; ``u``, ``v`` m/s; time s; ``f₀`` 1/s; ``β`` 1/(m·s))
and the EBM's ``T`` in °C. The forcing amplitude ``α`` is m/K. Latitudes in degrees.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from engines.fluid import ShallowWater, SWState, uniform_grid

from . import circulation as circ
from .albedo import EBMParams, present_day_climate
from .ebm import ClimateState

# --------------------------------------------------------------------------- #
# Coupler constants — the channel embedding + the forcing (calibrated, loose; plan §3).
# --------------------------------------------------------------------------- #
PHI_REF_DEG = 40.0          # °   — channel-centre latitude (midlatitudes; equatorward of the ice cliff)
CHANNEL_N_LR = 4.5          # —   — channel half-extent in deformation radii (brackets the baroclinic zone)
HEIGHT_PER_KELVIN = 28.0    # m/K — forcing amplitude η = α·θ. CALIBRATED so the jet is tens of m/s;
                            #       the jet *speed* scales with it (tuning), the latitude does not.
WINDOW_TAPER = 0.4          # —   — Tukey taper fraction → a C¹-periodic, zero-mean height target
THERMAL_RELAX_PERIODS = 3.0 # —   — τ_relax in inertial periods (the EBM-target relaxation time)
DRAG_PERIODS = 15.0         # —   — τ_drag in inertial periods; ≫ 1 so the steady jet is near-geostrophic
SPINUP_PERIODS = 35.0       # —   — inertial periods of forced spin-up cap (early-stops once the jet settles)
RELEASE_PERIODS = 10.0      # —   — inertial periods of the forcing-off release test
CONVERGE_TOL = 0.1          # m/s — steady when the zonal wind changes < this over a whole inertial period
                            #       (the jet creeps up monotonically on the drag timescale → ~99.9% settled)


def _tukey_window(n: int, taper: float) -> np.ndarray:
    """A Tukey (tapered-cosine) window of length ``n`` — flat interior, cosine ramps to 0 at both ends.

    ``taper`` is the fraction of the length spent in the two cosine ramps (0 → boxcar, 1 → Hann). The
    ramps reach value **and slope** zero at the edges, so multiplying a field by it yields a C¹-periodic
    field on the doubly-periodic grid (no seam discontinuity — the periodicity fix this module needs).
    """
    w = np.ones(n)
    if taper <= 0.0:
        return w
    edge = int(np.floor(taper * (n - 1) / 2.0))
    if edge < 1:
        return w
    k = np.arange(edge + 1)
    ramp = 0.5 * (1.0 + np.cos(np.pi * (2.0 * k / (taper * (n - 1)) - 1.0)))
    w[:edge + 1] = ramp
    w[-(edge + 1):] = ramp[::-1]
    return w


def channel_latitudes(grid, sw, phi_ref_deg: float = PHI_REF_DEG) -> np.ndarray:
    """Latitude (degrees) of each cell-centre row of the β-plane channel.

    The β-plane maps ``y`` linearly to latitude about the reference latitude: ``φ = φ_ref + (y − y_ref)/a``
    (radians), with ``a`` Earth's radius. ``y`` increases poleward (toward larger ``f``), so the row
    latitudes run equatorward-edge → poleward-edge.
    """
    y = grid.y_centers()
    return phi_ref_deg + np.degrees((y - sw.y_ref) / circ.R_EARTH)


def height_target(state: ClimateState, grid, sw, phi_ref_deg: float = PHI_REF_DEG,
                  alpha: float = HEIGHT_PER_KELVIN, taper: float = WINDOW_TAPER):
    """Build the doubly-periodic, zero-mean target height anomaly ``η_target(y)`` from the EBM ``T(φ)``.

    Samples the EBM equilibrium temperature at the channel latitudes, forms the anomaly
    ``θ(y) = T(φ(y)) − ⟨T⟩``, scales it to a height anomaly ``α·θ`` (warm → high), **windows** it to a
    C¹-periodic shape (the engine has no y-walls; see the module docstring), and removes the discrete
    mean so the relaxation conserves mass exactly. Returns ``(eta_field, eta_profile, phi)``:
    ``eta_field`` the full ``(ny, nx)`` target anomaly (broadcast across longitude — zonally symmetric),
    ``eta_profile`` its 1-D ``y`` profile (m), ``phi`` the channel latitudes (degrees).
    """
    phi = channel_latitudes(grid, sw, phi_ref_deg)
    T_chan = np.interp(phi, state.latitude_deg(), state.T)      # EBM T sampled at channel latitudes
    theta = T_chan - T_chan.mean()                              # warm equatorward (+), cold poleward (−)
    eta_profile = alpha * theta * _tukey_window(phi.size, taper)
    eta_profile = eta_profile - eta_profile.mean()              # discrete zero-mean → mass-conserving forcing
    eta_field = np.repeat(eta_profile[:, None], grid.nx, axis=1)
    return eta_field, eta_profile, phi


def gradient_peak_latitude(state: ClimateState, phi_lo: float, phi_hi: float) -> float:
    """Latitude (degrees) of the steepest EBM meridional temperature gradient within ``[phi_lo, phi_hi]``.

    Computed **independently of the flow** (straight off the EBM field) so the jet-latitude check is a
    genuine "the jet sits at the climate's gradient maximum" test, not a tautology. Restricted to the
    channel band (the smooth baroclinic zone, excluding the ice-line cliff).
    """
    phi = state.latitude_deg()
    dTdphi = np.gradient(state.T, phi)
    band = (phi >= phi_lo) & (phi <= phi_hi)
    idx = np.where(band)[0]
    return float(phi[idx[int(np.argmax(np.abs(dTdphi[idx])))]])


# --------------------------------------------------------------------------- #
# The coupled result — the banked Phase-4 artifact (plain arrays, the loose-coupling currency).
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class CoupledJet:
    """The emergent-jet record produced by the one-way coupler (plain arrays — no live solver crosses).

    ``phi`` channel latitudes (deg). Forced steady state: ``u_profile`` the zonal-mean zonal wind
    (m/s, +eastward/westerly), ``h_profile`` the zonal-mean thickness (m), ``eta_target`` the 1-D target
    anomaly (m), ``u_geostrophic`` the geostrophic estimate ``−(g/f)∂h/∂y`` (m/s). 2-D fields ``u2d``,
    ``v2d`` (m/s) and coordinates ``x``, ``y`` (m) for the map's circulation overlay. Scalars:
    ``jet_lat``/``jet_speed`` (the westerly maximum), ``gradient_peak_lat`` (the EBM gradient max,
    computed independently), ``core_balance_residual`` (the relative geostrophic imbalance at the jet
    core), ``L_R`` the deformation radius. Conservation diagnostics over the **forced** run
    (``times``/``mass``/``energy``/``enstrophy`` — relative drifts) and the **release** run
    (``*_release``, with ``u_profile_release`` to show the jet persists). ``converged``/``iterations``
    record the spin-up.
    """

    phi: np.ndarray
    u_profile: np.ndarray
    h_profile: np.ndarray
    eta_target: np.ndarray
    u_geostrophic: np.ndarray
    x: np.ndarray
    y: np.ndarray
    u2d: np.ndarray
    v2d: np.ndarray
    jet_lat: float
    jet_speed: float
    gradient_peak_lat: float
    core_balance_residual: float
    L_R: float
    times: np.ndarray
    mass: np.ndarray
    energy: np.ndarray
    enstrophy: np.ndarray
    times_release: np.ndarray
    mass_release: np.ndarray
    energy_release: np.ndarray
    enstrophy_release: np.ndarray
    u_profile_release: np.ndarray
    converged: bool
    iterations: int


def geostrophic_balance(sw, grid, state: SWState):
    """Zonal-mean geostrophic-balance diagnostic of a (near-zonal) state: ``(u_bar, h_bar, u_geo, f)``.

    Returns the zonal-mean zonal wind ``u_bar(y)``, thickness ``h_bar(y)``, the geostrophic estimate
    ``u_geo = −(g/f)·∂h_bar/∂y`` (central difference, periodic), and the Coriolis ``f(y)`` — all at the
    cell-centre rows (where both ``u`` and ``h`` live on the C-grid in ``y``), so the comparison is
    collocated. ``u_bar ≈ u_geo`` in the jet core is the Phase-4 analytic anchor.
    """
    u_bar = state.u.mean(axis=1)
    h_bar = state.h.mean(axis=1)
    dh_dy = (np.roll(h_bar, -1) - np.roll(h_bar, 1)) / (2.0 * grid.dy)
    f = sw.f0 + sw.beta * (grid.y_centers() - sw.y_ref)
    u_geo = -sw.g * dh_dy / f
    return u_bar, h_bar, u_geo, f


def _forcing_half(state: SWState, h_target: np.ndarray, decay_h: float, decay_uv: float) -> SWState:
    """One exact-exponential forcing half-step: thermal relaxation of ``h`` + Rayleigh drag on ``(u, v)``.

    ``h → h_target + (h − h_target)·decay_h`` (Newtonian relaxation; mass-neutral because ``h_target`` and
    ``h`` share the same area mean ``H``) and ``(u, v) → (u, v)·decay_uv`` (linear drag). Both are the
    *exact* solutions of their linear ODEs over the half-step — the EBM ``_radiation_half`` idiom.
    """
    return SWState(h=h_target + (state.h - h_target) * decay_h,
                   u=state.u * decay_uv, v=state.v * decay_uv)


def couple_jet(state: ClimateState | None = None, nx: int = 96, ny: int = 96,
               phi_ref_deg: float = PHI_REF_DEG, n_LR: float = CHANNEL_N_LR,
               alpha: float = HEIGHT_PER_KELVIN, taper: float = WINDOW_TAPER,
               spinup_periods: float = SPINUP_PERIODS, release_periods: float = RELEASE_PERIODS,
               tol: float = CONVERGE_TOL) -> CoupledJet:
    """Force the shallow-water engine with the EBM gradient → a steady geostrophic jet, then release it.

    Builds the channel and the EBM-derived height target, spins the flow up from rest by Strang-split
    thermal-relaxation + drag (stopping early when the zonal wind stops changing), then runs a
    forcing-off **release** phase on the bare engine. Returns the :class:`CoupledJet` artifact.

    Parameters
    ----------
    state : ClimateState | None
        The EBM equilibrium climate to force from; ``None`` → :func:`present_day_climate`.
    nx, ny : int
        Channel resolution.
    phi_ref_deg, n_LR, alpha, taper
        The channel embedding + forcing knobs (module-constant defaults).
    spinup_periods, release_periods : float
        Inertial periods of forced spin-up (capped; early-stops on convergence) and forcing-off release.
    tol : float
        Convergence tolerance on ``max|Δu_bar|`` (m/s) between iterations.
    """
    if state is None:
        state = present_day_climate(n_tau=0.05)

    f0 = circ.coriolis_f0(phi_ref_deg)
    beta = circ.coriolis_beta(phi_ref_deg)
    L_R = np.sqrt(circ.G_EARTH * circ.H_EQUIV) / f0
    L = n_LR * L_R
    grid = uniform_grid(L, L, nx, ny)
    sw = ShallowWater(grid, circ.G_EARTH, circ.H_EQUIV, f0=f0, beta=beta)

    eta_field, eta_profile, phi = height_target(state, grid, sw, phi_ref_deg, alpha, taper)
    h_target = sw.H + eta_field

    s = SWState(h=sw.H * np.ones((ny, nx)), u=np.zeros((ny, nx)), v=np.zeros((ny, nx)))
    inertial = 2.0 * np.pi / f0
    dt = sw.max_dt(s) * 0.5
    tau_relax = THERMAL_RELAX_PERIODS * inertial
    tau_drag = DRAG_PERIODS * inertial
    decay_h = np.exp(-0.5 * dt / tau_relax)
    decay_uv = np.exp(-0.5 * dt / tau_drag)
    n_max = int(spinup_periods * inertial / dt)

    m0, e0, z0 = sw.mass(s), sw.energy(s), sw.potential_enstrophy(s)
    times, mass, energy, enstrophy = [], [], [], []
    # Convergence is judged over a fixed *interval* (one inertial period), not step-to-step: the slow
    # relaxation makes per-step Δu tiny even far from steady, so a step-to-step tolerance would trip
    # immediately. Steady ⟺ the zonal wind stops changing over a whole inertial period.
    check_every = max(1, int(round(inertial / dt)))
    u_ref = s.u.mean(axis=1)
    converged = False
    t = 0.0
    it = 0
    for it in range(1, n_max + 1):
        s = _forcing_half(s, h_target, decay_h, decay_uv)
        s = sw.step(s, dt)
        s = _forcing_half(s, h_target, decay_h, decay_uv)
        t += dt
        times.append(t)
        mass.append(sw.mass(s) / m0 - 1.0)
        energy.append(sw.energy(s) / e0 - 1.0)
        enstrophy.append(sw.potential_enstrophy(s) / z0 - 1.0)
        if it % check_every == 0:
            u_now = s.u.mean(axis=1)
            if np.max(np.abs(u_now - u_ref)) < tol:
                converged = True
                break
            u_ref = u_now

    # -- release: forcing & drag OFF; the bare frozen engine must conserve and the jet must persist -- #
    mr0, er0, zr0 = sw.mass(s), sw.energy(s), sw.potential_enstrophy(s)
    n_rel = int(release_periods * inertial / dt)
    times_r, mass_r, energy_r, enstrophy_r = [], [], [], []
    tr = 0.0
    for _ in range(n_rel):
        s = sw.step(s, dt)
        tr += dt
        times_r.append(tr)
        mass_r.append(sw.mass(s) / mr0 - 1.0)
        energy_r.append(sw.energy(s) / er0 - 1.0)
        enstrophy_r.append(sw.potential_enstrophy(s) / zr0 - 1.0)
    u_profile_release = s.u.mean(axis=1)

    u_bar, h_bar, u_geo, _ = geostrophic_balance(sw, grid, s)
    ji = int(np.argmax(u_bar))
    return CoupledJet(
        phi=phi, u_profile=u_bar, h_profile=h_bar, eta_target=eta_profile, u_geostrophic=u_geo,
        x=grid.x_centers(), y=grid.y_centers(), u2d=s.u.copy(), v2d=s.v.copy(),
        jet_lat=float(phi[ji]), jet_speed=float(u_bar[ji]),
        gradient_peak_lat=gradient_peak_latitude(state, float(phi[0]), float(phi[-1])),
        core_balance_residual=float(abs(u_bar[ji] - u_geo[ji]) / abs(u_bar[ji])),
        L_R=float(L_R),
        times=np.array(times), mass=np.array(mass), energy=np.array(energy), enstrophy=np.array(enstrophy),
        times_release=np.array(times_r), mass_release=np.array(mass_r),
        energy_release=np.array(energy_r), enstrophy_release=np.array(enstrophy_r),
        u_profile_release=u_profile_release, converged=converged, iterations=it,
    )
