"""Emergent eddy heat flux on the barotropically-unstable jet — Planet rung 1, step 2, **Phase B**.

Phase A (:mod:`planet.transport`) built the two-way feedback **machinery** — the κ→D bridge and the
EBM re-equilibration — and drove it with a **synthetic** exactly-down-gradient flux, so the machinery
landed independent of the (messy) eddy sim. Phase B supplies the **real** thing: the meridional eddy
heat flux ``⟨v'θ'⟩(φ)`` diagnosed from a passive temperature tracer advected on the **released**,
barotropically-unstable Phase-4 jet. It is **genuinely emergent** — the eddies come from the resolved
barotropic instability (the step-0 probe established it exists: a v-perturbation grows ~200× then
saturates), with **no imposed stationary wave** and **no down-gradient closure assumed**. This is the
``flux_fn`` seam Phase A left open, now filled by physics rather than a stand-in.

What Phase B banks — and the honest scope of each leg
-----------------------------------------------------
* **(headline, DIRECTION banked) The eddy diffusivity is STATE-DEPENDENT — the loop is a real
  feedback, not cosmetic.** Across two climates with the height-forcing amplitude ``α`` held **fixed**,
  a *flatter* EBM gradient (a high-obliquity-like world, smaller ``|s₂|``) makes a **weaker jet** →
  **weaker eddies** → a **smaller** ``κ_eff``. So strengthening the mean state strengthens the
  transport: a right-signed negative feedback the two-way loop genuinely closes. The mechanism is the
  one the design predicts: in ``κ_eff = −∫F̄ dt / ∫θ̄_y dt`` the *gradient cancels*, so
  ``κ_eff ≈ v'·ℓ`` (the eddy **stirring rate**) and it tracks climate **only** through the jet — which
  is exactly why ``α`` must be held fixed (renormalizing the forcing to fix the jet speed would sever
  the one channel through which ``κ_eff`` tracks climate, and the test would go cosmetic by
  construction). This non-circularity is the leg that makes the loop *not* a re-labelling of a fixed D.
* **(MAGNITUDE — named, NOT banked) ``κ_eff ~ 10³ m²/s``, ~1000× below rung-0's ``2.2×10⁶``.** The
  *instantaneous* ``⟨v'θ'⟩`` is largely **reversible** — it oscillates in sign as the eddy meander
  sloshes — so the **irreversible fraction** (net ``∫F̄ dt`` over the time-integral of ``|F̄|``) is only
  ``~0.1``: the life-cycle integral cancels most of the sloshing and keeps that small down-gradient
  residual. (That ~10× reversibility loss is *one* contributor to the small κ, not the whole ~1000× —
  do not read it as a closed budget.) The value is **converged** in resolution but **not bankable**: it
  is suppressed by configuration choices that cannot be cleanly separated from the intrinsic physics —
  a single coherent seeded wavenumber mixes more reversibly than broadband turbulence, and the
  band-bulk estimator smears the jet-centred peak over the whole interior. So the magnitude is
  **window / forcing / configuration-tuned**, exactly as Phase 4's jet *speed* was — the *sign* and the
  *climate ordering* are what carry validation, not the number. (Chasing a bigger number is wasted
  effort: it is unbankable regardless.)
* **(the tight reduction — a FINDING, not a manufactured match) The barotropic eddy flux does NOT
  *tightly* reduce to the EBM diffusion operator at rung 1.** :func:`reduction_to_ebm_operator` tests —
  does **not** assume — whether the resolved flux-divergence ``−∂F̄/∂y`` has the *shape* of a smooth
  down-gradient diffusion ``(1/cosφ)∂/∂y[κ·cosφ·∂θ̄/∂y]`` (built from the band-bulk **scalar** κ + the
  smooth mean gradient, so it reuses neither the pointwise κ — *not circular* — nor an absolute scale,
  only the **normalised shape**). The result is a **partial correlation** (``~0.5–0.6`` for the
  strong-jet case, resolution-sensitive for the weak-jet one) — the gross down-gradient tendency is
  there, but the resolved divergence carries jet-localised structure the smooth operator misses, so it
  is **not a tight reduction**. And the comparison is *itself near the vacuous edge*: a uniform-κ
  diffusion of the
  **near-linear** midlatitude gradient produces almost no transport-divergence, so there is little
  signal to match. The honest reading: the tight reduction (flux-divergence = EBM operator) becomes
  **non-vacuous only at rung 3**, when a *strong baroclinic* flux drives a large divergence — and the
  **geometry correspondence it needs is already delivered** here
  (:func:`planet.transport.spherical_transport_tendency`: the spherical operator with its order-unity
  ``cos φ`` metric, anchored on the P₂ eigenvalue). "The reduction arrives at rung 3, with the geometry
  already correct" is stronger than a forced match on a ``~0`` field.

The diagnostic — a pure-release life-cycle integral (the settled fork)
----------------------------------------------------------------------
The jet is spun up **dry** by the Phase-4 forcing, then the tracer ``θ`` is initialised to the
**windowed EBM temperature profile** (the same C¹-periodic, interior-flat shape the coupler uses for
the height target, so the channel's mean meridional gradient lives in the window-flat interior) and a
small **deterministic** ``cos(kx)`` v-perturbation is added. The forcing is then switched **off**
(release mode) and ``(h, u, v, θ)`` advect on the bare engine: the instability grows, saturates, and
stirs the tracer. ``κ_eff`` is the **band-bulk** life-cycle integral
``−∫F̄ dt / ∫θ̄_y dt`` over the window-flat interior, integrated from a finite-amplitude onset through
the end of the run (the reversible sloshing then cancels in the time-integral). Initialising ``θ`` at
release — rather than relaxing it toward the target *during* spin-up — is equivalent for the
zonal-mean profile (a steady jet has ``v̄ ≈ 0``, so it does not advect the meridional mean) and avoids
doubling the tracer cost; relaxing ``θ`` *during release* was tried and **rejected** (it injects a
relaxation-timescale artefact that breaks the climate ordering).

Units — SI ([[shallow-water-source]], [[ebm-radiation-source]]); θ a °C temperature tracer, v in m/s,
the flux ``⟨v'θ'⟩`` in K·m/s, the diffusivity ``κ`` in m²/s, the EBM ``D`` in W m⁻² K⁻¹.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Optional

import numpy as np

from engines.fluid import ShallowWater, SWState, uniform_grid

from . import circulation as circ
from . import coupler
from . import transport as tr
from .albedo import EBMParams, present_day_climate
from .ebm import ClimateState

# --------------------------------------------------------------------------- #
# Phase-B run constants — the released life cycle (loose/tuned; see the module docstring).
# --------------------------------------------------------------------------- #
RELEASE_PERIODS = 60.0     # —   — inertial periods of forcing-off release (long enough that the
                           #       weak-jet climate also saturates within the run)
WINDOW_START = 10.0        # —   — integrate the life-cycle flux from this period to the end (past the
                           #       perturbation's finite-amplitude onset; κ is robust to this choice)
PERTURB_EPS = 0.5          # m/s — amplitude of the deterministic v-perturbation (≪ the jet)
PERTURB_WAVENUMBER = 2     # —   — zonal wavenumber of the cos(kx) seed (deterministic — reproducible)


@dataclass(frozen=True)
class EddyFrames:
    """Banked ``(h,u,v,θ)`` snapshots of the released eddy life cycle — the rung-A viz side-channel.

    **Diagnostic-pure:** produced only when ``eddy_life_cycle(..., n_frames>0)``; ``n_frames=0``
    leaves the κ diagnosis bit-for-bit unchanged (the inert-seam discipline applied to motion). Frames
    are snapshotted at even **time** thresholds over the **full** release ``[0, t_end]`` — the adaptive-
    dt span where growth→saturation, *the mechanism*, lives (not the κ window).

    ``times`` the release time (inertial periods) at each frame; ``h``/``u``/``v``/``theta`` the
    ``(n_frames, ny, nx)`` C-grid snapshots — ``h,u,v`` carried so the verification anchors recompute
    (∫hθ, eddy-KE), ``theta`` the stirred tracer the mechanism panel animates. ``x``/``y`` the
    cell-centre coordinates (m), ``phi`` the channel latitudes (deg), ``interior`` the window-flat band
    the flux is diagnosed over, ``cell_area`` the grid cell area (m²) the eddy-KE recompute needs,
    ``window_start`` the κ-window onset (periods) panel 2 marks.

    ``net_cum``/``thru_cum`` (``n_frames``) are the cumulative meridional-transport traces, integrated
    **per latitude first** over the full release: ``net_cum = Σ_interior|∫F̄ dt|`` and
    ``thru_cum = Σ_interior∫|F̄| dt`` (K·m, interior-band sums). The abs-after-time-integral kills
    *temporal* sloshing per latitude while the abs-before-spatial-sum keeps *spatial* structure from
    cancelling in — so ``thru_cum`` climbs (the swirls rage) while ``net_cum`` stays a small fraction,
    and their endpoint ratio **lands on** the (full-release) irreversible fraction by construction (the
    plan's "the real proof is numerical"). The banked, *windowed* :attr:`EddyFlux.irreversible_fraction`
    is the headline number; panel 2 marks ``window_start`` so the two read consistently.
    """

    times: np.ndarray
    h: np.ndarray
    u: np.ndarray
    v: np.ndarray
    theta: np.ndarray
    x: np.ndarray
    y: np.ndarray
    phi: np.ndarray
    interior: np.ndarray
    cell_area: float
    window_start: float
    net_cum: np.ndarray
    thru_cum: np.ndarray


@dataclass(frozen=True)
class EddyFlux:
    """One released eddy life cycle and its diagnosed effective diffusivity (plain arrays).

    ``phi`` channel latitudes (deg), ``y`` the meridional cell-centre coordinate (m), ``interior`` the
    window-flat mask κ is diagnosed over. ``jet_speed``/``jet_lat`` the spun-up jet (m/s, deg);
    ``rayleigh_kuo`` whether the necessary barotropic-instability condition (``β − u_yy`` sign change)
    is met. ``times``/``eddy_ke`` the eddy kinetic energy ``∫½h[(u−ū)²+(v−v̄)²]dA`` (J) over the
    release and ``saturation_period`` its peak time. ``theta_bar`` the (smooth, windowed) mean tracer
    profile (°C); ``F_int``/``G_int`` the life-cycle integrals ``∫F̄ dt`` (K·m) and ``∫θ̄_y dt`` (K)
    over the window; ``kappa_profile`` the pointwise ``κ(φ) = −F_int/G_int`` (m²/s, NaN at gradient
    zeros), ``kappa_bulk`` the band-bulk scalar (m²/s) and ``D_eff`` its EBM coefficient (W m⁻² K⁻¹).
    ``irreversible_fraction`` the net-over-throughput ratio ``|∫F̄ dt| / ∫|F̄| dt`` (~0.1 — how much of
    the largely-reversible oscillating flux survives as net down-gradient transport).
    """

    phi: np.ndarray
    y: np.ndarray
    interior: np.ndarray
    jet_speed: float
    jet_lat: float
    rayleigh_kuo: bool
    times: np.ndarray
    eddy_ke: np.ndarray
    saturation_period: float
    theta_bar: np.ndarray
    F_int: np.ndarray
    G_int: np.ndarray
    kappa_profile: np.ndarray
    kappa_bulk: float
    D_eff: float
    irreversible_fraction: float
    frames: Optional[EddyFrames] = None


def _eddy_kinetic_energy(s: SWState, cell_area: float) -> float:
    """``∫ ½ h[(u−ū)² + (v−v̄)²] dA`` (J) — the eddy KE (deviation from the zonal mean).

    Takes ``cell_area`` (not the solver) so it can be recomputed from a banked :class:`EddyFrames`
    snapshot — the frame-fidelity anchor (the banked ``(h,u,v)`` reproduce the eddy-KE series).
    """
    up = s.u - s.u.mean(axis=1, keepdims=True)
    vp = s.v - s.v.mean(axis=1, keepdims=True)
    u2c = 0.5 * (up ** 2 + np.roll(up, -1, axis=1) ** 2)
    v2c = 0.5 * (vp ** 2 + np.roll(vp, -1, axis=0) ** 2)
    return float(np.sum(0.5 * s.h * (u2c + v2c)) * cell_area)


def eddy_heat_flux(state: SWState) -> np.ndarray:
    """Zonal-mean meridional eddy heat flux ``F̄(y) = ⟨v'θ'⟩`` (K·m/s), with ``v`` averaged to centres.

    ``v`` lives on north–south faces and ``θ`` at cell centres; ``v`` is averaged to the centres
    (collocated with ``θ``) before forming the eddy covariance ``⟨v'θ'⟩ = ⟨vθ⟩ − v̄·θ̄`` (primes =
    deviations from the zonal mean). Requires ``state.tracer`` set.
    """
    if state.tracer is None:
        raise ValueError("state carries no tracer (state.tracer is None)")
    v_c = 0.5 * (state.v + np.roll(state.v, -1, axis=0))
    th = state.tracer
    return (v_c * th).mean(axis=1) - v_c.mean(axis=1) * th.mean(axis=1)


def eddy_life_cycle(climate: Optional[ClimateState] = None, params: Optional[EBMParams] = None, *,
                    nx: int = 80, ny: int = 80, phi_ref_deg: float = coupler.PHI_REF_DEG,
                    n_LR: float = coupler.CHANNEL_N_LR, alpha: float = coupler.HEIGHT_PER_KELVIN,
                    taper: float = coupler.WINDOW_TAPER, spinup_periods: float = coupler.SPINUP_PERIODS,
                    release_periods: float = RELEASE_PERIODS, window_start: float = WINDOW_START,
                    eps: float = PERTURB_EPS, wavenumber: int = PERTURB_WAVENUMBER,
                    n_frames: int = 0) -> EddyFlux:
    """Spin up the jet, release it, and diagnose the emergent eddy diffusivity ``κ_eff`` — the Phase-B spine.

    Builds the Phase-4 channel + EBM height target, spins the jet up from rest **dry** (Strang-split
    thermal-relaxation + drag, early-stopping on convergence — the coupler's own logic), initialises a
    passive tracer ``θ`` to the windowed EBM temperature profile, adds a deterministic ``cos(kx)``
    v-perturbation, then runs the bare engine **forcing-off** for ``release_periods``. The band-bulk
    life-cycle integral ``κ_bulk = −∫F̄ dt / ∫θ̄_y dt`` over the window-flat interior (from
    ``window_start`` to the end) is the emergent effective diffusivity; ``D_eff`` maps it through the
    physical κ→D bridge. Returns an :class:`EddyFlux`. **Hold ``alpha`` fixed across climates** for the
    non-circularity comparison (see the module docstring).

    Parameters
    ----------
    climate : ClimateState | None
        The EBM climate whose temperature gradient forces the jet and seeds the tracer; ``None`` →
        :func:`present_day_climate` (of ``params``).
    params : EBMParams | None
        EBM parameter bundle (defaults to Earth); only used to build the default ``climate``.
    nx, ny, phi_ref_deg, n_LR, alpha, taper, spinup_periods
        The Phase-4 channel embedding + forcing (coupler defaults; ``alpha`` held fixed across climates).
    release_periods, window_start, eps, wavenumber
        Release length, life-cycle integration start (periods), and the deterministic perturbation.
    n_frames : int
        ``>0`` banks an :class:`EddyFrames` viz side-channel (``(h,u,v,θ)`` snapshots + the cumulative
        transport traces) at ``n_frames`` even time thresholds over ``[0, t_end]``. **Diagnostic-pure:**
        ``0`` (the default) leaves the κ result bit-for-bit unchanged and returns ``frames=None``.
    """
    if params is None:
        params = EBMParams()
    if climate is None:
        climate = present_day_climate(params)

    f0 = circ.coriolis_f0(phi_ref_deg)
    beta = circ.coriolis_beta(phi_ref_deg)
    L_R = np.sqrt(circ.G_EARTH * circ.H_EQUIV) / f0
    L = n_LR * L_R
    grid = uniform_grid(L, L, nx, ny)
    sw = ShallowWater(grid, circ.G_EARTH, circ.H_EQUIV, f0=f0, beta=beta)
    eta_field, _eta_profile, phi = coupler.height_target(climate, grid, sw, phi_ref_deg, alpha, taper)
    h_target = sw.H + eta_field

    # -- spin up DRY: deliberately re-implements planet.coupler.couple_jet's spin-up (no tracer / no   #
    #    conservation recording) rather than sharing it — the two diverge here (this needs the bare     #
    #    jet, the coupler the full coupled state). DRIFT RISK: if coupler's forcing/convergence logic    #
    #    changes, this loop must be updated in lockstep (it calls coupler._forcing_half/_tukey_window/   #
    #    THERMAL_RELAX_PERIODS/DRAG_PERIODS/CONVERGE_TOL directly). Not factored into a shared helper:   #
    #    extracting one touches both delicate sim modules, guarded only by the slow Phase-B test. -- #
    s = SWState(h=sw.H * np.ones((ny, nx)), u=np.zeros((ny, nx)), v=np.zeros((ny, nx)))
    inertial = 2.0 * np.pi / f0
    dt = sw.max_dt(s) * 0.5
    decay_h = np.exp(-0.5 * dt / (coupler.THERMAL_RELAX_PERIODS * inertial))
    decay_uv = np.exp(-0.5 * dt / (coupler.DRAG_PERIODS * inertial))
    n_max = int(spinup_periods * inertial / dt)
    check_every = max(1, int(round(inertial / dt)))
    u_ref = s.u.mean(axis=1)
    for it in range(1, n_max + 1):
        s = coupler._forcing_half(s, h_target, decay_h, decay_uv)
        s = sw.step(s, dt)
        s = coupler._forcing_half(s, h_target, decay_h, decay_uv)
        if it % check_every == 0:
            u_now = s.u.mean(axis=1)
            if np.max(np.abs(u_now - u_ref)) < coupler.CONVERGE_TOL:
                break
            u_ref = u_now
    u_bar = s.u.mean(axis=1)
    ji = int(np.argmax(u_bar))
    # Rayleigh–Kuo necessary condition for barotropic instability (β − u_yy changes sign in the interior).
    uyy = (np.roll(u_bar, -1) - 2.0 * u_bar + np.roll(u_bar, 1)) / grid.dy ** 2
    rk = sw.beta - uyy
    rayleigh_kuo = bool(np.any(np.diff(np.sign(rk[2:-2])) != 0))

    # -- initialise the tracer = windowed EBM θ profile; deterministic v-perturbation on the jet -- #
    y = grid.y_centers()
    T_chan = np.interp(phi, climate.latitude_deg(), climate.T)
    theta_bar = T_chan.mean() + (T_chan - T_chan.mean()) * coupler._tukey_window(phi.size, taper)
    theta2d = np.repeat(theta_bar[:, None], nx, axis=1)
    X, Y = grid.center_mesh()
    k = 2.0 * np.pi * wavenumber / grid.Lx
    w = 6.0 * grid.dy
    s = SWState(h=s.h, u=s.u, v=s.v + eps * np.cos(k * X) * np.exp(-((Y - y[ji]) / w) ** 2),
                tracer=theta2d)
    interior = coupler._tukey_window(phi.size, taper) >= 1.0 - 1e-12

    # -- release: forcing OFF; advect (h,u,v,θ); integrate the life-cycle flux over the window -- #
    times, eke = [], []
    F_int = np.zeros(ny)
    Fabs_int = np.zeros(ny)
    G_int = np.zeros(ny)
    t = 0.0
    t_prev = 0.0
    t_end = release_periods * inertial
    t_window = window_start * inertial
    # Viz side-channel (rung A) — diagnostic-pure: built only when n_frames>0, so the κ accumulators
    # above are untouched and n_frames=0 is bit-for-bit. The panel-2 traces integrate the flux per
    # latitude FIRST over the FULL release (Fint_full/Fabs_full), so temporal sloshing cancels per-y
    # but spatial structure does not leak in (see EddyFrames). Frames snapshot at even time thresholds.
    frames = None
    if n_frames > 0:
        frame_times = np.linspace(0.0, t_end, n_frames)
        f_t, f_h, f_u, f_v, f_th, f_net, f_thru = [], [], [], [], [], [], []
        Fint_full = np.zeros(ny)
        Fabs_full = np.zeros(ny)
        next_f = 0
    while t < t_end:
        dt = sw.max_dt(s, safety=0.4)
        s = sw.step(s, dt)
        t += dt
        times.append(t / inertial)
        eke.append(_eddy_kinetic_energy(s, grid.cell_area))
        F = eddy_heat_flux(s) if (t >= t_window or n_frames > 0) else None
        if t >= t_window:
            thy = np.gradient(s.tracer.mean(axis=1), y)
            wdt = t - t_prev
            F_int += F * wdt
            Fabs_int += np.abs(F) * wdt
            G_int += thy * wdt
        if n_frames > 0:
            Fint_full += F * dt
            Fabs_full += np.abs(F) * dt
            while next_f < n_frames and t >= frame_times[next_f]:
                f_t.append(t / inertial)
                f_h.append(s.h.copy()); f_u.append(s.u.copy())
                f_v.append(s.v.copy()); f_th.append(s.tracer.copy())
                f_net.append(float(np.sum(np.abs(Fint_full[interior]))))
                f_thru.append(float(np.sum(Fabs_full[interior])))
                next_f += 1
        t_prev = t

    times = np.array(times)
    eke = np.array(eke)
    if n_frames > 0:
        frames = EddyFrames(
            times=np.array(f_t), h=np.array(f_h), u=np.array(f_u), v=np.array(f_v),
            theta=np.array(f_th), x=grid.x_centers(), y=y, phi=phi, interior=interior,
            cell_area=grid.cell_area, window_start=window_start,
            net_cum=np.array(f_net), thru_cum=np.array(f_thru),
        )
    kappa_profile = np.full(ny, np.nan)
    nz = G_int != 0.0
    kappa_profile[nz] = -F_int[nz] / G_int[nz]
    kappa_bulk = tr.bulk_diffusivity(F_int, G_int, interior)
    irr = float(np.sum(np.abs(F_int[interior])) / (np.sum(Fabs_int[interior]) + 1e-300))
    return EddyFlux(
        phi=phi, y=y, interior=interior,
        jet_speed=float(u_bar[ji]), jet_lat=float(phi[ji]), rayleigh_kuo=rayleigh_kuo,
        times=times, eddy_ke=eke, saturation_period=float(times[int(np.argmax(eke))]),
        theta_bar=theta_bar, F_int=F_int, G_int=G_int,
        kappa_profile=kappa_profile, kappa_bulk=kappa_bulk,
        D_eff=float(tr.kappa_to_ebm_D(kappa_bulk)), irreversible_fraction=irr,
        frames=frames,
    )


def reduction_to_ebm_operator(eddy: EddyFlux) -> dict:
    """Test (do **not** assume) whether the resolved flux-divergence reduces to the EBM operator's *shape*.

    Compares the **resolved** mean-tracer tendency ``−∂F̄/∂y`` (the divergence of the life-cycle eddy
    flux, the actual transport the eddies accomplished) against the EBM **down-gradient** prediction
    ``(1/cosφ)∂/∂y[κ_bulk·cosφ·∂θ̄/∂y]`` — built from the band-bulk **scalar** ``κ_bulk`` and the
    smooth mean gradient (so the test reuses neither the pointwise κ — *not circular* — nor an absolute
    scale: only the **normalised shape** is compared, via the correlation over the interior). A high
    correlation ⟹ the flux *is* smooth down-gradient diffusion; a **partial** one (measured ``~0.5–0.6``
    for the strong jet, resolution-sensitive for the weak one) ⟹ the gross down-gradient tendency is
    present but the resolved
    divergence is jet-localised and does **not** *tightly* reduce at rung 1 — the honest Phase-B
    finding. (Near the vacuous edge: the smooth operator on the near-linear midlatitude gradient is
    itself structure-poor; the tight reduction only becomes sharp at rung 3.) Returns
    ``shape_correlation`` plus the two interior profiles (``resolved``, ``predicted``) for inspection.
    """
    y, interior = eddy.y, eddy.interior
    resolved = -np.gradient(eddy.F_int, y)
    predicted = tr.spherical_transport_tendency(eddy.theta_bar, eddy.phi, y, eddy.kappa_bulk)
    r = resolved[interior] - resolved[interior].mean()
    p = predicted[interior] - predicted[interior].mean()
    denom = np.sqrt(np.sum(r ** 2) * np.sum(p ** 2))
    corr = float(np.sum(r * p) / denom) if denom > 0 else float("nan")
    return dict(shape_correlation=corr, resolved=resolved[interior], predicted=predicted[interior])


def close_loop(eddy: EddyFlux, params: Optional[EBMParams] = None) -> dict:
    """One feedback pass: route the **emergent** ``D_eff`` through the Phase-A bridge + re-equilibration.

    Re-equilibrates the EBM at the flow-diagnosed ``D_eff`` (the seam Phase A validated with a synthetic
    flux, here fed the *emergent* one) and reports the **direction** of the climate response — the
    Phase-A right-signed leg, now confirmed with an emergent κ. Because ``D_eff`` is ~1000× below rung-0
    (the named magnitude edge), the re-equilibrated climate is a **degenerate near-radiative-equilibrium
    state** and is **not banked** as a two-way climate; only the sign (weaker transport ⇒ steeper
    equator-to-pole contrast) is the result. Returns ``D_eff``, the contrasts before/after, the
    ``steeper`` direction flag, and whether the re-equilibration ``converged`` (a near-zero ``D`` with
    the ice-albedo feedback can fail to settle — checked, not assumed).
    """
    if params is None:
        params = EBMParams()
    before = present_day_climate(params)
    after = present_day_climate(replace(params, D=max(eddy.D_eff, 0.0)))
    c_before = float(before.T[0] - before.T[-1])
    c_after = float(after.T[0] - after.T[-1])
    return dict(D_eff=eddy.D_eff, contrast_before=c_before, contrast_after=c_after,
                steeper=c_after > c_before, converged=bool(after.converged))


def _main() -> None:  # pragma: no cover - manual inspection of the emergent diagnostics
    steep = eddy_life_cycle(present_day_climate(EBMParams(s2=-0.48)))
    flat = eddy_life_cycle(present_day_climate(EBMParams(s2=-0.32)))
    for label, e in (("steep s2=-0.48", steep), ("flat  s2=-0.32", flat)):
        red = reduction_to_ebm_operator(e)
        print(f"[{label}] jet={e.jet_speed:5.1f} m/s RK={e.rayleigh_kuo} sat@{e.saturation_period:.1f}P "
              f"kappa_bulk={e.kappa_bulk:.3e} D_eff={e.D_eff:.5f} irr_frac={e.irreversible_fraction:.2e} "
              f"reduce_corr={red['shape_correlation']:+.2f}")
    print(f"non-circularity: kappa_flat/kappa_steep = {flat.kappa_bulk/steep.kappa_bulk:.2f} "
          f"({'PASS flat<steep' if flat.kappa_bulk < steep.kappa_bulk else 'FAIL'})")
    print("close_loop(steep):", close_loop(steep))


if __name__ == "__main__":
    _main()
