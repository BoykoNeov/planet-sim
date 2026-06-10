"""Meridional eddy heat transport — closing the two-way EBM ⇄ circulation loop (Planet rung 1, step 2).

Phase 4 (:mod:`planet.coupler`) coupled the two shared engines **one way**: the EBM's
equilibrium temperature gradient forced the shallow-water engine and a jet emerged, but the flow
did **not** feed back on the climate. This module closes the loop — the **rung-1 two-way coupler**
(plan §5): the circulation's *resolved meridional heat flux* ``⟨v'θ'⟩`` is diagnosed as an
**effective diffusivity** ``D_eff`` and fed back into the EBM's transport, which re-equilibrates to
a new climate. That is the step on the GCM staircase where the **reduction-to-diffusive-EBM** and
**poleward-heat-transport** anchors — which the one-way dry coupler honestly *cannot* support —
finally live (plan §3–4).

The A/B split — an anchored spine driven by a (tuned) eddy engine
-----------------------------------------------------------------
The deliverable separates into two parts with very different validation status:

* **(A) the feedback *machinery* — where the tight anchor lives (this module, Phase A).** Given
  *any* meridional flux ``⟨v'θ'⟩(φ)``, the κ→D bridge + re-equilibration is a clean, physical,
  testable pipeline. Its anchor is **reduction-to-EBM**: the down-gradient closure
  ``⟨v'θ'⟩ = −D_eff·∂θ̄/∂y`` has the **same form** as the EBM's own transport term
  ``D·∂/∂x[(1−x²)∂T/∂x]``, so the two-way model with a constant flow-diagnosed ``D_eff`` *is* a
  rung-0 diffusive EBM with that ``D``. To land this independent of the (messy) eddy simulation,
  Phase A drives the machinery with a **synthetic, exactly down-gradient flux** (the Phase-4
  playbook — there a synthetic off-centre gradient proved the jet tracked the climate, not the
  channel; here a synthetic flux proves the feedback reduces to the EBM). :func:`two_way_pass`
  takes the flux through an injectable ``flux_fn`` — the seam where the Phase-B eddy simulation
  plugs in unchanged.
* **(B) the eddy *engine* — where the magnitude is a named tuned scope edge (Phase B, not built
  here).** A real ``⟨v'θ'⟩`` comes from advecting θ (relaxed toward the EBM target) on the
  **barotropically unstable** Phase-4 jet (the step-0 probe established the instability exists), so
  the transport is *emergent*. But its **magnitude** is window/forcing-tuned, so rung 1 banks *the
  loop closing* and *the reduction*, with the flux magnitude flagged — not the textbook ~5–6 PW
  (an eddy/baroclinic = rung-3 number).

The κ→D bridge (physical; only the flux it is fed is tuned)
----------------------------------------------------------
The EBM transport ``D·∂/∂x[(1−x²)∂T/∂x]`` (``D`` in W m⁻² K⁻¹) is, term-by-term, the spherical
diffusion operator of a temperature field with **physical diffusivity ``κ`` (m²/s)** and a column
heat capacity ``C_col``: in ``x = sin φ`` coordinates ``∇·(C_col κ ∇T) = (C_col κ / a²)
∂/∂x[(1−x²)∂T/∂x]``, so matching gives the bridge

    D = C_col · κ / a²        (a = Earth's radius).

The capacity that survives is the **advecting layer's** — the *atmosphere's* column heat capacity
``C_atm = c_p·p_s/g`` — **not** the mixed-layer ocean ``C`` the EBM carries only to set its
relaxation *timescale* (at equilibrium ``∂T/∂t = 0`` so the ocean ``C`` drops out of the steady
balance entirely; the transport coefficient ``D`` is an independent free parameter). The bridge is
therefore **physical and citable**, with the magnitude sanity check that rung-0's ``D = 0.555``
maps to ``κ ≈ 2.2×10⁶ m²/s`` — squarely the observed midlatitude eddy heat diffusivity
(~1–5×10⁶ m²/s). The tuning in rung 1 lives entirely in the *eddy flux fed through the bridge*, not
in the bridge.

Diagnosing the diffusivity — band-bulk, over the window-flat interior
---------------------------------------------------------------------
The down-gradient closure is ``F = −κ·g`` with ``F = ⟨v'θ'⟩`` and ``g = ∂θ̄/∂y``. The robust
estimator is the **least-squares constant** :func:`bulk_diffusivity`
``κ_bulk = −Σ(F·g)/Σ(g²)`` (well-defined where the pointwise ``κ = −F/g`` blows up at gradient
zeros), restricted to the **window-flat interior** of the periodic channel — the only band where
the meridional gradient is undistorted by the Tukey taper Phase 4 needed for periodicity (the
"window-construction-dependent" caveat applies doubly to a *flux*). ``κ_bulk > 0`` ⟺ the flux is
net **down-gradient** (warm → cold), the physical sign.

Scope edges (named)
-------------------
* **Headline feedback is uniform.** :func:`two_way_pass` feeds a *uniform* band-bulk ``D_eff``
  back across the whole EBM (cleanest; the exact limit in which the closure reduces to rung 0). A
  latitude-resolved ``D_eff(φ)`` informing only the midlatitude band — the array-``D`` EBM is
  built and ready (:class:`planet.ebm.EnergyBalanceModel` accepts a callable ``D(x)``) — is the
  Phase-B diagnostic, so the headline avoids a discontinuous-``D`` splice.
* **Band-limited information.** Only the midlatitude baroclinic zone (the channel) carries eddies;
  the tropics and poles keep the rung-0 ``D``.
* **Magnitude is forcing/window-tuned** (Phase B); the *bridge* and the *reduction* are not.

Units — SI ([[shallow-water-source]], [[ebm-radiation-source]])
---------------------------------------------------------------
``κ`` in m²/s; the EBM ``D`` in W m⁻² K⁻¹; the flux ``⟨v'θ'⟩`` in K·m/s (θ is a temperature tracer,
v in m/s); ``θ̄``/``T`` in °C; latitudes in degrees; the channel ``y`` in m.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Callable, Optional

import numpy as np

from engines.fluid import ShallowWater, uniform_grid

from . import circulation as circ
from . import coupler
from .albedo import EBMParams, present_day_climate
from .ebm import ClimateState

# --------------------------------------------------------------------------- #
# The κ→D bridge constant — the atmospheric column heat capacity.
# C_atm = c_p · p_s / g : the standard tropospheric heat capacity per unit area (~10⁷ J m⁻² K⁻¹),
# the textbook value (Hartmann, *Global Physical Climatology*; Pierrehumbert, *Principles of
# Planetary Climate*). It is the capacity of the *advecting* layer — the atmosphere — not the
# EBM's mixed-layer ocean (which only sets the relaxation timescale and cancels at equilibrium).
# --------------------------------------------------------------------------- #
CP_AIR = 1004.0            # J kg⁻¹ K⁻¹ — specific heat of dry air at constant pressure
P_SURFACE = 1.013e5        # Pa         — global-mean surface pressure
C_ATM = CP_AIR * P_SURFACE / circ.G_EARTH    # J m⁻² K⁻¹ — atmospheric column heat capacity ≈ 1.04e7


# --------------------------------------------------------------------------- #
# The κ ↔ D bridge (physical; cite C_ATM above).
# --------------------------------------------------------------------------- #
def kappa_to_ebm_D(kappa):
    """Physical eddy diffusivity ``κ`` (m²/s) → EBM transport coefficient ``D`` (W m⁻² K⁻¹).

    ``D = C_atm · κ / a²`` (``a`` = Earth's radius), the term-by-term match of the EBM transport
    operator to spherical temperature diffusion (see the module docstring). Physical and citable —
    in rung 1 only the *flux* fed through this bridge is tuned, not the bridge itself.
    """
    return np.asarray(kappa, dtype=float) * C_ATM / circ.R_EARTH ** 2


def ebm_D_to_kappa(D):
    """Inverse bridge: EBM ``D`` (W m⁻² K⁻¹) → physical eddy diffusivity ``κ`` (m²/s) ``= D·a²/C_atm``.

    The magnitude sanity check on the bridge: rung-0's ``D = 0.555`` maps to ``κ ≈ 2.2×10⁶ m²/s``,
    squarely the observed midlatitude eddy heat diffusivity (~1–5×10⁶ m²/s).
    """
    return np.asarray(D, dtype=float) * circ.R_EARTH ** 2 / C_ATM


# --------------------------------------------------------------------------- #
# Diagnosing the down-gradient eddy diffusivity from a meridional flux.
# --------------------------------------------------------------------------- #
def bulk_diffusivity(flux, dtheta_dy, mask: Optional[np.ndarray] = None) -> float:
    """Band-bulk down-gradient eddy diffusivity ``κ_bulk`` (m²/s) ``= −Σ(F·g)/Σ(g²)``.

    The least-squares constant ``κ`` in the closure ``F = −κ·g`` (``F = ⟨v'θ'⟩``, ``g = ∂θ̄/∂y``):
    the optimal single down-gradient diffusivity for the diagnosed flux, robust where the pointwise
    ``κ = −F/g`` blows up at gradient zeros. ``mask`` restricts the fit to the window-flat interior
    (the only honestly-diagnosable band). ``κ_bulk > 0`` ⟺ the flux is net **down-gradient**.
    Raises if the masked mean-square gradient is zero (``κ`` undefined).
    """
    F = np.asarray(flux, dtype=float)
    g = np.asarray(dtheta_dy, dtype=float)
    if mask is not None:
        F, g = F[mask], g[mask]
    denom = float(np.sum(g ** 2))
    if denom == 0.0:
        raise ValueError("zero mean-square gradient over the band — κ_bulk is undefined")
    return float(-np.sum(F * g) / denom)


def pointwise_diffusivity(flux, dtheta_dy) -> np.ndarray:
    """Local down-gradient diffusivity ``κ(φ) = −F/g`` (m²/s), ``NaN`` where ``g ≈ 0``.

    The latitude-resolved companion to :func:`bulk_diffusivity` — a diagnostic only (mask the
    ``NaN`` gradient-zeros before using it); the robust scalar feedback is the band-bulk value.
    """
    g = np.asarray(dtheta_dy, dtype=float)
    F = np.asarray(flux, dtype=float)
    out = np.full(g.shape, np.nan)
    nz = g != 0.0
    out[nz] = -F[nz] / g[nz]
    return out


# --------------------------------------------------------------------------- #
# Channel geometry — reuse the Phase-4 embedding so the flux lives on the jet's own grid.
# --------------------------------------------------------------------------- #
def channel_geometry(phi_ref_deg: float = coupler.PHI_REF_DEG, n_LR: float = coupler.CHANNEL_N_LR,
                     ny: int = 96, taper: float = coupler.WINDOW_TAPER):
    """The midlatitude β-plane channel's ``(phi, y, dy, interior)`` — reusing the Phase-4 embedding.

    Rebuilds the coupler's channel (:func:`planet.coupler.channel_latitudes`,
    :func:`~planet.coupler._tukey_window`) so the rung-1 flux is diagnosed on the **same grid** the
    one-way jet was forced on. Returns the channel-row latitudes ``phi`` (deg), the meridional
    cell-centre coordinate ``y`` (m) and spacing ``dy`` (m), and the boolean ``interior`` mask of
    the Tukey window's flat top (taper = 1) — the only band where ``∂θ̄/∂y`` is undistorted by the
    periodicity window, hence the only band the diffusivity is honestly diagnosed over.
    """
    f0 = circ.coriolis_f0(phi_ref_deg)
    L_R = np.sqrt(circ.G_EARTH * circ.H_EQUIV) / f0
    L = n_LR * L_R
    grid = uniform_grid(L, L, ny, ny)
    sw = ShallowWater(grid, circ.G_EARTH, circ.H_EQUIV, f0=f0, beta=circ.coriolis_beta(phi_ref_deg))
    phi = coupler.channel_latitudes(grid, sw, phi_ref_deg)
    window = coupler._tukey_window(phi.size, taper)
    interior = window >= 1.0 - 1e-12
    return phi, grid.y_centers(), grid.dy, interior


def diffusive_flux(theta, y, kappa: float) -> np.ndarray:
    """A purely down-gradient meridional flux ``F = −κ·∂θ̄/∂y`` (K·m/s) for a tracer profile ``θ̄(y)``.

    The **synthetic** Phase-A stand-in for the resolved eddy flux ``⟨v'θ'⟩`` the Phase-B eddy
    simulation will supply: an *exactly* down-gradient flux with a prescribed diffusivity ``κ``,
    used to drive and validate the feedback machinery independent of the (tuned) eddy sim.
    """
    return -float(kappa) * np.gradient(np.asarray(theta, dtype=float), np.asarray(y, dtype=float))


# --------------------------------------------------------------------------- #
# The two-way feedback pass — the banked rung-1 step-2 spine.
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class TwoWayResult:
    """One pass of the two-way feedback (plain arrays — the loose-coupling currency).

    ``phi`` channel latitudes (deg); ``theta``/``dtheta_dy`` the mean tracer profile (°C) and its
    meridional gradient (°C/m) on the channel ``y``-grid; ``flux`` the meridional eddy heat flux
    ``⟨v'θ'⟩`` (K·m/s) driving the feedback; ``interior`` the window-flat mask ``κ`` is diagnosed
    over; ``kappa_eff`` the band-bulk diffusivity (m²/s) and ``D_eff`` its EBM transport
    coefficient (W m⁻² K⁻¹); ``climate_before``/``climate_after`` the EBM climates before and after
    re-equilibrating at the flow-derived ``D_eff``; ``contrast_before``/``contrast_after`` the
    equator-to-pole temperature contrast (°C) — the right-signed-response diagnostic (more
    transport ⇒ smaller contrast).
    """

    phi: np.ndarray
    theta: np.ndarray
    dtheta_dy: np.ndarray
    flux: np.ndarray
    interior: np.ndarray
    kappa_eff: float
    D_eff: float
    climate_before: ClimateState
    climate_after: ClimateState
    contrast_before: float
    contrast_after: float


def equator_pole_contrast(state: ClimateState) -> float:
    """Equator-to-pole temperature contrast ``T(0°) − T(90°)`` (°C) — shrinks as transport strengthens."""
    return float(state.T[0] - state.T[-1])


def two_way_pass(climate: Optional[ClimateState] = None, params: Optional[EBMParams] = None,
                 flux_fn: Optional[Callable[[np.ndarray, np.ndarray], np.ndarray]] = None,
                 phi_ref_deg: float = coupler.PHI_REF_DEG, n_LR: float = coupler.CHANNEL_N_LR,
                 ny: int = 96, taper: float = coupler.WINDOW_TAPER) -> TwoWayResult:
    """One two-way feedback pass: EBM climate → channel tracer gradient → eddy flux → ``D_eff`` → EBM.

    The rung-1 step-2 spine. From an EBM ``climate`` (default present-day), samples its temperature
    onto the channel as the mean tracer ``θ̄(φ)``, obtains the meridional eddy heat flux
    ``⟨v'θ'⟩`` from ``flux_fn(theta, y) -> F`` — a **synthetic** down-gradient flux in Phase A
    (the seam where the Phase-B eddy simulation plugs in) — diagnoses the **band-bulk** down-gradient
    diffusivity over the window-flat interior, maps it to an EBM coefficient ``D_eff`` (the physical
    κ→D bridge), and **re-equilibrates the EBM at the uniform ``D_eff``**. Returns a
    :class:`TwoWayResult`.

    The **default** ``flux_fn`` is the down-gradient flux rung-0's *own* diffusivity implies
    (``F = −κ₀·∂θ̄/∂y``, ``κ₀ = ebm_D_to_kappa(params.D)``): the flux for which the two-way map has
    **rung-0 as a fixed point** — calling ``two_way_pass()`` recovers ``D_eff = params.D`` and the
    rung-0 climate, the reduction anchor's cleanest expression. The headline feedback is a *uniform*
    ``D_eff`` across the whole EBM (the exact limit in which the closure ``⟨v'θ'⟩ = −D_eff·∂θ̄/∂y``
    *is* a rung-0 diffusive EBM); a band-limited ``D_eff(φ)`` is the Phase-B diagnostic.

    Parameters
    ----------
    climate : ClimateState | None
        The EBM climate to read the tracer gradient from; ``None`` → :func:`present_day_climate`.
    params : EBMParams | None
        The EBM parameter bundle (defaults to Earth); ``D_eff`` re-equilibrates ``replace(params, D=…)``.
    flux_fn : callable ``(theta, y) -> F`` | None
        The meridional eddy heat flux ``⟨v'θ'⟩`` on the channel ``y``-grid; ``None`` → the rung-0
        diffusive flux (the fixed-point driver).
    phi_ref_deg, n_LR, ny, taper
        The channel embedding (the Phase-4 defaults).
    """
    if params is None:
        params = EBMParams()
    if climate is None:
        climate = present_day_climate(params)
    if flux_fn is None:
        kappa0 = float(ebm_D_to_kappa(params.D))
        flux_fn = lambda theta, y: diffusive_flux(theta, y, kappa0)

    phi, y, _dy, interior = channel_geometry(phi_ref_deg, n_LR, ny, taper)
    theta = np.interp(phi, climate.latitude_deg(), climate.T)      # EBM T sampled on the channel
    dtheta_dy = np.gradient(theta, y)
    flux = np.asarray(flux_fn(theta, y), dtype=float)
    kappa_eff = bulk_diffusivity(flux, dtheta_dy, interior)
    D_eff = float(kappa_to_ebm_D(kappa_eff))
    after = present_day_climate(replace(params, D=D_eff))
    return TwoWayResult(
        phi=phi, theta=theta, dtheta_dy=dtheta_dy, flux=flux, interior=interior,
        kappa_eff=kappa_eff, D_eff=D_eff,
        climate_before=climate, climate_after=after,
        contrast_before=equator_pole_contrast(climate),
        contrast_after=equator_pole_contrast(after),
    )
