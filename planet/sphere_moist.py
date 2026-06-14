"""Emergent ITCZ rain — the full-sphere moisture budget co-located with the energy-flux equator (rung 2.x).

Rung 2.x (:mod:`planet.sphere_ebm`) lifts the EBM to the full sphere and locates the ITCZ at the
**energy-flux equator** (EFE) — the zero of the atmospheric energy transport ``H(x)``. But its precip wire
(:func:`planet.sphere_ebm.itcz_informed_precip`) only **relocates a prescribed Gaussian band** to ``φ_EFE``:
a *dry* model painting a rain belt at the energetically-warmest latitude, not rainfall that falls out of a
water budget. This module closes that gap. It carries the rung-2 **column moisture budget**
(:mod:`planet.moist`) onto the full sphere and lets the ITCZ rain **emerge from a conserving ``P − E``
budget** whose convergence maximum **sits on the EFE** — *rained*, not *painted*.

The architectural upgrade is the headline: rain is now a conserving budget, not a band
-----------------------------------------------------------------------------------------
The atmospheric water budget in steady state is ``P − E = −∇·(moisture transport)`` (the moisture
convergence). On the full sphere it has two pieces, exactly as the hemisphere rung-2 budget does:

* the down-gradient **eddy** convergence ``(D/c_p)·∂/∂x[(1−x²) ∂q/∂x]`` with ``q = RH·q_sat(T)`` — the
  same diffusive operator the EBM transports heat with, ``L`` cancelled (:mod:`planet.moist`); and
* the mean **Hadley** convergence ``−∂F/∂x`` of an overturning moisture flux ``F`` — but now a **two-cell**
  circulation whose ascending branch is anchored at the **EFE** (:func:`hadley_streamfunction`), so the
  cell migrates with the energetics instead of being pinned to the equator.

Both are built in **conservative face-flux form** with insulated (Neumann 0) poles at **both** ends (the
full sphere has two real poles, no equatorial-symmetry boundary), so ``∫(P − E) dx = 0`` to **machine
precision** for *any* asymmetric cell — the ITCZ convergence is paid for, latitude by latitude, by the
descending dry belts. That is the genuinely-new "what": the migrated rain belt is a **budget**, not a
prescribed band, and the books balance exactly.

Two real nuggets — one positive falsifiable check, one clean negative result
----------------------------------------------------------------------------
The genuinely-*emergent* content of this rung is modest and honest: it is largely rung-2.x's **EFE
displacement** recombined with the rung-2 **Hadley convergence** into a conserving budget. Two findings
survive that are not given for free:

1. **Co-location of the NET ``P − E`` on the EFE — a falsifiable check (real-but-loose).** That the
   *Hadley* term converges at its own ascent is by-construction. That the **net** (eddy **+** Hadley)
   ``P − E`` peaks on the EFE is **not** guaranteed: the down-gradient eddy term *exports* moisture from
   the warm, moist EFE (it is backwards there, the rung-2 ITCZ trade), so the prescribed cell must **beat
   that export at the displaced latitude**. It does — by a comfortable ~2–3× margin for the calibrated
   ``HADLEY_STRENGTH`` (the eddy export is ~‑110 cm/yr at the EFE, the Hadley convergence ~+285) — so the
   net rain maximum lands on the EFE to within ~1°. True for the calibrated strength, checked, not assumed.

2. **The displaced-ITCZ intensification is GEOMETRIC, not emergent ``q(T)`` (a clean negative result).**
   When the ascent displaces while the subtropical descent edges stay pinned at ±``edge_deg`` (the physical
   cross-equatorial cell: the far cell widens, the near cell narrows), the ITCZ rain peak grows. **That
   growth is a prescribed-geometry artifact** — the narrowing near-side cell concentrates the convergence —
   **not** a moisture-content effect: replacing ``q(T)`` with a hemispherically-symmetric ``q`` leaves the
   peak unchanged, and letting the edges migrate with the ascent (constant cell width) removes the
   intensification entirely. So the warm-hemisphere ITCZ is **not** more intense *because* it is moister in
   this model. (The wet-NH / dry-SH **dipole** that accompanies the shift is likewise **displacement**-
   driven — it is present at full strength with a symmetric ``q`` — its *direction*, toward the warm
   hemisphere, is by-construction.) The peak intensity is therefore **not featured** as a result here; the
   only clean emergent ``q(T)`` signature is the **warming response** (the ITCZ convergence intensifies at
   the ~Clausius–Clapeyron rate under uniform warming — the rung-2 Hadley-fix nugget, re-confirmed on the
   full sphere, see :mod:`planet.moist`).

A SIBLING — ``moist.py`` / ``sphere_ebm.py`` untouched
------------------------------------------------------
Like every rung since 2.5, this is a **new model alongside** the existing ones, not an edit of them. It
reuses :mod:`planet.moist`'s thermodynamics (:func:`~planet.moist.specific_humidity`, the C–C ``q_sat``,
``HADLEY_STRENGTH``/``HADLEY_EDGE_DEG``, the cm/yr conversion) and re-derives the conservative operators on
the doubled grid. "Re-validate the hemisphere budget" is therefore a **cross-model reduction check**: at
``φ_EFE = 0`` the full-sphere operators reproduce :func:`planet.moist.moisture_convergence` /
:func:`planet.moist.moisture_budget` on the northern hemisphere to **machine precision** (a symmetric ``q``
makes the equatorial face-flux exactly zero, so the full-sphere stencil collapses to the hemisphere's
equatorial-symmetry boundary), the SW↔QG-bridge pattern.

Named scope edges (the walls this rung does NOT move)
-----------------------------------------------------
* **The cell is PRESCRIBED.** Its strength (``HADLEY_STRENGTH``) and shape are the named sub-grid wall;
  the *fully emergent* overturning (a streamfunction set by the energy transport through a gross-moist-
  stability closure, ``Ψ ∝ H/GMS``) is **rung 3+** — it double-counts the mean transport already lumped
  into the EBM's effective ``D`` and needs a GMS closure this column model lacks. Anchoring the ascent at
  the EFE is a *placement*, not a derivation (it is what makes "rain co-locates with the EFE" true).
* **The asymmetry is IMPOSED** (a cross-equatorial Q-flux / antisymmetric albedo, the rung-2.x knob), not
  an ocean.
* **The subtropical desert stays mislocated** — the hyper-peaked fixed-RH C–C ``q`` puts the dry belt
  equatorward of the canonical 25–35°, the same rung-2 wall the mean cell did not fix (and does not here).
* **``P − E``, not a full ``P``.** As in :mod:`planet.moist`, no honest zonal ``E`` keeps ``P ≥ 0``, so the
  budget reports the sign-meaningful convergence; the "rain belt" is the positive convergence, and rung-0
  :mod:`planet.precip` remains the precipitation *map*.

Validation triad (plan §3)
--------------------------
* **Tight (structural).** The cross-model reduction to the hemisphere :mod:`planet.moist` at ``φ_EFE = 0``
  (machine-exact); ``∫(P − E) dx = 0`` (machine-exact area-mean, conservative flux form, even for an
  asymmetric cell); a symmetric climate ⟹ a symmetric ``P − E`` peaking at the equator.
* **Real-but-loose (the unlock).** The net ``P − E`` peaks on the EFE (Hadley beats the eddy export at the
  displaced latitude — checked, ~2–3× margin, not assumed); the warming response intensifies at the ~C–C
  rate (transported from the Hadley fix).
* **Plumbing / named negatives.** Co-location *at all* is by-construction (ascent anchored at the EFE); the
  displaced-ITCZ peak intensification is **geometric** (pinned-edge cell narrowing), not emergent ``q``; the
  wet/dry dipole is displacement-driven with a by-construction direction.

Units — SI internally; ``P − E`` reported in **cm/yr**; ``x = sin φ`` on [−1, 1]; latitudes in degrees.
Sources (the ``[[…-source]]`` discipline): the energetic-ITCZ framework (Kang et al. 2008; Bischoff &
Schneider 2014) for the EFE/Hadley-ascent identification; Held & Soden 2006 for the rich-get-richer P−E
scaling; the diffusive moist EBM (Flannery 1984; Hwang & Frierson 2010) for the latent-diffusion closure —
all already pinned in :mod:`planet.moist` / :mod:`planet.sphere_ebm`. Extends [[planet-rung2x-itcz]],
[[planet-rung2-hadley-fix]], [[moist-ebm-source]].
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .ebm import D_TRANSPORT
from .moist import (
    HADLEY_EDGE_DEG, HADLEY_STRENGTH, RH_DEFAULT, _KGM2S_TO_CMYR,
    energy_constrained_rate, specific_humidity,
)
from .sphere_ebm import SphereClimate
from .transport import CP_AIR


# --------------------------------------------------------------------------- #
# Full-sphere conservative operators — both ends are real poles (Neumann 0).
# --------------------------------------------------------------------------- #
def _sphere_flux_divergence(field: np.ndarray, x: np.ndarray) -> np.ndarray:
    """The full-sphere spherical operator ``∂/∂x[(1 − x²) ∂field/∂x]`` in conservative flux form.

    Differences of interior **face** fluxes ``(1 − x_face²)·Δfield/Δx`` on the uniform ``x = sin φ`` grid,
    with **insulated (Neumann 0)** faces at **both** poles (no equatorial-symmetry boundary — the equator
    is interior on the full sphere). The discrete operator sums to **exactly zero** (``Σ div·Δx = 0``: the
    interior faces telescope and both boundary faces are zero), the machine-exact ``∫E = ∫P`` plumbing leg.
    The hemisphere counterpart :func:`planet.moist._spherical_flux_divergence` hard-codes the equatorial
    Neumann-0 at index 0; here both ends are poles. For a symmetric ``field`` the equatorial face-flux
    vanishes, so this stencil reduces to the hemisphere one on the NH (the cross-model reduction).
    """
    f = np.asarray(field, dtype=float)
    x = np.asarray(x, dtype=float)
    dx = x[1] - x[0]
    x_face = 0.5 * (x[:-1] + x[1:])                    # interior cell faces
    g = 1.0 - x_face ** 2                              # area weight (1 − x²) at faces (→ 0 at both poles)
    flux = g * (f[1:] - f[:-1]) / dx                   # interior-face down-gradient flux × (1 − x²)
    div = np.zeros_like(f)
    div[1:-1] = (flux[1:] - flux[:-1]) / dx            # interior cells: difference of bounding faces
    div[0] = flux[0] / dx                              # south pole: lower face is Neumann 0
    div[-1] = -flux[-1] / dx                           # north pole: upper face is Neumann 0
    return div


def _mean_flux_convergence(flux: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Convergence ``−∂F/∂x`` of a prescribed advective transport ``F(x)`` in conservative face form.

    The mean-circulation analogue of :func:`_sphere_flux_divergence` (full-sphere version of
    :func:`planet.moist._mean_flux_convergence`): ``F`` is a first-order transport, so the convergence is
    differences of interior **face** fluxes ``F_face = ½(F_i + F_{i+1})`` with the domain-boundary (polar)
    fluxes set to **zero**. Sums to **exactly zero** (``Σ conv·Δx = 0``) — moisture-mass conservation — for
    any ``F``, including an asymmetric NH/SH cell (both poles carry zero flux regardless).
    """
    f = np.asarray(flux, dtype=float)
    x = np.asarray(x, dtype=float)
    dx = x[1] - x[0]
    f_face = 0.5 * (f[:-1] + f[1:])                    # interior cell faces
    conv = np.zeros_like(f)
    conv[1:-1] = -(f_face[1:] - f_face[:-1]) / dx      # interior cells: difference of bounding faces
    conv[0] = -(f_face[0] - 0.0) / dx                  # south pole: boundary flux is 0
    conv[-1] = -(0.0 - f_face[-1]) / dx                # north pole: boundary flux is 0
    return conv


# --------------------------------------------------------------------------- #
# The two-cell Hadley overturning — ascent anchored at the EFE (migrates with it).
# --------------------------------------------------------------------------- #
def hadley_streamfunction(x: np.ndarray, phi_efe: float = 0.0,
                          edge_deg: float = HADLEY_EDGE_DEG) -> np.ndarray:
    """The signed two-cell overturning profile ``Ψ(x)`` (dimensionless) with its ascent at the EFE.

    Two Hadley cells meeting at the energy-flux equator ``x_efe = sin(φ_EFE)``: each low-level branch flows
    **toward** the ascent (equatorward, in the energetic sense), so ``Ψ`` is **positive south** of the EFE
    (northward low-level flow) and **negative north** of it (southward), vanishing at the EFE and at the
    fixed subtropical descent edges ``±edge_deg``. Each cell is the smooth cubic bump
    ``(27/4)·u·(1−u)²`` (``ψ(edge)=ψ′(edge)=0`` — a smooth merge to the extratropics; a finite ascent slope
    at the EFE) used by the hemisphere :func:`planet.moist.hadley_streamfunction`, here **reflected into two
    cells of asymmetric width**: with the descent edges pinned at fixed latitudes, a displaced ascent makes
    the **cross-equatorial cell wide** and the near cell narrow (the physical winter/summer-cell asymmetry).

    The moisture flux rides on this as ``F = strength·Ψ·q`` (:func:`hadley_moisture_convergence`). Reduces to
    the hemisphere single cell on the NH when ``φ_EFE = 0``. A **prescribed** kinematic profile — the mean
    circulation is imposed, and anchoring the ascent at the EFE is the *placement* that co-locates the rain
    with the energetics (the emergent cell = rung 3+; see the module docstring).
    """
    x = np.asarray(x, dtype=float)
    x_efe = math.sin(math.radians(float(phi_efe)))
    x_edge_n = math.sin(math.radians(float(edge_deg)))
    x_edge_s = math.sin(math.radians(-float(edge_deg)))
    w_n = max(x_edge_n - x_efe, 1e-12)                 # near/far widths set by the pinned descent edges
    w_s = max(x_efe - x_edge_s, 1e-12)
    nh = (x > x_efe) & (x < x_edge_n)
    sh = (x < x_efe) & (x > x_edge_s)
    psi = np.zeros_like(x)
    u_n = (x[nh] - x_efe) / w_n
    u_s = (x_efe - x[sh]) / w_s
    psi[nh] = -(27.0 / 4.0) * u_n * (1.0 - u_n) ** 2    # north of EFE: low-level flow southward (Ψ < 0)
    psi[sh] = +(27.0 / 4.0) * u_s * (1.0 - u_s) ** 2    # south of EFE: low-level flow northward (Ψ > 0)
    return psi


def eddy_convergence(climate: SphereClimate, RH: float = RH_DEFAULT,
                     D: float = D_TRANSPORT) -> np.ndarray:
    """The full-sphere down-gradient **eddy** moisture convergence ``P − E`` (cm/yr).

    ``(D/c_p)·∂/∂x[(1 − x²) ∂q/∂x]`` with ``q = RH·q_sat(T)`` (:func:`planet.moist.specific_humidity`),
    built with :func:`_sphere_flux_divergence` so ``∫ dx = 0`` exactly. The full-sphere twin of
    :func:`planet.moist.moisture_convergence`: positive (convergence) in the extratropics, **negative
    (export) at the warm tropics** — backwards at the ITCZ, the rung-2 trade, which the Hadley term repairs.
    ``D`` defaults to the rung-0 transport (``ClimateState`` does not carry its own — pass it for a
    non-default-``D`` world, the :mod:`planet.moist` "transport ``D``" scope edge). Scalar ``D`` only.
    """
    if callable(D):
        raise TypeError("eddy_convergence takes a scalar transport D; the array-D(x) EBM is a non-goal "
                        "for the column moisture diagnostic")
    q = specific_humidity(climate.T, RH)
    conv_mass = (float(D) / CP_AIR) * _sphere_flux_divergence(q, climate.x)   # kg m⁻² s⁻¹
    return conv_mass * _KGM2S_TO_CMYR                                         # cm/yr


def hadley_moisture_convergence(climate: SphereClimate, RH: float = RH_DEFAULT,
                                strength: float = HADLEY_STRENGTH,
                                edge_deg: float = HADLEY_EDGE_DEG) -> np.ndarray:
    """The two-cell mean **Hadley** moisture convergence ``P − E`` (cm/yr) — ascent on the EFE.

    The northward overturning moisture flux ``F(x) = strength·Ψ(x)·q(x)`` (the dry-upper-branch
    approximation ``Δq ≈ q_surface``) with the two-cell :func:`hadley_streamfunction` ``Ψ`` anchored at
    ``climate.phi_efe`` and ``q = RH·q_sat(T)``; its convergence ``−∂F/∂x`` (:func:`_mean_flux_convergence`,
    conservative ⟹ ``∫ = 0`` machine-exact) is **positive at the EFE** (the ascent gains moisture) and
    **negative under both descents** (the dry belts), vanishing poleward of the cell edges. The full-sphere,
    migrating-ascent twin of :func:`planet.moist.hadley_moisture_convergence`. ``strength`` is the named
    **prescribed wall**; see the module docstring on what is by-construction (the convergence structure, the
    displaced-peak intensification) versus emergent (only the ~C–C warming response).
    """
    q = specific_humidity(climate.T, RH)
    psi = hadley_streamfunction(climate.x, climate.phi_efe, edge_deg)
    flux = float(strength) * psi * q                                         # kg m⁻² s⁻¹
    return _mean_flux_convergence(flux, climate.x) * _KGM2S_TO_CMYR          # cm/yr


# --------------------------------------------------------------------------- #
# The full-sphere moisture budget — the frozen diagnostic (plain arrays).
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class SphereMoistureBudget:
    """The full-sphere column moisture budget for one :class:`~planet.sphere_ebm.SphereClimate` (plain arrays).

    ``phi`` latitudes (deg, −90 → +90); ``q`` the diagnostic specific humidity (kg/kg); ``p_minus_e`` the
    **net** moisture convergence ``P − E`` (cm/yr), with its ``p_minus_e_eddy`` and ``p_minus_e_hadley``
    components; ``phi_efe`` the energy-flux equator the cell is anchored on (deg). The banked diagnostics:
    ``rain_max_lat`` the latitude of the net ``P − E`` maximum (the emergent ITCZ — co-located with
    ``phi_efe`` to ~1°, the falsifiable check); ``net_p_minus_e`` the area integral ``∫(P − E) dx``
    (≈ 0 to machine precision); ``nh_p_minus_e`` / ``sh_p_minus_e`` the hemispheric integrals (the
    displacement-driven wet/dry dipole — equal and opposite, summing to ``net_p_minus_e``); ``hadley``
    whether the mean cell was added. **Honest reading (module docstring):** co-location is by-construction
    (ascent anchored at the EFE); the dipole *direction* is by-construction; the displaced-peak intensity is
    geometric, not emergent ``q``.
    """

    phi: np.ndarray
    q: np.ndarray
    p_minus_e: np.ndarray
    p_minus_e_eddy: np.ndarray
    p_minus_e_hadley: np.ndarray
    phi_efe: float
    rain_max_lat: float
    net_p_minus_e: float
    nh_p_minus_e: float
    sh_p_minus_e: float
    hadley: bool


def sphere_moisture_budget(climate: SphereClimate, RH: float = RH_DEFAULT,
                           D: float = D_TRANSPORT, hadley: bool = True,
                           strength: float = HADLEY_STRENGTH,
                           edge_deg: float = HADLEY_EDGE_DEG) -> SphereMoistureBudget:
    """Build the :class:`SphereMoistureBudget` — the emergent ``P − E`` and the co-location diagnostics.

    Reads the full-sphere :class:`~planet.sphere_ebm.SphereClimate` (its temperature, grid and ``phi_efe``)
    and returns the eddy + (opt-in, **default-on** here — the point of the rung) Hadley moisture budget. A
    **pure diagnostic** — it does not modify the climate. ``hadley=False`` gives the eddy-only budget (which
    is backwards at the ITCZ); ``hadley=True`` (default) adds the two-cell mean circulation anchored on the
    EFE so the net rain belt co-locates with the energetics. ``D`` is threaded to :func:`eddy_convergence`
    (defaults to rung-0 ``D_TRANSPORT``; pass the climate's own ``D`` for a non-default-``D`` world).

    The net ``P − E`` maximum (``rain_max_lat``) lands on ``phi_efe`` to ~1° because the prescribed cell
    beats the eddy export there; ``∫(P − E) dx`` is ~0 to machine precision (area-mean, the conservative
    flux form); the hemispheric integrals are an equal-and-opposite displacement-driven dipole.
    """
    phi = climate.latitude_deg()
    q = specific_humidity(climate.T, RH)
    eddy = eddy_convergence(climate, RH, D)
    hadley_term = (hadley_moisture_convergence(climate, RH, strength, edge_deg)
                   if hadley else np.zeros_like(eddy))
    pme = eddy + hadley_term

    # ∫(P − E) dx on the equal-area (uniform-Δx) grid is the area MEAN — the EBM's `total` convention under
    # which the conservative flux divergence sums to *exactly* zero (np.trapezoid would leave a boundary
    # residual). Hemispheric integrals are area-means × 2 (the half-domain length), so they sum to net × 2.
    net = float(np.mean(pme))
    nh = climate.x > 0.0
    nh_int = float(np.mean(np.where(nh, pme, 0.0)) * 2.0)
    sh_int = float(np.mean(np.where(~nh, pme, 0.0)) * 2.0)
    rain_max_lat = float(phi[int(np.argmax(pme))])

    return SphereMoistureBudget(
        phi=phi, q=q, p_minus_e=pme, p_minus_e_eddy=eddy, p_minus_e_hadley=hadley_term,
        phi_efe=float(climate.phi_efe), rain_max_lat=rain_max_lat,
        net_p_minus_e=net, nh_p_minus_e=nh_int, sh_p_minus_e=sh_int, hadley=bool(hadley),
    )
