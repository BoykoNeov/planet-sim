"""Circulation-informed precipitation — wiring the precip storm-track band to the emergent jet
(Planet rung 1, step 3).

Phase 2's precipitation (:mod:`planet.precip`) is a **prescribed** kinematic parameterization: a
fixed Gaussian band pattern × a global-mean Clausius–Clapeyron amplitude, with the midlatitude
**storm-track band hand-placed** at the cited :data:`planet.precip.MIDLAT_CENTER_DEG` (50°). Once the
rung-1 circulation exists (the one-way coupler's emergent jet, :mod:`planet.coupler`), the plan's
deep-end hook applies: make the pattern **circulation-informed — rain where the flow puts the storm
track — still without any moisture physics** (the §5 staircase, the rung-1/2 array-seam enhancement).
This module is that wiring. It centres the storm-track band on the **emergent jet latitude** instead
of the prescribed constant, so the band tracks the *dynamically-selected* circulation rather than a
hardcoded number. See [[precip-parameterization-source]], [[planet-phase4-coupler]], [[planet-plan]].

What this banks — and the honest scope of each leg (the trade, not a ranking)
-----------------------------------------------------------------------------
The work was **de-risked in two throwaway spikes before the module** (``outputs/rung1_circprecip*``,
gitignored — this project's discipline); the spikes set the headline and killed one anchor.

* **(banked, the seam + the reduction) The storm-track centre is now a parameter the circulation
  sets, and the rung-0 field is recovered exactly.** :func:`planet.precip.precip_pattern` takes a
  ``midlat_center_deg`` (default = the cited constant → rung-0 *bit-for-bit by construction*, the
  plumbing-not-an-independent-test honesty of :func:`planet.transport.two_way_pass`); this module
  feeds it the emergent ``jet_lat``. When the jet sits at the prescribed centre the circulation-informed
  field *is* the rung-0 field — the reduction.
* **(banked, the mechanism) The band tracks a flow-selected latitude — migration, shown via the
  synthetic-gradient playbook.** The decisive demonstration is the coupler's own non-circularity proof
  (an off-centre *synthetic* EBM gradient makes the jet follow it): the precip band then follows the
  jet, so the storm-track rain is anchored to the dynamically-selected circulation, not to the channel
  or a constant. Anchoring to ``jet_lat`` (not the EBM ``gradient_peak_lat``) is what makes this a
  **flow** response: the two nearly coincide at present-day (within ~2–3°), and the gap only opens for an
  off-centre gradient (the jet is pulled toward the channel — the circulation content).
* **(the rung-1 FINDING — named, NOT an accuracy gain) The dry circulation cannot *refine* the rain
  location at rung 1; it is a *trade*.** Two honest edges, both established in the spikes:
  - **It is a trade, not better.** Rung-0's 50° is **observation-calibrated** (the real storm track,
    [[precip-parameterization-source]]); the circulation-informed centre is the model's **own** jet at
    ~42°, ~8° **equatorward** — because the dry single-layer model's gradient peak (hence its jet) sits
    equatorward of Earth's observed storm track (the Phase-4 channel's known equatorward bias, which
    *excludes the ice cliff*). So the relocation **trades observational calibration for internal
    consistency**; it does **not** improve the band, and may land slightly worse against observations.
    This is why rung-0 :mod:`planet.precip` stays the **default** in the biome map and demos — the
    circulation-informed field is **opt-in**.
  - **For realistic knobs the band barely moves.** The EBM gradient peak (spike #1) is pinned at
    ~43–46° across realistic obliquity/CO₂/S₀ (the channel excludes the ice cliff, so the smooth midlat
    gradient maximum is robust) — it shifts only when a near-snowball cooling drags the ice cliff into
    the channel. So realistic **migration is mechanism-only**: the band tracking is real, but the demo
    is the *mechanism* (decisive only under a synthetic off-centre gradient), not a large realistic
    rain-belt shift.
* **(rejected anchor — a clean negative) "Rain where the flow converges" does NOT bank a band *shape*
  at rung 1.** The literal reading — centre the rain on the eddy heat-flux **convergence** ``−∂F̄/∂y``
  — was tested in spike #2 and **rejected**: on the released eddy life cycle the resolved convergence is
  **near-vacuous in the channel interior** (the near-linear midlatitude gradient diffuses to ~0
  divergence) with what little structure there is being a **window-taper edge dipole**, not a physical
  storm-track convergence. Anchoring a rain band to that would be building on a numerical artifact. This
  is the **same boundary** :mod:`planet.eddy_flux` already found: the eddy flux-divergence becomes
  non-vacuous only at **rung 3** (a strong baroclinic flux), where the live, *shape*-resolving precip
  refinement belongs — and the **geometry it needs is already delivered**
  (:func:`planet.transport.spherical_transport_tendency`). So rung 1 wires the *position* seam; the
  *shape* refinement is named for rung 3.

Scope edges (named)
-------------------
* **Midlatitude band only.** Only the storm-track centre is circulation-informed; the **ITCZ and
  subtropical** centres stay prescribed — the midlatitude β-plane channel does not represent the Hadley
  cell (the tropics/subtropics are out of its band).
* **Position, not amplitude.** The band's amplitude and width stay the prescribed Gaussian; the
  **wet-get-wetter** pattern amplification stays deferred (scaling amplitude by eddy/jet strength is
  un-anchored at rung 1 — the Phase-B ``κ`` *magnitude* is itself named-not-banked).
* **Large equatorward displacement merges the storm-track band into the ITCZ.** Moving the centre well
  equatorward (toward ~30°) pushes the Gaussian's tail over the ITCZ and shallows the subtropical
  trough; the band structure holds across the realistic-to-modest range (centre ≳ 36°) but degrades
  below it — so the band-tracking is asserted only where the trough survives.

Units — P in cm/yr (Whittaker's axis, the precip module's unit), latitudes in degrees.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import precip
from .coupler import CoupledJet
from .ebm import ClimateState


def storm_track_center(jet: CoupledJet) -> float:
    """The circulation-informed storm-track band centre (degrees) — the emergent ``|jet.jet_lat|``.

    The midlatitude precipitation band's centre, read from the one-way coupler's emergent jet
    (:class:`planet.coupler.CoupledJet`) rather than the prescribed
    :data:`planet.precip.MIDLAT_CENTER_DEG`. Magnitude (hemispherically symmetric pattern), so the
    northern channel's ``jet_lat`` sets both hemispheres' storm tracks. This is the one value the
    circulation contributes — anchored to the jet (a flow feature), the leg that makes the band
    *circulation*-informed rather than merely EBM-gradient-informed.
    """
    return abs(float(jet.jet_lat))


def circulation_informed_precip(state: ClimateState, jet: CoupledJet) -> np.ndarray:
    """Circulation-informed precipitation ``P(φ)`` (cm/yr): the storm-track band on the emergent jet.

    Like :func:`planet.precip.precip_field`, but the midlatitude band is centred on the emergent jet
    latitude (:func:`storm_track_center`) instead of the prescribed constant. Reduces to
    :func:`planet.precip.precip_field` **exactly** when the jet sits at
    :data:`planet.precip.MIDLAT_CENTER_DEG` (the reduction). The ITCZ/subtropics and the band
    amplitudes are unchanged (prescribed); only the storm-track *position* is circulation-informed.

    ``jet`` is the coupled-jet record for this climate, ``jet = planet.coupler.couple_jet(state)``
    (kept out of this function so it stays composable and fast — the coupler spin-up is the cost). The
    rung-0 :mod:`planet.precip` remains the default elsewhere; this is the **opt-in** trade (see the
    module docstring — it trades observational calibration for internal consistency, not an accuracy
    gain).
    """
    return precip.precipitation(state.latitude_deg(), state.global_mean_T,
                                midlat_center_deg=storm_track_center(jet))


@dataclass(frozen=True)
class Relocation:
    """The storm-track relocation — the rung-0 vs circulation-informed *trade*, made tangible (plain arrays).

    ``phi`` the latitudes (deg). ``center_rung0`` the prescribed (observation-calibrated) storm-track
    centre and ``center_circ`` the emergent jet's (the model's own); ``displacement`` their signed
    difference (``center_circ − center_rung0``, **negative = equatorward** — the dry channel's known
    bias). ``precip_rung0``/``precip_circ`` the two precip fields (cm/yr) on the climate's grid. The
    displacement is the headline number of the trade — *not* an error reduction.
    """

    phi: np.ndarray
    center_rung0: float
    center_circ: float
    displacement: float
    precip_rung0: np.ndarray
    precip_circ: np.ndarray


def relocate(state: ClimateState, jet: CoupledJet) -> Relocation:
    """Build the :class:`Relocation` — the rung-0 vs circulation-informed precip fields and the centres.

    Pairs the prescribed-centre field (:func:`planet.precip.precip_field`) with the circulation-informed
    one (:func:`circulation_informed_precip`) for the same climate, so the **trade** is explicit: the
    storm track moves from the observation-calibrated :data:`planet.precip.MIDLAT_CENTER_DEG` to the
    model's emergent ``jet_lat`` (~8° equatorward for present-day Earth). For the demo/notebook and the
    migration test; the signed ``displacement`` is the one quantity the relocation banks.
    """
    center_circ = storm_track_center(jet)
    return Relocation(
        phi=state.latitude_deg(),
        center_rung0=float(precip.MIDLAT_CENTER_DEG),
        center_circ=center_circ,
        displacement=center_circ - float(precip.MIDLAT_CENTER_DEG),
        precip_rung0=precip.precip_field(state),
        precip_circ=circulation_informed_precip(state, jet),
    )
