"""Diagnostic precipitation parameterization — the second input to the biome map (Planet Phase 2).

The Phase-2 payoff (the biome map, :mod:`planet.biomes`) needs two climate inputs: the
**temperature** ``T(φ)`` the EBM produces (:mod:`planet.ebm`) and a **precipitation**
field ``P(φ)``. This module supplies the latter. It is, deliberately and explicitly, a **prescribed
kinematic parameterization, NOT a simulated water cycle** (plan §3, the honesty flag): neither the
EBM (energy only) nor the later single-layer dry shallow-water engine (no thermodynamic moisture
variable) produces precipitation — a real hydrological cycle is moist thermodynamics, the GCM tar
pit (the rung-2 deferral on the §5 staircase). So ``P(φ)`` here **encodes the known observed
precip-by-latitude structure**, it does not derive it.

The two pieces — pattern (where) and amplitude (how much)
---------------------------------------------------------
``P(φ, T̄) = pattern(φ) · CC(T̄)``, a clean separation of the two physics:

* **pattern(φ)** — *where* it rains, set by the **general circulation** (Hadley/Ferrel cells), which
  is fixed geography in this annual-mean zonal-mean model: a wet equatorial **ITCZ**, dry
  **subtropics ~25–30°** (the world's great deserts, under the descending Hadley branch), wet
  **midlatitude ~50°** storm tracks, and dry **poles**. Built as a sum of Gaussian bands in latitude
  (the ITCZ peak + the two midlatitude peaks over a low polar/desert baseline); the subtropical
  minimum *emerges* as the trough between the ITCZ and midlatitude bands. Calibrated to the observed
  zonal-mean structure ([[precip-parameterization-source]]).
* **CC(T̄)** — *how much*, the thermodynamic amplitude. A warmer atmosphere holds more moisture, so
  the whole hydrological cycle intensifies. This is the **Clausius–Clapeyron** relation: saturation
  water-vapour content rises ≈ **7 %/K**, so ``CC(T̄) = exp(k·(T̄ − T_ref))`` with ``k = 0.07`` scales
  the pattern by the **global-mean** temperature relative to the present-day reference (where the
  pattern is calibrated, so ``CC = 1`` at present).

**Why global-mean T̄, not local T (a deliberate refinement).** The latitudinal *pattern* is set by
circulation, not by local warmth — the equator is wet because it is the ITCZ, not merely because it
is hot. Scaling the present-day pattern by *local* temperature against a fixed reference would wrongly
re-attribute that to local warmth and **over-amplify** the warm equator (a ×2.6 boost at 28 °C),
breaking the calibration. So the **circulation-set pattern is held fixed** and only the **global
moisture amplitude** responds to T̄ — the honest pattern/amplitude decomposition.

Validation triad (plan §3) — what is asserted tight vs loose
------------------------------------------------------------
* **Analytical/structural (tight).** The pattern has the right *qualitative band structure*: the
  equator (ITCZ) is the wettest, the subtropics (~25–30°) are a local dry minimum (drier than both
  the equator and the midlatitudes), and the poles are dry. These are partition-of-latitude facts,
  asserted firmly (:mod:`planet.tests.test_precip`).
* **Conservation — a *consistency* check, honestly weaker (named).** A prescribed precip field obeys
  no water-mass conservation law (there is no evaporation/condensation budget — that is the rung-2
  moist model). The honest leg is the **global-water budget *consistency***: the C–C-scaled global
  mean ``⟨P⟩`` moves **monotonically** with T̄ exactly as designed (warmer → wetter). It is a
  consistency check, not a conservation law — stated plainly, not dressed up as one.
* **Benchmark (loose).** Present-day Earth's major precip belts: a wet equator, the subtropical
  deserts, the midlatitude storm tracks. Cited ([[precip-parameterization-source]]); the exact band
  amplitudes/centres are calibration-dependent, asserted only in loose bands.

Non-circularity, named scope edge (plan §3)
-------------------------------------------
*Validated tight:* the band *structure* (ITCZ-wet / subtropics-dry / midlat-wet / poles-dry) and the
*monotone* global-water response. *Calibrated/flagged (loose):* the band amplitudes, centres, and the
7 %/K rate are cited, so the absolute precip values move only in loose bands.
*Scope edge, named:*
  - **No water cycle.** Prescribed, not derived — no moisture variable (rung-2).
  - **Band centres are fixed in v1; the storm-track centre is the rung-1 seam (now wired).** True band
    **migration** is a *circulation* response. The midlatitude storm-track centre is a parameter
    (``midlat_center_deg``) the **rung-1 circulation-informed** module :mod:`planet.circ_precip` sets
    from the **emergent jet latitude** — the band then tracks the dynamically-selected flow, not a
    constant (the array-seam the plan names). **Honest scope (see** :mod:`planet.circ_precip` **): a
    *trade* not an accuracy gain** (the dry single-layer jet sits ~8° equatorward of Earth's
    observation-calibrated 50°), and for *realistic* knobs the gradient/jet barely moves, so the
    present demo is the **migration mechanism** (decisive only under a synthetic off-centre gradient),
    not a large realistic shift. The **ITCZ/subtropical** centres stay fixed (the Hadley cell is out of
    the midlatitude channel's band); the uniform-``CC(T̄)`` **amplitude** is the rung-1
    *contrast-sharpening* seam, now built opt-in in :mod:`planet.moist` (the bullet below).
  - **Uniform amplitude here — "wet-get-wetter, dry-get-drier" is the opt-in rung-1 split.** ``CC(T̄)``
    in *this* module multiplies the *whole* pattern by one factor, so the dry subtropics get **wetter**
    under warming too — the opposite of observed. The **thermodynamic** pattern amplification (wet bands
    intensify faster than the mean while the deserts dry — Held & Soden 2006) is **BUILT** as the opt-in
    :func:`planet.moist.wet_get_wetter_precip_field`: the global *mean* scales at the slow
    energy-constrained rate, the *anomaly* (the wet−dry contrast) at the faster C–C rate. (The remaining
    deferral is the **dynamic** amplification — a *changing circulation* sharpening the pattern — which
    needs the circulation-set moisture convergence, rung 3+.) Named, distinct from the band-migration
    deferral above (that moves centres; this sharpens the wet/dry *contrast*).
  - **The 7 %/K is moisture-capacity, not the global precip rate.** Global-*mean* precipitation is
    **energy-constrained** to a slower ≈ 2–3 %/K (the atmosphere can only radiate away so much latent
    heating); the C–C 7 %/K is the *moisture-content* rate. v1 scales at the C–C rate as a moisture
    proxy and **does not** claim the energy-constrained global rate — that closure is the rung-2 moist
    energetics, named not modelled. So the "global-water budget" leg is a consistency check only.

Units — cm/yr (to match the Whittaker classifier's axes), latitude in degrees
-----------------------------------------------------------------------------
``P`` is **cm of precipitation per year** — the unit of Whittaker's diagram
(:mod:`planet.biomes`), chosen so the two modules share one unit with no conversion in the
path (the recurring units discipline). Latitude ``φ`` is in **degrees**; the EBM's area coordinate
``x = sin φ`` converts via ``φ = asin(x)`` (``ClimateState.latitude_deg``).
"""
from __future__ import annotations

import numpy as np

from .ebm import ClimateState

# --------------------------------------------------------------------------- #
# Pinned precipitation parameters ([[precip-parameterization-source]]).
# Calibrated to the observed annual zonal-mean precip-by-latitude structure (a wet ITCZ, dry
# subtropics, midlatitude storm tracks, dry poles) — cited and pinned at build, NOT from memory.
# An idealized teaching pattern: the band amplitudes are chosen to expose the desert/forest belts
# cleanly (the absolute values are loose / calibration-dependent, per the non-circularity split).
# --------------------------------------------------------------------------- #
P_BASELINE_CM = 20.0          # cm/yr — the dry polar/desert floor the bands ride on
ITCZ_AMP_CM = 215.0           # cm/yr — wet equatorial ITCZ peak (above baseline)
ITCZ_CENTER_DEG = 0.0         # °     — the annual-mean ITCZ sits ~on the equator
ITCZ_WIDTH_DEG = 12.0         # °     — Gaussian half-width of the equatorial rain belt
MIDLAT_AMP_CM = 75.0          # cm/yr — midlatitude storm-track peak (above baseline)
MIDLAT_CENTER_DEG = 50.0      # °     — the storm tracks sit near 50° (Ferrel-cell ascent)
MIDLAT_WIDTH_DEG = 15.0       # °     — Gaussian half-width of the storm-track belt

CC_RATE_PER_K = 0.07          # 1/K — Clausius–Clapeyron moisture-capacity rate (~7 %/K)
PRECIP_REF_TEMP_C = 15.0      # °C  — reference global-mean T where the pattern is calibrated (CC = 1);
                              #       ≈ present-day global-mean surface temperature (the EBM's 0-D anchor)


def precip_pattern(lat_deg: np.ndarray | float,
                   midlat_center_deg: float = MIDLAT_CENTER_DEG,
                   itcz_center_deg: float = ITCZ_CENTER_DEG) -> np.ndarray:
    """The circulation-set precipitation pattern ``pattern(φ)`` (cm/yr) at the reference climate.

    A sum of Gaussian latitude bands over a polar/desert baseline: the wet **ITCZ** peak (at
    ``itcz_center_deg``), plus the two **midlatitude** storm-track peaks at ``±midlat_center_deg``. The dry
    **subtropics** (~25–30°) are the *emergent trough* between the ITCZ and midlatitude bands; the
    **poles** decay to the baseline. The midlatitude bands are symmetric (``|φ|`` is used); the **ITCZ band
    uses the *signed* latitude**, so it can sit off the equator without breaking the midlatitude symmetry.
    This is *where* it rains; the thermodynamic amplitude is applied by :func:`clausius_clapeyron_factor`.

    Two band-centre seams, each defaulting to its cited constant so the rung-0 pattern is recovered
    **bit-for-bit by construction**:

    * ``midlat_center_deg`` — the **storm-track centre**, the **rung-1 circulation seam**. The rung-1
      module :mod:`planet.circ_precip` overrides it with the **emergent jet latitude** (a *trade*, not an
      accuracy gain — see that module).
    * ``itcz_center_deg`` — the **ITCZ centre**, the **rung-2.x energetic seam**. The full-sphere model
      :mod:`planet.sphere_ebm` overrides it with the **energy-flux equator** (:func:`planet.sphere_ebm.itcz_informed_precip`)
      so the rain belt migrates with the energetics. **Honest scope (see that module): a *dry* EBM
      relocates a *prescribed* band — not emergent rainfall.** The subtropical centres stay emergent
      troughs between the (now movable) bands.
    """
    lat = np.asarray(lat_deg, dtype=float)
    aphi = np.abs(lat)
    itcz = ITCZ_AMP_CM * np.exp(-(((lat - itcz_center_deg) / ITCZ_WIDTH_DEG) ** 2))      # signed (migratable)
    midlat = MIDLAT_AMP_CM * np.exp(-(((aphi - midlat_center_deg) / MIDLAT_WIDTH_DEG) ** 2))
    return P_BASELINE_CM + itcz + midlat


def clausius_clapeyron_factor(global_mean_T: float, ref_T: float = PRECIP_REF_TEMP_C,
                              rate: float = CC_RATE_PER_K) -> float:
    """The Clausius–Clapeyron moisture-amplitude factor ``CC(T̄) = exp(rate·(T̄ − T_ref))``.

    Scales the whole precipitation pattern by the **global-mean** temperature: a warmer atmosphere
    holds exponentially more moisture (≈ 7 %/K), so the hydrological cycle intensifies. ``CC = 1`` at
    the reference (present-day) global mean. This is the *moisture-capacity* rate — **not** the
    energy-constrained global-mean precip rate (~2–3 %/K, the named rung-2 scope edge); it is used as
    a moisture proxy, and the global-water leg checks only that ``⟨P⟩`` moves monotonically with T̄.
    """
    return float(np.exp(rate * (float(global_mean_T) - ref_T)))


def precipitation(lat_deg: np.ndarray | float, global_mean_T: float = PRECIP_REF_TEMP_C,
                  midlat_center_deg: float = MIDLAT_CENTER_DEG,
                  itcz_center_deg: float = ITCZ_CENTER_DEG) -> np.ndarray:
    """Diagnostic annual precipitation ``P(φ, T̄) = pattern(φ)·CC(T̄)`` (cm/yr).

    The circulation-set latitudinal pattern (:func:`precip_pattern`) times the global C–C moisture
    amplitude (:func:`clausius_clapeyron_factor`). Always non-negative (both factors are). This is the
    second input — beside the EBM temperature — the Whittaker classifier (:mod:`planet.biomes`)
    consumes. Pass the climate's ``global_mean_T`` to intensify the bands as the planet warms; the
    default reference leaves the pattern at its present-day calibration. ``midlat_center_deg`` (rung-1
    circulation) and ``itcz_center_deg`` (rung-2.x energetics) are the band-centre seams of
    :func:`precip_pattern`; their defaults recover the rung-0 field exactly.
    """
    return precip_pattern(lat_deg, midlat_center_deg, itcz_center_deg) * clausius_clapeyron_factor(global_mean_T)


def precip_field(state: ClimateState) -> np.ndarray:
    """Precipitation ``P(φ)`` (cm/yr) for an equilibrium :class:`~planet.ebm.ClimateState`.

    Convenience over :func:`precipitation`: reads the climate's latitudes (``state.latitude_deg()``)
    and its ``global_mean_T`` (the C–C amplitude knob), returning the precip field on the EBM's own
    grid — the array the biome map pairs with ``state.T``.
    """
    return precipitation(state.latitude_deg(), state.global_mean_T)
