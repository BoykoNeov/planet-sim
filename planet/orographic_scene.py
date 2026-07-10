"""Orographic scenes — placing the rain shadow on the sphere and into the biome map (Rung 5A.2, plan §12.5).

Rung 5A (:mod:`planet.orographic`) built the *engine*: the Smith & Barstad (2004) linear transfer
function that turns a 2-D terrain + a uniform cross-mountain wind into a windward-rain / lee-shadow
precipitation field, on a **regional Cartesian patch** with **exact analytic anchors**. What it left
**named and deferred** — the 5A.2 integration questions — was everything *around* that patch:

1. **Sphere placement.** The engine wants a patch with a uniform metre-spaced grid; the planet is a
   lat×lon globe in degrees. :func:`patch_spacings` supplies the tangent-plane metric
   ``dx = R·cos φ·Δλ``, ``dy = R·Δφ`` at the patch's reference latitude — exactly the tangent-plane
   the engine already documents as its honest scope (:mod:`planet.orographic`), now made explicit.
2. **Where the cross-mountain wind comes from.** On a globe whose emergent circulation
   (:mod:`planet.coupler`) is a **purely zonal** mid-latitude jet, the cross-mountain wind is the
   westerly sampled at the patch's latitude (:func:`wind_from_jet`). This is the honesty caveat carried
   forward from Rung 5A: the flow is **prescribed** (read off the zonal jet), **not emergent** — the
   linear theory does not solve the flow *over* the mountain, it is handed one.
3. **Units + combination with the biome map.** The engine's native unit is an *instantaneous* mm/hr
   condensation rate; the Phase-2 biome map (:mod:`planet.demo_biomes`) speaks *annual* cm/yr. The
   conversion is **not** a naive annualisation (``mm/hr × 8766 hr`` is an order of magnitude too big and
   swamps the biome classifier); it goes through :data:`OROGRAPHIC_HOURS_PER_YEAR`, an *effective annual
   duration of active orographic uplift* — a named, loose-magnitude calibration knob, in exactly the
   spirit of :mod:`planet.precip`'s loose band amplitudes.
4. **Serialization.** The regional scene is a :class:`~planet.planetmap.Grid` + a stack of
   :class:`~planet.planetmap.Layer` s, so it round-trips through the *grid-agnostic*
   :mod:`planet.planet_spec` schema for free — the interchange seam needs no new machinery.

The combination is *enhancement-only* — the honest scope (advisor-caught)
--------------------------------------------------------------------------
Smith & Barstad is a condensation-**from-forced-ascent** model: it produces windward **enhancement**
and a lee where the orographic bonus **decays to zero** (see the exact solution
:func:`planet.orographic.triangle_ridge_exact` — the lee term falls to 0 at ``x_c`` and is zero
beyond; it is *not* a strong negative that eats the ambient rain). So the combination here is

    ``P_total = P_zonal_baseline + P_orographic_bonus``          (enhancement-only)

and the rain shadow reads as **windward ≫ lee**, with the *immediate* lee left at (near) the zonal-mean
baseline — a dry shadow, never below baseline. (The *full* transfer function additionally produces a weak
**downstream secondary rain band** further out — driven by the **propagating-mode phase** (the
``1 − i m H_w`` factor), *real to the linear model*: verified by discriminator to **vanish at**
``H_w = 0`` and to **hold its physical downwind distance** when the domain is doubled, so it is neither
FFT wrap-around nor a pad artifact — see :mod:`planet.tests.test_orographic_scene`. It is *not* a trapped
lee wave (that needs a stratification profile this model does not carry). It is non-negative, so the
enhancement-only invariant ``P_total ≥ baseline`` still holds.) The biome payoff is therefore real and
directional on the
**windward** side (the extra rain shifts it toward forest / rain forest); a cell reads as a lee *desert*
only where the zonal baseline was already marginal.

*Background depletion* — the windward rainout drying the air so the lee baseline itself drops (the
mechanism behind the real Columbia-Basin desert *behind* the Cascades) — is an added moisture-budget
assumption *beyond* Smith & Barstad. It is built as the **opt-in Rung 5A.3** refinement
(:mod:`planet.orographic_depletion`, ``build_scene(..., deplete=True)``): a 1-D along-wind flux budget
that drains the advected column so the lee total drops *below* baseline. The default here stays
**enhancement-only** (``deplete=False``) — the honest 5A.2 combination — so the invariant above holds
unless depletion is explicitly turned on.

The demo range must sit under the westerlies (advisor-caught)
------------------------------------------------------------
The wind comes from the annual-mean **zonal** jet, which is zero outside the mid-latitude channel
(:func:`planet.coupler.couple_jet`). A mountain outside that band gets **zero cross-mountain wind →
zero orographic anomaly**. So a meaningful demo range must sit under the westerlies — the **Cascades**
(~47°N), the **southern Andes** (~40°S), the **Southern Alps** of New Zealand (~43°S). The **Western
Ghats** — a classic orographic wall — are driven by the *monsoon* (a seasonal SW flow) and cannot be
driven by the annual-mean zonal jet, so they are *not* a valid demo here (:data:`DEMO_RANGES`).

Resolution — a fine regional patch, not the 73-longitude globe (advisor-caught)
------------------------------------------------------------------------------
The biome globe is ``N_LON = 73`` (~5° ≈ ~400 km per cell at mid-latitude); a rain shadow is a
~20–100 km east–west feature that lives *inside one global cell*. Regridding the patch onto the coarse
globe would smear the shadow to nothing. So the scene is computed and displayed on its **own fine
regional grid** (its Δλ, Δφ in the tens of km), and the demo is a **regional inset**, not the native
globe. That resolution reality is a property of the deliverable, named here rather than hidden.

Validation triad (plan §3) — what is tight vs loose
---------------------------------------------------
* **Tight (exact):** the mm/hr ↔ cm/yr conversion round-trips exactly (:func:`mm_hr_to_cm_yr` /
  :func:`cm_yr_to_mm_hr`); the tangent-plane spacings are the closed-form metric.
* **Tight (structural):** the **placement/orientation** anchor — a *localized* hill
  (:func:`gaussian_hill`, compact in both axes) under a westerly casts its shadow **in longitude** (a
  windward/lee asymmetry with the peak displaced *upwind* in lon) and is ~**symmetric in latitude** (no
  N–S asymmetry, no lat displacement). That pins the patch mapping ``lon → x`` (the along-wind axis),
  ``lat → y``: a transposed placement would put the shadow in latitude. (Note the honest engine caveat:
  an *idealized infinite zonal ridge* would cast no shadow, but on a finite zero-padded patch a
  lon-uniform ridge is localised by the pad into a lon block that *does* respond — so the clean
  pad-safe null is a **compact hill's lat-symmetry**, not a zonal ridge's amplitude.) A mountain
  **outside the jet band** produces a strictly-zero anomaly; the enhancement-only combination leaves
  the **lee at the baseline** and lifts the **windward** above it.
* **Loose (magnitude):** the absolute cm/yr amplitude, set by :data:`OROGRAPHIC_HOURS_PER_YEAR` and the
  cited Smith & Barstad constants ([[smith-barstad-orographic-source]]) — a calibration band, as in
  :mod:`planet.precip`.

See [[planet-rung5a-orographic]]; the engine and its constants are :mod:`planet.orographic`.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from planet import biomes, orographic as og, orographic_depletion as ogd
from planet.planetmap import Grid, Layer, LayerKind, _biome_style

R_EARTH_M = 6.371e6              # m   — Earth's mean radius (the tangent-plane metric; matches circulation.R_EARTH)

# The one loose-magnitude calibration knob (§3 / advisor-caught). Smith & Barstad returns an
# *instantaneous* condensation rate in mm/hr; the biome map wants an *annual* accumulation in cm/yr.
# A naive annualisation (mm/hr × 8766 hr/yr) overstates the orographic bonus by ~an order of magnitude
# (a 4 mm/hr crest → ~3500 cm/yr, which swamps the Whittaker classifier). Instead we treat the S&B rate
# as active for an *effective annual duration* of forced orographic uplift — the fraction of the year
# the range actually sees saturated cross-flow at ~this rate. ~500 h ≈ 5.7 % of the year lands a
# ~4 mm/hr crest at a few-hundred cm/yr, a sane orographic surplus (cf. the Cascades' ~250–500 cm/yr).
# This is a calibration band, not a pinned constant — named, not hidden.
OROGRAPHIC_HOURS_PER_YEAR = 500.0    # h/yr — effective annual duration of active orographic uplift
_MM_PER_CM = 10.0                    # mm per cm (1 cm/yr = 10 mm/yr)


# Demo mountain ranges that sit under the annual-mean zonal westerlies (see the module docstring).
# (lat_center_deg, lon_center_deg) — the SH ranges use the hemispheric symmetry of the zonal climate.
DEMO_RANGES: dict[str, tuple[float, float]] = {
    "cascades": (47.0, -121.5),      # Pacific NW, USA — westerlies off the Pacific
    "southern_andes": (-40.0, -71.0),  # Patagonia — the roaring-forties westerlies
    "southern_alps_nz": (-43.5, 170.5),  # New Zealand — the Southern Alps in the westerly belt
}


# --------------------------------------------------------------------------- #
# Unit conversion — mm/hr (Smith & Barstad native) ↔ cm/yr (biome-map native).  Exact round-trip.
# --------------------------------------------------------------------------- #
def mm_hr_to_cm_yr(p_mm_hr, hours_per_year: float = OROGRAPHIC_HOURS_PER_YEAR):
    """Convert an instantaneous orographic rate (mm/hr) to an annual accumulation (cm/yr).

    ``cm/yr = mm/hr · hours_per_year / 10`` — the rate applied over the *effective* annual duration of
    active orographic uplift (:data:`OROGRAPHIC_HOURS_PER_YEAR`), then mm → cm. Not a naive annualisation
    over all 8766 hours (see the module docstring). Exactly inverted by :func:`cm_yr_to_mm_hr`.
    """
    return np.asarray(p_mm_hr, dtype=float) * hours_per_year / _MM_PER_CM


def cm_yr_to_mm_hr(p_cm_yr, hours_per_year: float = OROGRAPHIC_HOURS_PER_YEAR):
    """Inverse of :func:`mm_hr_to_cm_yr`: ``mm/hr = cm/yr · 10 / hours_per_year`` (exact round-trip)."""
    return np.asarray(p_cm_yr, dtype=float) * _MM_PER_CM / hours_per_year


# --------------------------------------------------------------------------- #
# Sphere placement — the tangent-plane metric that maps a lat×lon window to the engine's metre grid.
# --------------------------------------------------------------------------- #
def patch_spacings(lat_ref_deg: float, dlat_deg: float, dlon_deg: float) -> tuple[float, float]:
    """The metre grid spacings ``(dy, dx)`` for a lat×lon patch at reference latitude ``lat_ref_deg``.

    The linear theory lives on a tangent plane with a uniform metre grid. On the sphere, a step of
    ``Δφ`` degrees in latitude is ``dy = R·Δφ`` (in radians) northward, and a step of ``Δλ`` degrees in
    longitude is ``dx = R·cos φ·Δλ`` eastward — the cosine shrinking the zonal spacing toward the poles.
    Evaluating the cosine at a single ``lat_ref_deg`` is the tangent-plane approximation the engine
    already declares as its honest regional scope (:mod:`planet.orographic`); it is exact for a patch
    small compared with the planet's radius.
    """
    dy = R_EARTH_M * np.radians(dlat_deg)
    dx = R_EARTH_M * np.cos(np.radians(lat_ref_deg)) * np.radians(dlon_deg)
    return float(dy), float(dx)


# --------------------------------------------------------------------------- #
# The cross-mountain wind — read off the zonal jet (prescribed, not emergent; the Rung-5A caveat).
# --------------------------------------------------------------------------- #
def wind_from_jet(jet, lat_deg: float) -> tuple[float, float]:
    """The cross-mountain ``(speed, direction_deg)`` sampled from a coupled zonal jet at ``lat_deg``.

    The coupler's emergent circulation (:func:`planet.coupler.couple_jet`) is a **purely zonal**
    mid-latitude jet ``u(φ)`` living in a channel band (mirrored to both hemispheres by symmetry). The
    cross-mountain wind at a patch is that zonal wind sampled at the patch latitude: interpolated on
    ``|lat|`` (so a southern-hemisphere range reads the same symmetric jet), **zero outside the band**
    (a mountain off the westerlies gets no forcing — the ``left=right=0`` of
    :func:`planet.planetmap.circulation_layer`). A positive (eastward) wind is a **westerly**
    (meteorological ``direction = 270``); the sign flips to an easterly (90°) if the sampled wind is
    negative. This is **prescribed** flow — the honesty caveat carried from Rung 5A: the linear theory
    is *handed* the wind, it does not derive the flow over the mountain.
    """
    u = float(np.interp(abs(float(lat_deg)), jet.phi, jet.u_profile, left=0.0, right=0.0))
    direction_deg = og.DIRECTION_WESTERLY_DEG if u >= 0.0 else 90.0   # 270 = westerly (+x), 90 = easterly
    return abs(u), direction_deg


# --------------------------------------------------------------------------- #
# Regional grid + idealized terrain (the demo geometries; a real heightmap can be dropped in instead).
# --------------------------------------------------------------------------- #
def regional_grid(lat_center: float, lon_center: float, lat_span_deg: float, lon_span_deg: float,
                  n_lat: int, n_lon: int) -> tuple[np.ndarray, np.ndarray]:
    """A fine regional lat/lon axis pair centred on ``(lat_center, lon_center)`` — the patch grid (degrees).

    Spans ``lat_span_deg`` × ``lon_span_deg`` with ``n_lat`` × ``n_lon`` samples. The spans are small
    (a few degrees) and the sampling fine (Δ in the tens of km) so the rain shadow is resolved — the
    opposite of the coarse ``N_LON = 73`` globe, on which it would smear to nothing (module docstring).
    """
    lat = np.linspace(lat_center - lat_span_deg / 2, lat_center + lat_span_deg / 2, n_lat)
    lon = np.linspace(lon_center - lon_span_deg / 2, lon_center + lon_span_deg / 2, n_lon)
    return lat, lon


def meridional_ridge(lat_deg: np.ndarray, lon_deg: np.ndarray, lon_center: float | None = None,
                     amplitude_m: float = 2000.0, half_width_deg: float = 0.7) -> np.ndarray:
    """A **north–south** mountain ridge — a Gaussian in longitude, uniform in latitude → ``[lat, lon]`` (m).

    This is the orientation that **casts a shadow under a westerly**: the terrain varies along the wind
    (in longitude/x), so the forced ascent and its lee drying line up with the flow. ``half_width_deg``
    is the Gaussian σ in longitude; ``lon_center`` defaults to the grid centre. Pairs with
    :func:`zonal_ridge` (the null-orientation control).
    """
    lat_deg = np.asarray(lat_deg, dtype=float)
    lon_deg = np.asarray(lon_deg, dtype=float)
    lon0 = float(np.mean(lon_deg)) if lon_center is None else float(lon_center)
    profile = amplitude_m * np.exp(-((lon_deg - lon0) ** 2) / (2.0 * half_width_deg ** 2))   # in lon
    return np.repeat(profile[None, :], lat_deg.size, axis=0)                                  # [lat, lon]


def gaussian_hill(lat_deg: np.ndarray, lon_deg: np.ndarray, lat_center: float | None = None,
                  lon_center: float | None = None, amplitude_m: float = 2500.0,
                  half_width_deg: float = 0.5) -> np.ndarray:
    """A **localized** Gaussian hill — compact in both latitude and longitude → ``[lat, lon]`` (m).

    The **pad-safe placement anchor**: because it decays to ~0 well inside the patch in *both* axes, the
    engine's internal zero-padding adds no spurious structure. Under a westerly its shadow appears in
    **longitude** (peak displaced upwind, windward wetter than lee) and is ~symmetric in **latitude** —
    the clean discriminator that ``lon`` is the along-wind (x) axis and ``lat`` the y axis
    (:mod:`planet.tests.test_orographic_scene`). Preferred over a "zonal ridge null" because a
    lon-uniform ridge on a finite padded patch is localised by the pad into a responding block (see the
    module docstring's orientation note).
    """
    lat_deg = np.asarray(lat_deg, dtype=float)
    lon_deg = np.asarray(lon_deg, dtype=float)
    lat0 = float(np.mean(lat_deg)) if lat_center is None else float(lat_center)
    lon0 = float(np.mean(lon_deg)) if lon_center is None else float(lon_center)
    la, lo = np.meshgrid(lat_deg, lon_deg, indexing="ij")                                     # [lat, lon]
    return amplitude_m * np.exp(-(((la - lat0) ** 2 + (lo - lon0) ** 2)) / (2.0 * half_width_deg ** 2))


# --------------------------------------------------------------------------- #
# The regional scene — the assembled 5A.2 artifact.
# --------------------------------------------------------------------------- #
@dataclass(frozen=True, eq=False)
class OrographicScene:
    """A regional orographic-precipitation scene — the assembled Rung-5A.2 artifact (``eq=False``: arrays).

    A fine regional lat×lon patch carrying the zonal-mean climate sampled onto it, the orographic
    bonus computed on it, their enhancement-only sum, and the biomes that result — everything the demo
    figure and the serialization seam consume. All 2-D fields are indexed ``[lat, lon]``.

    Fields
    ------
    lat_deg, lon_deg : 1-D degree axes of the patch grid.
    elevation_m : the terrain (m) — now **live** (it drives ``orographic_precip_cm``), not the inert
        seam it is on the coarse globe.
    temperature_C : the zonal-mean EBM temperature sampled to the patch latitudes (°C), broadcast in lon
        — the climate underneath stays zonal-mean (the honest scope; only precipitation goes 2-D).
    baseline_precip_cm : the Phase-2 zonal-mean precip (cm/yr) sampled to the patch, broadcast in lon.
    orographic_precip_cm : the Smith & Barstad windward-enhancement bonus (cm/yr) on the patch.
    depletion_factor : ``g ∈ (0, 1]`` — the along-wind lee moisture depletion (Rung 5A.3,
        :mod:`planet.orographic_depletion`); ``g ≡ 1`` when depletion is off (the 5A.2 default). It
        multiplies the baseline, so the lee baseline drops **below** the zonal mean where ``g < 1``.
    precip_cm : ``depletion_factor·baseline + orographic`` — the total (cm/yr); enhancement-only reduces
        to ``baseline + orographic`` when ``g ≡ 1``.
    biome_codes / baseline_biome_codes : Whittaker biomes for the total vs the baseline precip — their
        difference is the *payoff* (where the mountain changes the map).
    wind_speed / wind_direction_deg : the prescribed cross-mountain wind used (m/s; meteorological deg).
    lat_ref_deg : the reference latitude for the tangent-plane metric and the Coriolis parameter.
    hours_per_year : the :data:`OROGRAPHIC_HOURS_PER_YEAR` calibration used for mm/hr → cm/yr.
    """

    lat_deg: np.ndarray
    lon_deg: np.ndarray
    elevation_m: np.ndarray
    temperature_C: np.ndarray
    baseline_precip_cm: np.ndarray
    orographic_precip_cm: np.ndarray
    depletion_factor: np.ndarray
    precip_cm: np.ndarray
    biome_codes: np.ndarray
    baseline_biome_codes: np.ndarray
    wind_speed: float
    wind_direction_deg: float
    lat_ref_deg: float
    hours_per_year: float

    @property
    def biome_changed_fraction(self) -> float:
        """Fraction of patch cells whose biome differs from the zonal-mean baseline — the payoff metric."""
        return float(np.mean(self.biome_codes != self.baseline_biome_codes))

    @property
    def lee_desert_fraction(self) -> float:
        """Fraction of cells whose total precip fell **below** the zonal-mean baseline — the 5A.3 payoff.

        Zero under the enhancement-only 5A.2 default (``g ≡ 1``); positive only when lee-side moisture
        depletion (:mod:`planet.orographic_depletion`) drains the baseline below the zonal mean — the
        real rain-shadow *desert* enhancement-only structurally cannot make.
        """
        return float(np.mean(self.precip_cm < self.baseline_precip_cm - 1e-9))


def build_scene(result, lat_deg: np.ndarray, lon_deg: np.ndarray, elevation_m: np.ndarray, *,
                jet=None, speed: float = og.U_REF_M_S, direction_deg: float = og.DIRECTION_WESTERLY_DEG,
                hours_per_year: float = OROGRAPHIC_HOURS_PER_YEAR,
                deplete: bool = False, pwv_in_mm: float = ogd.PWV_IN_MM,
                lat_ref_deg: float | None = None, **orographic_kwargs) -> OrographicScene:
    """Assemble an :class:`OrographicScene` — place the patch, source the wind, rain the shadow, re-map biomes.

    Parameters
    ----------
    result : BiomeResult
        The Phase-2 zonal-mean climate (:func:`planet.demo_biomes.compute`) — supplies the temperature
        and baseline precipitation that the orographic bonus is added *on top of*. Its climate is
        hemispherically symmetric, so the patch reads ``T``/``P`` by interpolating on ``|lat|`` against
        the hemisphere EBM grid (a southern-hemisphere range gets the mirror of the northern climate).
    lat_deg, lon_deg : 1-D degree axes of the (fine) patch grid — see :func:`regional_grid`.
    elevation_m : ``(n_lat, n_lon)`` terrain (m) on that grid — see :func:`meridional_ridge`.
    jet : CoupledJet | None
        If given, the cross-mountain wind is **sourced from the zonal jet** at the patch latitude
        (:func:`wind_from_jet`), overriding ``speed``/``direction_deg`` — the emergent-jet wiring. The
        flow is prescribed, not emergent (the Rung-5A honesty caveat).
    speed, direction_deg : the prescribed wind if no ``jet`` is given (default: the S&B reference
        westerly).
    hours_per_year : the mm/hr → cm/yr calibration (:data:`OROGRAPHIC_HOURS_PER_YEAR`).
    deplete : opt-in Rung 5A.3 lee moisture depletion (:mod:`planet.orographic_depletion`). Default
        ``False`` keeps the enhancement-only 5A.2 combination (``g ≡ 1``, lee left at the baseline). With
        ``True`` the windward rainout drains the advected column, so the lee baseline drops **below** the
        zonal mean — the real rain-shadow desert.
    pwv_in_mm : the incoming column precipitable water ``W₀`` for the depletion budget
        (:data:`planet.orographic_depletion.PWV_IN_MM`); only used when ``deplete=True``.
    lat_ref_deg : reference latitude for the tangent-plane metric + Coriolis ``f``; defaults to the
        patch's centre latitude.
    **orographic_kwargs : forwarded to :func:`planet.orographic.orographic_precip` (e.g. ``Cw``, ``Nm``).

    Returns the fully-populated :class:`OrographicScene`. With ``deplete=False`` this is the
    *enhancement-only* combination (windward enhancement, lee at the zonal baseline); with ``deplete=True``
    the lee is additionally drawn below baseline (see the module docstring's honest scope).
    """
    lat_deg = np.asarray(lat_deg, dtype=float)
    lon_deg = np.asarray(lon_deg, dtype=float)
    elevation_m = np.asarray(elevation_m, dtype=float)
    if elevation_m.shape != (lat_deg.size, lon_deg.size):
        raise ValueError(f"elevation must be {(lat_deg.size, lon_deg.size)} (lat×lon), got {elevation_m.shape}")

    ref = float(np.mean(lat_deg)) if lat_ref_deg is None else float(lat_ref_deg)
    if jet is not None:
        speed, direction_deg = wind_from_jet(jet, ref)

    # Sphere placement: the tangent-plane metre spacings (assumes a uniform degree grid).
    dlat = float(lat_deg[1] - lat_deg[0]) if lat_deg.size > 1 else 1.0
    dlon = float(lon_deg[1] - lon_deg[0]) if lon_deg.size > 1 else 1.0
    dy, dx = patch_spacings(ref, dlat, dlon)

    # The engine (mm/hr) → cm/yr. background_mm_hr = 0: this layer is the *anomaly* added to the baseline.
    p_mm_hr = og.orographic_precip(elevation_m, dx, dy, speed=speed, direction_deg=direction_deg,
                                   latitude_deg=ref, background_mm_hr=0.0, **orographic_kwargs)
    orographic_cm = mm_hr_to_cm_yr(p_mm_hr, hours_per_year)

    # Lee moisture depletion (Rung 5A.3, opt-in): the along-wind budget runs on the *instantaneous* mm/hr
    # rate (NOT the annualised cm/yr — the unit split), giving a dimensionless g that then multiplies the
    # annual baseline. g ≡ 1 (default) reduces to the 5A.2 enhancement-only combination exactly.
    if deplete:
        u0, _ = og.wind_components(speed, direction_deg)    # sign of u sets the downwind integration
        g = ogd.depletion_factor(p_mm_hr, dx, u0, pwv_in_mm=pwv_in_mm)
    else:
        g = np.ones_like(orographic_cm)

    # The zonal-mean climate sampled onto the patch (symmetric in |lat|), broadcast across longitude.
    state = result.state
    hemi_lat = state.latitude_deg()
    T_1d = np.interp(np.abs(lat_deg), hemi_lat, state.T)
    base_1d = np.interp(np.abs(lat_deg), hemi_lat, result.precip_cm)
    T_2d = np.repeat(T_1d[:, None], lon_deg.size, axis=1)
    baseline_cm = np.repeat(base_1d[:, None], lon_deg.size, axis=1)

    precip_cm = g * baseline_cm + orographic_cm             # depleted baseline + windward enhancement
    biome_codes = biomes.classify_field(T_2d, precip_cm)
    baseline_biome = biomes.classify_field(T_2d, baseline_cm)

    return OrographicScene(
        lat_deg=lat_deg, lon_deg=lon_deg, elevation_m=elevation_m, temperature_C=T_2d,
        baseline_precip_cm=baseline_cm, orographic_precip_cm=orographic_cm, depletion_factor=g,
        precip_cm=precip_cm, biome_codes=biome_codes, baseline_biome_codes=baseline_biome,
        wind_speed=float(speed), wind_direction_deg=float(direction_deg),
        lat_ref_deg=ref, hours_per_year=float(hours_per_year),
    )


# --------------------------------------------------------------------------- #
# Serialization seam — the scene is a Grid + Layers, so it rides the grid-agnostic planet_spec schema.
# --------------------------------------------------------------------------- #
def scene_to_view(scene: OrographicScene):
    """The regional scene as a :class:`~planet.planetmap.PlanetView` — the render + serialization seam.

    Registers the patch's layers on a regional :class:`~planet.planetmap.Grid`: **temperature**,
    the zonal-mean **precipitation** baseline, the **orographic precipitation** bonus, the enhancement
    **total precipitation**, the now-**live elevation** (``inert=False`` — on the patch it finally
    *drives* the rain, unlike the inert globe seam, §9.3), and the **biome** map. Because
    :mod:`planet.planet_spec` serializes any ``Grid`` + ``Layer`` stack, this round-trips through the v1
    interchange schema with no new machinery (build the seam, reuse it) — the 5A.2 serialization piece.

    **Re-runnability caveat (honest scope):** the round-trip is *array identity* of the rendered fields,
    not full re-computability. The spec's ``knobs`` are the EBM climate params only; the orographic
    *provenance* — the cross-mountain wind, ``hours_per_year``, and elevation-as-a-driver — is **not**
    carried, so a loaded scene can be *displayed and compared* but not *re-derived* from the spec alone
    (mirroring the §9.3 rule that array-identity is not a re-runnable climate).
    """
    from planet.planetmap import PlanetView   # local import: keep the module import-light

    grid = Grid(lat=scene.lat_deg, lon=scene.lon_deg)
    layers = (
        Layer("temperature", LayerKind.SCALAR_FIELD, scene.temperature_C, "°C",
              style={"colorscale": "RdBu_r"}, z_order=0),
        Layer("precipitation", LayerKind.SCALAR_FIELD, scene.baseline_precip_cm, "cm/yr",
              style={"colorscale": "GnBu"}, z_order=0),
        Layer("orographic precipitation", LayerKind.SCALAR_FIELD, scene.orographic_precip_cm, "cm/yr",
              style={"colorscale": "GnBu"}, z_order=0),
        Layer("total precipitation", LayerKind.SCALAR_FIELD, scene.precip_cm, "cm/yr",
              style={"colorscale": "GnBu"}, z_order=0),
        Layer("elevation", LayerKind.SCALAR_FIELD, scene.elevation_m, "m",
              style={"colorscale": "Earth"}, z_order=-1, inert=False),   # live here: it drives the rain
        Layer("biome", LayerKind.SCALAR_FIELD, scene.biome_codes.astype(int),
              "Whittaker biome code", style=_biome_style(), z_order=1),
    )
    return PlanetView(grid=grid, layers=layers)


# --------------------------------------------------------------------------- #
# Demo — a rain shadow behind a range that sits under the westerly jet.
# --------------------------------------------------------------------------- #
def demo_scene(range_name: str = "cascades", *, n_lat: int = 41, n_lon: int = 161,
               lat_span_deg: float = 6.0, lon_span_deg: float = 6.0,
               amplitude_m: float = 2500.0, use_jet: bool = True,
               deplete: bool = False) -> OrographicScene:
    """Build a demo :class:`OrographicScene`: a meridional ridge under the westerlies for a named range.

    Solves the present-day zonal-mean climate, (optionally) the emergent coupled jet, places a fine
    patch on the named :data:`DEMO_RANGES` entry with a north–south ridge, and returns the scene. With
    ``use_jet`` the cross-mountain wind is read off the emergent jet at that latitude; otherwise the S&B
    reference westerly is used (a fallback that does not need the shallow-water spin-up). With
    ``deplete`` the Rung-5A.3 lee moisture budget is turned on, drawing the lee below baseline (the real
    rain-shadow desert); the default is the enhancement-only 5A.2 combination.
    """
    from planet import demo_biomes

    if range_name not in DEMO_RANGES:
        raise KeyError(f"unknown range {range_name!r}; have {sorted(DEMO_RANGES)}")
    lat_c, lon_c = DEMO_RANGES[range_name]

    result = demo_biomes.compute()
    lat, lon = regional_grid(lat_c, lon_c, lat_span_deg, lon_span_deg, n_lat, n_lon)
    elevation = meridional_ridge(lat, lon, lon_center=lon_c, amplitude_m=amplitude_m)

    jet = None
    if use_jet:
        from planet import coupler
        jet = coupler.couple_jet(result.state)

    return build_scene(result, lat, lon, elevation, jet=jet, lat_ref_deg=lat_c, deplete=deplete)
