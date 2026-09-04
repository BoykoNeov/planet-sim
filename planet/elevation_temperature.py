"""Elevation → temperature: the mountain is *cold*, not only wet (Rung 5A.4, plan §12.5).

Rung 5A woke the dormant elevation seam on the **precipitation** side: terrain forces air upward, the
air rains (:mod:`planet.orographic`), and the lee dries (:mod:`planet.orographic_depletion`). But the
**temperature** underneath stayed the zonal-mean climate broadcast across longitude — a 2500 m crest and
the valley beside it were handed the *same* temperature. That is the last piece of the §12.5 "cheap
tier" line (*elevation → a lapse-rate map diagnostic*), and it is what this module supplies: a
**diagnostic** surface-temperature correction

    ``T_surface(φ, λ) = T_zonal(φ) − ∫₀^z Γ dz'``

so that the terrain finally cools its own air. It is the second half of what a mountain does to a map,
and on the biome side it is the bigger half: the crest of the demo ridge drops ~16 K and re-classifies
from temperate forest to **alpine tundra** — the vertical biome zonation of Whittaker's own Great Smoky
Mountains transect, now emerging on a planet whose climate was computed, not drawn.

Two lapse rates — the pinned constant (default) and the emergent moist adiabat (opt-in)
--------------------------------------------------------------------------------------
* **Constant ``Γ`` (default).** :data:`planet.radiation.LAPSE_RATE` = 6.5 K/km — the convective-adjustment
  lapse rate rung 4 already pins and uses, so this path *reuses* an in-repo constant rather than inventing
  one. The correction is then the closed form ``T − Γ·z``.
* **The moist adiabat (``moist=True``).** ``Γ_m(T, p)`` is rung 4's own
  :func:`planet.radiation.moist_adiabatic_lapse_rate` — steep (→ ``g/c_p`` ≈ 9.8 K/km) in cold dry air,
  flattened toward ~3.5 K/km where latent-heat release is large. Integrating it up the column with a
  hydrostatic pressure makes the cooling **emergent and temperature-dependent** instead of prescribed:
  the *same* mountain cools its summit by different amounts in different climates.

The emergent rate was the bet, and the bet came back a NEGATIVE — the constant stays the default
--------------------------------------------------------------------------------------------------
The repository's habit is that making a prescribed number *emergent* retires a wall (rung 4 retired the
fixed ``B``). The obvious version of that here is to retire the 6.5 K/km constant in favour of the moist
adiabat. **It does not survive its own benchmark**, and that negative is the substance of this rung:

* **At mid-latitudes it merely confirms the constant.** At the demo range (47°, surface ≈ 6.6 °C) the
  moist adiabat integrated over 2500 m gives an **effective 6.31 K/km** — within ~3 % of the 6.5 it was
  meant to replace. The constant is revealed as a *mid-latitude* calibration, which is exactly what the
  standard atmosphere is. No wall retired; a pin confirmed from an independent direction.
* **Away from mid-latitudes it is worse, not better.** The one observational check available is the
  **freezing level** — the altitude of the 0 °C isotherm, which is measured. The deep-tropical level is the
  highest on the planet at **≈ 5 km** (Harris, Bowman & Shin 2000, a 20-year NCEP + TRMM climatology;
  :data:`OBSERVED_TROPICAL_FREEZING_LEVEL_M`). The 6.5 K/km constant puts it at **4.38 km — just *below*
  the band**, close but not inside it; the moist adiabat, whose effective rate there is only ≈ 3.7 K/km,
  puts it at **7.09 km — ~45 % above**. The verdict rests on that **ordering** (one is close, the other is
  far), *not* on the constant landing in the band — said this way deliberately, because the band is narrow
  enough for the difference to matter.
* **Why, physically.** ``Γ_m`` is the lapse rate of a *saturated ascending parcel*. The **environmental**
  profile of the tropical atmosphere is not that: the mean column is unsaturated, and it is the average of
  narrow saturated updraughts with broad dry subsidence. So the parcel adiabat systematically
  *under-cools* the real mean column where the air is warm and moist.
* **And at the other end it fails in the opposite direction.** ``Γ_m`` steepens toward ``g/c_p`` as the air
  cools, so this model says the same 2500 m mountain cools its summit **~9 K in the tropics and ~22 K at
  85°** (a factor ~2.4). The observed high-latitude lower troposphere is strongly **stably stratified and
  frequently inverted** — low-level inversions appear in over 95 % of Eurasian-Arctic winter soundings and
  are typically surface-based (Serreze, Kahl & Schnell 1992) — so a real polar mountain does *not* see
  8.8 K/km; if anything the near-surface air gets *warmer* with height. A saturated adiabat is the wrong
  idealisation where the air is neither saturated nor convecting. The model's latitude contrast is therefore
  **anti-correlated with the real one at both ends**, and is reported as what the *idealisation* says, not
  as a property of the planet.

So: ``moist=False`` (the pinned constant) stays the **default** everywhere, the emergent path ships
**opt-in as a diagnostic of what the parcel idealisation says**, and its failure is named rather than
hidden. The crest biome is the same under either rate, so the rung's **payoff is robust to the choice** —
what the choice changes is a claim this module declines to make. Same treatment rung 4 gave its
global-mean lapse-rate feedback (a single tropical column cannot carry the global mean) and rung 5B.4
gave obliquity (the axis was rejected after being tried, not after being assumed).

Validation triad (plan §3) — what is tight vs loose
---------------------------------------------------
* **Tight (exact).** The constant-``Γ`` path is the closed form ``T − Γz`` to machine precision, and the
  *integrator* reproduces it exactly when handed a constant ``lapse_rate_fn`` — a **reduction of the
  emergent path to the prescribed one** that makes the two structurally the same code. Flat terrain
  (``z ≡ 0``) returns the input temperature unchanged; the whole rung is **default-off**
  (:func:`planet.orographic_scene.build_scene` keeps ``lapse=False``), so the 5A.2/5A.3 scene is
  bit-for-bit what it was.
* **Tight (convergent).** The moist column is a **Heun (predictor–corrector) march in height** →
  second-order in ``dz``: halving the step quarters the error (the §4.3 "convergence check when there is
  a step" leg). It calls rung 4's ``Γ_m`` *directly*, so it cannot drift from the radiation column, and
  its small-``z`` limit is ``T₀ − Γ_m(T₀, p₀)·z`` with an O(z²) residual.
* **Conservation → a *consistency* leg, honestly weaker (named).** Cooling the surface without re-solving
  the EBM breaks the TOA budget over the patch, so there is **no** energy law to close here — the same
  honesty :mod:`planet.biomes` states for the classifier. The substitute is an exact **identity**: for the
  constant rate, the area-weighted patch-mean cooling equals ``Γ × ⟨z⟩`` exactly, and the biome partition
  still tiles the patch (area fractions sum to 1).
* **Loose (calibrated / benchmark) — and it is the leg that decided the default.** The **freezing-level
  altitude** ``z(T = 0 °C)`` against the observed tropical ~4.5–5 km: the constant rate lands just below
  the band, the moist adiabat overshoots it by ~45 %. The **ordering** is the whole basis for keeping
  ``moist=False`` as the default; it is a *loose* benchmark, asserted as a band and an ordering, never to
  a number.

Named scope edges — what a higher rung would derive
---------------------------------------------------
* **Diagnostic, one-way.** The cooled surface does **not** feed back: it does not re-solve the EBM, does
  not change the albedo (a cold crest should grow snow — rung 0's ice-albedo step applied *per cell* would
  be the obvious next slice), and does not re-run the orographic model. Smith & Barstad's uptake
  coefficient ``C_w`` stays evaluated at the *upstream sea-level* condition, which is defensible (that is
  where S&B evaluates it) but is the clearest "a higher rung would derive this".
* **The classifier's cold bands are precipitation-independent** (:mod:`planet.biomes`: "a very wet
  sub-zero climate is still called boreal here"). The windward crest is *wet and cold*, so on exactly the
  cells this rung and rung 5A both act on, the Whittaker rule **ignores the orographic rain** and the two
  effects partly degenerate. Cited, not rediscovered — it is why the crest reads tundra regardless of how
  much it rains.
* **No cold-trap floor.** The march is a tropospheric one; :data:`MAX_ELEVATION_M` refuses terrain above
  the troposphere rather than silently clamping at :data:`planet.radiation.STRATOSPHERE_T`.
* **Sea-level reduction.** The EBM temperature is treated as a *sea-level* temperature. On a real planet a
  zonal mean over land already contains the continents' elevation, so applying the full ``Γz`` on top
  slightly double-counts; on this planet the zonal mean is genuinely elevation-free (there is no terrain in
  rung 0), so the treatment is exact *for this model* and approximate for Earth.

Units: ``T`` in **°C** (the EBM/Whittaker unit) on the public surface, kelvin only inside the march;
``z`` in **m**; ``Γ`` in **K m⁻¹** (rung 4's unit).

Sources: the constant and ``Γ_m`` are rung 4's ([[planet-rung4-radiation]]); the moist-adiabat formula and
its constants are pinned there. See [[planet-rung5a-orographic]].
"""
from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np

from planet.radiation import (
    GRAVITY,
    LAPSE_RATE,
    P_SURFACE,
    R_DRY_AIR,
    moist_adiabatic_lapse_rate,
)

KELVIN_0C = 273.15                  # K — °C ↔ K offset (public API is °C, the march is kelvin)
MAX_ELEVATION_M = 12_000.0          # m — terrain above this leaves the troposphere; refused, not clamped
DEFAULT_STEPS = 64                  # Heun steps per column — converged to <1e-3 K over 3 km (test-pinned)

# The one observational benchmark this module is scored against (loose, asserted as a band + an ordering).
# The deep tropics carry the planet's HIGHEST freezing level, cited at ~5000 m — Harris, Bowman & Shin
# (2000) J. Climate 13, 4137 (a 20-yr NCEP-reanalysis climatology of the 0 °C isotherm cross-checked
# against TRMM radar brightband heights), corroborated by Bradley et al. (2009) GRL 36, L17701. The band
# below is a generous read of "the tropics", its warm edge being the pinned number. It is what rejects the
# emergent moist adiabat as the default (which puts the level above 7 km) and keeps the pinned constant,
# which lands ~4.4 km — just BELOW the band, close but not inside it. The verdict is the ORDERING; do not
# upgrade it to "in band". Pinned in [[smith-barstad-orographic-source]].
OBSERVED_TROPICAL_FREEZING_LEVEL_M = (4500.0, 5000.0)

# The moist march's lapse-rate callable: rung 4's own Γ_m(T, p), used directly so it cannot drift.
MOIST_LAPSE_RATE_FN: Callable[..., np.ndarray] = moist_adiabatic_lapse_rate


def constant_lapse_rate_fn(lapse_rate: float = LAPSE_RATE) -> Callable[..., np.ndarray]:
    """A ``Γ(T, p)`` callable that ignores its arguments and returns ``lapse_rate`` everywhere.

    Handing this to :func:`column_temperature` makes the integrator reproduce the closed form
    ``T₀ − Γ·z`` to machine precision — the **reduction of the emergent path to the prescribed one**
    (the module's tight structural anchor: the two paths are the same code, not two implementations).
    """
    def _gamma(T_kelvin, p_pa):
        return np.full(np.shape(np.asarray(T_kelvin, dtype=float)), float(lapse_rate))
    return _gamma


def surface_temperature(T_sea_C, elevation_m, lapse_rate: float = LAPSE_RATE):
    """The terrain-cooled surface temperature (°C) under a **constant** lapse rate — the closed form.

    ``T = T_sea − Γ·z``, broadcast over ``T_sea_C`` and ``elevation_m``. This is the default path: ``Γ``
    is rung 4's pinned :data:`planet.radiation.LAPSE_RATE` (6.5 K/km), not a constant invented here.
    """
    T = np.asarray(T_sea_C, dtype=float)
    z = _checked_elevation(elevation_m)
    return T - float(lapse_rate) * z


def column_temperature(T_sea_C, elevation_m, *, lapse_rate_fn=None, n_steps: int = DEFAULT_STEPS,
                       p_surface: float = P_SURFACE):
    """Integrate a lapse rate up from the surface to each cell's own elevation → temperature (°C).

    A vectorised **Heun (predictor–corrector) march in height**, one column per cell, each with its own
    step ``dz = z/n_steps`` so every cell lands exactly on its own summit:

        ``dT/dz = −Γ(T, p)``,   ``dp/dz = −p·g/(R_d·T)``   (hydrostatic, layer-mean ``T``)

    ``lapse_rate_fn`` defaults to rung 4's saturated moist adiabat
    (:func:`planet.radiation.moist_adiabatic_lapse_rate`), which makes the cooling *emergent* — steep in
    cold air, flat in warm moist air. Passing :func:`constant_lapse_rate_fn` recovers
    :func:`surface_temperature` **exactly** (the reduction anchor). Second order in ``dz``.

    The pressure is integrated alongside because ``Γ_m`` depends on it; note that the *whole* z→p
    relation here is hydrostatic in the marched ``T``, so unlike
    :func:`planet.radiation.moist_adiabat_temperature` (which maps z→p through a **fixed** isothermal
    scale height) this path does not lean on any prescribed constant — the circularity is genuinely
    absent, not merely small.
    """
    gamma = MOIST_LAPSE_RATE_FN if lapse_rate_fn is None else lapse_rate_fn
    if n_steps < 1:
        raise ValueError(f"n_steps must be >= 1, got {n_steps}")

    T0 = np.asarray(T_sea_C, dtype=float)
    z = _checked_elevation(elevation_m)
    T, z = np.broadcast_arrays(T0, z)
    T = T.astype(float) + KELVIN_0C                      # march in kelvin
    p = np.full(T.shape, float(p_surface))
    dz = z / float(n_steps)

    for _ in range(n_steps):
        g1 = np.asarray(gamma(T, p), dtype=float)
        T_pred = T - g1 * dz
        p_next = p * np.exp(-GRAVITY * dz / (R_DRY_AIR * 0.5 * (T + T_pred)))    # hydrostatic over the layer
        g2 = np.asarray(gamma(T_pred, p_next), dtype=float)
        T = T - 0.5 * (g1 + g2) * dz
        p = p_next
    return T - KELVIN_0C


def terrain_temperature(T_sea_C, elevation_m, *, moist: bool = False, lapse_rate: float = LAPSE_RATE,
                        n_steps: int = DEFAULT_STEPS):
    """The one entry point the scene wires: terrain-cooled surface temperature (°C).

    ``moist=False`` (the default) is the closed-form constant-``Γ`` correction
    (:func:`surface_temperature`); ``moist=True`` integrates rung 4's emergent moist adiabat
    (:func:`column_temperature`). Both reduce to the input temperature where ``z = 0``.
    """
    if moist:
        return column_temperature(T_sea_C, elevation_m, n_steps=n_steps)
    return surface_temperature(T_sea_C, elevation_m, lapse_rate=lapse_rate)


def effective_lapse_rate(T_sea_C, elevation_m, *, moist: bool = True, lapse_rate: float = LAPSE_RATE,
                         n_steps: int = DEFAULT_STEPS):
    """The **column-averaged** rate actually realised over the terrain: ``(T_sea − T_summit)/z`` (K m⁻¹).

    For the constant path this is trivially ``lapse_rate``; for the moist path it is the diagnostic that
    lands the module's headline — at the mid-latitude demo range it comes out ≈ 6.3 K/km, within ~3 % of
    the 6.5 constant, while the tropics give ≈ 3.7 and the high latitudes ≈ 8.8. ``NaN`` where ``z = 0``
    (no column to average over).
    """
    z = _checked_elevation(elevation_m)
    T_top = terrain_temperature(T_sea_C, z, moist=moist, lapse_rate=lapse_rate, n_steps=n_steps)
    T0, z = np.broadcast_arrays(np.asarray(T_sea_C, dtype=float), z)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(z > 0.0, (T0 - T_top) / np.where(z > 0.0, z, 1.0), np.nan)


def freezing_level(T_sea_C, *, moist: bool = False, lapse_rate: float = LAPSE_RATE,
                   n_steps: int = DEFAULT_STEPS, z_max: float = MAX_ELEVATION_M, n_z: int = 121):
    """Altitude (m) of the **0 °C isotherm** above a surface at ``T_sea_C`` — the freezing level.

    The physical read of the lapse rate: how high you must climb before the air freezes. Under the
    constant rate it is the closed form ``T_sea/Γ``; under the moist adiabat the (strictly decreasing)
    profile is marched on an ``n_z``-point ladder and inverted by interpolation. Returns ``0`` where the
    surface is already at or below freezing, and is capped at ``z_max``.

    **Loose (benchmark), not a snow line.** The permanent snow line sits *above* the annual 0 °C level
    (summer melt, not the isotherm, sets it), so this is compared to observation only as a band and an
    ordering (:data:`OBSERVED_TROPICAL_FREEZING_LEVEL_M`).
    """
    T0 = np.asarray(T_sea_C, dtype=float)
    if not moist:
        return np.clip(np.where(T0 > 0.0, T0 / float(lapse_rate), 0.0), 0.0, z_max)

    ladder = np.linspace(0.0, z_max, n_z)                                       # (n_z,)
    profile = column_temperature(T0[..., None], ladder, n_steps=n_steps)        # (..., n_z), falls with z
    flat = np.atleast_2d(profile.reshape(-1, n_z))
    out = np.empty(flat.shape[0], dtype=float)
    for i, row in enumerate(flat):
        out[i] = 0.0 if row[0] <= 0.0 else float(np.interp(0.0, row[::-1], ladder[::-1]))
    return out.reshape(T0.shape)


@dataclass(frozen=True, eq=False)
class LapseDiagnostic:
    """The two lapse rates set side by side along a latitude sweep (``eq=False``: arrays).

    The object the demo figure plots and the rung's verdict is read off: for each latitude's sea-level
    temperature, the **effective** rate the moist adiabat realises over ``elevation_m``, the pinned
    constant it was meant to replace, and the **freezing level** each one implies — the leg on which the
    constant wins (:data:`OBSERVED_TROPICAL_FREEZING_LEVEL_M`).
    """

    latitude_deg: np.ndarray
    sea_level_T_C: np.ndarray
    gamma_constant: np.ndarray        # K m⁻¹ — flat by construction (the pinned rung-4 rate)
    gamma_moist: np.ndarray           # K m⁻¹ — the column-average the moist adiabat realises over `elevation_m`
    freezing_constant_m: np.ndarray   # m — 0 °C isotherm under the constant rate
    freezing_moist_m: np.ndarray      # m — 0 °C isotherm under the moist adiabat
    elevation_m: float


def lapse_diagnostic(latitude_deg, sea_level_T_C, *, elevation_m: float = 2500.0,
                     lapse_rate: float = LAPSE_RATE, n_steps: int = DEFAULT_STEPS) -> LapseDiagnostic:
    """Build a :class:`LapseDiagnostic` for a latitude sweep — the constant vs the emergent adiabat."""
    lat = np.asarray(latitude_deg, dtype=float)
    T = np.asarray(sea_level_T_C, dtype=float)
    return LapseDiagnostic(
        latitude_deg=lat,
        sea_level_T_C=T,
        gamma_constant=np.full(T.shape, float(lapse_rate)),
        gamma_moist=np.asarray(effective_lapse_rate(T, elevation_m, moist=True, n_steps=n_steps)),
        freezing_constant_m=np.asarray(freezing_level(T, lapse_rate=lapse_rate)),
        freezing_moist_m=np.asarray(freezing_level(T, moist=True, n_steps=n_steps)),
        elevation_m=float(elevation_m),
    )


def _checked_elevation(elevation_m) -> np.ndarray:
    """Elevation as a float array, refusing terrain that leaves the troposphere (named, not clamped)."""
    z = np.asarray(elevation_m, dtype=float)
    if np.any(z > MAX_ELEVATION_M):
        raise ValueError(
            f"elevation up to {float(np.max(z)):.0f} m exceeds MAX_ELEVATION_M={MAX_ELEVATION_M:.0f} m: "
            "the march is a tropospheric one (no cold-trap floor) — refused rather than silently clamped"
        )
    if np.any(z < 0.0):
        raise ValueError("elevation must be >= 0 (heights above the sea-level reference temperature)")
    return z
