"""The Whittaker biome classifier — climate → habitability (Planet Phase 2, the consequence payoff).

This is the project's **dramatic end-to-end win** (plan §3, banked early): two climate inputs — the
EBM temperature ``T(φ)`` (:mod:`planet.ebm`) and the diagnostic precipitation ``P(φ)``
(:mod:`planet.precip`) — are mapped to a **biome** at every latitude, producing the bands
of life that **migrate as the planetary knobs turn**. It is the planet analogue of Steel's
*microstructure* and Chip's *device*: knobs → climate → **habitability**, no lookup table.

The classifier is an **original, total partition** of the ``(T, P)`` plane — not an embedded copy of
Whittaker's figure (the **Irvin precedent**, Microchip Phase 1a: a copyrighted *graphical* diagram is
reproduced by an *independent* computation calibrated to it and benchmarked loosely, never digitized
and redistributed). The boundary thresholds are **read off** Whittaker's diagram and pinned as cited
constants ([[whittaker-biome-source]]); the partition is then a deterministic rule that covers the
*whole* plane (no gaps — so the map tiles the planet and biome area fractions sum to 1, the §3
consistency leg).

The partition's structure mirrors the diagram's two physics
-----------------------------------------------------------
* **Cold biomes are temperature-limited** → the **tundra** (T < ``TUNDRA_MAX_C``) and **boreal forest**
  (``TUNDRA_MAX_C`` ≤ T < ``BOREAL_MAX_C``) boundaries are *vertical* (precip-independent): nothing
  grows tall where it is cold enough, regardless of rain.
* **Warm biomes are moisture-limited, and the boundaries slope** → above ``BOREAL_MAX_C`` the
  desert → semi-arid → forest → rain-forest transition happens at **higher precipitation as it gets
  warmer** (warm air evaporates more, so a forest needs more rain). The three precip thresholds are
  therefore **linear in T** (sloped, not axis-aligned — the diagonal Whittaker boundaries). The same
  three precip bins carry different biome *names* in the **temperate** band (``BOREAL_MAX_C`` ≤ T <
  ``TROPICAL_MIN_C``) and the **tropical** band (T ≥ ``TROPICAL_MIN_C``) — the only thing the warm
  T-boundary does is relabel, so it introduces no precip-threshold jump.

Validation triad (plan §3) — what is asserted tight vs loose
------------------------------------------------------------
* **Analytical (tight).** The classifier is an **exact, deterministic partition**: a uniform-climate
  planet maps to a single biome, and chosen ``(T, P)`` **probe points** land in their textbook biome.
  The probe points are **independent canonical facts** (tropical rain forest = warm + wet, tundra =
  cold, subtropical desert = warm + dry, …), *not* points drawn inside lines this module invented —
  the non-circularity guard (the boundaries are calibrated to Whittaker; the probes are facts).
* **Conservation — a *consistency*/partition check, honestly weaker (named).** A classifier has no
  energy or mass law. The honest leg: every ``(T, P)`` maps to **exactly one** biome (totality), the
  biome map **tiles the planet with no unclassified gaps**, and the **area fractions sum to 1**
  (x = sin φ is the equal-area coordinate, so a biome's area fraction is its share of the cells). A
  consistency check, not a conservation law — stated plainly.
* **Benchmark (loose).** Present-day Earth reproduces the **observed major biome bands** in order
  equator → pole: tropical rain forest, then savanna/desert near 15–30°, temperate forest/grassland
  in midlatitudes, boreal forest, then tundra toward the poles ([[whittaker-biome-source]]). The
  *absolute* band latitudes depend on the calibrated precip param, so they are asserted in loose bands.

Non-circularity, named scope edge (plan §3)
-------------------------------------------
*Validated tight:* the partition's totality/determinism and the present-Earth band **ordering**.
*Calibrated/flagged (loose):* the threshold values + slopes (cited to Whittaker/Ricklefs), so absolute
biome latitudes move only in loose bands. *Scope edge:* **Whittaker** (annual ``T, P``) is used, not
**Köppen** — Köppen needs *seasonal/monthly* precip the annual-mean v1 does not produce (named). No
continentality / orography (zonal-mean planet); precip is prescribed (:mod:`planet.precip`);
no dynamic vegetation or carbon feedback. The cold-limited bands are precip-independent (a
simplification: a very wet sub-zero climate is still called boreal here — outside the planet's actual
trajectory).

Units — T in °C, P in cm/yr (Whittaker's axes), shared with the precip module
-----------------------------------------------------------------------------
``T`` in **°C** (the EBM's unit), ``P`` in **cm/yr** (Whittaker's unit and the precip module's).
"""
from __future__ import annotations

from enum import IntEnum

import numpy as np


class Biome(IntEnum):
    """The nine Whittaker biomes (cited [[whittaker-biome-source]]), as integer codes for the map.

    Ordered coldest → warmest then dry → wet within a thermal band, so the code value reads roughly as
    "polar → tropical, arid → humid" — convenient for an ordered colour map. The integer *value* is
    the biome field's pixel value (a scalar layer the map paints); :data:`BIOME_NAMES` / :data:`BIOME_COLORS`
    give the label and colour.
    """

    TUNDRA = 0
    BOREAL_FOREST = 1
    TEMPERATE_GRASSLAND_DESERT = 2
    WOODLAND_SHRUBLAND = 3
    TEMPERATE_SEASONAL_FOREST = 4
    TEMPERATE_RAIN_FOREST = 5
    SUBTROPICAL_DESERT = 6
    TROPICAL_SEASONAL_FOREST_SAVANNA = 7
    TROPICAL_RAIN_FOREST = 8


BIOME_NAMES: dict[int, str] = {
    Biome.TUNDRA: "Tundra",
    Biome.BOREAL_FOREST: "Boreal forest",
    Biome.TEMPERATE_GRASSLAND_DESERT: "Temperate grassland / desert",
    Biome.WOODLAND_SHRUBLAND: "Woodland / shrubland",
    Biome.TEMPERATE_SEASONAL_FOREST: "Temperate seasonal forest",
    Biome.TEMPERATE_RAIN_FOREST: "Temperate rain forest",
    Biome.SUBTROPICAL_DESERT: "Subtropical desert",
    Biome.TROPICAL_SEASONAL_FOREST_SAVANNA: "Tropical seasonal forest / savanna",
    Biome.TROPICAL_RAIN_FOREST: "Tropical rain forest",
}

# Stable, roughly ecology-conventional colours (cold/pale → warm/green); the render layer's palette.
BIOME_COLORS: dict[int, str] = {
    Biome.TUNDRA: "#c7d4d8",
    Biome.BOREAL_FOREST: "#6f9b7a",
    Biome.TEMPERATE_GRASSLAND_DESERT: "#d9c789",
    Biome.WOODLAND_SHRUBLAND: "#bfae5c",
    Biome.TEMPERATE_SEASONAL_FOREST: "#5fa364",
    Biome.TEMPERATE_RAIN_FOREST: "#2f7d5a",
    Biome.SUBTROPICAL_DESERT: "#e3c16f",
    Biome.TROPICAL_SEASONAL_FOREST_SAVANNA: "#a7c34f",
    Biome.TROPICAL_RAIN_FOREST: "#1f6b3b",
}

# --------------------------------------------------------------------------- #
# Pinned Whittaker boundary thresholds ([[whittaker-biome-source]]).
# Read off Whittaker (1975) / Ricklefs (2008) Fig 5.5 and pinned as cited constants — NOT a digitized
# copy of the figure (the Irvin precedent). Temperatures (°C) are the cold-limited vertical cuts; the
# precip thresholds (cm/yr) are LINEAR IN T (the diagonal, moisture-limited boundaries): a forest needs
# more rain when it is warmer. Values are loose / calibration-dependent (the non-circularity split).
# --------------------------------------------------------------------------- #
TUNDRA_MAX_C = -5.0            # °C — below this it is tundra (cold-limited), any precipitation
BOREAL_MAX_C = 3.0            # °C — [TUNDRA_MAX_C, this) is boreal forest (cold-limited)
TROPICAL_MIN_C = 20.0        # °C — at/above this the warm bins take tropical names (else temperate)

# Sloped precip thresholds p(T) = intercept + slope·T  (cm/yr). p_arid < p_semiarid < p_humid.
P_ARID_INTERCEPT, P_ARID_SLOPE = 25.0, 0.8          # desert ↔ semi-arid
P_SEMIARID_INTERCEPT, P_SEMIARID_SLOPE = 50.0, 2.0  # semi-arid ↔ sub-humid (forest)
P_HUMID_INTERCEPT, P_HUMID_SLOPE = 150.0, 2.0       # sub-humid ↔ humid (rain forest)


def _precip_thresholds(T: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The three sloped precip thresholds ``(p_arid, p_semiarid, p_humid)`` (cm/yr) at temperature ``T``.

    Linear in T (the diagonal Whittaker boundaries): the rain needed to cross from desert → forest →
    rain forest rises with warmth. Strictly ordered ``p_arid < p_semiarid < p_humid`` over the model's
    temperature range, so they bin precipitation into arid / semi-arid / sub-humid / humid.
    """
    p_arid = P_ARID_INTERCEPT + P_ARID_SLOPE * T
    p_semiarid = P_SEMIARID_INTERCEPT + P_SEMIARID_SLOPE * T
    p_humid = P_HUMID_INTERCEPT + P_HUMID_SLOPE * T
    return p_arid, p_semiarid, p_humid


def classify_field(T, P) -> np.ndarray:
    """Classify temperature ``T`` (°C) and precipitation ``P`` (cm/yr) arrays into :class:`Biome` codes.

    Vectorized total partition (the §3 consistency leg: every cell gets exactly one biome, no gaps).
    Cold cells are temperature-limited (tundra / boreal); warm cells are binned by the sloped precip
    thresholds, with the temperate/tropical name set by ``TROPICAL_MIN_C``. Returns an integer array
    (the biome codes) broadcast to the shape of ``T``/``P`` — the scalar layer the biome map paints.
    """
    T = np.asarray(T, dtype=float)
    P = np.asarray(P, dtype=float)
    T, P = np.broadcast_arrays(T, P)
    p_arid, p_semiarid, p_humid = _precip_thresholds(T)
    tropical = T >= TROPICAL_MIN_C

    # First-match priority list (np.select): cold cuts, then tropical bins, then temperate bins.
    conditions = [
        T < TUNDRA_MAX_C,                                   # tundra (cold-limited)
        T < BOREAL_MAX_C,                                   # boreal forest (cold-limited)
        tropical & (P < p_arid),                            # subtropical desert
        tropical & (P < p_semiarid),                        # tropical savanna (semi-arid)
        tropical & (P < p_humid),                           # tropical savanna (sub-humid) — merged name
        tropical,                                           # tropical rain forest (humid)
        P < p_arid,                                         # temperate grassland / desert
        P < p_semiarid,                                     # woodland / shrubland
        P < p_humid,                                        # temperate seasonal forest
    ]
    choices = [
        int(Biome.TUNDRA),
        int(Biome.BOREAL_FOREST),
        int(Biome.SUBTROPICAL_DESERT),
        int(Biome.TROPICAL_SEASONAL_FOREST_SAVANNA),
        int(Biome.TROPICAL_SEASONAL_FOREST_SAVANNA),
        int(Biome.TROPICAL_RAIN_FOREST),
        int(Biome.TEMPERATE_GRASSLAND_DESERT),
        int(Biome.WOODLAND_SHRUBLAND),
        int(Biome.TEMPERATE_SEASONAL_FOREST),
    ]
    # default = the only remaining case: temperate band, P ≥ p_humid → temperate rain forest.
    return np.select(conditions, choices, default=int(Biome.TEMPERATE_RAIN_FOREST)).astype(int)


def classify(T: float, P: float) -> Biome:
    """Classify a single ``(T, P)`` point (°C, cm/yr) into its :class:`Biome` — the scalar partition."""
    return Biome(int(classify_field(np.asarray(T), np.asarray(P))))


def biome_area_fractions(codes: np.ndarray) -> dict[Biome, float]:
    """Area fraction of each present biome from a field of codes on the equal-area ``x = sin φ`` grid.

    Because ``x = sin φ`` makes every cell equal area on the sphere, a biome's planetary area fraction
    is simply its share of the grid cells. Returns ``{Biome: fraction}`` for the biomes present; the
    fractions **sum to 1** (the §3 consistency leg — the map tiles the planet, no gaps). Pass a 1-D
    latitudinal field (the v1 zonal-mean planet) or any equal-area field.
    """
    codes = np.asarray(codes, dtype=int).ravel()
    n = codes.size
    values, counts = np.unique(codes, return_counts=True)
    return {Biome(int(v)): float(c) / n for v, c in zip(values, counts)}
