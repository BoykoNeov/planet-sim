"""Planet Phase-2 validation: the Whittaker biome classifier (plan §3, the biome leg).

The classifier half of the biome-map triad (planet-earth-system.md). What is asserted **tight** is the
partition's *totality/determinism* (every (T,P) → exactly one biome, the map tiles the planet, area
fractions sum to 1) and the present-day **band ordering** (equator → pole: rain forest → savanna →
grassland/desert → temperate forest → boreal → tundra). The probe points are **independent canonical
facts** (rain forest = warm+wet, tundra = cold, desert = warm+dry), so matching them is not circular
with boundaries this module drew. What is **loose** is the calibrated threshold values/slopes (cited
[[whittaker-biome-source]]) — so the *absolute* biome latitudes move only in loose bands.
"""
import numpy as np
import pytest

from planet import biomes, precip
from planet.albedo import present_day_climate
from planet.biomes import Biome, classify, classify_field


# --------------------------------------------------------------------------- #
# Pinned thresholds ([[whittaker-biome-source]]) + ordering invariant
# --------------------------------------------------------------------------- #
def test_pinned_thresholds_and_their_ordering():
    assert biomes.TUNDRA_MAX_C == pytest.approx(-5.0)
    assert biomes.BOREAL_MAX_C == pytest.approx(3.0)
    assert biomes.TROPICAL_MIN_C == pytest.approx(20.0)
    # the three precip thresholds are strictly ordered across the model's whole temperature range
    T = np.linspace(-20.0, 35.0, 200)
    p_arid, p_semiarid, p_humid = biomes._precip_thresholds(T)
    assert np.all(p_arid < p_semiarid) and np.all(p_semiarid < p_humid)


# --------------------------------------------------------------------------- #
# Conservation = a *consistency* check: the partition is total and tiles the planet
# --------------------------------------------------------------------------- #
def test_classify_is_total_over_the_TP_plane():
    # every (T,P) point maps to exactly one VALID biome — no gaps, no errors (the §3 consistency leg)
    T = np.linspace(-40.0, 45.0, 120)
    P = np.linspace(0.0, 500.0, 120)
    TT, PP = np.meshgrid(T, P)
    codes = classify_field(TT, PP)
    assert codes.shape == TT.shape
    valid = {int(b) for b in Biome}
    assert set(np.unique(codes)).issubset(valid)               # only valid codes, every cell assigned


def test_area_fractions_sum_to_one():
    rng = np.random.default_rng(0)
    codes = classify_field(rng.uniform(-30, 40, 500), rng.uniform(0, 400, 500))
    fr = biomes.biome_area_fractions(codes)
    assert sum(fr.values()) == pytest.approx(1.0)
    # and on the actual planet (equal-area x = sin φ grid)
    st = present_day_climate(n_tau=0.02)
    fr_planet = biomes.biome_area_fractions(classify_field(st.T, precip.precip_field(st)))
    assert sum(fr_planet.values()) == pytest.approx(1.0)


def test_classify_scalar_matches_field():
    for T, P in [(28.0, 230.0), (18.0, 30.0), (-12.0, 30.0), (0.0, 50.0), (10.0, 250.0)]:
        assert classify(T, P) == Biome(int(classify_field(T, P)))


# --------------------------------------------------------------------------- #
# Analytical (tight): canonical probe points land in their textbook biome
# --------------------------------------------------------------------------- #
def test_probe_points_land_in_textbook_biomes():
    # independent canonical facts (not points drawn inside this module's own lines)
    probes = {
        (26.0, 350.0): Biome.TROPICAL_RAIN_FOREST,             # warm + very wet
        (25.0, 10.0): Biome.SUBTROPICAL_DESERT,                # warm + dry
        (24.0, 120.0): Biome.TROPICAL_SEASONAL_FOREST_SAVANNA, # warm + intermediate
        (-12.0, 30.0): Biome.TUNDRA,                           # cold (any precip)
        (0.0, 50.0): Biome.BOREAL_FOREST,                      # cool
        (12.0, 120.0): Biome.TEMPERATE_SEASONAL_FOREST,        # temperate + moist
        (10.0, 250.0): Biome.TEMPERATE_RAIN_FOREST,            # cool + very wet
        (15.0, 55.0): Biome.WOODLAND_SHRUBLAND,                # temperate + semi-arid
        (10.0, 25.0): Biome.TEMPERATE_GRASSLAND_DESERT,        # temperate + arid
    }
    for (T, P), want in probes.items():
        assert classify(T, P) == want, f"({T}°C, {P} cm) → {classify(T, P).name}, expected {want.name}"


def test_uniform_climate_planet_is_a_single_biome():
    # the analytical limit: a constant (T, P) planet maps everywhere to ONE biome
    codes = classify_field(np.full(50, 26.0), np.full(50, 300.0))
    assert np.all(codes == int(Biome.TROPICAL_RAIN_FOREST))


# --------------------------------------------------------------------------- #
# The diagonal (sloped) warm boundary + the vertical cold boundary
# --------------------------------------------------------------------------- #
def test_warm_boundary_slopes_a_fixed_precip_is_forest_when_cool_savanna_when_warm():
    # THE diagonal Whittaker structure: a forest needs more rain when warmer. The same P that supports
    # temperate forest at a cool temperature is only savanna/semi-arid at a hot one.
    P = 80.0
    assert classify(8.0, P) == Biome.TEMPERATE_SEASONAL_FOREST            # cool → forest
    assert classify(28.0, P) == Biome.TROPICAL_SEASONAL_FOREST_SAVANNA    # hot → savanna (drier-relative)


def test_cold_biomes_are_temperature_limited_independent_of_precip():
    for P in (5.0, 80.0, 300.0):
        assert classify(-10.0, P) == Biome.TUNDRA                         # below TUNDRA_MAX_C → tundra
        assert classify(0.0, P) == Biome.BOREAL_FOREST                    # [-5, 3) → boreal


# --------------------------------------------------------------------------- #
# Benchmark (loose): the present-day Earth band ordering, equator → pole
# --------------------------------------------------------------------------- #
def test_present_day_band_ordering():
    st = present_day_climate(n_tau=0.02)
    lat = st.latitude_deg()
    codes = classify_field(st.T, precip.precip_field(st))

    def biome_at(target_deg):
        return Biome(int(codes[int(np.argmin(np.abs(lat - target_deg)))]))

    assert biome_at(0) == Biome.TROPICAL_RAIN_FOREST                      # wet warm equator
    assert biome_at(15) == Biome.TROPICAL_SEASONAL_FOREST_SAVANNA         # savanna belt
    assert biome_at(30) in (Biome.SUBTROPICAL_DESERT,
                            Biome.TEMPERATE_GRASSLAND_DESERT)             # the subtropical deserts/grassland
    assert biome_at(45) == Biome.TEMPERATE_SEASONAL_FOREST                # midlatitude forest
    assert biome_at(60) == Biome.BOREAL_FOREST                            # boreal
    assert biome_at(75) == Biome.TUNDRA and biome_at(90) == Biome.TUNDRA  # polar tundra
    # a desert/grassland band genuinely appears in the subtropics (the great deserts at ~15–35°)
    sub = (lat >= 15.0) & (lat <= 35.0)
    desert_codes = {int(Biome.SUBTROPICAL_DESERT), int(Biome.TEMPERATE_GRASSLAND_DESERT)}
    assert desert_codes & set(np.unique(codes[sub]))
