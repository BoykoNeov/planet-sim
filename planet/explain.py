"""Plain-language *what changed + why* for a climate what-if — rule-based, keyed to deltas.

The **single source** of the explanation prose shown on every interactive surface. The notebook's
live what-if cells call :func:`explain` on each re-run; the browser interactive page bakes the same
strings into its precomputed grid (``planet.interactive``). One implementation, so the two surfaces
can never drift.

It is deliberately **rule-based, not generative**: every clause is templated and keyed to the sign
and magnitude of a *computed* delta, so the prose can never claim more than the validated model
produced — the same anti-overclaim discipline the rest of the simulator holds (nothing is fit in
prose; every number traces to the model). Two depths are offered, per request:

* :attr:`Explanation.oneline` — a single causal-chain sentence (knob → mechanism → result).
* :attr:`Explanation.paragraph` — the fuller mechanism: every active knob, the temperature and
  gradient response, the ice line, the biome migration, and any honest caveat (the Snowball is
  path-dependent; an ice-free planet's bands have stopped moving).

NumPy-only, no viz/notebook deps — it imports on a bare core install and is unit-tested directly.
"""
from __future__ import annotations

from dataclasses import dataclass

from planet.albedo import A_OLR, D_TRANSPORT, S0_EARTH
from planet.biomes import Biome
from planet.exoplanet import T_SUN
from planet.obliquity import OBLIQUITY_EARTH
from planet.ocean import OCEAN_FRACTION_EARTH

# Ice-line sentinels (verified against the model: hothouse pins 90.0, a frozen planet pins 0.0).
_ICE_FREE_LAT = 89.5     # ice line at/above this ≈ ice-free
_SNOWBALL_LAT = 1.0      # ice line at/below this ≈ frozen to the equator (Snowball)


@dataclass(frozen=True)
class Knobs:
    """The knobs every interactive surface exposes; defaults are present-day Earth (the baseline)."""

    S0: float = S0_EARTH
    A: float = A_OLR
    D: float = D_TRANSPORT
    obliquity_deg: float = OBLIQUITY_EARTH
    ocean_fraction: float = OCEAN_FRACTION_EARTH
    T_star: float = T_SUN
    size: float = 1.0


@dataclass(frozen=True)
class Diagnostics:
    """The handful of climate/habitability numbers the prose reasons over (all read from the model)."""

    global_mean_T: float    # °C
    ice_line_lat: float     # degrees; ``>= _ICE_FREE_LAT`` ice-free, ``<= _SNOWBALL_LAT`` Snowball
    rainforest_pct: float   # % of the planet that is tropical rain forest
    tundra_pct: float       # % tundra
    desert_pct: float       # % desert (subtropical desert + temperate grassland/desert)

    @property
    def ice_free(self) -> bool:
        return self.ice_line_lat >= _ICE_FREE_LAT

    @property
    def snowball(self) -> bool:
        return self.ice_line_lat <= _SNOWBALL_LAT


def diagnose(result) -> Diagnostics:
    """Pull the reasoning numbers out of a :class:`planet.demo_biomes.BiomeResult`."""
    def pct(*biomes: Biome) -> float:
        return 100.0 * sum(result.fraction(b) for b in biomes)

    return Diagnostics(
        global_mean_T=float(result.state.global_mean_T),
        ice_line_lat=float(result.state.ice_line_lat),
        rainforest_pct=pct(Biome.TROPICAL_RAIN_FOREST),
        tundra_pct=pct(Biome.TUNDRA),
        desert_pct=pct(Biome.SUBTROPICAL_DESERT, Biome.TEMPERATE_GRASSLAND_DESERT),
    )


@dataclass(frozen=True)
class Explanation:
    """The rendered prose at two depths plus a short title (see the module docstring)."""

    headline: str       # a short label for the panel title
    oneline: str        # one causal-chain sentence
    paragraph: str      # the fuller mechanism (3–6 sentences)


# --- the rule base ---------------------------------------------------------------------------- #
# Each active knob contributes (cause-clause, gradient-effect or None). ``cause`` is phrased as
# "<doing X> <mechanism>"; clauses are stitched into both depths. Tolerances keep a knob that
# only nudged off its default (slider rounding) from generating spurious prose.
_TOL = {"S0": 0.5, "A": 0.5, "D": 0.005, "obliquity_deg": 0.05, "ocean_fraction": 0.02,
        "T_star": 1.0, "size": 0.005}


def _causes(knobs: Knobs, base: Knobs) -> tuple[list[str], list[str]]:
    """Return (cause clauses, gradient-effect clauses) for every knob that moved off its baseline."""
    causes: list[str] = []
    gradient: list[str] = []

    dA = base.A - knobs.A                                   # +dA  == more greenhouse (lower A = warmer)
    if abs(dA) > _TOL["A"]:
        if dA > 0:
            causes.append(f"adding greenhouse gas (a {dA:.0f} W/m² cut in how readily heat escapes to "
                          "space) traps outgoing infrared")
        else:
            causes.append(f"stripping out greenhouse gas (a {-dA:.0f} W/m² rise) lets more infrared "
                          "escape to space")

    if abs(knobs.S0 - base.S0) > _TOL["S0"]:
        pct = 100.0 * (knobs.S0 - base.S0) / base.S0
        if pct > 0:
            causes.append(f"a brighter sun delivers {pct:.0f}% more sunlight")
        else:
            causes.append(f"a dimmer sun delivers {-pct:.0f}% less sunlight")

    if abs(knobs.D - base.D) > _TOL["D"]:
        if knobs.D > base.D:
            causes.append("stronger heat transport carries more warmth from equator to pole")
            gradient.append("flattening the equator-to-pole temperature gradient (poles warm, "
                            "tropics cool)")
        else:
            causes.append("weaker heat transport leaves more warmth stuck in the tropics")
            gradient.append("steepening the equator-to-pole gradient (colder poles, hotter tropics)")

    if abs(knobs.obliquity_deg - base.obliquity_deg) > _TOL["obliquity_deg"]:
        if knobs.obliquity_deg > base.obliquity_deg:
            causes.append(f"a larger axial tilt ({knobs.obliquity_deg:.0f}°) sends more yearly-average "
                          "sunlight to the poles")
            gradient.append("flattening the gradient (the poles warm)")
        else:
            causes.append(f"a smaller axial tilt ({knobs.obliquity_deg:.0f}°) starves the poles of "
                          "sunlight")
            gradient.append("steepening the gradient (the poles cool)")

    if abs(knobs.ocean_fraction - base.ocean_fraction) > _TOL["ocean_fraction"]:
        pct = round(100.0 * knobs.ocean_fraction)
        if knobs.ocean_fraction > base.ocean_fraction:
            causes.append(f"a wetter world — {pct}% ocean, more sea and less land — darkens the "
                          "surface, so it soaks up more sunlight")
            gradient.append("flattening the gradient (the added ocean carries more heat poleward)")
        else:
            causes.append(f"a drier world — {pct}% ocean, more land and less sea — brightens the "
                          "surface, so it reflects more sunlight")
            gradient.append("steepening the gradient (less ocean to carry heat poleward)")

    if abs(knobs.T_star - base.T_star) > _TOL["T_star"]:
        if knobs.T_star > base.T_star:
            causes.append("a hotter, bluer star shifts its light out of the near-infrared, so snow and "
                          "ice reflect more of it (a stronger ice-albedo)")
        else:
            causes.append("a cooler, redder star puts more of its light in the near-infrared, which "
                          "snow and ice absorb better (a weaker ice-albedo)")

    if abs(knobs.size - base.size) > _TOL["size"]:
        if knobs.size > base.size:
            causes.append(f"a larger planet ({knobs.size:.2f} R⊕) spreads the same heat transport over "
                          "more area (D ∝ 1/size²), so less warmth reaches the poles")
            gradient.append("steepening the gradient")
        else:
            causes.append(f"a smaller planet ({knobs.size:.2f} R⊕) concentrates the heat transport "
                          "(D ∝ 1/size²), so more warmth reaches the poles")
            gradient.append("flattening the gradient")

    return causes, gradient


def _join(clauses: list[str]) -> str:
    """Oxford-comma join of clauses, capitalising the first letter of the result."""
    if not clauses:
        return ""
    if len(clauses) == 1:
        s = clauses[0]
    elif len(clauses) == 2:
        s = f"{clauses[0]} and {clauses[1]}"
    else:
        s = ", ".join(clauses[:-1]) + f", and {clauses[-1]}"
    return s[0].upper() + s[1:]


def _warm_phrase(dT: float) -> str:
    """A finite temperature-response clause graded by magnitude (the verb carries the sign)."""
    a = abs(dT)
    if a < 0.25:
        return "the global mean temperature barely changes"
    adverb = "slightly " if a < 1.5 else "" if a < 4 else "sharply " if a < 8 else "dramatically "
    verb = "warms" if dT > 0 else "cools"
    return f"the planet {adverb}{verb} by {a:.1f} °C"


def _ice_phrase(base: Diagnostics, now: Diagnostics) -> str:
    """A finite ice-line / regime clause comparing now to the baseline."""
    if now.snowball:
        return "it freezes over entirely — a Snowball, ice all the way to the equator"
    if now.ice_free:
        return "the last of the polar ice melts away — an ice-free hothouse"
    d = now.ice_line_lat - base.ice_line_lat
    if abs(d) < 0.6:
        return f"the polar ice cap holds near {now.ice_line_lat:.0f}°"
    direction = "retreats poleward" if d > 0 else "advances toward the equator"
    return (f"the ice line {direction}, from {base.ice_line_lat:.0f}° to {now.ice_line_lat:.0f}°")


def _biome_sentence(base: Diagnostics, now: Diagnostics) -> str:
    """How the habitability bands moved (the §3 payoff), keyed to the rain-forest/tundra deltas."""
    if now.snowball:
        return ("Almost nothing is habitable: the biome bands have collapsed under planet-wide ice.")
    dr = now.rainforest_pct - base.rainforest_pct
    dt = now.tundra_pct - base.tundra_pct
    if abs(dr) < 0.5 and abs(dt) < 0.5:
        return (f"The biome bands barely move: tropical rain forest still covers {now.rainforest_pct:.0f}% "
                f"of the planet and tundra {now.tundra_pct:.0f}%.")
    shift = "expand poleward" if dr > 0 or dt < 0 else "retreat toward the equator"
    return (f"The warm biome bands {shift}: tropical rain forest now covers {now.rainforest_pct:.0f}% of "
            f"the planet (was {base.rainforest_pct:.0f}%) and cold tundra {now.tundra_pct:.0f}% "
            f"(was {base.tundra_pct:.0f}%).")


def _caveat(now: Diagnostics) -> str:
    """An honest footnote where the physics has a wrinkle the static number hides."""
    if now.snowball:
        return ("Getting back out is harder than getting in: a frozen planet reflects so much sunlight "
                "that re-brightening the sun to today's value won't melt it — the climate is "
                "path-dependent (hysteresis). The notebook's §2 runs both branches live.")
    if now.ice_free:
        return ("With no ice left, the ice-albedo feedback is switched off, so further warming now moves "
                "the bands far more gently than it did while there was still ice to melt.")
    return ""


def _obliquity_cooling_note(knobs: Knobs, baseline: Knobs, now: Diagnostics) -> str:
    """Ice-free + more tilt: the geometric cooling the ice-albedo warming normally hides.

    With no cap left to melt, the ice-albedo amplifier is off, and raising the tilt only
    *redistributes* the same total sunlight (the global mean is ``S₀/4`` at every obliquity) onto the
    model's intrinsically brighter high latitudes (the fixed ``a₂`` poleward-albedo term), so a little
    less is absorbed and the planet cools — exactly the ``−(S₀/20)·a₂·Δs₂/B`` cross-term. The sign is
    guaranteed *only* in this regime, so the clause fires only when the world is ice-free **and** the
    tilt was raised above the baseline; with a cap still present, more tilt instead *warms* (it melts
    that ice — the opposite sign), so the note must not appear there.
    """
    if not now.ice_free:
        return ""
    if knobs.obliquity_deg <= baseline.obliquity_deg + _TOL["obliquity_deg"]:
        return ""
    return ("There is a subtler twist: tilt only redistributes the planet's sunlight — its total "
            "intake is fixed — and an ice-free planet's poles are slightly more reflective than its "
            "tropics, so steering more sunlight poleward has it absorb a little less. With the "
            "ice-melting bonus exhausted, that is why adding still more tilt now cools the planet "
            "slightly instead of warming it.")


def _ocean_caveat(knobs: Knobs, baseline: Knobs) -> str:
    """The two honest limitations of the ocean knob, fired only when the ocean fraction has moved.

    Both are real ceilings of an *equilibrium annual-mean* EBM whose precipitation is a prescribed
    pattern (see :mod:`planet.ocean`): the sea's heat-capacity role drops out of the steady state, and
    the rain bands do not migrate with the water. Surfaced so the prose never lets the knob imply more
    than the model carries.
    """
    if abs(knobs.ocean_fraction - baseline.ocean_fraction) <= _TOL["ocean_fraction"]:
        return ""
    return ("Two things this steady-state view leaves out: the ocean's real superpower — its huge heat "
            "capacity, which softens the seasons — does not show here, because this is the equilibrium "
            "annual-mean climate; and the rain bands do not move with the water (precipitation follows a "
            "fixed latitude pattern), so a wetter world rains a little more everywhere rather than "
            "somewhere new.")


def explain(knobs: Knobs, base_diag: Diagnostics, now_diag: Diagnostics,
            baseline: Knobs = Knobs()) -> Explanation:
    """Render the *what changed + why* prose for moving from ``baseline`` knobs to ``knobs``.

    ``base_diag`` is the baseline climate's diagnostics, ``now_diag`` the current one — both produced
    by the model (this function never recomputes anything; it only narrates the deltas).
    """
    causes, gradient = _causes(knobs, baseline)
    dT = now_diag.global_mean_T - base_diag.global_mean_T

    if not causes:
        return Explanation(
            headline="Present-day Earth — the baseline",
            oneline=("This is present-day Earth: global mean "
                     f"{now_diag.global_mean_T:.1f} °C, a polar ice cap near "
                     f"{now_diag.ice_line_lat:.0f}°. Move a knob to see what changes."),
            paragraph=("Nothing has moved off its present-day value yet — this is the baseline every "
                       "what-if below is measured against. The planet sits at a global-mean "
                       f"{now_diag.global_mean_T:.1f} °C with the ice line near {now_diag.ice_line_lat:.0f}°; "
                       f"tropical rain forest covers {now_diag.rainforest_pct:.0f}% of it and tundra "
                       f"{now_diag.tundra_pct:.0f}%. Brighten or dim the sun, add greenhouse gas, tilt "
                       "the axis or cover more of it in ocean, and the chain from knob to climate to "
                       "biomes plays out here."))

    cause_text = _join(causes)
    warm = _warm_phrase(dT)
    ice = _ice_phrase(base_diag, now_diag)

    # --- one-line causal chain --- #
    oneline = f"{cause_text} — {warm}, and {ice}."

    # --- headline --- #
    if now_diag.snowball:
        headline = "Snowball — the planet has frozen over"
    elif now_diag.ice_free:
        headline = f"Ice-free hothouse — {now_diag.global_mean_T:.0f} °C, no polar cap"
    elif abs(dT) < 0.25:
        headline = "A different planet, a familiar temperature"
    else:
        word = "Warmer" if dT > 0 else "Cooler"
        headline = f"{word} world — {dT:+.1f} °C, ice line at {now_diag.ice_line_lat:.0f}°"

    # --- paragraph --- #
    parts = [f"{cause_text}."]
    grad = f", {_join(gradient).lower()}" if gradient else ""
    parts.append(f"In response, {warm}{grad}, and {ice}.")
    parts.append(_biome_sentence(base_diag, now_diag))
    caveat = _caveat(now_diag)
    if caveat:
        parts.append(caveat)
    obl_note = _obliquity_cooling_note(knobs, baseline, now_diag)
    if obl_note:
        parts.append(obl_note)
    ocean_note = _ocean_caveat(knobs, baseline)
    if ocean_note:
        parts.append(ocean_note)
    paragraph = " ".join(parts)

    return Explanation(headline=headline, oneline=oneline, paragraph=paragraph)


# --- the Snowball (cold) branch: a DELIBERATE second prose category --------------------------- #
# explain() above narrates the WARM branch as a function of knob *deltas*. The Snowball branch is a
# different kind of statement: the SAME knobs, but a world that *began frozen* — the bistability /
# path-dependence, which a knob-delta rule base structurally cannot express (no knob moved). So this
# is intentionally a SEPARATE generator, not a drift from explain() — do not "unify" them. The two
# stay consistent by tracing the same demo_snowball numbers (the freeze and re-melt folds below),
# and the warm-branch _caveat() already cross-references this branch. Kept numpy-free like the module.
_SNOWBALL_FREEZE_S0 = 1252      # ≈ demo_snowball.compute().freeze_S0 — the warm branch's fold (dim → freeze)
_SNOWBALL_MELT_S0 = 1831        # ≈ demo_snowball.compute().melt_S0   — the Snowball's fold (bright → re-melt)


def snowball_branch_explain(knobs: Knobs, diag: Diagnostics,
                            warm_is_snowball: bool = False) -> Explanation:
    """Narrate the cold (Snowball) branch — the same knobs as :func:`explain`, but a world that *began
    frozen*. The warm and cold branches together *are* the bistability; this names the path-dependence
    the warm-branch prose only footnotes. ``warm_is_snowball`` flags the dimmed corner where even a
    warm start freezes (the hysteresis loop has closed), so the prose never claims a temperate twin
    that does not exist there. ``knobs`` is accepted for symmetry with :func:`explain` (the branch's
    character is the same across tilt/ocean — a frozen planet ignores them — so only ``diag`` is read).
    """
    Tbar = diag.global_mean_T
    if warm_is_snowball:
        headline = "Snowball — and now the only climate"
        oneline = (f"Started frozen, this world is a Snowball at {Tbar:.0f} °C — but here the Sun is so "
                   f"dim that even a warm start freezes: below about {_SNOWBALL_FREEZE_S0:.0f} W/m² the "
                   "bistability is gone and the Snowball is the only stable climate.")
        twin = ("At this dim a Sun the two branches have merged — a warm start no longer survives, so "
                "there is just one climate left, the frozen one.")
    else:
        headline = "Snowball — frozen by its own history"
        oneline = (f"Started frozen, this world stays a Snowball at the same Sun — ice reflects most of "
                   f"the sunlight, so it holds itself at {Tbar:.0f} °C, frozen pole to pole, while a "
                   "world that started warm sits temperate at the very same settings.")
        twin = ("Turn no knob — just change how the world began. A planet that started warm settles "
                "temperate; one that started frozen stays a Snowball: two stable climates for one Sun "
                "(path-dependence, or hysteresis).")
    paragraph = (
        f"{twin} A white planet reflects most of the sunlight that reaches it, which is what holds it "
        f"frozen at a global-mean {Tbar:.0f} °C. Its tilt and ocean barely matter — a frozen surface is "
        "uniformly reflective whatever lies beneath — so only the Sun and the greenhouse nudge its "
        "(still sub-zero) temperature. Climbing back out is far harder than falling in: you would have "
        f"to brighten the Sun by about a third, past ~{_SNOWBALL_MELT_S0:.0f} W/m², before a Snowball "
        "would melt. The notebook's §2 traces the full loop and its catastrophic jumps live.")
    return Explanation(headline=headline, oneline=oneline, paragraph=paragraph)
