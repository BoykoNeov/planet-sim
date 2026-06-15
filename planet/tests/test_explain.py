"""Guards for the rule-based what-if explainer (:mod:`planet.explain`).

The explainer is the single source of the "what changed + why" prose shown in the notebook's live
cells and baked into the browser what-if. These tests pin the two things that matter for that role:
the prose tracks the *sign* of the computed change (it never says "warming" when the model cooled),
and both depths are produced — without asserting exact wording (which is free to be tuned).
"""
from __future__ import annotations

from planet.albedo import A_OLR, EBMParams, S0_EARTH
from planet.demo_biomes import compute
from planet.explain import Diagnostics, Knobs, diagnose, explain
from planet.obliquity import OBLIQUITY_EARTH

_BASE = diagnose(compute(EBMParams()))

# A synthetic present-day-ish baseline for the prose-only tests below (explain() never recomputes —
# it narrates the deltas — so the regime can be dialled directly via Diagnostics).
_BASE_DIAG = Diagnostics(global_mean_T=14.7, ice_line_lat=73.0, rainforest_pct=10.0,
                         tundra_pct=8.0, desert_pct=15.0)


def _diag(Tbar: float, ice_line_lat: float) -> Diagnostics:
    return Diagnostics(global_mean_T=Tbar, ice_line_lat=ice_line_lat,
                       rainforest_pct=20.0, tundra_pct=2.0, desert_pct=12.0)


def _explain(**knob_kw):
    """Run the model for the given knob overrides and explain it against present-day Earth."""
    knobs = Knobs(**knob_kw)
    result = compute(EBMParams(S0=knobs.S0, A=knobs.A, D=knobs.D))
    return knobs, explain(knobs, _BASE, diagnose(result))


def test_diagnose_reads_the_model():
    d = diagnose(compute(EBMParams()))
    assert 13.0 < d.global_mean_T < 16.0           # present-day ~14.7 °C
    assert 70.0 < d.ice_line_lat < 75.0            # finite polar cap ~73°
    assert 0.0 <= d.rainforest_pct <= 100.0
    assert not d.ice_free and not d.snowball


def test_baseline_says_baseline():
    ex = explain(Knobs(), _BASE, _BASE)
    assert "baseline" in ex.headline.lower()
    assert "present-day earth" in ex.oneline.lower()
    assert ex.paragraph and ex.oneline


def test_both_depths_present_and_distinct():
    _, ex = _explain(A=A_OLR - 8)
    assert ex.headline and ex.oneline and ex.paragraph
    assert ex.paragraph != ex.oneline              # the paragraph is the *fuller* mechanism
    assert len(ex.paragraph) > len(ex.oneline)


def test_warming_reads_as_warming():
    """Lower A (more greenhouse) warms — the prose must say so, and name the mechanism."""
    _, ex = _explain(A=A_OLR - 8)
    text = (ex.oneline + " " + ex.paragraph).lower()
    assert "warm" in text and "cool" not in ex.oneline.lower()
    assert "greenhouse" in text and "infrared" in text


def test_cooling_reads_as_cooling():
    """A dimmer sun (still above the Snowball cliff) cools without freezing over."""
    knobs, ex = _explain(S0=1300.0)
    assert "cool" in ex.oneline.lower() and "warm" not in ex.oneline.lower()
    assert "sun" in ex.oneline.lower()


def test_snowball_names_the_hysteresis():
    """Dim far enough → Snowball, and the caveat must flag the path-dependence (the lookup can't show it)."""
    _, ex = _explain(S0=1235.0)
    assert "snowball" in ex.headline.lower()
    assert "hysteresis" in ex.paragraph.lower() or "path-dependent" in ex.paragraph.lower()


def test_hothouse_flags_ice_free():
    _, ex = _explain(A=A_OLR - 20)
    assert "ice-free" in ex.headline.lower() or "hothouse" in ex.headline.lower()


def test_multiple_knobs_each_appear():
    """Two knobs moved → both causes are narrated (not just the dominant one)."""
    _, ex = _explain(S0=1330.0, A=A_OLR - 4)
    text = ex.paragraph.lower()
    assert "sun" in text and "greenhouse" in text


def test_obliquity_cooling_note_fires_ice_free_and_tilt_up():
    """Ice-free + more tilt → the paragraph carries the geometric-cooling twist (and only there)."""
    ex = explain(Knobs(obliquity_deg=45.0), _BASE_DIAG, _diag(20.0, 90.0))
    para = ex.paragraph.lower()
    assert "redistributes" in para and "cools the planet slightly" in para
    assert "more reflective than its tropics" in para


def test_obliquity_cooling_note_absent_when_ice_present():
    """With a cap still present, more tilt *warms* (melts ice) — the cooling note must NOT appear."""
    ex = explain(Knobs(obliquity_deg=45.0), _BASE_DIAG, _diag(5.0, 55.0))
    assert "redistributes" not in ex.paragraph.lower()


def test_obliquity_cooling_note_absent_when_tilt_not_raised():
    """Ice-free but reached by another knob (tilt at the baseline) → no tilt-cooling claim."""
    ex = explain(Knobs(obliquity_deg=OBLIQUITY_EARTH, A=A_OLR - 20), _BASE_DIAG, _diag(35.0, 90.0))
    assert "redistributes" not in ex.paragraph.lower()


def test_obliquity_cooling_note_absent_when_tilt_lowered():
    """A *smaller* tilt starves the poles (a steeper gradient) — not the redistribution-cooling case."""
    ex = explain(Knobs(obliquity_deg=6.0), _BASE_DIAG, _diag(35.0, 90.0))
    assert "redistributes" not in ex.paragraph.lower()


def test_deterministic():
    a = explain(Knobs(A=A_OLR - 6), _BASE, diagnose(compute(EBMParams(A=A_OLR - 6))))
    b = explain(Knobs(A=A_OLR - 6), _BASE, diagnose(compute(EBMParams(A=A_OLR - 6))))
    assert a == b
