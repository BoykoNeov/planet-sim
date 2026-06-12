# 0002 — Visualization & UX strategy

Status: Accepted — 2026-06-08 (status note appended 2026-06-09 — see end:
`viz/` is not yet built; rule-of-three holds by *count* but not by *substance*)
Scope: Program-level invariant; inherited by every per-project plan.

## Context

Education and experimentation are two of the program's three core targets
(ARCHITECTURE.md §1), yet "interactivity" was named as a goal with no mechanism:
nothing in the architecture said how a simulator is *shown* or *explored*, and
the §10 plan template had no visualization component. This ADR sets that
doctrine.

Forces in play:

- Engines are headless and data-oriented (arrays / numeric records out) behind
  their stable data contracts (ADR 0001). Visualization must not contaminate that.
- Correctness is established by the validation triad on the *numbers*
  (ARCHITECTURE.md §7) — not by how a plot looks.
- The product teaches mechanisms and runs what-ifs on a laptop; it is not a
  maximum-fidelity rendering engine.
- Built by an LLM agent under a fixed context window across the portfolio: the
  viz layer must be agent-fluent and must not need loading when working on a
  solver.
- Compute is light by design (scope ceilings, ADR 0001), so a
  slider → re-run → re-plot loop is responsive without special engineering.

## Decision

**1. Separate compute from render.** Engines never import a plotting library.
Visualization consumes the same plain data the validation tests and any compiled
reimplementation consume — the array-out contract (ADR 0001) serves all three.

**2. Visualization is never in the correctness path.** A figure is a consumer of
already-validated data, never evidence of validity. Test the numbers (the triad),
then draw them. This guards against the "looks plausible ⇒ correct" trap the
whole validation methodology exists to prevent.

**3. A shared `viz/` toolkit, peer to `engines/`, governed by the same rules.**
Reusable primitives — line/series, 2-D field/heatmap, time-animation,
parameter-sweep comparison grid, annotated overlay — promoted from project-local
to shared only by rule-of-three (ARCHITECTURE.md §6). The same primitives feed
every rendering target.

**4. Progressive enhancement** (mirrors the bank-an-artifact-then-layer phasing):

- *Floor (universal):* a matplotlib static figure per phase — the banked
  artifact; testable against golden/numeric; zero deployment.
- *Interactive:* notebooks (matplotlib + ipywidgets) for teaching narrative,
  and/or a thin Streamlit/Gradio app for slider-driven what-ifs.
- *Selective deep-end:* Plotly / web / WebGL only where a specific sim's payoff
  demands it (3-D galaxy, planet maps). Opt-in, behind the same data boundary.

**5. Visualize the mechanism, not just the output.** The target is "teach real
mechanisms, not black boxes," so views are designed to reveal *why* — e.g.,
animating a cooling path across the TTT C-curve so the learner sees why a fast
quench misses the nose and lands in martensite — rather than plotting a bare
result.

**The "both" consistency guard.** "Support both notebooks and a web app" means
the **toolkit** targets both (they share the primitives and the headless engine)
and the matplotlib floor is universal — *not* that every sim must ship both. Each
project builds the interactive surface its pedagogy needs. Steel, as the flagship,
ships both a notebook and a thin Streamlit app as the demonstrator; later sims
pick what fits. This keeps the doctrine lean-by-default (ARCHITECTURE.md §6).

## Consequences

- `+` Engines stay testable, headless, and small in context; the viz layer never
  loads when working on a solver.
- `+` One set of primitives serves static, interactive, and web targets — DRY by
  rule-of-three, not premature abstraction.
- `+` The array-out contract (ADR 0001) now demonstrably serves three consumers
  (tests, compiled reimpl, viz) — evidence the boundary is right.
- `−` Two interactive targets (notebook + app) is more surface than one; mitigated
  by making the floor universal and the interactive surface per-need, not
  mandatory-both.
- `−` Animations / web views can drift from the headless data if discipline slips;
  the rule "viz consumes validated data only" must hold.

## Alternatives considered

- **Pure notebooks (no web app)** — great teaching narrative, but sharing means
  "reader runs the notebook"; weaker for a shareable interactive what-if. Kept as
  a target, not the sole one.
- **A unified web/JS frontend (D3 / three.js) for all sims** — max polish, but a
  two-language split against the agent-context-economy doctrine and a large
  maintenance surface. Reserved for selective deep-end cases only.
- **Native GUI (Qt / Tk)** — worse for sharing and education than web/notebook;
  rejected.
- **Plotly / Bokeh everywhere as the floor** — heavier and less agent-fluent than
  matplotlib for the universal static floor; used selectively at the deep end
  instead.

## Status note — 2026-06-09: `viz/` not yet built; the rule-of-three substance check

Three projects now each carry a project-local `plots.py` (steel, chip, planet), so
Decision 3's rule-of-three trigger is satisfiable **by count**. A read of all three
bodies (the "visual engine" session, 2026-06-09) found it is **not** satisfied **by
substance**: they share **conventions and styling** — near-identical doctrine
docstrings, color-constant blocks, the `axvline/axhline + marker + annotate-arrow`
idiom, grid/legend defaults — but **almost no copy-pasted code**. Of the named
primitives, the 2-D field/heatmap and the stacked bars are effectively *steel-only*
today, and **time-animation exists nowhere** in the repo. Promoting on the strength of
the count would build the thin styling shim Decision 3 itself warns against ("DRY by
rule-of-three, **not premature abstraction**").

**So `viz/` stays unbuilt, deliberately.** A *thin* extraction (the doctrine docstring +
axis/marker/styling helpers) is the leaned-to next step — **not yet done**; the
genuinely additive primitives (a real 2-D field/heatmap and a **time-animation**) are
better built *with* their first consumer than extracted from thin air. The Planet
**interactive-map layer registry** (ADR 0004 #1) — which exercises exactly the 2-D-field
and annotated-overlay primitives — is the awaited **third consumer** that will *earn* the
promotion to `viz/`. This note is recorded so a future session does not read the three
`plots.py` as a standing mandate to stand up `viz/` on the count alone.

## Status note — 2026-06-12: honest-by-disclosure for the showcase layer

A scoped refinement to Decision 2 (*a figure is never evidence of validity*), recorded for the
planet-sim visualization "Rung C" (the *Perpetual-Ocean*-style particle globe; build plan §9.5,
decided with the user 2026-06-12). It does **not** weaken the doctrine — it follows from it.

Most planet-sim figures are **honest-by-construction**: the geometry cannot misrepresent the model
(the eddy globe is one true latitude band, never a 360° wrap). A pure **showcase** renderer, whose
purpose is *reach and delivery, not teaching*, is granted a second mode — **honest-by-disclosure**:
it **may** present a view the model does not literally compute (a global-looking flow, particles that
*imply* "currents carrying heat" that the ~90 %-reversible flux does not), **provided a visible
on-screen disclaimer documents the departure.** This is licensed *precisely because* a figure is never
in the correctness path (Decision 2): an illustration cannot fabricate validity, since validity lives
only in the tested numbers.

The carve-out is **narrow and asymmetric**, and one guardrail is non-negotiable:

- **Physics-fidelity verification relaxes** for the showcase only — approximate rendering is fine, no
  byte-golden, no numerical-transport proof. The engine/science layer keeps the full validation-triad
  discipline untouched; so do the honest-by-construction figures.
- **Documentation verification tightens.** Because the disclaimer **is the entire license**, it is the
  one thing **machine-checked**: a structural test asserts the artifact carries the honesty caption,
  on-screen and legible to a casual viewer. **Principle: the documentation is machine-checked even when
  the physics is not.**

So Decision 2 stands — *test the numbers, then draw them; a figure proves nothing.* This note adds only
that a showcase figure may *illustrate beyond* the numbers when it *says so on its face*, and that the
saying-so is itself tested. Scope: this is applied at planet-sim Rung C (§9.5); it is recorded here, in
the program-level viz ADR, so a future session reads it as the sanctioned form rather than an erosion.
