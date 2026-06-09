# 0004 — Interactive deep-end views: the layer registry & state interchange

Status: Accepted — 2026-06-09
Scope: Program-level; **extends ADR 0002** (visualization & UX). First instance:
the Planet interactive map (`docs/plans/planet-earth-system.md` §9).

## Context

ADR 0002 §4 reserves a *selective deep-end* (Plotly / web / WebGL) for the sims
whose payoff is inherently interactive and spatial — "planet maps," a 3-D galaxy.
It fixed the **floor** (a universal matplotlib static figure) and the **boundary**
(viz consumes already-validated arrays, never a live solver — §2/§3), but it left
three things unspecified, because no project had reached the deep end yet. Planet
is the first that does, and building its interactive map surfaced exactly those
gaps:

1. **A deep-end view accretes across a project's phases.** Planet's map grows
   temperature → ice line → biomes → precipitation → circulation across Phases
   1–4, and further up the §5 GCM staircase. A view re-written each phase is the
   cross-cutting-change hazard ARCHITECTURE.md §6 warns against, in miniature.
2. **The interaction model is not fixed across compute tiers.** At rung-0 compute
   (laptop-seconds) a knob can drive an instant recompute-and-repaint; as compute
   climbs the staircase (3-D, many timesteps) that live loop becomes impossible —
   but nothing said so, inviting a future session to over-engineer to preserve it.
3. **Users want the world to round-trip.** Take a simulated planet *out* (share,
   inspect) and bring an externally-authored one *in* (a future geography-editing
   app paints elevation/coastlines, the model imports it). That needs a stable
   interchange contract, not an ad-hoc array dump.

## Decision

**1. Deep-end views are a *layer registry*.** A view is an ordered stack of
self-contained layers, each `(name, kind, data-array, style, z-order)` with
`kind ∈ {scalar field, vector/line overlay, annotation}`. The renderer is generic
over `kind`; a phase contributes to a view by **registering layers**, never by
editing the renderer. This is precisely ADR 0002 §3's "2-D field/heatmap" +
"annotated overlay" primitives instantiated for one concrete surface. The registry
stays **project-local** until a third consumer (after the project-local
field/animation code already living in steel/chip/planet `plots.py`) promotes the
primitives to the shared `viz/` toolkit by rule-of-three (ARCHITECTURE.md §6). The
registry does **not** pre-empt that promotion — it is the consumer that will
eventually *earn* it.

**2. The interaction model is a function of the compute tier; the renderer is
not.** Because the renderer only ever consumes arrays (the ADR 0001/0002
boundary), it is **invariant across the whole staircase** — only the *trigger*
changes:

> **live knob → instant remap** (light compute) **→ set parameters → launch a run
> → view the result** (heavy compute).

This is named as an explicit staircase consequence (planet plan §5) so a future
session does not try to keep the live loop alive past where compute allows. The
slider is a *driver of compute*; the map is a *consumer of arrays*; they are
separate concerns, and only the former is tier-dependent.

**3. State interchange: pin a *schema*, not a file format.** The portable artifact
is a documented **state spec** —

- the **grid geometry**,
- **explicit units** (self-describing, so an external consumer cannot misread it —
  this program is unit-obsessive by §7 discipline),
- the **layer list** (the registry of Decision 1 *is* the export manifest — one
  structure serialized, not a second one invented),
- the **parameter / knob values**, and
- a **`schema_version`** for forward/backward compatibility (round-trip "at any
  point of development" needs versioning).

The **encoding** is chosen per *consumer*, behind that schema:

- *lean default* — a JSON manifest + numeric arrays (`.npz`);
- *editable-geometry layers* (heightmaps — elevation, bathymetry, land/ocean
  masks) interchange as **16-bit grayscale images** (PNG): the native currency of
  paint/terrain tools and web canvases, which is the round-trip the editing-app use
  case actually needs (8-bit is too coarse for elevation — noted);
- *domain-standard* encodings (e.g. **NetCDF** for gridded geophysical fields) are
  documented future options behind the same schema.

The tension this resolves: the climate-domain-standard format (NetCDF) and the
"easy for a web/JS editor to read" goal pull in **opposite** directions (NetCDF in
a browser needs WASM or a server). Resolving by *consumer* — and pinning the schema,
not the encoding — keeps both reachable without betting the contract on either.

**4. Round-trip identity is a *real* correctness property — test it.** Unlike a
rendered view (whose only test is an execution smoke-test, ADR 0002 §2),
`import(export(state)) == state` (array identity within the schema) is a genuine
invariant and gets a genuine test. **Importing an externally-authored state and
running a model that *responds* to it is a separate, staircase-gated capability:**
until the consuming physics exists, an imported geometry is **inert** — carried,
displayed, and round-tripped, but not yet changing the climate. That honesty
boundary is named, not blurred (it is the same "preplanned, not built" discipline
as the geography seam in the planet plan).

## Consequences

- `+` A deep-end view grows by adding layers, not by editing a renderer — the
  per-phase accretion the deep end demands is cheap and blast-radius-localized.
- `+` The renderer survives the entire compute staircase unchanged; only the
  separately-owned compute trigger changes tier to tier.
- `+` A stable, versioned, self-describing interchange contract lets an external
  editor round-trip without coupling to internals — and gives the deep end one
  place with a *real* test amid the smoke-tests.
- `+` Reinforces the ADR 0001/0002 array boundary as serving a **fourth** consumer
  (an external interchange format), after tests, compiled reimplementations, and
  the render layer — more evidence the boundary is drawn in the right place.
- `−` A schema with a version field and several encodings is more surface than a
  single dump; mitigated by building only the lean encoding in v1 and documenting
  the rest behind the schema (build the seam, not the machinery).
- `−` "Inert on import until the physics climbs" can disappoint ("I painted a
  mountain and the climate did not change"); mitigated by stating it plainly in the
  UI/docs and on the staircase, where the consuming rung is named.

## Alternatives considered

- **A monolithic per-phase view function** — simplest at Phase 1, but every later
  phase edits it (the cross-cutting hazard). Rejected for the layer registry, which
  is the same primitives ADR 0002 already names.
- **Promote field/animation to shared `viz/` now** — premature: rule-of-three is
  satisfied by *count* (three `plots.py`) but not by *substance* (they share
  conventions and styling, not copy-pasted code). The layer registry is the third
  *consumer* that earns the promotion later; pre-promoting is the "premature
  abstraction" ADR 0002 §3 explicitly guards against.
- **Pick NetCDF as the interchange format** — domain-standard but browser-hostile,
  contradicting the actual interchange consumer (a future web editor). Rejected for
  pin-the-schema + per-consumer encoding.
- **Keep the live-slider loop at every tier** — impossible once a single run is
  minutes-to-hours. Rejected for the tier-dependent interaction model.
