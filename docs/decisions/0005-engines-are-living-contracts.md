# 0005 — Engines are living contracts (not frozen)

Status: Accepted — 2026-06-10 (user direction)
Scope: Program-level invariant. **Supersedes the freeze clauses of ADR 0001 and ADR 0003**
(see those ADRs' superseded-in-part notes). Keeps ADR 0001's array-boundary interface contract
and the validation-triad discipline intact.

## Context

The program's original doctrine **froze** each shared engine behind its one-page `CONTRACT.md`
once its validation suite passed: "reuse only frozen modules" (the invariant-5 framing), and
"changing the frozen surface ⇒ a new ADR + re-running the seal." That fit a *build-once,
reuse-many* shape — `engines/diffusion` sealed in Steel and reused unchanged by Microchip and
Planet.

But Planet's **GCM staircase** (`docs/plans/planet-earth-system.md` §5) makes the explicit
growth axis *extending `engines/fluid` at every rung*: an advected tracer (rung 1), vertical
layers (rung 3), the sphere (rung 5). The `engines/fluid` contract itself anticipates this — the
"GCM-climb seam" declares the stacked-field + tracer slots as planned extensions. Under the
freeze doctrine, **every rung would pay a per-change ADR + re-seal ceremony** — friction with no
payoff, because:

- the **validation suite** already guarantees correctness (and is run in full on any engine edit
  by the cross-cutting gate, ADR 0003), and
- an **additive** extension (e.g. advecting the tracer only when one is present) leaves the prior
  behaviour *bit-for-bit identical*, which a regression test proves directly.

The trigger was rung 1 needing to advect the long-declared `engines/fluid` tracer slot. The
freeze ceremony added nothing the suite + gate do not already provide.

## Decision

1. **Engines are living, versioned contracts — not frozen artifacts.** A consumer may **extend
   an engine directly** when it needs to, provided:
   - (a) the existing suite stays **green** — and where the change is additive, the prior path is
     **bit-for-bit unchanged** (proven by a regression test);
   - (b) the **new surface gets its own tests** (the engine's validation triad grows with it);
   - (c) the change is recorded in the engine `CONTRACT.md` **Changelog**.
2. **The guardrail is the test suite + the gate, not a freeze.** Editing a shared engine remains
   the cross-cutting **full-repo-gate** case (ADR 0003) and runs the **import-drift guard** —
   that is what catches a regression or an undeclared dependency.
3. **No per-change ADR.** Extending an engine needs no ADR. (A genuinely *breaking* change —
   removing or renaming a public surface a consumer relies on — still warrants updating those
   consumers in the same change; that is ordinary engineering, not a freeze ceremony.)

**Explicitly NOT superseded (still in force):**
- **ADR 0001's array-boundary interface contract** — only plain arrays / numeric records cross
  the per-step boundary (the compiled-core / extension seam). A *stable interface* is not a
  *frozen implementation*.
- **The validation-triad discipline** (ARCHITECTURE.md §7) — *more* central now: the suite **is**
  the contract.
- **ADR 0003's gate mechanics** — unchanged; only its freeze-based *justification* was reworded
  ("frozen" → "untouched": an engine's tests can only break when the engine is edited).

## Consequences

- `+` The GCM climb (and any future engine growth) proceeds **without per-rung ceremony**; the
  extension seam the contracts already advertised is usable as intended.
- `+` **One source of truth** for "is the engine correct": its suite. No freeze status to keep in
  sync with reality.
- `+` **Honest:** the engines were never truly immutable (they were built extension-ready by
  design); the doctrine now matches the practice.
- `−` Loses the freeze's blunt "an engine cannot change under you" guarantee. *Mitigation:*
  additive changes are bit-for-bit (a regression test proves it), the full gate runs on every
  engine edit, and the Changelog records surface changes. A consumer pins behaviour by its
  **tests**, not by a freeze.
- `−` Slightly more responsibility on the editor: keep the suite green and add tests for the new
  surface (the freeze previously *forbade* the edit; it is now allowed, but must be validated).

## Alternatives considered

- **Keep the freeze; pay the per-change ADR + re-seal each rung.** Rejected (user, 2026-06-10):
  pure friction on a project whose design *is* iterative engine growth; the ceremony duplicates
  what the suite + gate already enforce.
- **Un-freeze only `engines/fluid`** (the one rung 1 touches). Rejected: the freeze was a
  program-level invariant; a split doctrine (one engine frozen, one not) is incoherent, and the
  rationale applies to both engines equally.
- **Drop `CONTRACT.md` entirely.** Rejected: the one-page contract is still the unit of context
  downstream code loads, and the home of the stable interface + the guaranteed invariants + the
  Changelog — value independent of the freeze.
