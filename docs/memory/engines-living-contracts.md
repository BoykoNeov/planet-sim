---
name: engines-living-contracts
description: "Engines are NO LONGER frozen (ADR 0005, 2026-06-10) — living/versioned contracts extended directly; supersedes the \"freeze-before-reuse\" framing in older memories"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 88f3e2c8-5355-42ce-bbb4-2fe776581b65
---

**The freeze-before-reuse doctrine was DROPPED on 2026-06-10 (user direction).** `docs/decisions/0005-engines-are-living-contracts.md` makes `engines/*` **living, versioned contracts** — extend an engine **directly** when a consumer needs it; no per-change ADR, no "re-seal". The guardrails that replace the freeze: (a) the engine's validation suite stays green (additive changes bit-for-bit, proven by a regression test), (b) the new surface gets its own tests, (c) record it in the `CONTRACT.md` **Changelog** (each contract now has one). Editing a shared engine still triggers the **full-repo gate + import-drift guard** (ADR 0003, mechanics unchanged).

**Why:** Planet's GCM staircase *is* iterative engine growth (rung 1 advected tracer → rung 3 layers → rung 5 sphere); the freeze ceremony was pure friction duplicating what the suite + gate already enforce. Triggered by rung 1 needing to advect `engines/fluid`'s long-declared `tracer` slot.

**How to apply:** Do NOT propose "new ADR + re-run the seal" to change an engine — just extend + test + changelog. **Older memories that invoke "freeze before reuse" / "reuse only frozen modules" / "frozen behind CONTRACT.md" ([[planet-plan]], [[planet-phase3-engine]], [[planet-phase4-coupler]]) describe the pre-0005 doctrine — historically accurate, now superseded for go-forward work.** KEPT intact: ADR 0001's array-boundary interface (plain arrays cross the per-step boundary) + the validation-triad discipline. Status headers flipped FROZEN→LIVING in both `engines/{fluid,diffusion}/CONTRACT.md`. Part of starting [[planet-rung1-two-way-coupler]].
