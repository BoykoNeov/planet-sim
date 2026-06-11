---
name: memory-vendored-in-repo
description: auto-memory files live in the repo at docs/memory/, junctioned to the global Claude path
metadata:
  type: project
---

As of 2026-06-11 the Claude auto-memory files for planet-sim are vendored into
the repo at `docs/memory/` (committed) instead of living only in the per-machine
global Claude store. The global path
`C:\Users\boiko\.claude\projects\M--claud-projects-planet-sim\memory` is a
Windows **directory junction** pointing at `M:\claud_projects\planet-sim\docs\memory`,
so harness recall and future memory writes still resolve through the hardcoded
global path but the bytes land in the repo.

**Why:** the global store is per-machine and not version-controlled; vendoring it
makes the memories durable, reviewable, and travel with the code. The junction
keeps the harness working without any path reconfiguration.

**How to apply:** on a fresh clone / new machine the junction won't exist (git
stores the files, not the link), so the harness would create an empty global
memory dir and recall would miss this history. To restore: move/merge any global
files into `docs/memory/`, delete the empty global dir, then recreate the link —
`New-Item -ItemType Junction -Path "<global memory path>" -Target "<repo>\docs\memory"`
(no admin needed; works across the C:→M: volumes). Keep editing memories the
normal way — writes flow through the junction into the repo; commit them per
[[always-push-commits]].
