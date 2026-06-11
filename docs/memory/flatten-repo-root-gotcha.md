---
name: flatten-repo-root-gotcha
description: "Post-flatten parents[N] off-by-one + CI editable-install masks bare-checkout test failures"
metadata: 
  node_type: memory
  type: project
  originSessionId: 8fbc0ac7-489a-4716-ba8a-72a7aa11fee1
---

The standalone-flatten (commit 438a222, 2026-06-10) removed one directory level
(monorepo `BigSim/planet/...` → standalone `planet-sim/...`) but did NOT decrement
every `Path(__file__).resolve().parents[N]` filesystem-root computation. Fixed in
two passes: `fd6305c` caught only `test_planet_notebook.py`; the remaining nine were
fixed 2026-06-10 — tests at `planet/tests/` need `parents[2]`, demos + `planetmap.py`
at `planet/` need `parents[1]` (each builds its figure/HTML bank path or a subprocess
`cwd` from `_REPO_ROOT`).

**Why it stayed hidden:** the full-gate CI does `pip install -e ".[...]"` (full-gate.yml:57),
so `import planet` resolves from site-packages regardless of `cwd` — the wrong index
passes in CI. The two `*_stays_headless` subprocess tests only go RED in a **bare,
non-installed checkout**. So "167 green" in the README is true only WITH the editable
install.

**How to apply:** (1) When verifying the suite, run it in the bare checkout condition
(no `pip install -e`) to catch cwd/root bugs CI masks. (2) Any NEW file that computes a
filesystem root via `parents[N]` is layout-depth-sensitive — count from the file's actual
depth, don't copy-paste the index from a file at a different depth. A guard asserting
`_REPO_ROOT` contains the `planet/` package would kill this class permanently (offered,
not yet built). See [[planet-plan]].
