---
name: teaching-surfaces-resynced
description: notebook + README + interactive re-synced to rungs 1–4 (2026-06-15); the generated-HTML / big-notebook / DOCS_FIGURE gotchas
metadata:
  type: project
---

**Teaching surfaces brought up to date with the built rungs, 2026-06-15.** Audit finding: the
interactive what-if was current (4 knobs), but the **notebook + `planet/README.md` described built
rungs (1–4) as "still ahead"**. User chose *prose-sync + grow notebook*.

What changed:
- **`planet/planet.ipynb`** gained a new group **"§8 — Up the staircase"** (inserted before the
  provenance cell): four showcase sections — **rung 2** wet-get-wetter (`demo_wet_get_wetter`),
  **rung 2.5** polar amplification (`moist_ebm.polar_amplification`, inline figure), **rung 3**
  baroclinic QG turbulence (`demo_baroclinic_qg`), **rung 4** spectral log law
  (`demo_spectral_band`). Cheap demos run **live** + embed the banked PNG; the QG turbulence run is
  ~minutes so it **embeds-only**. Stale forward-prose de-staled: §6 "what's next" now reports rungs
  1–4 built; §3's two "named gaps" (band migration, energy-constrained rate) point at the modules
  that closed them (`circ_precip`/`sphere_ebm`, `moist`).
- **`planet/README.md`** rung log extended Rung 1-step-3 → **rungs 2 / 2.5 / 2.x / 3 / 4** + a
  "teaching surfaces re-synced" bullet.
- **Interactive header** gained the ocean knob.

Non-obvious gotchas (worth keeping):
- **`docs/interactive/index.html` is GENERATED** by `interactive.write_app()` and **byte-pinned** by
  the slow `test_committed_page_is_up_to_date` (committed == fresh full-grid build). So edit `_BODY`
  in `interactive.py` and **regenerate** (`python -m planet.interactive`, ~minutes, recomputes the
  S0×CO2×tilt×ocean grid) — never hand-edit the HTML, or that test goes red.
- **`planet.ipynb` is too big for the Read tool** (banked image outputs ~27k tokens); edit it by
  loading/patching the JSON in a script (markdown-source edits don't touch outputs). New code-cell
  outputs were banked by a **surgical nbclient run of only the new cells** + output injection, so
  existing banked outputs stayed byte-identical.
- New cells use **`demo.DOCS_FIGURE` (absolute)** for `Image()`, not a relative `../docs/figures`
  path — the notebook smoke test (`test_planet_notebook`) runs the kernel with cwd = **repo root**,
  so a relative figure path would break.
- `demo_wet_get_wetter` self-labels "rung-1 evidence" but lives in `moist.py` (rung 2); §8.1 carries
  a one-line note reconciling the printed label.

[[interactive-what-if]] [[planet-rung3-qg-built]] [[planet-rung25-mse-diffusion]] [[planet-rung1-two-way-coupler]]
