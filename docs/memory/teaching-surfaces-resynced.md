---
name: teaching-surfaces-resynced
description: "notebook §8 re-synced to rungs 1–4 (2026-06-15) then EXTENDED to 2.x/5A/5B (2026-07-10, now 7 sections); the generated-HTML / big-notebook / DOCS_FIGURE gotchas"
metadata: 
  node_type: memory
  type: project
  originSessionId: 44ec96b4-30cb-4ab9-b968-c8dc0d64f3af
---

**UPDATE 2026-07-10 — §8 extended from four showcase sections to SEVEN** (commit `feat(pedagogy)`):
added **§8.3 rung 2.x** (energetic ITCZ + radiation-limit; `demo_sphere_itcz` + `demo_itcz_radiation_limit`,
both live — fast), **§8.6 rung 5A** (orographic rain shadow; NEW `demo_orographic.py` module + banked
`planet-orographic.png` + `test_demo_orographic.py` — embed-only in the notebook, the jet spin-up is ~2 min),
**§8.7 rung 5B** (continentality→ice→2-D map; `demo_seasonal` live, `demo_seasonal_ice`/`demo_seasonal_map`
figures embed-only). Renumber shifted old rung-3/4 to §8.4/§8.5 (fix the two README `§8.N` cross-refs too).
Gotchas that held: patch the JSON via a script (notebook too big for Read/NotebookEdit); bank new-cell
outputs with a **surgical nbclient run of ONLY the new self-contained code cells** (`json.dump(indent=1,
ensure_ascii=False)` reproduced the original formatting → diff stayed +248/−8, no reformat churn); the
figure smoke test must render into a **tmp_path** (monkeypatch DOCS_FIGURE/OUTPUT_FIGURE) so the fast
`use_jet=False` render never clobbers the committed emergent-jet figure. **Always re-run
`test_planet_notebook.py` on the FINAL committed state** (advisor-caught: an intermediate bank had passed,
the re-banked jet-figure state had not). NB: the README *rung log* was missing dedicated 5B.1+/5B.2 bullets
(the 5B.1 bullet still framed 5B.2 as future) — folded a note into the new pedagogy bullet, did not backfill.

---
Original 2026-06-15 sync (rungs 1–4) below.

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
