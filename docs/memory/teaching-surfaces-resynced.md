---
name: teaching-surfaces-resynced
description: "notebook §8 re-synced to rungs 1–4 (2026-06-15), EXTENDED to 2.x/5A/5B (2026-07-10), then §8.8/§8.9/§8.10 (2026-09-02/04) → TEN sections, plus the 2026-09-04 predict-then-check + missions pass; the generated-HTML / big-notebook / surgical-re-bank / DOCS_FIGURE gotchas"
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

**Extended again 2026-09-02 — §8.8 added** (cells 51–53: markdown + a LIVE code cell running
`bifurcation.equilibrium_curve(EBMParams(n_cells=720))` (0.1 s; output banked *surgically* by running the cell's
source in a subprocess and pasting stdout as a `stream` output — no full nbclient re-execution) + an embed-only
cell for `planet-bifurcation.png`, `planet-seasonal-ice-map.png`, `planet-seasonal-ice-map.gif`). §8 heading
de-numbered ("the research rungs beyond rung 0", eight sections §8.1–§8.8); the "Where the numbers come from"
table gained two rows. The notebook's original file ends WITHOUT a trailing newline — strip it after
`json.dump` to keep the diff honest. `docs/index.html` regenerated for the two new catalogue rows
(`bifurcation`, `seasonal_ice_map` + its `interactive=` globe).
- **Notebook-hang data point (2026-09-02):** `test_planet_notebook` failed ONCE with the documented "Ubuntu" signature
  (the §2 snowball cell timing out) when run under `pytest -n 3` alongside the slow 2-D-march tests on a 4-core Linux
  container, then passed clean solo. Correlates with CPU contention, not content; run the notebook test solo (or the full
  gate with fewer workers) before blaming a cell. Root cause still open.

**Extended again 2026-09-04 — §8.9 added** (cells 54–56: markdown + a LIVE code cell running rung 5B.4's
`seasonal_sici.annual_mean_curve` + a 4-point `hysteresis_loop` at 360 cells / 90 steps — **~9 s**, the
slowest live cell in §8, so keep the grid at 360 if it is ever re-tuned — plus an embed-only cell for
`planet-seasonal-sici.png`). Output banked **surgically** again (run the cell's source in a subprocess,
paste stdout as a `stream` output; `execution_count` = max + 1 = 21). §8 now has **nine** sections
§8.1–§8.9. `docs/index.html` regenerated for the new `seasonal_sici` catalogue row. Two gotchas re-learned:
`json.dumps(nb, indent=1, ensure_ascii=False)` round-trips this notebook **byte-identically** (verify that
before editing — it is the cheapest proof the rewrite is diff-clean), and the file still ends with **no**
trailing newline. Also: `plots.py` is **not** imported by most of the suite, so a syntax error there passes
the fast lane and surfaces only in the notebook test — `python -c "import planet.plots"` after every edit.


**Extended again 2026-09-04 — §8.10 added** (cells 57–59: markdown + a LIVE code cell running rung 5A.4's
`demo_alpine_biomes.compute(use_jet=False)` — **~1.4 s**, cheap enough to run live unlike §8.6's jet-sourced
5A cell — plus an embed-only cell for `planet-alpine-biomes.png`). §8 now has **ten** sections §8.1–§8.10.
Same surgical bank (subprocess → `stream` output, `execution_count` = 22). Two things to carry:
- the subprocess needs `env=dict(os.environ, PYTHONIOENCODING="utf-8")` or its stdout comes back cp1252 on
  Windows and the capture dies on `→` / `°`;
- the §8.10 markdown deliberately explains the rung's **negative** (deriving the 6.5 °C/km lapse rate
  reproduces it at mid-latitudes and loses to it in the tropics) in plain words — the notebook is where the
  "what did NOT work" half is most valuable, and it needs no jargon to land.


**Extended again 2026-09-04 — the pedagogy pass (markdown-only, NO re-bank).** Six markdown cells added
(five "🔮 Predict first" boxes before the live sliders in §§1–5, one "Four missions" cell at the end of §7):
buckets A and B of [[pedagogy-novice-intermediate]]. The cheap-lane recipe, for the next markdown-only edit:
- **No code cells ⇒ no output banking and no `docs/index.html` regeneration** — the two most expensive steps
  of every previous §8 extension simply do not apply. `git diff --stat` came back **+222/−0** on the notebook.
- **Insert by CONTENT, never by index** (advisor-caught): every insertion shifts the later cells, so the
  patch script iterates the cell list and matches each anchor's lead text (the `# Live:` comment, the `## 8 —`
  heading), then asserts the final count is exactly `old + 6` and that the §8 / provenance headings still
  parse. Script kept at `M:\claud_projects	emp\pedagogy-ab\patch_notebook.py`.
- The `json.dumps(nb, indent=1, ensure_ascii=False)` byte-identity check is now an **assert at the top of the
  patch script**, not a manual pre-step; the no-trailing-newline rule still holds.
- **Cells 54–56 in HEAD (the §8.9 block) carry no `id` field** — `nbformat.validate` warns about it. Pre-existing,
  not introduced here; left alone rather than churning the diff.
- `test_planet_notebook.py` re-run **solo** on the final state: `pytest planet/tests/test_planet_notebook.py -n0`
  (the addopts carry `-n auto`, so `-p no:xdist` is a *parse error* — `-n0` is the in-process form).
