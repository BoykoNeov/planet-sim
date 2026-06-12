---
name: planet-viz-animation-rungs
description: "Animated eddy-flow visualization, 3 rungs A/B/C: rung A (matplotlib two-panel) + shared (h,u,v,θ) frame side-channel BUILT 2026-06-11; rung B (Plotly-globe anim, planet/eddy_globe.py) BUILT 2026-06-12; rung C build-approach DECIDED 2026-06-12, AMENDED same day to a three.js/WebGL true-3D-sphere particle flow-globe (was original Canvas2D) + honest-by-disclosure carve-out (NOT built); honesty edges = channel-not-globe + ~90% reversible"
metadata: 
  node_type: memory
  type: project
  originSessionId: 8df81636-13d1-4e22-a947-a3050d3b8286
---

**Decided 2026-06-11 (user): build all three visualization rungs A→B→C** to animate the
emergent eddy life cycle (`eddy_flux.eddy_life_cycle` — the only time-varying,
longitudinally-structured 2-D flow the project produces). Recorded in plan **§9.5** + the
§10 running log. **Rung A + the shared frame side-channel BUILT 2026-06-11; Rung B BUILT 2026-06-12;
Rung C build-approach DECIDED 2026-06-12 (plan-only, NOT built — see the Rung C bullet below).** These are
**visualization** rungs A/B/C — do NOT conflate with the §5 GCM staircase rungs 0–6.

**BUILT (2026-06-11) — rung A + the shared prerequisite:** `n_frames` side-channel on
`eddy_flux.eddy_life_cycle` (opt-in `EddyFrames` snapshotting `(h,u,v,θ)` + cumulative
transport traces over the full release), **diagnostic-pure** (`n_frames=0` bit-for-bit, test
asserts `==` on `kappa_bulk`/`F_int`/`G_int`/`jet_speed`/`eddy_ke`); `plots.eddy_life_animation`
(two-panel `FuncAnimation` — θ-stir+eddy-quiver beside the throughput-vs-net transport budget,
fixed colour range + fixed quiver scale); `demo_eddy_life.py` banks `docs/figures/planet-eddy-life.gif`
(Pillow, `dpi=90`/48 frames to bound size). **Advisor caught the one design error:** the net trace
must integrate the flux **per latitude FIRST** then `Σ|·|` — `mean_interior(signed F̄)` leaks
*spatial* cancellation into a *temporal*-reversibility curve (makes flux look MORE reversible than
it is = overclaim). Per-latitude-first makes the endpoint ratio land on `irreversible_fraction` by
construction; a marked `window_start` line reconciles full-release curve vs windowed banked number.
Tests green (3 new + fast lane 208). §9.4 rule-of-three NOT yet acted on (still project-local in
`plots.py`; promotion-to-`viz/` is a future call). **Prose pass 2026-06-12:** the GIF got the
novice→intermediate relabel+caption sweep — de-jargoned legend (`back-and-forth stirring (throughput)` /
`net heat moved poleward (κ residual)`, the old "swirls rage"/"small residual" + raw formulas gone from
the legend), a plain-language `fig.text` caption (β-plane band / reversible-% / κ residual, with the
formulas restored as **mathtext** `$\bar{F}$` and glossed in words), the constrained-layout engine
`rect=(0,0.17,1,0.83)` reserving the bottom band, %s off `irreversible_fraction` (banked GIF irr≈0.081 →
92/8, agrees with the globe); pinned by a caption assertion. This **completes the owed sweep** (planetmap
biome + coupler globes + the A·B·Δ triptych also done same day) — see [[viz-prose-novice-intermediate]].

**Why:** user wants the NASA *Perpetual Ocean* / Ventusky "flowing" look (named as *broad*
references — weather, not climate). The rungs rise in cost, *fall* in pedagogical return,
*rise* in overclaim risk, so each is deliberate.

**How to apply:**
- **Shared prerequisite for ALL three:** bank `(h,u,v,θ)` frames from the `eddy_flux`
  release loop, which currently discards them. Add an opt-in `n_frames` snapshotting into an
  `EddyFrames` side-channel — include `h` (verification anchors need it), snapshot on even
  *time* thresholds over the **full** release `[0,t_end]` (`dt` is adaptive; full span = the
  growth→saturation mechanism, NOT just the κ window), and keep it **diagnostic-pure**
  (`n_frames=0` ⇒ κ result bit-for-bit unchanged — the inert-seam discipline).
- **Rung A** = matplotlib `FuncAnimation`, in-repo, the repo's FIRST time-animation primitive
  (`plots.eddy_life_animation` + `demo_eddy_life.py`) → finally trips §9.4's rule-of-three
  promotion-to-`viz/` trigger (steel/chip + this = 3rd consumer). **Two-panel** (θ swirls
  beside running cumulative `∫F̄dt` vs `|F̄|` throughput) so the **~90%-reversible** flux
  finding (`irreversible_fraction ~0.1`) is *visible* — a plain stirring movie would
  CONTRADICT the module's headline. GIF/Pillow = CI-safe default, MP4/ffmpeg optional.
  Frame-fidelity tests: `∫hθ` machine-exact, `eddy_ke` recomputed-from-frame matches series,
  `n_frames=0`-vs-N bit-for-bit. **Build A first** (honest by construction), judge B-vs-C
  after seeing real frames move.
- **Rung B** = animate the existing Plotly globe — **BUILT 2026-06-12** (`planet/eddy_globe.py`:
  `eddy_globe_figure` + `save_eddy_globe_html`; `demo_eddy_globe.py` banks
  `docs/figures/planet-eddy-globe.html` ~3.8 MB; `test_eddy_globe.py` 4 tests). **No new stack**
  (reuses `[webviz]` Plotly + `planetmap._sphere_xyz`; plotly-free at import — lazy, subprocess
  headless guard). **Both honesty edges carried GEOMETRICALLY:** the doubly-periodic channel maps to a
  **bounded ~55° longitude sector at TRUE physical width** `Δlon=Lx/(a·cosφ_c)` (Earth radius `a`
  recovered from the frames' own linear β-plane `(y,φ)` metric — no new constant/schema field), a
  **single NH band** (measured 18.9°–61.1° lat × 55.1° lon) on a bare base sphere — **NOT** 360°-wrapped,
  **NOT** mirrored to the SH (the `circulation_layer` two-hemisphere broadcast is valid only for a
  *zonal-mean* jet; the eddy field is the project's ONLY longitudinally-structured field → does NOT
  transfer — this was my design error the advisor overturned pre-build); + the Rung-A flux-budget panel
  beside the globe (throughput-rages/net-small) keeps ~90%-reversible on screen, cursor + `κ diagnosed`
  line tie the two. **Advisor's pre-push gate:** "builds ≠ animates" — frames update only the band's
  `surfacecolor` (lean HTML, `x/y/z` static-merged), a known Plotly-3D soft spot → `cmin/cmax` re-stated
  per frame (no colour-autoscale) + `redraw=True` on play/slider (force gl3d repaint) + a structural
  test pinning each frame's `surfacecolor`/`traces`/that-the-band-changes. Browser play-through is the
  ONE thing not headlessly self-verifiable (no kaleido) → handed to the user to eyeball.
  **Post-eyeball polish (2026-06-12):** user confirmed the globe animates (a/b/c) but flagged the
  flux panel + overlapping text. Fixed: the "κ diagnosed" annotation's hardcoded `xref="x2"/yref="y2"`
  spawned a **phantom overlaid axis** (the xy subplot is on `x`/`y` — a 3-D scene consumes no cartesian
  axis number), producing the double x-axis + collapsed ~33.2M y-range — re-anchored via `row/col`;
  dropped the redundant "Throughput rages" annotation; **normalized the panel to a fraction of total
  throughput** (throughput→1.0, net plateaus at ~`irr`) so it reads 0..1 not ~3e7 K·m; dropped the
  floating `t=N` marker; reseated the legend. Also removed hover spike lines per the new standing
  preference → [[no-pointer-spike-lines]]. **Second polish pass (2026-06-12):** legend moved to the
  panel's upper-RIGHT corner (still overlapped the rising curves on the left; near-opaque bg masks the
  throughput endpoint it now sits on); + **richer novice→intermediate prose** = relabel-don't-append
  (de-jargoned legend/subtitle, technical term parenthetical) + ONE plain-language caption (β-plane band
  / "reversible" + the % / κ residual) seated below the controls — the Rung-B figure is the **exemplar**
  for the new standing rule [[viz-prose-novice-intermediate]] (sweep of planetmap globes + triptych +
  GIFs still owed). **Third pass (2026-06-12):** user confirmed legend OK but caption overlapped the
  slider's "… periods" labels + wanted it bigger + the formulas BACK (explained, not stripped) → caption
  reseated low (`y=-0.40`, `b=360`/`height=860`), font 11→14, formulas `Σ∫|F̄|dt`/`Σ|∫F̄dt|` restored
  with word-glosses + curve-coloured `<span>` tints — these refinements folded into
  [[viz-prose-novice-intermediate]].
- **Rung C** = the showcase — **build approach DECIDED 2026-06-12 (user), AMENDED the same day to a true
  3-D sphere, PLAN-ONLY (not built this session).** User reframed it as a **general-purpose
  flow-on-a-globe renderer** aimed at *one day* visualizing a full **GCM/ESM** field; renders lesser models
  (today: the one eddy band) meanwhile. **First locked as an *original Canvas2D orthographic* globe, then
  amended (user) to a `three.js` / WebGL *perspective* globe — particles streaming on a REAL rotatable
  sphere (PerspectiveCamera + orbit control, back-face occlusion).** LOCKED: (1) **Tech = three.js / WebGL
  3-D sphere; particle technique stays ORIGINAL** (our own advection updating a `BufferGeometry`, CPU-side
  to start); three.js is the only **vendored** piece — the scene/camera framework, not a particle lib we
  copy. **§6 REVERSES from sidestepped → owed:** three.js is a *library* (vendored, not reimplemented) so
  its MIT licence **requires a `NOTICE`/attribution file** (currently absent → a NAMED build deliverable);
  **vendor three.js INLINE, not CDN** (CDN rejected — breaks offline self-containment; inline keeps the
  `interactive.py` file:// / deterministic / golden-able property). Original Canvas2D demotes to a lighter
  §6-free **fallback**. (2) **Renderer-agnostic data contract UNCHANGED** = a generic
  **vector-field-on-a-globe** layer: grid + per-cell `(u,v)` + optional scalar (colour) + optional frames +
  **coverage-extent** + **provenance/honesty label** — NOTHING about projection/particles (renderer-side).
  **This is the WIN: the contract was built to absorb exactly this swap → the amendment touches ONLY the
  renderer** (contract, §9.3, honesty carve-out, ADR note all unchanged = proof the boundary was right).
  **Coverage-extent carries the band-vs-globe truth into the data** (label-the-band today,
  illustrate-the-globe-with-disclosure tomorrow). Joins the §9.1 layer registry / planet-spec as a new
  layer TYPE. (3) **GPU-advection seam SHIFTED DOWN a level** — WebGL is the v1 renderer now (not the future
  upgrade), so the seam is CPU `BufferGeometry` advection v1 → **GPU ping-pong** (webgl-wind technique, frag
  shader) upgrade *within* WebGL; same **trigger:** frame-rate < ~30 fps at GCM resolution. (4) The one
  **policy change = honest-by-disclosure** (see below), UNCHANGED by the amendment. **B vs C now both 3-D
  globes → ROLES separate them** (not dimensionality): B = faithful **scalar field** on the one true band
  (honest-by-construction); C = **immersive particle-streaming** showcase (honest-by-disclosure).
  Deliverables (future build): `planet/flow_globe.py` (generic renderer, three.js vendored inline) +
  `demo_eddy_particles.py` → `docs/figures/planet-eddy-particles.html` + a **`NOTICE`/attribution file** +
  catalog entry + drift-guarded site regen + structural & **disclaimer-presence** tests (disclaimer now a
  **DOM overlay** over the WebGL canvas → easier to machine-check than canvas-drawn text; mirror
  `test_eddy_globe.py` minus byte-golden, plus the caption assertion). Doctrine recorded in plan §9.5 +
  §9.3 + ADR 0002 status note.

**Two honesty edges carried through all three (hardest to preserve at C):** (1) the flow is a
**doubly-periodic midlatitude β-plane band patch, NOT a global field** (same edge as
`planetmap.circulation_layer`) → on a globe it is one honest latitude band; (2) the
instantaneous flux is **~90% reversible** → particles stream but mostly slosh, net transport
is only the small κ residual. The prettier the rung, the more the medium implies global net
transport the model lacks. **Edges carried two WAYS (the 2026-06-12 carve-out):** A/B are
**honest-by-construction** — the geometry cannot lie (B keeps the flux indicator + a true band label,
never a 360° wrap); **C, the showcase, is honest-by-disclosure** — it MAY illustrate beyond the model
(global-looking flow, "currents carry heat" though flux is reversible) **IF a VISIBLE on-screen
disclaimer documents the departure** (user 2026-06-12: illustrate freely *"if documented"*). The carve-out
is **narrow + asymmetric**: physics-fidelity verification RELAXES for C only (approximate OK, no
byte-golden — a figure was never in the correctness path, ADR 0002 #2), but **documentation verification
TIGHTENS** — a test machine-checks the disclaimer is present (the disclaimer IS the entire license).
Science layer + A/B untouched. Verification (answers the user's screenshot concern): rung A is a
fixed-camera flat field (frames directly comparable) + numerical proof (`∫hθ` machine-exact,
θ bounded, `eddy_ke` saturates, cumulative flux → diagnosed `kappa_bulk`), not the eye.

Related: [[planet-rung1-two-way-coupler]] (eddy_flux source of the frames),
[[planet-interactive-map-design]] (the Plotly globe rung B reuses + layer registry),
[[planet-rung2-scoped]].
