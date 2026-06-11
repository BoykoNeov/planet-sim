---
name: planet-viz-animation-rungs
description: "Animated eddy-flow visualization decided as 3 rungs A/B/C (all to be built, not yet); shared prerequisite = bank (h,u,v,θ) frames from eddy_flux; honesty edges = channel-not-globe + ~90% reversible"
metadata: 
  node_type: memory
  type: project
  originSessionId: 8df81636-13d1-4e22-a947-a3050d3b8286
---

**Decided 2026-06-11 (user): build all three visualization rungs A→B→C** to animate the
emergent eddy life cycle (`eddy_flux.eddy_life_cycle` — the only time-varying,
longitudinally-structured 2-D flow the project produces). Recorded in plan **§9.5** + the
§10 running log. **NOT built yet** — this is the decision + the design, not code. These are
**visualization** rungs A/B/C — do NOT conflate with the §5 GCM staircase rungs 0–6.

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
- **Rung B** = animate the existing Plotly globe (`planetmap`) — globe view, no new stack.
- **Rung C** = WebGL particle globe; vendor **mapbox/webgl-wind (ISC, GPU particles)** +
  **cambecc/earth (MIT, orthographic projection)** (licenses verified). webgl-wind is
  flat/equirectangular → needs a three.js sphere or a port of cambecc's projection. New
  JS/WebGL stack; reach/delivery not new teaching; §6 attribution diligence applies.

**Two honesty edges carried through all three (hardest to preserve at C):** (1) the flow is a
**doubly-periodic midlatitude β-plane band patch, NOT a global field** (same edge as
`planetmap.circulation_layer`) → on a globe it is one honest latitude band; (2) the
instantaneous flux is **~90% reversible** → particles stream but mostly slosh, net transport
is only the small κ residual. The prettier the rung, the more the medium implies global net
transport the model lacks → **B/C must keep a flux indicator + band label** (§9.3 inert-honesty
applied to motion). Verification (answers the user's screenshot concern): rung A is a
fixed-camera flat field (frames directly comparable) + numerical proof (`∫hθ` machine-exact,
θ bounded, `eddy_ke` saturates, cumulative flux → diagnosed `kappa_bulk`), not the eye.

Related: [[planet-rung1-two-way-coupler]] (eddy_flux source of the frames),
[[planet-interactive-map-design]] (the Plotly globe rung B reuses + layer registry),
[[planet-rung2-scoped]].
