---
name: planet-rung5a-orographic
description: "Rung 5A (first step off the zonal mean) BUILT 2026-07-09 (planet/orographic.py = Smith & Barstad linear orographic precip) + 5A.2 integration BUILT 2026-07-10 (planet/orographic_scene.py: tangent-plane patch + jet-sourced wind + mm/hr→cm/yr + enhancement-only biome re-map + serialization + demo figure) + 5A.3 lee-depletion BUILT 2026-07-10 (planet/orographic_depletion.py: opt-in along-wind moisture budget draws the lee BELOW baseline = the real desert; drying-ratio-calibrated) — PAYOFF: the mountain changes the biome map AND casts a real rain-shadow desert; + 5A.4 elevation->temperature BUILT 2026-09-04 (planet/elevation_temperature.py, build_scene(lapse=True)): the crest turns ALPINE TUNDRA (16 K colder, 42% of the patch re-classified by the cooling alone) — and the bet to make the 6.5 K/km lapse rate EMERGENT came back a NEGATIVE (the moist adiabat reproduces the constant at midlatitudes and loses to it on the tropical freezing level), so the constant stays the default"
metadata:
  node_type: memory
  type: project
  originSessionId: fdb095c9-cdad-43a3-ac05-2e4a17bc3da3
---

**Rung 5A = the first laptop-scale step toward the "north star"** (true 2-D longitudinal geography →
regional climate, orographic precip, rain shadows; plan §12.5 / §5). Literal rung 5 (a full idealized
GCM) is infeasible-tier, so it is climbed as a **spike-first sub-ladder of reduced models** — the same
decomposition the gas-giant sketch uses ([[gas-giant-feasibility]]). Continues [[planet-plan]];
follows the diagnostic-precip precedent [[planet-phase2]] one rung out.

**BUILT 2026-07-09** — `planet/orographic.py` + `planet/tests/test_orographic.py` (10 anchors green,
fast lane). **Smith & Barstad (2004) linear theory of orographic precipitation:** a *diagnostic* on a
**prescribed uniform** background wind over a **2-D terrain** → a wavenumber-space transfer function
(one FFT, laptop-trivial) → windward rain + a lee **rain shadow**. Numerics + constants pinned in
[[smith-barstad-orographic-source]].

**What it banks (the payoff):**
- **Wakes the dormant elevation seam.** Since v1 `planet_spec` carried an `elevation` layer tagged
  `inert=True` (carried/displayed/round-tripped, climate-inert; §9.3). This is the module the plan
  named as the one that finally makes elevation *do something*.
- **Tight anchor = convergence to the closed-form triangle-ridge solution** — but honestly scoped
  (advisor-caught): in that limit `(1 − i·m·H_w) → 1`, so `m` and its `sgn(σ)` branch DROP OUT — the
  triangle anchor pins only the *reduced* transfer function (`C_w` scaling, upslope `iσ`, fallout
  `τ_f`), NOT the branch or `H_w`. **The branch is guarded by exactly one test — the rain-shadow
  direction** (windward wetter than lee, peak upwind); empirically, flipping the branch sign reddens
  *only* that test (shadow flips to the windward side). The wind-reversal **mirror** (machine precision,
  1.8e-15) is a *reflection self-consistency* check — invariant under a branch flip, so NOT a branch
  test. Plus: upslope limit `C_w·max(0,U·∇h)` to the analytic Gaussian gradient; flat-ground null.

**The honesty flag (a TRADE, not a win — do not overclaim):** it makes the **precipitation** 2-D, NOT
the engine. The **temperature** climate underneath stays **zonal-mean** (the EBM). So Rung 5A does
*not* claim the engine has left the zonal mean — only the rain shadow has. Regional Cartesian patch,
uniform prescribed wind; cross-mountain flow is **prescribed, not emergent**.

**Advisor-caught traps (all now handled):** (1) the upslope limit is `H_w=τ_c=τ_f=0`, NOT `U→∞`
(large U drives precip to *zero* via the τ factors); (2) the `sgn(σ)` branch is the #1 bug — its ONLY
validator is the **rain-shadow direction** test (NOT the triangle anchor, where the branch drops out;
NOT wind-reversal, which is branch-flip-invariant); (3) zero the `σ=0` locus explicitly; (4) zero-pad
the FFT so lee drying doesn't wrap. Also: compare the **truncated** model to the exact solution
(untruncated has a −C anti-rain spike at the downwind kink).

**5A.2 — BUILT 2026-07-10** (`planet/orographic_scene.py` + `test_orographic_scene.py` (12 anchors) +
`plots.orographic_scene_figure`). All five named pieces shipped: **sphere placement** (`patch_spacings`
= tangent-plane metric `dx=R·cosφ·Δλ`, `dy=R·Δφ`); **cross-mountain wind from the emergent zonal jet**
(`wind_from_jet` — the westerly sampled at the patch latitude via [[planet-phase4-coupler]], **prescribed
not emergent** = the Rung-5A caveat, zero outside the westerly band); **mm/hr→cm/yr** via a named
loose-magnitude knob `OROGRAPHIC_HOURS_PER_YEAR=500` (an *effective annual uplift duration* — NOT a naive
×8766-h annualisation, which overstates ~10× and swamps the classifier); **enhancement-only combination**
`P_total = zonal_baseline + orographic_bonus` + biome re-classification ([[planet-phase2]]); **serialization
for free** (the regional scene is a `Grid`+`Layer` stack → rides the grid-agnostic [[planet-interactive-map-design]]
schema, round-trip green). **THE PAYOFF: the mountain finally changes the biome map** (~40% of a
Cascades-scale patch re-classified, windward → temperate rain forest — elevation is no longer inert).

**Advisor-caught honesty flags (all handled):** (1) S&B is *enhancement* physics — windward wet, *dry
immediate lee at* baseline; it does NOT model **background depletion** (windward rainout drying the lee
*below* baseline = the real Columbia-Basin desert) → named+deferred as a future **5A.3 moisture budget**.
(2) The full model has a weak **downstream secondary rain band** from the **propagating-mode phase** —
VERIFIED real to the model (not FFT wrap / pad artifact) by a discriminator: it vanishes at `H_w=0` AND
holds its physical downwind distance under domain-doubling (locked as a test); it is *not* a trapped lee
wave. ≥0 so the enhancement-only invariant `P_total≥baseline` holds; the reduced-limit "clean decay to
baseline" is only the immediate lee. (3) The "zonal ridge casts no shadow" idealisation **fails on the finite zero-padded
engine** (the pad localises a lon-uniform ridge into a responding block) → the pad-safe placement anchor
is a **compact hill's latitude-symmetry** (shadow in lon, ~symmetric in lat → pins `lon→x`), not a zonal
ridge's amplitude. (4) Coarse `N_LON=73` globe (~400 km/cell) would smear the ~20–100 km shadow to nothing
→ the scene is a **fine regional patch** (Δ in tens of km); the demo is a regional inset, not the globe.
(5) Demo range must sit under the **annual-mean zonal westerlies** (Cascades/Andes/NZ); the Western Ghats
(monsoon-driven) can't be driven by the zonal jet.

**5A.3 — BUILT 2026-07-10** (`planet/orographic_depletion.py` + `test_orographic_depletion.py` (9
anchors) + `build_scene(..., deplete=True)` + a depletion-aware demo figure). The one piece 5A.2 named
and deferred: the real Columbia-Basin desert is **not** "no orographic bonus" but a lee baseline drawn
**below** the zonal mean — the windward rainout drains the passing air — which enhancement-only 5A.2
structurally cannot make. Built as an **opt-in** (default `deplete=False` keeps 5A.2 exactly; the
`moist.py` opt-in precedent).

**The model = a 1-D along-wind moisture-flux budget** (the wind is zonal → each lat row is an independent
streamline): `d(U·W)/dx = P_base·(1−g) − P_oro`, for the dimensionless depletion factor `g = W/W₀`.
Combination `P_total = g·P_base + P_oro` → the lee drops below baseline where `g < 1`.

**Advisor's load-bearing correction (do not re-litigate):** the **refill** term `P_base·(1−g)` and the
**depletion** term `−P_oro` are a *forced package* from the **same** premise (`S = P_base`, evaporation
resupplies at baseline) — you cannot keep bonus-only depletion and separately "choose" to drop refill.
The no-refill form `g = 1 − (1/U·W₀)∫P_oro` (cumulative, DOWNWIND) is the budget's **L→∞ limit**, honest
**because** the *derived* evaporative-refill length `L = U·W₀/P_base` (~16 000 km) ≫ the ~450 km patch
(verified numerically + a test). Frame it as the large-L limit, NOT "we chose no refill." The converse is
the honest cost: the modelled desert does not relax back within the window (a real one does, over ~L).

**Two traps the advisor flagged, both handled:** (1) **unit split** — the budget runs on the
*instantaneous* mm/hr S&B rate (physical U, physical `PWV_IN_MM`) → a dimensionless `g` → applied to the
*annual* cm/yr baseline; `OROGRAPHIC_HOURS_PER_YEAR` must NOT enter the budget. (2) **integration
direction = the new `sgn(σ)`** — the cumulative integral runs downwind so the drain lands in the **lee,
not windward**; the sole guard is a "depletion-in-the-lee" test (reverse the wind → the desert flips
sides), the mirror of Rung 5A's rain-shadow-direction test.

**Anchors:** tight/exact = conservation (water off the flux = orographic water rained → `DR = 1 − g_lee`
identically) + reduction (`g ≡ 1` recovers 5A.2 bit-for-bit); tight/structural = monotone-`g`-downwind +
depletion-in-lee; **directional payoff** = lee total below baseline; loose = the incoming column
`PWV_IN_MM ≈ 30 mm` calibrated so the demo **drying ratio** `DR ≈ 0.47` sits in the cited observed band
(~0.3–0.5, Sierra/Cascades 0.48±0.02; Roe 2005, Smith & Evans 2007 — pinned in
[[smith-barstad-orographic-source]]; first-pass cites were fabricated, advisor+WebSearch-corrected).
**Payoff numbers:** Cascades lee ~90 → ~55 cm/yr, ~⅓ of the patch turns lee-desert, reclassified
41 % → 56 % (windward temperate rain forest, lee woodland/shrubland). **Honest scope (named, not fixed):**
per-streamline (zonal wind only); the S&B bonus is not itself depleted (simplest-first); no on-patch
refill (the L→∞ cost).

**5A.4 — BUILT 2026-09-04** (`planet/elevation_temperature.py` + `test_elevation_temperature.py`
(15 anchors) + `build_scene(..., lapse=True)` + 7 new scene cases + `demo_alpine_biomes.py` +
`plots.alpine_biomes_figure` + catalogue `alpine` **and** `orographic` — the 5A demo had never been
catalogued, a real drift against the ARCHITECTURE §4.5 recipe, fixed here under a new "Geography —
mountains" section). The last unstruck line of §12.5's *cheap tier*: 5A/5A.2/5A.3 woke the elevation seam
on the **rain** side only — the temperature stayed the zonal mean broadcast across longitude, so a 2500 m
crest and its valley got the **same** number. Now each cell is cooled by its own terrain height before the
Whittaker classifier runs. **Diagnostic and one-way** (no EBM re-solve, no snow/albedo feedback on the cold
crest, no S&B re-run; `C_w` stays at its upstream sea-level value).

**THE PAYOFF:** the Cascades crest cools **16.2 K** (6.58 → −9.67 °C), crosses both Whittaker cold
thresholds and turns **woodland/shrubland → tundra**; the cooling **alone** (measured against the *same*
rainfall at the uncooled temperature — the new `sea_level_biome_codes` control) re-classifies **42 %** of
the patch; total vs the zonal-mean map 56 % (rain only) → 60 % (rain + cooling).

**THE BET → A NEGATIVE (the more interesting half; the rung-4 "retire the prescribed number" habit fails
here).** The target was `radiation.LAPSE_RATE = 6.5 K/km`, to be replaced by rung 4's **own** `Γ_m(T,p)`
integrated up a hydrostatic column. It loses twice: (i) at the demo latitude the emergent effective rate is
**6.31 K/km** — it *reproduces* the constant (~3 %) instead of retiring it, revealing the textbook number as
a **mid-latitude calibration**; (ii) on the one observation available — the **freezing level** (0 °C isotherm,
the deep tropics' being the planet's highest at **≈ 5 km**; Harris, Bowman & Shin 2000, pinned in
[[smith-barstad-orographic-source]]) — the constant lands **4.38 km, just BELOW the band** and the moist
adiabat **7.09 km, ~45 % above**, because a saturated **parcel** adiabat is not the **environmental** lapse
rate of an unsaturated mean column. **Say "just below", not "in band"** — the verdict is the *ordering*
(one close, one far), and the tests assert exactly that. Its predicted latitude contrast (same mountain
cooling ~9 K in the tropics, ~22 K at 85°) is **anti-correlated with reality at both ends** (polar surface
inversions: >95 % of Eurasian-Arctic winter soundings, Serreze, Kahl & Schnell 1992). ⟹ `moist=False` is the **default**,
the emergent path is **opt-in as a diagnostic**, the contrast is reported as the idealisation's and **not
banked as the planet's** — the [[planet-rung5b4-seasonal-sici]] pattern (an axis rejected after being tried).

**Anchors:** tight/exact = the closed form `T − Γz`; **the integrator reduces to it exactly under a constant
`Γ` callable** (so emergent and prescribed are the *same code*, not two implementations); `Γ_m` is rung 4's
function **by identity** ([[planet-rung4-radiation]]); two exact nulls on the scene (default-off, and
flat terrain with the correction on). Tight/convergent = Heun march **second order in the height step**
(ratios ~4.0; shipped `n_steps=64` good to 9e-5 K over 3 km) + the small-`z` limit `T₀ − Γ_m(T₀,p₀)·z` with
an O(z²) residual. **No conservation leg exists** (cooling a surface without re-solving the EBM breaks the
TOA budget — the [[planet-phase2]] classifier precedent, stated not papered over); the substitute is the
exact identity `⟨T_sea − T⟩ = Γ·⟨z⟩` + the partition still tiling.

**Traps + named scope:** (1) **the z→p circularity** the advisor flagged — `radiation.moist_adiabat_temperature`
maps height to pressure through a *fixed* isothermal scale height, so integrating there would leave the
"emergent" rate leaning on a prescribed constant; this module integrates pressure **hydrostatically in the
marched `T`** instead, so the circularity is *absent*, not merely small. (2) **The classifier's cold bands
are precipitation-independent**, so on the crest (wet *and* cold) the rain shadow can no longer change the
answer — 5A and 5A.4 genuinely **degenerate** there; that is why the total only moves 56 → 60 %, and it is
pinned as a test rather than rediscovered later. (3) Terrain above 12 km is **refused, not clamped** (a
tropospheric march, no cold-trap floor). (4) The EBM temperature is treated as a **sea-level** temperature —
exact for this model (rung 0 has no terrain), approximate for Earth.

Then Rung 5B/5C = the 2-D EBM engine step (needs seasonality for continentality) / Charney–Eliassen
stationary waves.
