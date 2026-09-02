---
name: planet-bifurcation-sici
description: Rung 0+ BUILT 2026-09-02 — the COMPLETE equilibrium diagram of the ice-albedo EBM (`bifurcation.py`, the inverse solve, dt-free) — every branch incl. unstable, BOTH folds (Snowball + the small-ice-cap instability θ_c≈10°), the slope-stability theorem checked by marching, Legendre-mode anchor at 2nd order, the relax O(Δt) bias quantified + retired, present day 0.18 % below the SICI cliff; SICI no longer deferred
metadata:
  type: project
---

**Rung 0+ — the complete equilibrium diagram + the small-ice-cap instability · BUILT 2026-09-02**
(`planet/bifurcation.py`, `test_bifurcation.py` (16 fast + 1 slow), `demo_bifurcation.py` +
`test_demo_bifurcation.py`, `plots.bifurcation_figure` → `docs/figures/planet-bifurcation.png`; catalogue
key `bifurcation`; notebook §8.8). Closes the **"SICI named-deferred"** line of [[planet-rung5b-seasonal]]
(5B.1+) — and closes it on rung 0, where it belongs (the annual-mean EBM is where North found it).

**The method — the INVERSE problem, and it is exact.** Phase 1's continuation sweep ([[planet-phase1]])
can only sit on *stable* climates, and every point rides the Strang-split relaxation whose fixed point has
an **O(Δt) bias in the profile shape** (backward-Euler transport substep) — the reason `demo_snowball` runs
`PRESENT_N_TAU=0.01`. Instead: *prescribe the ice line `x_s`* → the step albedo is a known field → the EBM
is **linear** → one tridiagonal solve on the engine-pinned operator (`EnergyBalanceModel._transport_tridiag`)
gives `u(x)` per unit S₀ → `T = S₀u − A/B`, `T(x_s)=T_f` ⟹ **`S₀(x_s) = (T_f + A/B)/u(x_s)`**. Sweep `x_s`
over every grid face = **every equilibrium with an ice line, stable AND unstable**, as one curve (North
1975's S-curve, off the *same* discrete operator the marcher uses). Ice-free (`x_s=1`) and Snowball
(`x_s=0`) branches cap the ends. 0.02 s at 180 cells.

**Stability = the slope (Cahalan & North 1979 slope-stability theorem): stable ⟺ `dS₀/dx_s > 0`.** Read
off the curve, then **CHECKED by marching** (`relax_from_curve`): ±1 K nudge on a stable segment returns
(<1.5°), on the unstable mid-branch departs to Snowball (cold) / a warm branch (warm), on the unstable
small cap → ice-free (warm) / the stable finite cap (cold). Two traps found + pinned: (1) **critical
slowing near a fold** — the per-step `tol` on the relaxation declares convergence early, so use
`tol≤1e-11`, big `max_iter`, and read the ice line away from the folds; (2) the marcher's ice edge is
**cell-quantized** (a whole cell flips albedo) vs the curve's face-interpolated edge → they differ by ≤ one
cell's latitude width and the gap **halves per grid doubling** (−1.4/−0.75/−0.41/−0.05° at n=60/120/240/480)
— a Δx effect, NOT Δt. So the "relax → exact as dt→0" reduction is stated to within a cell.

**Anchors banked (tight):** FV curve → North's **even-Legendre-mode** solution (`legendre_equilibrium_curve`,
`T=ΣTₙPₙ`, transport diagonal `−n(n+1)D`, albedo step by exact piecewise Gauss–Legendre; North 1975 = the
`n_modes=2` truncation, ~1 % off) at **2nd order** (exact face 7e-4→1.9e-4→4.7e-5 rel. at n=45/90/180;
harmonic 6e-3→1e-3→2e-4); **slope theorem by marching** (above); the Phase-1 **sweep jumps within one
sweep step of the exact folds** (freeze 1252 vs fold 1259; melt 1831 vs branch-end 1835; step 15) and every
finite-cap point it visits has a stable curve twin; net TOA exact on every point (conservation);
precedence trap caught: `a * b @ P` parses as `(a*b)@P` — parenthesize the matmul.

**The findings (loose/calibrated on `A,B,D,α,T_f`; the structure is the bank):**
- **TWO folds at Earth params:** the Snowball fold at **33.0°, S₀=1259 (−7.8 %)** — Phase 1's freeze — and
  the **small-ice-cap fold at 79.1° (θ_c = 10.9°, converged on a 720-cell grid; 10.3° on 180), S₀=1367.3
  (+0.16 %)**. The **finite-cap window** (the only suns holding a polar cap) is **1259…1367 W/m²**, ~8 %
  wide, and **present day sits 2.1–2.5 W/m² (0.16–0.18 %) below its top** — a brightening of that size
  loses the cap in a jump (the SICI, North 1984).
- **Five equilibria at today's sun:** ice-free (stable, S₀ ≥ 1359 — coexists with the finite cap between
  1359 and 1367), unstable small cap (84.8°), **stable finite cap 74.8°** (Earth's; benchmark ~70°),
  unstable mid-branch separator (13.8°), Snowball (stable, S₀ ≤ 1835).
- **θ_c vs D** (`critical_cap_sweep`, 720 cells — the uniform-`x` grid is ~2°/cell at 80° on 180 cells, the
  fold converges only to ~1° there; 720 → 0.1°): **≈10° floor for weak transport** (10.06 / 10.13 at D=0.2 /
  0.3, grid-converged, both face modes), **growing past D≈0.4** (10.9 at 0.555, 12.2 at 0.8, 14.0 at 1.0,
  18.7 at 1.3) while the **window narrows** (448 → 109 → 16 → 1 W/m²) and **vanishes at D≈1.4** (folds
  merge: an efficiently-mixed planet is ice-free or Snowball only — North). NOT a `√(D/B)` law: the ratio
  θ_c/√(D/B) drifts 0.79→0.34; report numbers, not the scaling.
- **The retired bias:** relaxed present-day ice line 64.8 / 69.3 / 71.2 / 72.8 / 73.3° at n_tau 0.5 / 0.2 /
  0.1 / 0.05 / 0.02 vs exact 74.8° — the O(Δt) bias, then the ~1° cell floor. The exact diagram needs no
  step. (`ebm.py`/`albedo.py` untouched — a sibling, ADR 0005.)

**Named edges:** uniform scalar `D` only for the Legendre anchor (a callable `D(x)` is not mode-diagonal —
the FV curve accepts it); the theorem is for this EBM class (diffusive, step albedo, linear OLR); the
seasonal SICI (Huang & Bowman 1992 — does the seasonal marcher's cap vanish by a jump?) is the natural
next sweep on [[planet-rung5b-seasonal]]'s 5B.1+ marcher, named not built. Sources: North 1975 JAS 32;
Cahalan & North 1979 JAS 36 1178; North 1984 JAS 41 3390. → [[ebm-radiation-source]], [[planet-phase1]].
