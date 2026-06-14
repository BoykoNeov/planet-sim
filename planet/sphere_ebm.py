"""The pole-to-pole energy-balance model and the **energetic ITCZ** (Planet rung 2.x).

Rung 0's EBM (:mod:`planet.ebm`) is solved on a **single hemisphere** — ``x = sin φ ∈ [0, 1]`` with a
hemispheric-**symmetry** boundary at the equator (``ebm.py`` line: *"the equator (x=0) end by hemispheric
symmetry"*). That symmetry is load-bearing and also limiting: an Intertropical Convergence Zone that
**migrates off the equator breaks equatorial symmetry by definition**, so the hemisphere grid *cannot*
represent it. This module lifts rung 0 to the **full sphere** ``x ∈ [−1, 1]`` (two real poles, the equator
an interior point) and adds the **energetic ITCZ**: the latitude where the atmosphere's emergent meridional
energy transport crosses zero — the **energy-flux equator** (EFE) — which migrates toward the warmer
hemisphere under an interhemispheric energy imbalance (Kang et al. 2008; Bischoff & Schneider 2014;
Schneider, Bischoff & Haug 2014).

A SIBLING model — ``ebm.py`` is untouched (the every-rung-is-a-sibling discipline)
-----------------------------------------------------------------------------------
Like :mod:`planet.moist_ebm` (rung 2.5) and :mod:`planet.baroclinic_qg` (rung 3), this is a **new model
alongside rung 0, not an edit of it**. It reuses :mod:`engines.diffusion` exactly as ``ebm.py`` does
(``grid_from_edges(linspace(-1, 1, n+1))``, the scaled transport coefficient ``(D/C)(1−x²)`` vanishing at
**both** poles, ``Neumann(0)`` at each). "Re-validate the protected Phase-1 climate" therefore becomes a
**cross-model reduction check** (:func:`SphereEBM` reproduces the hemisphere ``ebm.py`` climate to the
relaxation tolerance under symmetric forcing — see :mod:`planet.tests.test_sphere_ebm`), exactly the
SW↔QG rigid-lid bridge pattern, not a modification of the spine.

The energy-flux equator — and what it is (and is NOT)
----------------------------------------------------
The northward atmospheric energy transport across latitude ``x`` is
``H(x) = −2π a² · D (1 − x²) ∂T/∂x`` (W; :meth:`SphereEBM.atmospheric_transport`); the **EFE** is the zero
of ``H`` nearest the equator (:meth:`SphereEBM.energy_flux_equator`). In the modern energetic framework the
ITCZ sits at the EFE, so this latitude is **identified with the ITCZ**. Three honesty edges are baked into
that identification and named so the rung is not over-read:

* **This is a DRY EBM.** There is no moisture, no convection, no convergence — in a diffusive EBM the EFE
  is simply the latitude of the **temperature maximum** (``∂T/∂x = 0``). Calling it "the ITCZ" is an
  *identification* by appeal to external moist theory, not an emergent rain belt. Wiring it into the
  prescribed precipitation band (:func:`itcz_informed_precip`) **relocates** that band to the
  energetically-warmest latitude; it does **not** make rainfall emergent. (The *moisture-convergence* sign
  in the deep tropics is a separate fix — the eddy-only :func:`planet.moist.moisture_convergence` stays
  backwards there by default; the opt-in :func:`planet.moist.hadley_moisture_convergence` adds the mean
  Hadley cell that flips it. This rung is about ITCZ *position*, that one about the deep-tropical *sign*.)
* **The asymmetry is IMPOSED, not emergent.** The interhemispheric imbalance that moves the EFE is a knob
  — a prescribed cross-equatorial energy flux ``Q(x)`` (``∫Q dx = 0``) or an antisymmetric albedo — not
  the output of an ocean model (the coupler's synthetic-gradient precedent, :mod:`planet.circ_precip`).
* **The shift DIRECTION is by-construction; only the MAGNITUDE is a test.** That the ITCZ moves toward the
  warm hemisphere is guaranteed by the sign of the response, not evidence of anything (the same
  "guaranteed result" the QG rung carried). The non-vacuous question is the *rate*.

The sensitivity is a CLOSED FORM of the calibrated D — banked at that altitude
------------------------------------------------------------------------------
For a small antisymmetric perturbation about the symmetric mean state, the EFE shift ``δ`` and the
cross-equatorial transport ``AHT_eq = H(0)`` depend on the perturbation **only** through the equatorial
temperature-gradient anomaly, so their ratio is the forcing-independent

    δ / AHT_eq = 1 / (2π a² · D · T̄ₓₓ(0))          (:func:`itcz_sensitivity_closed`)

— a **near-algebraic consequence of the already-calibrated transport ``D`` and the mean-state curvature
``T̄ₓₓ(0)``**, *not* an emergent prediction. The forcing-independence (a cross-equatorial Q-flux and an
albedo asymmetry give the **same** number) is this linear-operator identity, not robustness. Evaluated:
**≈ −6.3 deg/PW** (constant-albedo, the splitting-free :meth:`SphereEBM.steady_linear` value matching the
closed form) and **≈ −4.9 deg/PW** (the present-day ice climate, a steeper equatorial curvature) — the
**same order** as the observed ``~3 deg/PW`` (Donohoe et al. 2013) but a **factor ~1.5–2 high**, and
``∝ 1/D``. So what is banked is *"the ITCZ sensitivity is a closed-form consequence of the calibrated D and
the mean-state curvature, of the observed order"* — a low-degrees-of-freedom corroboration that ``D`` is
realistic, **not** "the EBM predicts the observed ITCZ migration."

The operator-splitting gotcha (spike-caught)
--------------------------------------------
The EFE and its sensitivity must be read off a **converged** temperature profile. The default Strang-split
relaxation at ``n_tau = 0.5`` (``ebm.py``'s default step) carries an **O(Δt) steady-state splitting error
in the shape** (first-order — the backward-Euler transport substep, split against the exact radiation
half-step; the global *mean* stays exact) that *steepens* the equatorial curvature (so ``ebm.py`` anchors
its own North check with ``method="direct"``, never the relaxation). :meth:`steady_linear` (the dt-free
solve) is splitting-free
and is the path the tight sensitivity uses; the ice (nonlinear) climate is converged with a small
``n_tau``. A naive ``n_tau = 0.5`` reading inflates the curvature and spuriously *lowers* ``|deg/PW|``
toward ``−3.8`` — a numerical artifact, not the physical number.

Validation triad (plan §3)
--------------------------
* **Tight (analytic/structural).** The North (1975) two-mode on the full sphere via :meth:`steady_linear`
  (constant albedo) at ~2nd order in Δx; the **reduction** to the hemisphere ``ebm.py`` climate (1e-9 under
  symmetric forcing); the closed form ``δ/AHT = 1/(2π a² D T̄ₓₓ(0))`` reproduced by the engine to ~1 %; and
  EFE ``= 0`` exactly for symmetric forcing.
* **Real-but-loose (the unlock, banked at the lower altitude).** The ITCZ sensitivity ``≈ −5 deg/PW`` (ice)
  / ``−6.3`` (no-ice) — same order as observed ``~3``, a factor ~1.5–2 high; ``∝ 1/D``; direction
  by-construction.
* **Plumbing.** Symmetric forcing ⟹ EFE ``= 0``; the precip wiring reduces to the rung-0 ITCZ centre
  bit-for-bit when ``φ_EFE = 0``.

Units — SI, climlab-conventional (W m⁻², °C, x = sin φ on [−1, 1]); transport in PW (10¹⁵ W).
Sources (cited at build, the ``[[…-source]]`` discipline): Kang, Held, Frierson & Zhao 2008 (hemispheric
forcing → ITCZ shift); Bischoff & Schneider 2014 (*J. Climate*, the energetic ITCZ); Schneider, Bischoff &
Haug 2014 (*Nature*); Donohoe et al. 2013 (the ~3 deg/PW cross-equatorial sensitivity); North 1975 (the
two-mode, already pinned in :mod:`planet.ebm`). Extends [[ebm-radiation-source]], [[precip-parameterization-source]].
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np
from scipy.linalg import solve_banded

from engines.diffusion import Diffusion1D, grid_from_edges, Neumann
from . import precip
from .circulation import R_EARTH
from .ebm import (
    A_OLR, B_OLR, D_TRANSPORT, T_FREEZE, WATER_DEPTH, RHO_WATER, CW_WATER,
    legendre_P2, two_mode_solution,
)

AREA_FACTOR = 2.0 * math.pi * R_EARTH ** 2     # dA = 2π a² dx  (x = sin φ); H = AREA_FACTOR·flux (W)
PW = 1.0e15                                    # W — petawatt


# --------------------------------------------------------------------------- #
# The frozen result — plain arrays (the loose-coupling currency).
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class SphereClimate:
    """A pole-to-pole equilibrium climate plus its energetic-ITCZ diagnostics (plain arrays/scalars).

    ``x`` the cell-center area coordinates ``sin φ`` on [−1, 1]; ``T`` the equilibrium profile (°C);
    ``global_mean_T`` the area mean ``∫T dx / 2`` (°C). ``phi_efe`` the energy-flux-equator latitude
    (degrees, the ITCZ identification — 0 for a symmetric climate, signed toward the warmer hemisphere
    under an imbalance); ``aht_eq`` the cross-equatorial northward atmospheric energy transport ``H(0)``
    (PW). ``converged``/``iterations`` record the relaxation.
    """

    x: np.ndarray
    T: np.ndarray
    global_mean_T: float
    phi_efe: float
    aht_eq: float
    converged: bool
    iterations: int

    def latitude_deg(self) -> np.ndarray:
        """Grid latitudes ``φ = asin(x)`` in degrees (south pole −90° → north pole +90°)."""
        return np.degrees(np.arcsin(np.clip(self.x, -1.0, 1.0)))


# --------------------------------------------------------------------------- #
# The closed-form sensitivity (the banked altitude — a function of D and curvature).
# --------------------------------------------------------------------------- #
def itcz_sensitivity_closed(Txx0: float, D: float = D_TRANSPORT) -> float:
    """The closed-form ITCZ sensitivity ``δ/AHT_eq = 1/(2π a² D T̄ₓₓ(0))`` in **deg / PW**.

    Derived in the module docstring: in a dry EBM the EFE is the temperature maximum, so for a small
    antisymmetric perturbation the shift and the cross-equatorial transport both depend on the response
    only through the equatorial gradient anomaly, leaving this forcing-independent ratio. ``T̄ₓₓ(0)`` is the
    equatorial curvature of the symmetric base state (°C, negative — the temperature maximum), ``D`` the
    EBM transport coefficient. The radians→degrees and W→PW conversions are folded in. **This is the
    banked quantity's altitude** — a property of the calibrated ``D`` and the mean state, not a prediction.
    """
    return (180.0 / math.pi) * PW / (AREA_FACTOR * float(D) * float(Txx0))


# --------------------------------------------------------------------------- #
# The pole-to-pole EBM (x in [-1, 1]) — engines/diffusion reused, ebm.py untouched.
# --------------------------------------------------------------------------- #
class SphereEBM:
    """Full-sphere energy-balance model on ``x = sin φ ∈ [−1, 1]``; see the module docstring.

    Mirrors :class:`planet.ebm.EnergyBalanceModel` (Strang-split relaxation + a direct linear solve) on the
    doubled grid with **two** insulated poles and the equator interior. ``D`` is a uniform scalar transport
    coefficient (the array-``D(x)`` rung-1 feedback is out of scope for this diagnostic). The forcing is
    injected as an ``absorbed(x, T) → W m⁻²`` callable (forcing-agnostic, like the parent); an optional
    steady source ``Q(x)`` (the imposed cross-equatorial flux, ``∫Q dx ≈ 0``) is the interhemispheric knob.
    """

    def __init__(self, A: float = A_OLR, B: float = B_OLR, D: float = D_TRANSPORT,
                 T_freeze: float = T_FREEZE, water_depth: float = WATER_DEPTH, n_cells: int = 360):
        if B <= 0.0:
            raise ValueError(f"B (OLR slope) must be positive for a stable relaxation, got {B}")
        if float(D) < 0.0:
            raise ValueError(f"D (transport) must be non-negative, got {D}")
        self.A, self.B, self.D = float(A), float(B), float(D)
        self.T_freeze = float(T_freeze)
        self.n_cells = int(n_cells)
        self.C = RHO_WATER * CW_WATER * float(water_depth)       # J m⁻² K⁻¹ (timescale only)
        self.tau_rad = self.C / self.B
        # Full sphere: x = sin φ on [-1, 1], insulated at BOTH poles, equator interior.
        self.grid = grid_from_edges(np.linspace(-1.0, 1.0, self.n_cells + 1))
        self.x = self.grid.centers
        self._Dcells = (self.D / self.C) * (1.0 - self.x ** 2)   # vanishes at x = ±1
        self.solver = Diffusion1D(self.grid, self._Dcells, Neumann(0.0), Neumann(0.0))

    # -- means ------------------------------------------------------------- #
    def global_mean(self, T: np.ndarray) -> float:
        """Area-mean temperature ``∫T dx / 2`` (°C) — the engine ``total`` over the length-2 domain."""
        return float(self.solver.total(T) / self.grid.length)

    # -- the Strang-split relaxation (general / nonlinear path) ------------- #
    def equilibrate(self, absorbed_fn, T_init, Q: Optional[np.ndarray] = None,
                    n_tau: float = 0.5, tol: float = 1e-10, max_iter: int = 200000) -> SphereClimate:
        """Relax to the steady climate by half-radiation / full-transport / half-radiation stepping.

        ``Q`` is an optional steady energy source (W m⁻²; the imposed interhemispheric flux). **Note the
        splitting-error gotcha (module docstring):** the EFE/sensitivity wants a *converged* profile — use a
        **small ``n_tau``** here (or :meth:`steady_linear` in the constant-albedo case); the default
        ``n_tau = 0.5`` carries an O(Δt) curvature-steepening (shape) error — first-order, from the
        backward-Euler transport split against the exact radiation half-step (the mean stays exact).
        """
        T = np.full(self.x.shape, float(T_init)) if np.isscalar(T_init) else np.array(T_init, float)
        Qarr = np.zeros_like(self.x) if Q is None else np.asarray(Q, dtype=float)
        dt = n_tau * self.tau_rad
        decay = math.exp(-0.5 * dt * self.B / self.C)

        def rad_half(T):
            T_eq = (absorbed_fn(self.x, T) + Qarr - self.A) / self.B
            return T_eq + (T - T_eq) * decay

        converged, it = False, 0
        for it in range(1, max_iter + 1):
            T_old = T
            T = rad_half(T)
            T = self.solver.step(T, dt)
            T = rad_half(T)
            if np.max(np.abs(T - T_old)) < tol:
                converged = True
                break
        return self._result(T, converged, it)

    # -- the direct linear steady solve (the splitting-free / tight path) --- #
    def _transport_tridiag(self):
        """Tridiagonals of ``L_T = D·d/dx[(1−x²)d/dx]`` — assembled exactly as the engine assembles it."""
        dx = self.grid.widths
        Dc = self._Dcells * self.C                               # unscale: (D/C)(1−x²)·C = D(1−x²)
        Dface = 2.0 * Dc[:-1] * Dc[1:] / (Dc[:-1] + Dc[1:])      # harmonic-mean interior faces
        Tt = Dface / np.diff(self.grid.centers)
        n = self.grid.n
        sub, diag, sup = np.zeros(n), np.zeros(n), np.zeros(n)
        sup[:-1] += Tt / dx[:-1]; diag[:-1] += -Tt / dx[:-1]
        sub[1:] += Tt / dx[1:];   diag[1:] += -Tt / dx[1:]
        return sub, diag, sup

    def steady_linear(self, absorbed_fn, Q: Optional[np.ndarray] = None) -> SphereClimate:
        """Direct dt-free steady solve of the **linear** EBM — splitting-free; constant-albedo only.

        Solves ``(L_T − B·I) T = A − S(x)(1−α) − Q`` in one tridiagonal solve (no time-stepping ⟹ **no
        splitting error**), reproducing the North two-mode to the engine's spatial order and giving the
        clean (no-ice) EFE sensitivity. **Raises** if ``absorbed_fn`` is state-dependent (an ice feedback
        must go through :meth:`equilibrate` with a small ``n_tau``).
        """
        a_cold = absorbed_fn(self.x, np.full(self.x.shape, -100.0))
        a_warm = absorbed_fn(self.x, np.full(self.x.shape, 100.0))
        if not np.allclose(a_cold, a_warm):
            raise ValueError("steady_linear requires a state-independent (constant-albedo) absorbed field; "
                             "use equilibrate(..., n_tau=small) for the ice-albedo feedback")
        Qarr = np.zeros_like(self.x) if Q is None else np.asarray(Q, dtype=float)
        sub, diag, sup = self._transport_tridiag()
        diag = diag - self.B
        rhs = self.A - np.asarray(a_cold, dtype=float) - Qarr
        n = self.grid.n
        ab = np.zeros((3, n))
        ab[0, 1:] = sup[:-1]; ab[1, :] = diag; ab[2, :-1] = sub[1:]
        T = solve_banded((1, 1), ab, rhs)
        return self._result(T, True, 0)

    # -- energetic-ITCZ diagnostics ---------------------------------------- #
    def atmospheric_transport(self, T: np.ndarray) -> np.ndarray:
        """Northward atmospheric energy transport ``H(x) = −2π a² D (1−x²) ∂T/∂x`` (PW).

        The EBM's diffusive heat transport, in petawatts (the realistic ~few-PW poleward transport). Its
        zero nearest the equator is the energy-flux equator / ITCZ; its peak is the poleward transport
        maximum. ``∂T/∂x`` by centered differences on the uniform grid.
        """
        dTdx = np.gradient(np.asarray(T, dtype=float), self.x)
        return -AREA_FACTOR * self.D * (1.0 - self.x ** 2) * dTdx / PW

    def energy_flux_equator(self, T: np.ndarray) -> tuple[float, float]:
        """``(φ_EFE, AHT_eq)``: the ITCZ latitude (deg, zero of ``H`` nearest the equator) and ``H(0)`` (PW).

        The EFE migrates toward the warmer hemisphere under an interhemispheric imbalance (signed φ); the
        cross-equatorial transport ``AHT_eq = H(0)`` is the energetic driver the ``deg/PW`` sensitivity is
        measured against. A symmetric climate gives ``φ_EFE = 0`` exactly (``H`` is odd ⟹ ``H(0) = 0``).
        """
        H = self.atmospheric_transport(T)
        x = self.x
        band = np.abs(x) < 0.6                                   # tropical bracket
        xb, Hb = x[band], H[band]
        crossings = np.where(np.diff(np.sign(Hb)) != 0)[0]
        aht_eq = float(np.interp(0.0, x, H))
        if crossings.size == 0:
            return float("nan"), aht_eq
        i = crossings[np.argmin(np.abs(xb[crossings]))]          # crossing nearest the equator
        x0 = xb[i] - Hb[i] * (xb[i + 1] - xb[i]) / (Hb[i + 1] - Hb[i])
        return float(np.degrees(np.arcsin(np.clip(x0, -1.0, 1.0)))), aht_eq

    def equatorial_curvature(self, T: np.ndarray, band: float = 0.15) -> float:
        """``T̄ₓₓ(0)`` (°C) — equatorial curvature from a degree-2 fit on ``|x| < band`` (the closed form's input)."""
        m = np.abs(self.x) < band
        return float(2.0 * np.polyfit(self.x[m], np.asarray(T, dtype=float)[m], 2)[0])

    def _result(self, T: np.ndarray, converged: bool, it: int) -> SphereClimate:
        phi_efe, aht_eq = self.energy_flux_equator(T)
        return SphereClimate(x=self.x, T=T, global_mean_T=self.global_mean(T),
                             phi_efe=phi_efe, aht_eq=aht_eq, converged=converged, iterations=it)

    # -- the ITCZ-migration sensitivity (deg/PW) --------------------------- #
    def itcz_sensitivity(self, absorbed_fn, q_amps=(0.0, 1.0, 2.0, 4.0), linear: bool = True,
                         n_tau: float = 0.02) -> tuple[float, float]:
        """``(slope_measured, slope_closed)`` in deg/PW — the EFE migration rate and its closed-form value.

        Sweeps a cross-equatorial Q-flux ``Q(x) = q·x`` (odd, ``∫Q dx = 0``) over ``q_amps`` (W m⁻²),
        fits ``φ_EFE`` against ``AHT_eq``, and returns that slope beside the closed form
        :func:`itcz_sensitivity_closed` evaluated at the symmetric base state's curvature — they must agree
        (the tight check). ``linear=True`` uses the splitting-free :meth:`steady_linear` (constant albedo);
        ``linear=False`` converges the ice climate with the given small ``n_tau`` (the splitting-error
        gotcha — slow). Sign is negative: the ITCZ shifts toward the warm hemisphere.
        """
        ahts, phis = [], []
        for q in q_amps:
            Q = q * self.x
            c = (self.steady_linear(absorbed_fn, Q=Q) if linear
                 else self.equilibrate(absorbed_fn, T_init=15.0, Q=Q, n_tau=n_tau, tol=1e-12,
                                       max_iter=400000))
            phis.append(c.phi_efe); ahts.append(c.aht_eq)
            if q == 0.0:
                Txx0 = self.equatorial_curvature(c.T)
        slope = float(np.polyfit(ahts, phis, 1)[0])
        return slope, itcz_sensitivity_closed(Txx0, self.D)


# --------------------------------------------------------------------------- #
# Convenience: a present-day full-sphere climate (reuses the rung-0 ice-albedo forcing).
# --------------------------------------------------------------------------- #
def present_day_sphere_climate(Q: Optional[np.ndarray] = None, n_cells: int = 360,
                               n_tau: float = 0.02, ic_equator: float = 30.0,
                               ic_pole: float = -30.0, **kw) -> SphereClimate:
    """The present-day full-sphere climate (ice-albedo feedback), optionally with an imposed Q-flux.

    Uses the rung-0 ice-albedo forcing (:func:`planet.albedo.absorbed_shortwave`) on the full sphere, from
    an Earth-like start (warm equator, frozen poles). Defaults to a **converged** ``n_tau = 0.02`` (the
    splitting-error gotcha — the EFE wants a converged profile). With ``Q = None`` the climate is symmetric
    (``φ_EFE = 0``); pass an odd ``Q`` (e.g. ``q·model.x``) to drive the ITCZ off the equator.
    """
    from .albedo import absorbed_shortwave                       # local import: albedo imports ebm
    m = SphereEBM(n_cells=n_cells, **kw)
    T_init = ic_equator + (ic_pole - ic_equator) * np.abs(m.x)   # warm equator → cold (iced) both poles
    return m.equilibrate(lambda x, T: absorbed_shortwave(x, T), T_init, Q=Q, n_tau=n_tau)


# --------------------------------------------------------------------------- #
# The precip wiring — the "moist precip pattern" half (opt-in; ITCZ band ← φ_EFE).
# --------------------------------------------------------------------------- #
def itcz_informed_precip(climate: SphereClimate) -> np.ndarray:
    """Full-sphere precipitation ``P(φ)`` (cm/yr) with the **ITCZ band centred on the energy-flux equator**.

    The rung-0 pattern (:func:`planet.precip.precip_pattern`) with its ITCZ band shifted from the equator to
    ``climate.phi_efe`` (the new ``itcz_center_deg`` seam), times the rung-0 global C–C amplitude. Reduces
    to the symmetric rung-0 field **bit-for-bit when ``φ_EFE = 0``** (the plumbing reduction). **Honest
    scope (module docstring):** this *relocates a prescribed band* to the energetically-warmest latitude in
    a **dry** model — it is **not** emergent rainfall, and the rung-0 :mod:`planet.precip` stays the default
    everywhere else (opt-in, like :mod:`planet.circ_precip`).
    """
    return precip.precipitation(climate.latitude_deg(), climate.global_mean_T,
                                itcz_center_deg=climate.phi_efe)
