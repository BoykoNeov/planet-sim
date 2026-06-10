"""The latitudinal energy-balance model — the planet's climate spine (Planet Phase 1).

Project #3 (the capstone) opens where Steel and Chip both started: **on the frozen
diffusion spine** (:mod:`engines.diffusion`). Latitudinal heat transport on a sphere
*is* a diffusion equation, so the EBM reuses the sealed parabolic solver a **third**
time — now as a planet's pole-to-equator heat transport. The governing 1-D
annual-mean energy-balance model (Budyko 1969 / Sellers 1969 / North 1975) is

    C ∂T/∂t = D ∂/∂x[(1−x²) ∂T/∂x] + S(x)(1−α(x,T)) − (A + B·T),    x = sin φ

with φ latitude, ``x = sin φ`` the area coordinate (equal Δx = equal area on the
sphere — that is *why* x=sinφ is the natural coordinate, and why the global mean is a
plain ``∫₀¹ T dx``). ``D·∂/∂x[(1−x²)∂T/∂x]`` is meridional diffusive heat transport,
``S(x)(1−α)`` the absorbed shortwave, and ``A + B·T`` the linearized outgoing longwave
(the OLR *parameterization* — the named scope edge). The headline payoff (Phase 1's
banked artifact) is the **Snowball-Earth hysteresis** the ice-albedo feedback produces
(:mod:`planet.albedo`); this module is the engine-coupled foundation it stands on.

How the frozen engine is reused — transport in the engine, radiation split around it
-------------------------------------------------------------------------------------
The frozen solver advances ``∂u/∂t = ∂/∂x(D_eng ∂u/∂x) + S(x,t)``. The **transport**
maps on exactly: in heat mode ``u = T`` and the spatially-varying array diffusivity is

    D_eng(x) = (D / C)·(1 − x²)        (W m⁻² K⁻¹ → per-step diffusivity; vanishes at the pole)

with **Neumann(0)** at both ends — the equator (x=0) end by hemispheric symmetry, the
pole (x=1) end by no-flux through a pole. That insulated pair is the engine's exact
conservation BC, so the transport **only redistributes** heat: ``∫T dx`` is conserved
structurally (the engine's frozen no-flux invariant).

The **radiation cannot be the engine's source.** ``−(A+B·T)`` is *linear in the state*
and ``S(x)(1−α(T))`` is a *nonlinear pointwise function of T* (the albedo threshold),
while the engine's ``S`` is only ``S(x,t)`` — not ``S(u)``. So radiation is composed
*around* the engine by **Strang operator splitting**, the **identical idiom Steel's
Jominy Phase-2a** used to graft its state-dependent lateral sink ``−h(T−T_air)`` onto
the same solver (:mod:`projects.steel.jominy`). Each step is: a half-step of the local
radiation ODE, a full implicit transport step (the frozen solver, untouched), a second
radiation half-step. The local ODE ``dT/dt = [S(1−α) − A − B·T]/C`` has, with α frozen
at the substep's T, the **exact analytic** solution

    T ← T_eq + (T − T_eq)·exp(−½ Δt·B/C),     T_eq = (S(1−α) − A)/B,   τ_rad = C/B

— exact for the linear ``−B·T`` relaxation (as Jominy's lateral decay was), so the
composition inherits the engine's unconditional stability with a 2nd-order splitting
error in Δt. The **albedo threshold makes that pointwise step a nonlinear relaxation
(α re-evaluated each substep) — which is precisely what creates the multiple equilibria**
the Snowball demo rides. The reuse is therefore of *both* the frozen engine and the
frozen splitting pattern — a stronger reuse than "the EBM is the engine."

Because equilibria are the target, ``C`` (heat capacity) sets only the relaxation
*timescale*, not the steady state: at ``∂T/∂t = 0`` it cancels (``test_ebm`` asserts the
equilibrium is independent of ``water_depth``). It is kept physical (a mixed-layer ocean,
``C = ρ_w c_w · depth``) so the time-stepping is honest and ``dt`` is expressed in ``τ_rad``.

Three interchangeable steady-state modes (two orthogonal knobs)
--------------------------------------------------------------
Reaching the steady climate carries two independent, user-selectable choices — useful both
for accuracy/speed trade-offs and as a mutual cross-check web for validation:

* **Face-diffusivity representation** (``face=`` on the model). The transport coefficient
  ``(1−x²)`` *vanishes* at the pole, where the engine's harmonic-mean face averaging is biased
  ~25 % low (harmonic mean is exact for a *discontinuous* D, the layered-media case it was
  designed for, but not for a smoothly-vanishing one). ``"harmonic"`` feeds the plain
  cell-centered ``(1−x²)`` and accepts that bias (a named ~0.1 °C polar floor on the North
  check); ``"exact"`` instead feeds cell values **pre-distorted** so the engine's harmonic
  mean *reproduces the true face coefficient* (:func:`cell_diffusivity_for_exact_faces`),
  removing the floor (North → ~0.01 °C). Both feed the *same* relaxation and the *same* frozen
  engine — the engine is never modified.
* **Steady-state method** (``method=`` on :meth:`EnergyBalanceModel.equilibrium`).
  ``"relax"`` is the **Strang-split relaxation** above — the general path, the *only* one that
  handles the nonlinear ice-albedo feedback (so the Snowball sweep uses it). ``"direct"`` is a
  **dt-free linear solve** (:meth:`steady_linear`) of ``(L_T − B·I)T = A − S(1−α)`` — exact (no
  splitting error), but valid **only for a state-independent (constant-albedo) absorbed field**,
  so it serves as the fast splitting-error-free *reference* the no-feedback North check is
  cross-validated against (it raises if handed the ice feedback). Its reconstructed transport
  operator ``L_T`` is **pinned to the frozen engine** by a test (the engine's transport ``step``
  equals solving ``(I − dt·L_T/C)``), so the "direct" path cannot silently drift from the engine.

The default — ``face="harmonic"``, ``method="relax"`` — is the simple, general, snowball-capable
combination; the others are validation/accuracy alternates.

Validation triad (plan §3) — what is asserted tight vs loose
------------------------------------------------------------
* **Analytical limit (tight).** (a) The **0-D** global-mean equilibrium
  ``T̄ = [(S₀/4)(1−ᾱ) − A]/B`` (:func:`equilibrium_temperature_0d`) — the transport is
  mean-preserving, so the relaxed mean obeys the *discrete* energy balance exactly (net-TOA
  ~machine, the conservation leg) and matches the *continuous* ``T̄`` to the grid's O(1/n²)
  quadrature limit (point-sampled insolation does not integrate ``P₂`` to exactly 0). (b) The **North (1975)
  two-mode** solution ``T(x) = T₀ + T₂·P₂(x)`` (:func:`two_mode_solution`), exact because
  ``d/dx[(1−x²) dP₂/dx] = −6 P₂`` (Legendre's equation) makes the transport diagonal in
  the P₂ mode — reproduced by the FV engine at **~2nd order** in Δx (the contract's
  convergence invariant). This validates that transport + linear radiation are assembled
  correctly, in the no-feedback limit where the splitting is exact.
* **Conservation (tight).** At equilibrium the **global energy balance**
  ``⟨S(1−α)⟩ = A + B·⟨T⟩`` (absorbed solar = OLR, area-mean ``⟨·⟩ = ∫₀¹·dx``) holds to
  machine precision; the diffusive transport conserves ``∫T dx`` structurally (the frozen
  no-flux invariant, re-confirmed for the Neumann(0)/Neumann(0) pair).
* **Benchmark (loose).** climlab's EBM — present-day ice line ~70°, the Snowball threshold,
  the hysteresis width (:mod:`planet.climate_reference`, the pycalphad pattern:
  a frozen reference table keeps the triad green without the ``[climate]`` extra).

Non-circularity, named scope edge (plan §3)
-------------------------------------------
*Validated tight:* the structural reuse (frozen transport + the analytic two-mode it must
reproduce) and the global-balance conservation. *Calibrated/flagged:* the radiation/albedo
constants ``A, B, D, α, Tf`` are the climlab/North/Budyko defaults ([[ebm-radiation-source]],
pinned at build) — so the *exact* threshold numbers are calibration-dependent, asserted only
in loose bands (the way Steel's 1045 knee and Chip's contrast curve were).
*Scope edge:* the **linear OLR ``A+B·T``** is a parameterization accurate only near the
present climate (deep-snowball / hot states need real radiative transfer — the rung-4
deferral, named not modeled); v1 is **annual-mean** (no seasonal cycle / obliquity) and
**zonal-mean** (no land/ocean contrast, no orography). T is in **°C** throughout (the
climlab convention; ``A+B·T`` and ``Tf`` are defined for °C).

Units — SI, climlab-conventional (W m⁻², °C, x = sin φ dimensionless)
---------------------------------------------------------------------
S0, S(x), A in **W m⁻²**; B, D in **W m⁻² K⁻¹**; T, Tf in **°C**; C in **J m⁻² K⁻¹**;
x = sin φ dimensionless on [0, 1] (equator → pole). Latitudes reported in **degrees**.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.linalg import solve_banded

from engines.diffusion import Diffusion1D, uniform_grid, Neumann

# --------------------------------------------------------------------------- #
# Pinned radiation / transport / albedo constants ([[ebm-radiation-source]]).
# The climlab `EBM` defaults (North 1975 / Budyko 1969 lineage); load-bearing for every
# climate number, so cited and pinned at build — NOT carried from memory.
# --------------------------------------------------------------------------- #
S0_EARTH = 1365.2          # W m⁻² — present-day solar constant (climlab `const.S0`)
S2_INSOLATION = -0.48      # —      — P₂ coefficient of the annual-mean insolation (more sun at the equator)
A_OLR = 210.0              # W m⁻²  — OLR offset  (OLR = A + B·T)
B_OLR = 2.0                # W m⁻² K⁻¹ — OLR slope (the radiative-damping / longwave feedback)
D_TRANSPORT = 0.555        # W m⁻² K⁻¹ — meridional diffusive transport coefficient
T_FREEZE = -10.0           # °C     — ice-line isotherm (surface freezes / ice albedo below this)
ALBEDO_A0 = 0.30           # —      — ice-free planetary albedo, base value      α_open = a0 + a2·P₂(x)
ALBEDO_A2 = 0.078          # —      — ice-free planetary albedo, P₂ (poleward) term
ALBEDO_ICE = 0.62          # —      — ice/snow albedo where T < Tf
WATER_DEPTH = 10.0         # m      — mixed-layer ocean depth → heat capacity (timescale only)

RHO_WATER = 1.0e3          # kg m⁻³  — for C = ρ_w c_w · depth
CW_WATER = 4.181e3         # J kg⁻¹ K⁻¹


def legendre_P2(x: np.ndarray | float) -> np.ndarray | float:
    """The second Legendre polynomial ``P₂(x) = ½(3x² − 1)``.

    The single mode the annual-mean insolation, the ice-free albedo, and the North
    analytic solution are all written in. ``P₂(0) = −½`` (equator), ``P₂(1) = +1`` (pole),
    and ``∫₀¹ P₂ dx = 0`` (so the P₂ structure carries no global-mean — the global mean is
    the P₀ term). It is the transport operator's eigenfunction: ``d/dx[(1−x²)dP₂/dx] = −6 P₂``.
    """
    x = np.asarray(x, dtype=float)
    return 0.5 * (3.0 * x * x - 1.0)


def insolation(x: np.ndarray | float, S0: float = S0_EARTH, s2: float = S2_INSOLATION) -> np.ndarray:
    """Annual-mean absorbed-free insolation ``S(x) = (S₀/4)(1 + s₂·P₂(x))`` (W m⁻²).

    The factor ``1/4`` is the disk/sphere geometric mean; the ``s₂·P₂`` term concentrates
    sunlight at the equator (``s₂ < 0`` → ``S`` larger where ``P₂ < 0``). Its area mean is
    ``∫₀¹ S dx = S₀/4`` (since ``∫₀¹ P₂ dx = 0``), so the global-mean forcing is the textbook
    ``S₀/4``. This is *incident* flux; the absorbed flux multiplies by the coalbedo ``1−α``
    (supplied by :mod:`planet.albedo`).
    """
    return (S0 / 4.0) * (1.0 + s2 * legendre_P2(x))


def equilibrium_temperature_0d(S0: float = S0_EARTH, albedo: float = ALBEDO_A0,
                               A: float = A_OLR, B: float = B_OLR) -> float:
    """0-D global-mean equilibrium temperature ``T̄ = [(S₀/4)(1−ᾱ) − A]/B`` (°C).

    The exact zero-dimensional anchor: globally, absorbed solar ``(S₀/4)(1−ᾱ)`` balances
    OLR ``A + B·T̄``. For present-day Earth (``S₀ = 1365.2``, ``ᾱ = 0.30``) this is
    ``≈ 14.4 °C`` — the global-mean surface temperature. The 1-D model's relaxed *mean*
    reproduces this to the grid's O(1/n²) quadrature limit (the transport is mean-preserving;
    the *discrete* energy balance / net-TOA is machine-exact), independent of the spatial
    structure — the tight analytic leg's first half.
    """
    return ((S0 / 4.0) * (1.0 - albedo) - A) / B


def two_mode_solution(x: np.ndarray | float, S0: float = S0_EARTH, albedo: float = ALBEDO_A0,
                      A: float = A_OLR, B: float = B_OLR, D: float = D_TRANSPORT,
                      s2: float = S2_INSOLATION) -> np.ndarray:
    """North (1975) exact two-mode steady solution ``T(x) = T₀ + T₂·P₂(x)`` (°C), constant albedo.

    With a **constant** albedo (no ice feedback) and the two-mode insolation, the steady EBM
    has an exact closed form because the transport operator is diagonal in the Legendre modes
    (``d/dx[(1−x²)dP₂/dx] = −6 P₂``): matching the P₀ and P₂ terms gives

        T₀ = [(S₀/4)(1−α) − A]/B            (= the 0-D global mean),
        T₂ = (S₀/4)(1−α)·s₂ / (6D + B).

    The transport spreads the equator-to-pole contrast: larger ``D`` (more efficient
    transport) shrinks ``|T₂|`` (a flatter planet). This is the **headline analytic check** —
    the FV engine must reproduce it at ~2nd order in Δx (the spatial half of the tight leg);
    valid only in the no-feedback limit where the splitting is exact.
    """
    T0 = equilibrium_temperature_0d(S0, albedo, A, B)
    T2 = (S0 / 4.0) * (1.0 - albedo) * s2 / (6.0 * D + B)
    return T0 + T2 * legendre_P2(x)


def _harmonic_mean(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Harmonic mean ``2ab/(a+b)`` — the engine's interior-face averaging (assumed positive a, b)."""
    return 2.0 * a * b / (a + b)


def cell_diffusivity_for_exact_faces(grid, coeff) -> np.ndarray:
    """Cell-centered D whose engine harmonic-mean *faces* reproduce ``coeff(x)`` exactly (mode ``"exact"``).

    The frozen engine builds each interior face diffusivity as the harmonic mean of the two
    adjacent **cell** values; for a smoothly-varying coefficient that is ~2nd-order accurate in
    the interior but O(1)-biased where the coefficient nearly vanishes (the pole, ``(1−x²)→0``).
    This inverts the harmonic mean so the engine sees the **true face values**: writing
    ``rₖ = 1/D_cellₖ``, the constraint ``harmonic(D_k, D_{k+1}) = coeff(x_face,k)`` is the
    bidiagonal recurrence ``rₖ₊₁ = 2/coeff(x_face,k) − rₖ`` (seeded from the equator cell). The
    result is a *cell* array the engine consumes unchanged — the engine is not touched — chosen
    so its harmonic-mean faces equal ``coeff`` at every interior face (``test_ebm`` asserts this).
    Raises if the construction yields a non-positive diffusivity (the coefficient varies too
    sharply for the grid) — the positivity guard.
    """
    edges = np.asarray(grid.edges, dtype=float)
    centers = np.asarray(grid.centers, dtype=float)
    g_face = np.asarray(coeff(edges[1:-1]), dtype=float)        # true coefficient at interior faces
    if np.any(g_face <= 0.0):
        raise ValueError("exact-face construction needs a positive coefficient at every interior face")
    r = np.empty(centers.size)
    r[0] = 1.0 / float(coeff(centers[0]))
    for k in range(centers.size - 1):
        r[k + 1] = 2.0 / g_face[k] - r[k]
    if np.any(r <= 0.0):
        raise ValueError("exact-face construction produced a non-positive cell diffusivity; "
                         "the coefficient varies too sharply for this grid")
    return 1.0 / r


def ice_line_latitude(x: np.ndarray, T: np.ndarray, T_freeze: float = T_FREEZE) -> float:
    """Latitude (degrees) of the ice line — where ``T`` crosses ``T_freeze`` — from the field.

    The edge of the polar ice cap: the lowest latitude poleward of which ``T < T_freeze``.
    Found by interpolating the (equator-warm, pole-cold) profile for the crossing
    ``x_ice = sin φ_ice`` and returning ``φ_ice = asin(x_ice)`` in degrees. Two limits are
    handled as the cap's degenerate ends: an **ice-free** planet (``T > Tf`` everywhere) →
    ``90°`` (the cap has shrunk to the pole), a **Snowball** (``T < Tf`` everywhere) → ``0°``
    (ice all the way to the equator). The diagnostic the hysteresis loop tracks.
    """
    x = np.asarray(x, dtype=float)
    T = np.asarray(T, dtype=float)
    if np.all(T > T_freeze):
        return 90.0                                  # ice-free: cap shrunk to the pole
    if np.all(T <= T_freeze):
        return 0.0                                   # Snowball: ice to the equator
    # T decreases with x; reverse to ascending T for np.interp, read the crossing x.
    x_ice = float(np.interp(T_freeze, T[::-1], x[::-1]))
    x_ice = min(1.0, max(0.0, x_ice))
    return math.degrees(math.asin(x_ice))


@dataclass(frozen=True)
class ClimateState:
    """An equilibrium climate: the temperature field plus its banked diagnostics.

    ``x`` are the cell-center area coordinates ``sin φ`` on [0, 1]; ``T`` the equilibrium
    profile (°C, equator → pole). ``global_mean_T`` is ``∫₀¹ T dx`` (°C — the area mean, since
    Δx = equal area); ``ice_line_lat`` the ice-line latitude (degrees); ``net_toa`` the
    residual top-of-atmosphere imbalance ``⟨S(1−α)⟩ − A − B⟨T⟩`` (W m⁻², ~0 at equilibrium —
    the conservation diagnostic). ``converged`` / ``iterations`` record the relaxation. Plain
    arrays/scalars — the loose-coupling currency Phase 2's biomes consume.
    """

    x: np.ndarray
    T: np.ndarray
    global_mean_T: float
    ice_line_lat: float
    net_toa: float
    converged: bool
    iterations: int

    def latitude_deg(self) -> np.ndarray:
        """The grid latitudes ``φ = asin(x)`` in degrees (equator 0° → pole 90°)."""
        return np.degrees(np.arcsin(np.clip(self.x, 0.0, 1.0)))


class EnergyBalanceModel:
    """The transport + radiation EBM machinery: frozen-engine diffusion, Strang-split radiation.

    Holds the OLR (``A``, ``B``) and transport (``D``) constants, the ice-line isotherm ``Tf``
    (for the diagnostic), and the heat capacity ``C`` (from ``water_depth``). Builds the frozen
    :class:`~engines.diffusion.Diffusion1D` **once** in heat mode — array diffusivity
    ``D_eng(x) = (D/C)(1−x²)``, insulated (Neumann 0) at both ends — and reuses it across a
    whole continuation sweep (only the radiation forcing changes with S₀). The radiation is
    **injected** as an ``absorbed(x, T) → W m⁻²`` callable (the absorbed shortwave
    ``S(x)(1−α)``), so this machinery is *forcing-agnostic*: the no-feedback North check feeds
    a constant-albedo callable, the Snowball demo feeds the ice-albedo one
    (:mod:`planet.albedo`). That mirrors the engine's own "machinery here, physical
    constants in the consumer" boundary.
    """

    def __init__(self, A: float = A_OLR, B: float = B_OLR, D: float = D_TRANSPORT,
                 T_freeze: float = T_FREEZE, water_depth: float = WATER_DEPTH,
                 n_cells: int = 180, face: str = "harmonic"):
        if B <= 0.0:
            raise ValueError(f"B (OLR slope) must be positive for a stable relaxation, got {B}")
        if D < 0.0:
            raise ValueError(f"D (transport) must be non-negative, got {D}")
        if face not in ("harmonic", "exact"):
            raise ValueError(f"face must be 'harmonic' or 'exact', got {face!r}")
        self.A = float(A)
        self.B = float(B)
        self.D = float(D)
        self.T_freeze = float(T_freeze)
        self.water_depth = float(water_depth)
        self.n_cells = int(n_cells)
        self.face = face
        self.C = RHO_WATER * CW_WATER * self.water_depth     # J m⁻² K⁻¹ (timescale only)
        self.tau_rad = self.C / self.B                       # s — radiative relaxation time
        # Build the frozen heat-transport solver ONCE: x = sin φ on [0, 1], insulated poles.
        self.grid = uniform_grid(1.0, self.n_cells)
        self.x = self.grid.centers                           # cell-center sin φ
        # The (scaled) transport coefficient (D/C)(1−x²), vanishing at the pole. In "harmonic"
        # mode the plain cell-centered values are handed to the engine (harmonic-mean faces, the
        # ~0.1 °C polar bias); in "exact" mode they are pre-distorted so the engine's harmonic
        # mean reproduces the true face coefficient (no polar bias). Either way the *engine* is
        # untouched — only the cell array it is constructed with differs.
        coeff = lambda x: (self.D / self.C) * (1.0 - np.asarray(x, dtype=float) ** 2)
        if face == "harmonic":
            self._Dcells = coeff(self.x)
        else:
            self._Dcells = cell_diffusivity_for_exact_faces(self.grid, coeff)
        self.solver = Diffusion1D(self.grid, self._Dcells, Neumann(0.0), Neumann(0.0))

    def global_mean(self, T: np.ndarray) -> float:
        """Area-mean temperature ``∫₀¹ T dx`` (°C) — the frozen engine's ``total`` (Δx = equal area)."""
        return float(self.solver.total(T))                   # L = 1, so total = area mean

    def _radiation_half(self, T: np.ndarray, absorbed: np.ndarray, dt: float) -> np.ndarray:
        """Analytic Strang half-step of the local radiation ODE (exact for the frozen-α linear sink).

        ``dT/dt = (absorbed − A − B·T)/C`` integrates exactly over Δt/2 to
        ``T ← T_eq + (T − T_eq)·exp(−½ Δt B/C)`` with ``T_eq = (absorbed − A)/B`` — the Jominy
        ``lateral_half`` pattern (there the sink was ``−h(T−T_air)``; here ``−B(T − T_eq)``).
        ``absorbed`` is frozen at the substep's T by the caller (the albedo nonlinearity).
        """
        T_eq = (absorbed - self.A) / self.B
        decay = math.exp(-0.5 * dt * self.B / self.C)
        return T_eq + (T - T_eq) * decay

    def equilibrate(self, absorbed_fn, T_init, n_tau: float = 0.5,
                    tol: float = 1e-9, max_iter: int = 20000) -> ClimateState:
        """Relax to the steady climate by Strang-split stepping; return its :class:`ClimateState`.

        Each iteration is half-radiation / full-transport / half-radiation, with the absorbed
        shortwave ``absorbed_fn(x, T)`` re-evaluated (α re-frozen) at the start of each
        radiation half-step — the nonlinear pointwise relaxation that, with the ice feedback,
        admits multiple equilibria (which one is reached depends on ``T_init``: that path
        dependence *is* the hysteresis). Iterates with ``dt = n_tau·τ_rad`` until the field
        stops changing (``max|ΔT| < tol``); both sub-operators are unconditionally stable
        (implicit transport; exact-exponential radiation), so the relaxation cannot blow up.

        Parameters
        ----------
        absorbed_fn : callable ``(x, T) -> ndarray``
            The absorbed shortwave ``S(x)(1−α(x,T))`` (W m⁻²) — encapsulates S₀ and the albedo.
        T_init : array | float
            Initial field (°C). A scalar broadcasts to a uniform start. The warm/cold start
            selects the branch under the ice feedback.
        n_tau, tol, max_iter
            Step size as a multiple of ``τ_rad`` (gentle by default so the ice line tracks
            smoothly through thresholds), convergence tolerance (°C), and the iteration cap.
        """
        T = np.array(np.broadcast_to(np.asarray(T_init, dtype=float), self.x.shape), dtype=float)
        dt = n_tau * self.tau_rad
        converged = False
        it = 0
        for it in range(1, max_iter + 1):
            T_old = T
            T = self._radiation_half(T, absorbed_fn(self.x, T), dt)
            T = self.solver.step(T, dt)
            T = self._radiation_half(T, absorbed_fn(self.x, T), dt)
            if np.max(np.abs(T - T_old)) < tol:
                converged = True
                break
        Tbar = self.global_mean(T)
        net_toa = float(np.mean(absorbed_fn(self.x, T)) - self.A - self.B * Tbar)
        return ClimateState(
            x=self.x, T=T, global_mean_T=Tbar,
            ice_line_lat=ice_line_latitude(self.x, T, self.T_freeze),
            net_toa=net_toa, converged=converged, iterations=it,
        )

    # -- the direct linear steady solve (mode C — constant-albedo reference) -------- #
    def _transport_tridiag(self):
        """Tridiagonals of the unscaled transport operator ``L_T = D·d/dx[(1−x²)d/dx]``.

        Assembled **exactly as the frozen engine assembles its operator** — harmonic-mean faces
        of the model's cell diffusivity (plain or exact-face, per ``self.face``), divided by the
        center spacing, with insulated (Neumann 0) ends carrying no exterior-face term. Equals
        ``C × (the engine's transport rate operator)``; ``test_ebm`` pins the two together (the
        engine's transport ``step`` reproduces ``(I − dt·L_T/C)``), so this reconstruction cannot
        drift from the engine.
        """
        dx = self.grid.widths
        Dc = self._Dcells * self.C                          # unscale (D/C)(1−x²)·C = D(1−x²)
        Dface = _harmonic_mean(Dc[:-1], Dc[1:])
        Tt = Dface / np.diff(self.grid.centers)             # interior-face transmissibilities
        n = self.grid.n
        sub = np.zeros(n)
        diag = np.zeros(n)
        sup = np.zeros(n)
        sup[:-1] += Tt / dx[:-1]
        diag[:-1] += -Tt / dx[:-1]
        sub[1:] += Tt / dx[1:]
        diag[1:] += -Tt / dx[1:]
        return sub, diag, sup

    def steady_linear(self, absorbed_fn) -> ClimateState:
        """Direct (dt-free) steady solve of the **linear** EBM — the constant-albedo reference (mode C).

        Solves ``(L_T − B·I) T = A − S(x)(1−α)`` in one tridiagonal solve (no time-stepping, hence
        no operator-splitting error), reproducing the North two-mode profile to the engine's spatial
        order. Valid **only** when ``absorbed_fn(x, T)`` is state-independent (no ice feedback): it
        probes the field at a cold and a warm uniform T and **raises** if they differ — the ice
        feedback must go through :meth:`equilibrate` (``method="relax"``).
        """
        a_cold = absorbed_fn(self.x, np.full(self.x.shape, -100.0))
        a_warm = absorbed_fn(self.x, np.full(self.x.shape, 100.0))
        if not np.allclose(a_cold, a_warm):
            raise ValueError("steady_linear requires a state-independent (constant-albedo) absorbed "
                             "field; use method='relax' (equilibrate) for the ice-albedo feedback")
        absorbed = np.asarray(a_cold, dtype=float)
        sub, diag, sup = self._transport_tridiag()
        diag = diag - self.B                                # (L_T − B·I): diagonally dominant, invertible
        rhs = self.A - absorbed                             # = A − S(1−α)
        n = self.grid.n
        ab = np.zeros((3, n))
        ab[0, 1:] = sup[:-1]
        ab[1, :] = diag
        ab[2, :-1] = sub[1:]
        T = solve_banded((1, 1), ab, rhs)
        Tbar = self.global_mean(T)
        net_toa = float(np.mean(absorbed_fn(self.x, T)) - self.A - self.B * Tbar)
        return ClimateState(
            x=self.x, T=T, global_mean_T=Tbar,
            ice_line_lat=ice_line_latitude(self.x, T, self.T_freeze),
            net_toa=net_toa, converged=True, iterations=0,
        )

    def equilibrium(self, absorbed_fn, T_init=None, method: str = "relax", **kw) -> ClimateState:
        """Steady climate by the selected method — the interchangeable entry point.

        ``method="relax"`` (default) runs the Strang-split relaxation (:meth:`equilibrate`, the
        general / nonlinear-capable path — needs ``T_init``); ``method="direct"`` runs the
        dt-free linear solve (:meth:`steady_linear`, constant-albedo only). Combined with the
        model's ``face`` knob, the two give the A/B/C trio (harmonic-relax / exact-relax / direct).
        """
        if method == "relax":
            if T_init is None:
                raise ValueError("method='relax' needs an initial field T_init")
            return self.equilibrate(absorbed_fn, T_init, **kw)
        if method == "direct":
            if kw:
                raise ValueError(f"method='direct' takes no relaxation kwargs, got {sorted(kw)}")
            return self.steady_linear(absorbed_fn)
        raise ValueError(f"method must be 'relax' or 'direct', got {method!r}")
