"""The **seasonal** zonal energy-balance model — waking heat capacity, casting continentality (rung 5B.1).

Every EBM up this staircase — rung 0 (:mod:`planet.ebm`), the full sphere (:mod:`planet.sphere_ebm`),
the moist diffuser (:mod:`planet.moist_ebm`) — solves for an **equilibrium**: ``∂T/∂t = 0``. At
equilibrium the heat capacity ``C`` **cancels** (``ebm.py``: *"C sets only the relaxation timescale, not
the steady state"*; ``test_ebm`` asserts the equilibrium is independent of ``water_depth``). That is
exactly why those models are **annual-mean** and why an imported ocean-depth / land-mask is *inert*
(§9.3): with ``C`` gone, a shallow-``C`` land column and a deep-``C`` ocean column at the same latitude
reach the **identical** temperature. **Continentality — the large seasonal temperature range over
continental interiors — is precisely zero without a seasonal cycle.**

This module builds that seasonal cycle. It marches the *same* diffusive transport + linear radiation
**forward in time** under **periodic** insolation ``S(x, t)`` (the axial-tilt seasons), to a converged
**annual limit cycle** rather than to a dead steady state. Now ``C`` is **load-bearing** — it sets the
seasonal amplitude and the thermal-lag phase — and resolving **two heat-capacity tiles per latitude**
(a small-``C`` land tile, a large-``C`` ocean tile) makes the land tile swing hard and prompt while the
ocean tile swings little and lags: **continentality, emergent from the ``C`` contrast alone.** This is
the 1-D-in-latitude reduction (North & Coakley 1979) of the 2-D seasonal-and-continents EBM (North,
Mengel & Short 1983); the full ``T(φ, λ, t)`` land–sea *map* is rung 5B.2.

A SIBLING model — the spine EBMs are untouched (the every-rung-is-a-sibling discipline)
----------------------------------------------------------------------------------------
Like :mod:`planet.sphere_ebm` and :mod:`planet.moist_ebm`, this is a **new model beside rung 0, not an
edit of it**. It reuses :mod:`engines.diffusion` exactly as they do — the full sphere ``x = sin φ ∈
[−1, 1]`` (the seasonal cycle is hemispherically **anti**-symmetric — NH summer *is* SH winter — so the
hemisphere-symmetry grid of ``ebm.py`` **cannot** carry it; both poles insulated ``Neumann(0)``, the
equator an interior point), the scaled transport coefficient ``(D/C)(1−x²)`` vanishing at each pole. The
seasonal forcing **reuses the pinned daily-insolation kernel** (:func:`planet.obliquity.daily_mean_insolation`)
— the same formula the obliquity knob averages over a year, here kept on the time axis.

Two solvers — a marcher (the engine-reuse method) and a spectral reference (the tight anchor)
---------------------------------------------------------------------------------------------
* **Time-marcher** (:meth:`SeasonalEBM.march`) — the method. Strang-split stepping
  (half-radiation / full-transport / half-radiation), the ``ebm.py`` idiom, but **periodically forced**
  and run **to a limit cycle** (year-over-year change ``< tol``) instead of to equilibrium. The transport
  substep advances the **zonal-mean** ``T̄ = f_L T_L + f_O T_O`` through the engine with an *effective
  transport heat capacity* ``C_a = (f_L/C_L + f_O/C_O)⁻¹``; the resulting ``ΔT̄`` is handed back to each
  tile as an **energy flux** ``ΔT_i = (C_a/C_i)·ΔT̄``, which is **exactly consistent** with the engine's
  ``T̄`` step (``f_L·(C_a/C_L) + f_O·(C_a/C_O) = 1``) and conserves column energy (uniform ``f_L`` ⟹
  constant ``C_a`` ⟹ the engine's ``∫T̄dx`` no-flux invariant *is* energy conservation). This is the
  path that later carries the nonlinear ice-albedo feedback.
* **Spectral solve** (:meth:`SeasonalEBM.spectral`) — the tight reference, and the reason the payoff is
  *analytic*, not merely *converged*. Because the radiation is **linear** (fixed albedo), the limit cycle
  has a closed form in the frequency domain: Fourier-transform the forcing in time and solve, per temporal
  harmonic ``n`` (angular frequency ``ω_n``), one **complex** banded system

      (I − G_n·L_T)·T̄̂_n = G_n·F̂_n,      G_n = Σ_i f_i /(iω_n C_i + B),   F = S(1−α) − A,

  with ``L_T`` the engine's transport operator (reconstructed exactly as the engine assembles it, so it
  cannot drift). Per-tile spectra follow from ``T̂_{i,n} = (L_T T̄̂_n + F̂_n)/(iω_n C_i + B)``; an inverse
  transform returns ``T_i(x, t)`` with **no time-stepping error**. Every *tight* anchor below is blind to
  time-accuracy in the backward-Euler transport substep (they either turn transport off, or time-average
  the cycle — ``⟨∂T/∂t⟩ = 0`` for *any* periodic state); the spectral solve is what pins the seasonal
  **amplitude and phase**, and cross-validating the marcher against it is what proves the marcher's
  backward-Euler transport is not quietly damping the midlatitude swing.

  The ``n = 0`` harmonic is special and beautiful: ``G_0 = Σ_i f_i / B = 1/B``, so its system collapses
  to ``(B·I − L_T)·T̄̂_0 = F̂_0`` — **exactly** the annual-mean EBM (:meth:`planet.sphere_ebm.SphereEBM.steady_linear`)
  forced by the year-average insolation ``⟨S⟩``. So the reduction-to-the-annual-mean-parent check is not
  bolted on — it *is* the DC component of the seasonal solution, and it forces ``T̂_{L,0} = T̂_{O,0}``
  (``C``-independent): **the annual means of land and ocean are identical; continentality lives entirely
  in the seasonal amplitude.**

Validation triad (plan §3)
--------------------------
* **Tight (analytic/structural).** (a) The **0-D slab** limit — transport off (``D = 0``), a pure
  sinusoidal forcing ``F_0 + F_1 cos ωt`` → amplitude ``F_1/√(B² + ω²C²)`` and phase lag
  ``arctan(ωC/B)`` (:func:`slab_amplitude_lag`), reproduced by both solvers (pins the time integration and
  the ``C``-dependence — the mechanism of continentality). (b) The **reduction**: the spectral ``n = 0``
  harmonic equals the annual-mean :class:`~planet.sphere_ebm.SphereEBM` steady solve to machine precision,
  and ``⟨T_L⟩ = ⟨T_O⟩`` to the same. (c) **Marcher ≡ spectral** amplitude/phase to the stepping tolerance
  (the anti-damping cross-check). (d) **Hemispheric antisymmetry**: symmetric geometry ⟹
  ``T(x, t) = T(−x, t + ½ yr)``.
* **Conservation (tight).** Over one converged year the global-and-annual-mean net TOA
  ``⟨⟨S(1−α) − A − B T⟩⟩ ≈ 0`` (the transport conserves column energy each step).
* **Benchmark (loose / calibrated).** Continentality — land seasonal **range ≫** ocean, ocean **lags
  more** — is *banked in direction*; the *magnitude* rides the calibrated heat capacities
  (land ≈ the atmospheric column ``c_p p_s/g``; ocean ``+`` a ~50 m seasonal mixed layer,
  [[seasonal-ebm-source]]) and is reported only in the observed ballpark, the way rung 0's ice line was.

Named scope edges
-----------------
Fixed (ice-free) albedo — so the reduction is *exact* and the spectral solve *applies*; the seasonal
ice-albedo feedback (and its small-ice-cap instability) is the marcher's future, not built here. **Same
albedo on both tiles** — so continentality is unambiguously a **heat-capacity** effect (land/ocean
albedo contrast is an available knob, deliberately off). **Uniform land fraction** ``f_L`` (constant
``C_a``, exact energy conservation) — the latitude-varying land mask and the true ``T(φ,λ)`` map are
rung 5B.2. Transport ``D`` is the **atmospheric** heat transport shared over both tiles (the ocean's own
circulation is elsewhere on the staircase); a well-mixed-atmosphere closure, named.

Units — SI, climlab-conventional (W m⁻², °C, x = sin φ ∈ [−1, 1]); ``C`` in J m⁻² K⁻¹, time in seconds,
reported in days/months; amplitude = half the peak-to-peak range (K), phase lag in days.
Sources ([[seasonal-ebm-source]], extends [[ebm-radiation-source]] / [[obliquity-insolation-source]]):
North & Coakley 1979 (*J. Atmos. Sci.* — the seasonal EBM, land/ocean thermal inertia); North, Mengel &
Short 1983 (*JGR* — resolving the seasons and the continents); Hartmann *Global Physical Climatology*
(the atmospheric-column and mixed-layer heat capacities).
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np
from scipy.linalg import solve_banded

from engines.diffusion import Diffusion1D, grid_from_edges, Neumann
from .ebm import (
    A_OLR, B_OLR, D_TRANSPORT, S0_EARTH, RHO_WATER, CW_WATER,
    ALBEDO_A0, ALBEDO_A2, legendre_P2,
)
from .obliquity import daily_mean_insolation, OBLIQUITY_EARTH

# --------------------------------------------------------------------------- #
# Pinned time / heat-capacity constants ([[seasonal-ebm-source]]).
# The heat capacities are BUILT, not memorized: land = the atmospheric column heat capacity
# c_p·p_s/g (a textbook constant); ocean = that + a mixed layer ρ_w c_w h_ml. The mixed-layer /
# soil depths are the calibrated (loose) knobs; everything else is a physical constant.
# --------------------------------------------------------------------------- #
SECONDS_PER_YEAR = 365.25 * 86400.0          # s — the orbital period (circular-orbit year)
OMEGA_YEAR = 2.0 * math.pi / SECONDS_PER_YEAR  # rad s⁻¹ — the annual angular frequency

CP_AIR = 1004.0            # J kg⁻¹ K⁻¹ — specific heat of dry air at constant pressure (Hartmann)
P_SURF = 1.013e5           # Pa        — mean surface pressure
G_GRAV = 9.81              # m s⁻²     — gravitational acceleration
C_ATMOS = CP_AIR * P_SURF / G_GRAV           # J m⁻² K⁻¹ — the atmospheric column heat capacity (~1.04e7)

OCEAN_MIXED_DEPTH = 50.0   # m — seasonal ocean mixed-layer depth (Hartmann; the calibrated ocean knob)
LAND_SOIL_DEPTH = 2.0      # m water-equiv — the seasonal soil active layer (thermal penetration ~2 m; the calibrated land knob)
LAND_FRACTION = 0.30       # — global land fraction (Earth ≈ 0.29); uniform per the named scope edge


def land_heat_capacity(soil_depth: float = LAND_SOIL_DEPTH) -> float:
    """Land-column effective heat capacity ``C_L = C_atmos + ρ_w c_w·soil_depth`` (J m⁻² K⁻¹).

    Land responds to the seasons with ~the **atmospheric column** heat capacity (``c_p p_s/g``) plus a
    thin conducting soil skin — a *small* reservoir, so the land tile heats and cools fast and far. The
    soil depth is the calibrated land knob ([[seasonal-ebm-source]]); the atmospheric term is a textbook
    constant.
    """
    return C_ATMOS + RHO_WATER * CW_WATER * float(soil_depth)


def ocean_heat_capacity(mixed_depth: float = OCEAN_MIXED_DEPTH) -> float:
    """Ocean-column effective heat capacity ``C_O = C_atmos + ρ_w c_w·mixed_depth`` (J m⁻² K⁻¹).

    The seasonal ocean mixed layer (~50 m) is a *large* reservoir — ~15–20× the land column — so the
    ocean tile barely warms over the year and lags the sun by ~2 months. The mixed-layer depth is the
    calibrated ocean knob ([[seasonal-ebm-source]]).
    """
    return C_ATMOS + RHO_WATER * CW_WATER * float(mixed_depth)


def slab_amplitude_lag(C: float, F1: float, B: float = B_OLR,
                       omega: float = OMEGA_YEAR) -> tuple[float, float]:
    """0-D slab analytic response to ``F_0 + F_1 cos ωt``: ``(amplitude K, phase-lag days)``.

    A single column with **no transport** obeys ``C dT/dt = F_0 + F_1 cos ωt − A − B T``; its periodic
    solution is ``T = T̄ + [F_1/√(B²+ω²C²)]·cos(ωt − ψ)`` with ``tan ψ = ωC/B``. So the seasonal
    **amplitude** is ``F_1/√(B²+ω²C²)`` (K) and the temperature maximum **lags** the forcing maximum by
    ``ψ/ω`` seconds. A larger ``C`` (ocean) shrinks the amplitude and grows the lag toward the ¼-year
    (``ψ→π/2``) quadrature limit; a smaller ``C`` (land) gives a big, prompt swing. This is the exact
    mechanism of continentality — the tight anchor both solvers must reproduce with ``D = 0``.
    """
    denom = math.hypot(B, omega * C)
    amplitude = float(F1) / denom
    lag_days = math.atan2(omega * C, B) / omega / 86400.0
    return amplitude, lag_days


@dataclass(frozen=True)
class SeasonalClimate:
    """A converged annual limit cycle: land + ocean seasonal fields on ``(x, day)`` plus diagnostics.

    ``x`` the cell-center ``sin φ`` on [−1, 1]; ``days`` the time samples over one year (days, endpoint
    excluded). ``T_land`` / ``T_ocean`` the tile temperatures ``[n_x, n_day]`` (°C); ``T_mean`` the
    zonal-mean ``f_L T_land + f_O T_ocean``. ``land_fraction`` the ``f_L`` used. ``method`` is
    ``"march"`` or ``"spectral"``; ``converged`` / ``years`` record the marcher's spin-up (``years = 0``
    for the spectral solve). Plain arrays — the loose-coupling currency.
    """

    x: np.ndarray
    days: np.ndarray
    T_land: np.ndarray
    T_ocean: np.ndarray
    T_mean: np.ndarray
    land_fraction: float
    method: str
    converged: bool
    years: int

    def latitude_deg(self) -> np.ndarray:
        """Grid latitudes ``φ = asin(x)`` in degrees (south pole −90° → north pole +90°)."""
        return np.degrees(np.arcsin(np.clip(self.x, -1.0, 1.0)))

    def annual_mean(self, field: str = "mean") -> np.ndarray:
        """Year-average temperature profile (°C) of ``"land"``/``"ocean"``/``"mean"`` (the reduction)."""
        arr = {"land": self.T_land, "ocean": self.T_ocean, "mean": self.T_mean}[field]
        return arr.mean(axis=1)

    def amplitude(self, field: str = "mean") -> np.ndarray:
        """Seasonal amplitude ``(max − min)/2`` (K) per latitude — half the peak-to-peak swing."""
        arr = {"land": self.T_land, "ocean": self.T_ocean, "mean": self.T_mean}[field]
        return 0.5 * (arr.max(axis=1) - arr.min(axis=1))


def _time_of_max(series: np.ndarray, days: np.ndarray) -> float:
    """Day-of-year of a periodic series' maximum, parabola-refined around the peak sample.

    ``argmax`` gives the coarse peak; a 3-point parabolic fit through the neighbouring samples (wrapping
    the year) refines it to sub-sample resolution — enough to read the thermal-lag phase against the
    analytic slab value. Assumes uniform spacing (the seasonal ``days`` grid).
    """
    n = series.size
    k = int(np.argmax(series))
    y0, y1, y2 = series[(k - 1) % n], series[k], series[(k + 1) % n]
    denom = y0 - 2.0 * y1 + y2
    shift = 0.5 * (y0 - y2) / denom if denom != 0.0 else 0.0
    dt_day = days[1] - days[0]
    return float((days[k] + shift * dt_day) % (n * dt_day))


class SeasonalEBM:
    """The seasonal two-tile (land/ocean) EBM on ``x = sin φ ∈ [−1, 1]`` — see the module docstring.

    Holds the OLR (``A``, ``B``) and transport (``D``) constants, the land fraction ``f_L`` and the tile
    heat capacities ``C_L``/``C_O``, the tilt ``ε`` and solar constant ``S₀``, and the year discretization
    (``n_steps`` samples). Builds the transport operator ``L_T`` once (as the engine assembles it) and a
    :class:`~engines.diffusion.Diffusion1D` for the marcher's implicit transport substep (heat capacity
    ``C_a``). The absorbed forcing is the seasonal insolation × a fixed co-albedo; a caller may inject a
    custom ``absorbed[x, t]`` array (the slab anchor supplies a synthetic sinusoid).
    """

    def __init__(self, A: float = A_OLR, B: float = B_OLR, D: float = D_TRANSPORT,
                 land_fraction: float = LAND_FRACTION,
                 land_soil_depth: float = LAND_SOIL_DEPTH,
                 ocean_mixed_depth: float = OCEAN_MIXED_DEPTH,
                 obliquity_deg: float = OBLIQUITY_EARTH, S0: float = S0_EARTH,
                 n_cells: int = 180, n_steps: int = 360):
        if B <= 0.0:
            raise ValueError(f"B (OLR slope) must be positive, got {B}")
        if float(D) < 0.0:
            raise ValueError(f"D (transport) must be non-negative, got {D}")
        if not 0.0 <= float(land_fraction) <= 1.0:
            raise ValueError(f"land_fraction must be in [0, 1], got {land_fraction}")
        self.A, self.B, self.D = float(A), float(B), float(D)
        self.f_land = float(land_fraction)
        self.f_ocean = 1.0 - self.f_land
        self.C_land = land_heat_capacity(land_soil_depth)
        self.C_ocean = ocean_heat_capacity(ocean_mixed_depth)
        # The effective transport heat capacity that makes the per-tile energy-flux redistribution
        # EXACTLY reproduce the engine's zonal-mean step (docstring): 1/C_a = f_L/C_L + f_O/C_O.
        self.C_a = 1.0 / (self.f_land / self.C_land + self.f_ocean / self.C_ocean)
        self.obliquity_deg = float(obliquity_deg)
        self.S0 = float(S0)
        self.n_cells = int(n_cells)
        self.n_steps = int(n_steps)

        # Full sphere: x = sin φ on [-1, 1], insulated at BOTH poles, equator interior.
        self.grid = grid_from_edges(np.linspace(-1.0, 1.0, self.n_cells + 1))
        self.x = self.grid.centers
        self.phi = np.arcsin(np.clip(self.x, -1.0, 1.0))
        # The engine for the marcher's transport substep: (D/C_a)(1−x²) ⟹ C_a ∂T̄/∂t = L_T T̄.
        self._Dcells = (self.D / self.C_a) * (1.0 - self.x ** 2)
        # D = 0 is the transport-off slab limit (the 0-D anchor): the engine's harmonic-mean faces would
        # be 0/0, so skip the transport substep entirely rather than build a degenerate solver.
        self.solver = (Diffusion1D(self.grid, self._Dcells, Neumann(0.0), Neumann(0.0))
                       if self.D > 0.0 else None)
        # L_T = D·d/dx[(1−x²)d/dx], assembled exactly as the engine assembles it (harmonic-mean faces).
        self._LT = self._transport_tridiag()

        # Year discretization (endpoint excluded so day 0 and day 365 aren't double-counted).
        self.dt = SECONDS_PER_YEAR / self.n_steps
        self.day_seconds = np.arange(self.n_steps) * self.dt
        self.days = self.day_seconds / 86400.0
        # Phase the year so day 0 ≈ 1 January (NH winter solstice): δ(0) = −ε, rising to +ε at mid-year,
        # so the figure's "month of year" axis reads as a familiar Jan-start calendar (NH summer near
        # month 6–7). A uniform time-shift — it changes no amplitude, lag, or anchor, only the labels.
        self._sin_delta = math.sin(math.radians(self.obliquity_deg)) * (
            -np.cos(2.0 * math.pi * self.day_seconds / SECONDS_PER_YEAR))

    # -- forcing ----------------------------------------------------------- #
    def coalbedo(self, albedo: Optional[np.ndarray | float] = None) -> np.ndarray:
        """Co-albedo ``1 − α(x)`` per latitude (default: the ice-free EBM ``α = a0 + a2 P₂``, same on both tiles)."""
        if albedo is None:
            alb = ALBEDO_A0 + ALBEDO_A2 * legendre_P2(self.x)
        else:
            alb = np.broadcast_to(np.asarray(albedo, dtype=float), self.x.shape)
        return 1.0 - alb

    def insolation_series(self) -> np.ndarray:
        """Seasonal insolation ``S(x, t) = (S₀/π)·Q_kernel(φ, δ(t))`` (W m⁻²), shape ``[n_x, n_steps]``.

        The pinned daily-insolation kernel (:func:`planet.obliquity.daily_mean_insolation`) evaluated at
        the declination ``δ(t) = arcsin(sin ε·sin(2π t/yr))`` of each day, restored to absolute flux by the
        dropped ``S₀/π`` factor (global-annual mean ``= S₀/4``). This is the *incident* flux; multiply by
        the co-albedo for the absorbed forcing.
        """
        kernel = daily_mean_insolation(self.phi, self._sin_delta)        # [n_x, n_steps], relative
        return (self.S0 / math.pi) * kernel

    def absorbed_series(self, albedo: Optional[np.ndarray | float] = None) -> np.ndarray:
        """Absorbed shortwave ``S(x, t)(1 − α(x))`` (W m⁻²), shape ``[n_x, n_steps]`` — the forcing both solvers take."""
        return self.insolation_series() * self.coalbedo(albedo)[:, None]

    # -- the transport operator (reconstructed to match the engine) -------- #
    def _transport_tridiag(self):
        """Tridiagonals ``(sub, diag, sup)`` of ``L_T = D·d/dx[(1−x²)d/dx]`` (W m⁻² K⁻¹), engine-exact.

        Harmonic-mean interior faces of ``D(1−x²)`` over the center spacing, insulated (Neumann 0) ends —
        the identical assembly :meth:`planet.sphere_ebm.SphereEBM._transport_tridiag` uses, so ``L_T``
        cannot drift from the engine that marches it.
        """
        dx = self.grid.widths
        Dc = self.D * (1.0 - self.x ** 2)
        s = Dc[:-1] + Dc[1:]
        Dface = np.divide(2.0 * Dc[:-1] * Dc[1:], s, out=np.zeros_like(s), where=s > 0.0)
        Tt = Dface / np.diff(self.grid.centers)
        n = self.grid.n
        sub, diag, sup = np.zeros(n), np.zeros(n), np.zeros(n)
        sup[:-1] += Tt / dx[:-1]; diag[:-1] += -Tt / dx[:-1]
        sub[1:] += Tt / dx[1:];   diag[1:] += -Tt / dx[1:]
        return sub, diag, sup

    def _apply_LT(self, v: np.ndarray) -> np.ndarray:
        """``L_T · v`` (the transport convergence, W m⁻²) for a real or complex profile ``v``."""
        sub, diag, sup = self._LT
        out = diag * v
        out[:-1] += sup[:-1] * v[1:]
        out[1:] += sub[1:] * v[:-1]
        return out

    # -- the spectral limit cycle (the tight, splitting-free reference) ---- #
    def spectral(self, albedo: Optional[np.ndarray | float] = None,
                 absorbed: Optional[np.ndarray] = None) -> SeasonalClimate:
        """Exact frequency-domain limit cycle (no time-stepping error) — the tight amplitude/phase anchor.

        Solves, per temporal harmonic ``n``, the complex banded system ``(I − G_n L_T) T̄̂_n = G_n F̂_n``
        with ``F = absorbed − A`` and ``G_n = Σ_i f_i/(iω_n C_i + B)``; recovers the tile spectra and
        inverse-transforms. The ``n = 0`` harmonic is the annual-mean EBM (docstring). A caller may pass a
        precomputed ``absorbed[x, t]`` (e.g. the slab anchor's synthetic sinusoid); otherwise it is built
        from the seasonal insolation × co-albedo. Requires state-independent (fixed-albedo) forcing.
        """
        F = (self.absorbed_series(albedo) if absorbed is None else np.asarray(absorbed, float)) - self.A
        Fhat = np.fft.rfft(F, axis=1)                              # [n_x, n_freq], complex
        freqs = np.fft.rfftfreq(self.n_steps, d=self.dt)           # cycles s⁻¹
        omegas = 2.0 * math.pi * freqs
        n = self.grid.n
        sub, diag, sup = self._LT
        TL_hat = np.zeros_like(Fhat)
        TO_hat = np.zeros_like(Fhat)
        for k, omega in enumerate(omegas):
            gL = 1.0 / (1j * omega * self.C_land + self.B)
            gO = 1.0 / (1j * omega * self.C_ocean + self.B)
            Gk = self.f_land * gL + self.f_ocean * gO
            # (I − G_k L_T): banded complex, ab rows [sup; diag; sub].
            ab = np.zeros((3, n), dtype=complex)
            ab[0, 1:] = -Gk * sup[:-1]
            ab[1, :] = 1.0 - Gk * diag
            ab[2, :-1] = -Gk * sub[1:]
            Tbar_k = solve_banded((1, 1), ab, Gk * Fhat[:, k])
            forcing_k = self._apply_LT(Tbar_k) + Fhat[:, k]        # L_T T̄̂ + F̂
            TL_hat[:, k] = forcing_k * gL
            TO_hat[:, k] = forcing_k * gO
        T_land = np.fft.irfft(TL_hat, n=self.n_steps, axis=1)
        T_ocean = np.fft.irfft(TO_hat, n=self.n_steps, axis=1)
        T_mean = self.f_land * T_land + self.f_ocean * T_ocean
        return SeasonalClimate(self.x, self.days, T_land, T_ocean, T_mean,
                               self.f_land, "spectral", True, 0)

    # -- the time-marcher (the engine-reuse method) ------------------------ #
    def march(self, albedo: Optional[np.ndarray | float] = None,
              absorbed: Optional[np.ndarray] = None, T_init: Optional[float] = None,
              tol: float = 1e-6, max_years: int = 80) -> SeasonalClimate:
        """March the Strang-split model to a converged annual limit cycle; return the last year's fields.

        Each substep is half-radiation (per tile, analytic exact for the frozen-α linear sink) /
        full-transport (the engine, on ``T̄`` with ``C_a``, redistributed to the tiles as an energy flux) /
        half-radiation. Runs whole years, comparing each year's day-0 state to the previous year's; stops
        when ``max|ΔT| < tol`` (K) — the limit cycle. Seeded from a uniform ``T_init`` (or 15 °C); the
        limit cycle is unique for the linear system, so the seed only sets the spin-up length (ocean
        ``τ = C_O/B`` is a few years — hence the tens-of-years cap).
        """
        A = (self.absorbed_series(albedo) if absorbed is None else np.asarray(absorbed, float))
        if T_init is None:
            # Seed each tile at its per-latitude annual-mean radiative equilibrium ⟨S(1−α)−A⟩/B: the
            # correct annual mean (independent of the spectral solve), so only the seasonal anomaly — not
            # the slow global-mean offset — has to spin up, cutting the ocean τ=C_O/B transient.
            TL = (A.mean(axis=1) - self.A) / self.B
            TO = TL.copy()
        else:
            TL = np.full(self.x.shape, float(T_init))
            TO = np.full(self.x.shape, float(T_init))
        decayL = math.exp(-0.5 * self.dt * self.B / self.C_land)
        decayO = math.exp(-0.5 * self.dt * self.B / self.C_ocean)

        def rad_half(TL, TO, absorbed_t):
            TeqL = (absorbed_t - self.A) / self.B
            TeqO = TeqL
            return TeqL + (TL - TeqL) * decayL, TeqO + (TO - TeqO) * decayO

        converged, year = False, 0
        # Store the running year so the *last* year is what we return.
        TL_year = np.zeros((self.n_cells, self.n_steps))
        TO_year = np.zeros((self.n_cells, self.n_steps))
        for year in range(1, max_years + 1):
            TL_ref, TO_ref = TL.copy(), TO.copy()
            for s in range(self.n_steps):
                TL_year[:, s] = TL
                TO_year[:, s] = TO
                a_now = A[:, s]
                a_next = A[:, (s + 1) % self.n_steps]
                TL, TO = rad_half(TL, TO, a_now)                  # half radiation @ t
                if self.solver is not None:                       # transport (skipped in the D=0 slab limit)
                    Tbar = self.f_land * TL + self.f_ocean * TO
                    dTbar = self.solver.step(Tbar, self.dt) - Tbar  # implicit transport on the mean
                    TL = TL + (self.C_a / self.C_land) * dTbar    # energy-flux redistribution
                    TO = TO + (self.C_a / self.C_ocean) * dTbar
                TL, TO = rad_half(TL, TO, a_next)                 # half radiation @ t+dt
            if max(np.max(np.abs(TL - TL_ref)), np.max(np.abs(TO - TO_ref))) < tol:
                converged = True
                break
        T_mean = self.f_land * TL_year + self.f_ocean * TO_year
        return SeasonalClimate(self.x, self.days, TL_year, TO_year, T_mean,
                               self.f_land, "march", converged, year)

    # -- diagnostics ------------------------------------------------------- #
    def phase_lag_days(self, T_series: np.ndarray) -> np.ndarray:
        """Thermal lag (days) of each latitude's temperature max behind its insolation max.

        The seasonal signature of heat capacity: ``time_of_max(T) − time_of_max(S)`` (mod year) per
        latitude, both parabola-refined (:func:`_time_of_max`). Read the ocean tile's midlatitude value
        against the ~2-month observed SST lag. Non-sinusoidal polar-night latitudes are meaningless here
        (the max-based read is only physical where the cycle is roughly single-peaked) — keep to
        ``|φ| ≲ 60°``.
        """
        S = self.insolation_series()
        lags = np.empty(self.n_cells)
        year_days = self.n_steps * (self.days[1] - self.days[0])
        for i in range(self.n_cells):
            lag = _time_of_max(T_series[i], self.days) - _time_of_max(S[i], self.days)
            lags[i] = (lag + 0.5 * year_days) % year_days - 0.5 * year_days
        return lags

    def nearest_index(self, lat_deg: float) -> int:
        """Grid index nearest latitude ``lat_deg`` (degrees) — for reading a single band (e.g. 45°)."""
        return int(np.argmin(np.abs(self.latitude_deg() - float(lat_deg))))

    def latitude_deg(self) -> np.ndarray:
        """Grid latitudes ``φ = asin(x)`` in degrees (south −90° → north +90°)."""
        return np.degrees(np.arcsin(np.clip(self.x, -1.0, 1.0)))
