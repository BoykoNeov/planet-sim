"""The **two-dimensional** seasonal energy-balance model — ``T(φ, λ, t)`` over land and sea (rung 5B.2).

Rung 5B.1 (:mod:`planet.seasonal`) woke the heat capacity ``C`` by marching the seasons, and cast
**continentality** — the large seasonal temperature range of continental interiors — from the land/ocean
``C`` contrast. But it carried the continents only as a **zonal-mean caricature**: two heat-capacity
*tiles* per latitude (a land tile, an ocean tile) sharing one meridional transport, giving a single
number per latitude rather than a *map*. "*2-D EBM*" is a term of art for **two spatial dimensions**
(latitude × longitude; North, Mengel & Short 1983, *"…resolving the seasons and the continents"*), and
this module builds that: a **single** temperature field ``T(φ, λ, t)`` on the sphere, with a
spatially-varying heat capacity ``C(φ, λ)`` set by a **real land mask**, marched under the axial-tilt
seasons to a converged annual limit cycle. Where the mask is land, ``C`` is small and the point swings
hard and prompt; where it is ocean, ``C`` is large and the point barely moves and lags — and now the
transport *diffuses heat between neighbouring longitudes*, so a continental **interior** swings more than
its **coasts** (which the adjacent ocean moderates). Continentality becomes a **map**, the NMS83 result.

The equation — a single field, a mask-set heat capacity, the sphere Laplacian
-----------------------------------------------------------------------------
On ``x = sin φ ∈ [−1, 1]`` and longitude ``λ ∈ [0, 2π)`` (both cell-centered), the model is

    C(φ, λ) ∂T/∂t = D[ ∂ₓ((1 − x²) ∂ₓT) + (1/(1 − x²)) ∂²T/∂λ² ] + S(x, t)(1 − α) − (A + B T).

The **meridional** operator ``Lₓ = D ∂ₓ((1 − x²) ∂ₓ·)`` is *exactly* rung 5B.1 / :mod:`planet.sphere_ebm`'s
1-D transport (the ``(1 − x²)`` weight vanishing at both poles). The **zonal** operator
``L_λ = (D/(1 − x²)) ∂²/∂λ²`` is a *periodic* 1-D diffusion along each latitude circle, its coefficient
growing toward the poles — the **meridian-convergence** wrinkle: near a pole all longitudes crowd onto
one point, so zonal mixing there is near-instantaneous. The cell-centered grid keeps every ``1 − x²`` off
the singular endpoints (no cell sits *at* a pole), so the polar coefficient is large-but-finite, and the
unconditionally-stable backward-Euler sweep simply relaxes each polar ring toward the zonal mean it should
physically have. Insolation and albedo are **longitude-independent** here (same named scope as 5B.1: fixed
ice-free albedo, so continentality is a *pure* heat-capacity effect) — every drop of zonal structure comes
from the mask.

The method — ADI operator-splitting REUSES the 1-D tridiagonal engine (a bounded step, not a rewrite)
-----------------------------------------------------------------------------------------------------
This is the plan's named method. Each time step is split into one-dimensional implicit **sweeps**, so the
2-D solve is a sequence of 1-D tridiagonal solves — the engine spine (:mod:`engines.diffusion`) reused one
axis at a time, never a 2-D rewrite. The Strang-ordered step (matching 5B.1's split, so the reductions
below are exact) is **½ radiation → meridional sweep → zonal sweep → ½ radiation**:

* **radiation** — analytic per cell (the linear sink is exact): ``T ← T_eq + (T − T_eq)·exp(−½ΔtB/C)``
  with ``T_eq(x, t) = (S(1 − α) − A)/B``. Each cell relaxes on *its own* ``C(φ, λ)`` timescale — the
  seat of continentality.
* **meridional sweep** — for every longitude column, a backward-Euler solve of ``C ∂T/∂t = Lₓ T``:
  ``(diag(C) − Δt Lₓ) T* = C·Tⁿ``. ``Lₓ`` is reconstructed from the engine's *exact* harmonic-mean-face
  assembly (:attr:`planet.seasonal.SeasonalEBM._LT`, itself the engine's), and the columns are solved
  together by a batched Thomas sweep. The per-cell ``C`` sits on the diagonal — the varying-``C`` wrinkle
  the engine's uniform-``C`` ``step()`` cannot swallow, so we reconstruct the operator exactly as
  :mod:`planet.sphere_ebm` does for its dt-free solve and add the ``C/Δt`` diagonal ourselves.
* **zonal sweep** — for every latitude row, the **periodic** counterpart ``(diag(C) − Δt L_λ) T* = C·Tⁿ``,
  a *cyclic* tridiagonal (the wrap-around couples the last longitude back to the first) solved by a batched
  Sherman–Morrison correction of the same Thomas sweep. This is the one genuinely new numerical object; it
  is pinned machine-tight by the circulant-eigenmode anchor (``cos mλ`` is an exact eigenvector).

Each backward-Euler sweep conserves the column energy ``∫ C T dA`` exactly (the finite-volume no-flux /
periodic operators sum to zero, and the uniform area element ``dA = a² dx dλ`` in these coordinates makes
the area weighting a constant that factors out) — so the split marches to the correct global balance.

Two reductions make this honest (both machine-tight)
----------------------------------------------------
* **Down to rung 5B.1.** A **zonally-uniform** mask (all-land, all-ocean, or latitude bands) has no
  longitude structure, so the field stays zonally flat, the zonal sweep is the identity, and the model
  collapses **bit-for-bit** to the 5B.1 single-field marcher (all-land ≡ ``land_fraction=1``, all-ocean ≡
  ``0``). 5B.1's marcher is itself validated against its exact spectral solve — so this transitively
  inherits the anti-damping guarantee for the meridional + time integration, and the new zonal sweep is
  pinned separately by its eigenmode anchor.
* **Down to the annual-mean parent — the NMS83 headline.** Average the equation over a converged year:
  ``⟨C ∂T/∂t⟩ = 0`` (periodic), so ``⟨T⟩`` solves the *annual-mean* EBM with **longitude-independent**
  forcing ``⟨S⟩(1 − α) − A`` — whose solution is therefore **zonally uniform** and equals the 1-D parent
  :meth:`planet.sphere_ebm.SphereEBM.steady_linear`, *for any mask*. This is 5B.1's ``⟨T_L⟩ = ⟨T_O⟩``
  generalized to the map: **the land/sea contrast lives entirely in the seasonal amplitude — the annual
  mean is blind to the mask.** So the continentality map is a map of *seasonal range*, laid over an
  annual-mean climate that is as smooth as the aquaplanet's.

Named scope edges (carried from 5B.1)
-------------------------------------
Fixed ice-free albedo, identical over land and sea (continentality is *pure* heat capacity — the
land/sea albedo contrast is a deliberately-off knob; the seasonal ice-albedo feedback is the marcher's
future). No zonal insolation or advective (wind-driven) asymmetry — this is **diffusive** continentality
(interior extremes + coastal moderation, symmetric about a continent), not a downwind maritime/continental
tilt (that needs the atmosphere's mean winds, elsewhere on the staircase). The land mask is prescribed
geography, not an emergent coastline.

Units — SI, climlab-conventional (W m⁻², °C, ``x = sin φ ∈ [−1, 1]``, ``λ`` in radians / reported degrees);
``C`` in J m⁻² K⁻¹, time in seconds, seasonal range/amplitude in K, lag in days. Sources
([[seasonal-ebm-source]], extends [[ebm-radiation-source]] / [[obliquity-insolation-source]]):
North, Mengel & Short 1983 (*JGR* 88, C11 — the 2-D seasons-and-continents EBM); North & Coakley 1979
(the seasonal land/ocean thermal-inertia model 5B.1 built). See [[planet-rung5b-seasonal]].
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional

import numpy as np

from .seasonal import SeasonalEBM, LAND_FRACTION, LAND_SOIL_DEPTH, OCEAN_MIXED_DEPTH
from .ebm import A_OLR, B_OLR, D_TRANSPORT, S0_EARTH
from .obliquity import OBLIQUITY_EARTH


# --------------------------------------------------------------------------- #
# Batched 1-D implicit solvers — the ADI sweeps (reuse the engine's assembly).
# Each solves MANY independent tridiagonal systems that share off-diagonals but
# carry a per-system main diagonal (the spatially-varying heat capacity).
# --------------------------------------------------------------------------- #
def _thomas_columns(lower: np.ndarray, diag: np.ndarray, upper: np.ndarray,
                    rhs: np.ndarray) -> np.ndarray:
    """Solve a tridiagonal system down **axis 0** for every column, batched (the Thomas algorithm).

    ``lower``/``upper`` are length-``n`` 1-D arrays shared by all columns (``lower[0]`` and ``upper[-1]``
    unused); ``diag`` and ``rhs`` are ``[n, m]`` (a distinct main diagonal per column — the meridional
    sweep's per-column heat capacity). Returns the ``[n, m]`` solution. Pure-NumPy forward-elimination /
    back-substitution, vectorized over the ``m`` columns — the meridional ADI sweep, one tridiagonal solve
    per longitude done at once.
    """
    n = diag.shape[0]
    cp = np.empty_like(diag)                       # modified super-diagonal
    dp = np.empty_like(rhs)                        # modified rhs
    cp[0] = upper[0] / diag[0]
    dp[0] = rhs[0] / diag[0]
    for i in range(1, n):
        denom = diag[i] - lower[i] * cp[i - 1]
        cp[i] = upper[i] / denom
        dp[i] = (rhs[i] - lower[i] * dp[i - 1]) / denom
    x = np.empty_like(rhs)
    x[-1] = dp[-1]
    for i in range(n - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]
    return x


def _cyclic_thomas_rows(offdiag: np.ndarray, diag: np.ndarray, rhs: np.ndarray) -> np.ndarray:
    """Solve a **periodic** (cyclic) tridiagonal system along **axis 1** for every row, batched.

    Each latitude row's zonal operator is a *constant-coefficient* periodic Laplacian: the sub-, super-,
    and both wrap-around corner entries all equal ``offdiag[row]`` (a per-row scalar, ``[n_rows]``), while
    the main ``diag`` (``[n_rows, m]``) carries the per-cell heat capacity. Sherman–Morrison (Numerical
    Recipes' ``cyclic``) reduces the cyclic solve to two ordinary tridiagonal solves plus a rank-1 update,
    here batched over rows. The constant vector is in the operator's null space, so a flat row is returned
    unchanged — which is what makes the zonal sweep the identity on zonally-uniform fields.
    """
    m = diag.shape[1]
    e = offdiag[:, None]                                   # [n_rows, 1] broadcast over longitude
    lower = np.broadcast_to(e, diag.shape)
    upper = np.broadcast_to(e, diag.shape)
    gamma = -diag[:, 0]                                    # NR's γ = −b₀ (avoids cancellation)
    bmod = diag.copy()
    bmod[:, 0] = diag[:, 0] - gamma
    bmod[:, -1] = diag[:, -1] - e[:, 0] * e[:, 0] / gamma
    # Two tridiagonal solves (non-cyclic), along axis 1 → transpose into the axis-0 Thomas kernel.
    def tri(rhs_):
        return _thomas_columns(lower.T, bmod.T, upper.T, rhs_.T).T
    u = np.zeros_like(diag)
    u[:, 0] = gamma
    u[:, -1] = e[:, 0]
    y = tri(rhs)
    z = tri(u)
    # v = (1, 0, …, 0, e/γ);  fact = (v·y)/(1 + v·z)
    vy = y[:, 0] + (e[:, 0] / gamma) * y[:, -1]
    vz = z[:, 0] + (e[:, 0] / gamma) * z[:, -1]
    fact = vy / (1.0 + vz)
    return y - fact[:, None] * z


# --------------------------------------------------------------------------- #
# Land masks — prescribed geography on the (x, λ) grid.
# --------------------------------------------------------------------------- #
def uniform_mask(n_cells: int, n_lon: int, land: bool = True) -> np.ndarray:
    """An all-land (``land=True``) or all-ocean mask, ``[n_cells, n_lon]`` — the reduction end-members."""
    return np.full((n_cells, n_lon), bool(land))


def zonal_band_mask(x: np.ndarray, n_lon: int, bands_deg) -> np.ndarray:
    """A **zonally-uniform** land mask: land wherever ``|φ|`` (or ``φ``) falls in one of ``bands_deg``.

    ``bands_deg`` is a list of ``(lat_min, lat_max)`` pairs (degrees, signed latitude). Every longitude at
    a given latitude is the same — so the mask has no longitude structure and the model reduces to a 1-D
    seasonal EBM with a latitude-varying ``C`` (a reduction anchor, and the closest 2-D grid can come to
    5B.1's zonal caricature).
    """
    lat = np.degrees(np.arcsin(np.clip(x, -1.0, 1.0)))
    is_land = np.zeros(x.shape, dtype=bool)
    for lo, hi in bands_deg:
        is_land |= (lat >= lo) & (lat <= hi)
    return np.broadcast_to(is_land[:, None], (x.size, n_lon)).copy()


def box_mask(x: np.ndarray, lon_rad: np.ndarray, lat_deg, lon_deg) -> np.ndarray:
    """A rectangular continent: land where latitude ∈ ``lat_deg`` **and** longitude ∈ ``lon_deg``.

    ``lat_deg`` / ``lon_deg`` are ``(min, max)`` pairs in degrees (longitude on ``[0, 360)``). The
    idealized single-continent demonstrator — a clean block whose interior shows the extreme seasonal range
    and whose edges are moderated by the surrounding ocean.
    """
    lat = np.degrees(np.arcsin(np.clip(x, -1.0, 1.0)))
    lon = np.degrees(lon_rad) % 360.0
    lat_ok = (lat >= lat_deg[0]) & (lat <= lat_deg[1])
    lon_ok = (lon >= lon_deg[0]) & (lon <= lon_deg[1])
    return lat_ok[:, None] & lon_ok[None, :]


def earthlike_mask(x: np.ndarray, lon_rad: np.ndarray) -> np.ndarray:
    """A coarse idealized-Earth land mask (a few blocky continents) — the NMS83-style *map* demonstrator.

    A deliberately schematic arrangement (a broad Eurasia-like NH landmass with a deep interior, a narrow
    meridional Americas strip, an equatorial Africa block, and an Australia block), enough to show the
    map's headline: interiors of the wide continent reach Siberian-style seasonal extremes while narrow or
    coastal land stays maritime. Not a real coastline — the mask is a plain boolean array, so a real one is
    a drop-in replacement.
    """
    m = np.zeros((x.size, lon_rad.size), dtype=bool)
    m |= box_mask(x, lon_rad, (8.0, 78.0), (35.0, 145.0))     # broad Eurasia (deep interior)
    m |= box_mask(x, lon_rad, (-38.0, 12.0), (10.0, 42.0))    # Africa (equatorial)
    m |= box_mask(x, lon_rad, (-55.0, 62.0), (270.0, 310.0))  # the Americas (a narrow meridional strip)
    m |= box_mask(x, lon_rad, (-38.0, -12.0), (200.0, 240.0)) # Australia
    return m


# --------------------------------------------------------------------------- #
# The frozen result — plain arrays (the loose-coupling currency).
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class SeasonalMapClimate:
    """A converged 2-D annual limit cycle ``T(x, λ, day)`` plus its continentality diagnostics.

    ``x`` cell-center ``sin φ`` on [−1, 1]; ``lon`` cell-center longitude (radians); ``days`` the time
    samples over one year. ``T`` the field ``[n_x, n_lon, n_day]`` (°C); ``land_mask`` the ``[n_x, n_lon]``
    boolean geography. ``method`` is ``"march"``; ``converged``/``years`` record the spin-up. Plain arrays.
    """

    x: np.ndarray
    lon: np.ndarray
    days: np.ndarray
    T: np.ndarray
    land_mask: np.ndarray
    method: str
    converged: bool
    years: int

    def latitude_deg(self) -> np.ndarray:
        """Grid latitudes ``φ = asin(x)`` in degrees (south −90° → north +90°)."""
        return np.degrees(np.arcsin(np.clip(self.x, -1.0, 1.0)))

    def longitude_deg(self) -> np.ndarray:
        """Grid longitudes in degrees on ``[0, 360)``."""
        return np.degrees(self.lon) % 360.0

    def annual_mean(self) -> np.ndarray:
        """Year-average temperature map ``⟨T⟩(x, λ)`` (°C) — the reduction (should be ~zonally uniform)."""
        return self.T.mean(axis=2)

    def seasonal_range(self) -> np.ndarray:
        """Peak-to-peak seasonal range ``max_t T − min_t T`` (K) per grid point — the continentality map."""
        return self.T.max(axis=2) - self.T.min(axis=2)

    def amplitude(self) -> np.ndarray:
        """Seasonal amplitude ``(max − min)/2`` (K) per grid point — half the peak-to-peak range."""
        return 0.5 * self.seasonal_range()


# --------------------------------------------------------------------------- #
# The 2-D seasonal EBM.
# --------------------------------------------------------------------------- #
class SeasonalMapEBM:
    """The 2-D (lat × lon) seasonal EBM on ``x = sin φ ∈ [−1, 1]`` × ``λ ∈ [0, 2π)`` — see the module docstring.

    Wraps a rung-5B.1 :class:`~planet.seasonal.SeasonalEBM` as the **source of truth** for the latitude
    grid, the meridional transport operator ``Lₓ``, the seasonal insolation ``S(x, t)``, the co-albedo, and
    the time discretization — so the reduction to 5B.1 is *bit-identical*, not merely close. It adds the
    longitude axis, the mask-set heat-capacity field ``C(x, λ)``, and the ADI sweeps. The mask is a
    ``[n_cells, n_lon]`` boolean array (``True`` = land); build one with :func:`box_mask`,
    :func:`zonal_band_mask`, :func:`uniform_mask`, or :func:`earthlike_mask`.
    """

    def __init__(self, land_mask: Optional[np.ndarray] = None, A: float = A_OLR, B: float = B_OLR,
                 D: float = D_TRANSPORT, land_soil_depth: float = LAND_SOIL_DEPTH,
                 ocean_mixed_depth: float = OCEAN_MIXED_DEPTH, obliquity_deg: float = OBLIQUITY_EARTH,
                 S0: float = S0_EARTH, n_cells: int = 90, n_lon: int = 128, n_steps: int = 360):
        # The 5B.1 model carries the grid, operator, insolation, co-albedo, dt — reused verbatim (the
        # land_fraction it needs is irrelevant here: we take its per-tile C_land/C_ocean, not its C_a).
        self.zonal = SeasonalEBM(A=A, B=B, D=D, land_fraction=LAND_FRACTION,
                                 land_soil_depth=land_soil_depth, ocean_mixed_depth=ocean_mixed_depth,
                                 obliquity_deg=obliquity_deg, S0=S0, n_cells=n_cells, n_steps=n_steps)
        self.A, self.B, self.D = self.zonal.A, self.zonal.B, self.zonal.D
        self.C_land, self.C_ocean = self.zonal.C_land, self.zonal.C_ocean
        self.n_cells, self.n_steps = int(n_cells), int(n_steps)
        self.n_lon = int(n_lon)

        self.x = self.zonal.x
        self.grid = self.zonal.grid
        self.dt = self.zonal.dt
        self.days = self.zonal.days
        # Longitude: cell centers on [0, 2π), uniform, periodic (no cell at the seam).
        self.dlon = 2.0 * math.pi / self.n_lon
        self.lon = (np.arange(self.n_lon) + 0.5) * self.dlon

        if land_mask is None:
            land_mask = uniform_mask(self.n_cells, self.n_lon, land=True)
        self.land_mask = np.asarray(land_mask, dtype=bool)
        if self.land_mask.shape != (self.n_cells, self.n_lon):
            raise ValueError(f"land_mask must be shape {(self.n_cells, self.n_lon)}, "
                             f"got {self.land_mask.shape}")
        # C(x, λ): small over land, large over ocean — the whole seat of continentality.
        self.C = np.where(self.land_mask, self.C_land, self.C_ocean)   # [n_cells, n_lon]

        # Meridional operator Lₓ (W m⁻² K⁻¹), the engine's exact harmonic-mean-face tridiagonal (reused).
        self._LT_sub, self._LT_diag, self._LT_sup = self.zonal._LT
        # Zonal coefficient a(x) = D/((1−x²)·Δλ²) (W m⁻² K⁻¹ per unit T) — the per-row off-diagonal of L_λ.
        # 1−x² never hits 0 (cell-centered), but grows the polar rows' zonal coupling (meridian convergence).
        self._zonal_a = self.D / ((1.0 - self.x ** 2) * self.dlon ** 2)   # [n_cells]

    # -- forcing (longitude-independent — reused from 5B.1) ---------------- #
    def absorbed_series(self, albedo=None) -> np.ndarray:
        """Absorbed shortwave ``S(x, t)(1 − α(x))`` (W m⁻²), ``[n_x, n_steps]`` — same at every longitude."""
        return self.zonal.absorbed_series(albedo)

    def _radiative_equilibrium(self, absorbed: np.ndarray) -> np.ndarray:
        """Per-latitude annual-mean radiative-equilibrium seed ``⟨S(1−α)−A⟩/B`` (°C), broadcast over λ."""
        Teq = (absorbed.mean(axis=1) - self.A) / self.B          # [n_x]
        return np.broadcast_to(Teq[:, None], (self.n_cells, self.n_lon)).copy()

    # -- the ADI sweeps ---------------------------------------------------- #
    def _meridional_sweep(self, T: np.ndarray) -> np.ndarray:
        """Backward-Euler ``(diag(C) − Δt Lₓ) T* = C·T`` down each longitude column (batched Thomas)."""
        lower = -self.dt * self._LT_sub
        upper = -self.dt * self._LT_sup
        diag = self.C - self.dt * self._LT_diag[:, None]          # [n_cells, n_lon]
        return _thomas_columns(lower, diag, upper, self.C * T)

    def _zonal_sweep(self, T: np.ndarray) -> np.ndarray:
        """Backward-Euler ``(diag(C) − Δt L_λ) T* = C·T`` around each latitude row (batched cyclic Thomas)."""
        offdiag = -self.dt * self._zonal_a                       # [n_cells], per-row constant off-diagonal
        diag = self.C + 2.0 * self.dt * self._zonal_a[:, None]   # [n_cells, n_lon]
        return _cyclic_thomas_rows(offdiag, diag, self.C * T)

    # -- the time-marcher (the engine-reuse method) ------------------------ #
    def march(self, albedo=None, absorbed: Optional[np.ndarray] = None, T_init: Optional[float] = None,
              tol: float = 1e-6, max_years: int = 60) -> SeasonalMapClimate:
        """March the split model to a converged annual limit cycle; return the last year's field.

        Each step is ½ radiation / meridional sweep / zonal sweep / ½ radiation (the module docstring's
        split). Runs whole years, comparing each year's day-0 state to the previous year's; stops when
        ``max|ΔT| < tol`` (K). Seeded from each latitude's annual-mean radiative equilibrium (so only the
        seasonal anomaly, not the slow global-mean offset, has to spin up — the ocean ``τ = C_O/B`` is a
        few years). ``absorbed[x, t]`` may be injected (the slab / synthetic-forcing anchors).
        """
        absorbed = (self.absorbed_series(albedo) if absorbed is None
                    else np.asarray(absorbed, dtype=float))       # [n_x, n_steps]
        if T_init is None:
            T = self._radiative_equilibrium(absorbed)
        else:
            T = np.full((self.n_cells, self.n_lon), float(T_init))
        # Per-cell half-step radiation decay factor (its own C at every point).
        decay = np.exp(-0.5 * self.dt * self.B / self.C)          # [n_cells, n_lon]

        def rad_half(T, absorbed_t):
            Teq = ((absorbed_t - self.A) / self.B)[:, None]       # [n_x, 1] broadcast over λ
            return Teq + (T - Teq) * decay

        converged, year = False, 0
        T_year = np.zeros((self.n_cells, self.n_lon, self.n_steps))
        for year in range(1, max_years + 1):
            T_ref = T.copy()
            for s in range(self.n_steps):
                T_year[:, :, s] = T
                a_now = absorbed[:, s]
                a_next = absorbed[:, (s + 1) % self.n_steps]
                T = rad_half(T, a_now)                            # ½ radiation @ t
                T = self._meridional_sweep(T)                     # implicit meridional transport
                T = self._zonal_sweep(T)                          # implicit periodic zonal transport
                T = rad_half(T, a_next)                           # ½ radiation @ t+dt
            if np.max(np.abs(T - T_ref)) < tol:
                converged = True
                break
        return SeasonalMapClimate(self.x, self.lon, self.days, T_year, self.land_mask,
                                  "march", converged, year)

    # -- diagnostics ------------------------------------------------------- #
    def latitude_deg(self) -> np.ndarray:
        """Grid latitudes ``φ = asin(x)`` in degrees (south −90° → north +90°)."""
        return np.degrees(np.arcsin(np.clip(self.x, -1.0, 1.0)))

    def longitude_deg(self) -> np.ndarray:
        """Grid longitudes in degrees on ``[0, 360)``."""
        return np.degrees(self.lon) % 360.0

    def nearest_index(self, lat_deg: float) -> int:
        """Latitude-row index nearest ``lat_deg`` (degrees)."""
        return int(np.argmin(np.abs(self.latitude_deg() - float(lat_deg))))

    def nearest_lon_index(self, lon_deg: float) -> int:
        """Longitude-column index nearest ``lon_deg`` (degrees, [0, 360))."""
        return int(np.argmin(np.abs(self.longitude_deg() - (float(lon_deg) % 360.0))))
