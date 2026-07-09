"""Orographic precipitation — the linear mountain-wave rain-shadow model (Rung 5A, plan §12.5).

This is the **first step off the zonal mean** (the "north star", plan §5/§12.5): a *2-D
longitude×latitude* precipitation field that depends on where the **mountains** are, not just on
latitude. It is the physics that makes the windward Cascades wet and the Columbia Basin behind them a
desert, the Western Ghats a monsoon wall and the Deccan a rain shadow — **regional** climate, set by
**geography**, that the zonal-mean :mod:`planet.precip` cannot express (it rains the same at every
longitude of a given latitude).

Concretely it implements the **Smith & Barstad (2004) Linear Theory of orographic precipitation** — a
*diagnostic* on a **prescribed** background wind and a 2-D terrain. Air forced up a mountain cools,
condenses, and rains on the windward slope; the descending, dried air on the lee slope leaves a rain
**shadow**. Smith & Barstad linearize this whole moist mountain-wave-plus-microphysics chain, so the
precipitation field is a single **wavenumber-space transfer function** applied to the terrain — one
FFT, laptop-trivial, with **exact analytic limits** to anchor it against.

What this is, and is honestly NOT (the honesty flag, plan §3/§9.3)
------------------------------------------------------------------
This module **wakes the dormant elevation seam.** Since v1 the geography spec has carried an
``elevation`` layer tagged ``inert=True`` — imported, displayed, round-tripped, but climate-inert
(:mod:`planet.planet_spec`, §9.3). This is the module the plan named as the one that finally makes
elevation *do something*.

It is, deliberately and explicitly:

* **A diagnostic, not an engine step — a *trade*, not a win.** Like :mod:`planet.precip` before it,
  this is a **prescribed kinematic parameterization**, not a simulated water cycle: it *maps* a given
  wind + terrain to a rain pattern, it does not *derive* the wind or close a moisture budget. It makes
  the *precipitation* 2-D; the *temperature* climate underneath stays zonal-mean (the EBM). So it does
  **not** claim the engine has left the zonal mean — it claims the rain shadow has (the Phase-2
  diagnostic-precip precedent, one rung further out).
* **On a regional Cartesian patch with a *prescribed, uniform* background wind.** The linear theory
  lives on a tangent plane with a constant cross-mountain wind — that is where its analytic anchors are
  exact. Where the wind vector *comes from* on a sphere whose emergent jet is purely **zonal**
  (:mod:`planet.coupler`), and how a patch is placed under a real mountain range, are the **5A.2**
  integration questions (named, deferred). Cross-mountain flow is **prescribed here, not emergent** —
  the honesty caveat carried into that step.

Validation triad (plan §3) — what is asserted tight vs loose
------------------------------------------------------------
* **Analytical (tight, *exact*).** For a **triangle ridge** in the reduced limit ``H_w = τ_c = 0``
  (no airborne-water advection), the model has a **closed-form** solution
  (:func:`triangle_ridge_exact`): a windward exponential rise to the crest and an exponential decay
  into the **lee** cut off at ``x_c`` — the rain shadow, analytically. The FFT model converges to it
  as the grid refines (:mod:`planet.tests.test_orographic`). This one anchor validates the transfer
  function, the vertical-wavenumber **branch**, and the windward-wet / lee-dry structure at once.
* **Structural (tight).** The **upslope limit** ``H_w = τ_c = τ_f = 0`` recovers the classic upslope
  model ``P = C_w·max(0, U·∂h/∂x)`` to machine precision; reversing the wind **mirrors** the pattern;
  flat terrain gives zero orographic anomaly; precipitation is non-negative.
* **Benchmark / magnitude (loose).** The pinned constants (``C_w``, ``H_w``, ``N_m``, ``τ_c``, ``τ_f``,
  ``U``) are the cited Smith & Barstad values ([[smith-barstad-orographic-source]]); absolute
  mm/hr amplitudes move only in loose bands, exactly as :mod:`planet.precip`'s band amplitudes do.

Units — mm/hr (the linear-theory native unit), distances in metres, wind in m/s
-------------------------------------------------------------------------------
The transfer function returns a condensation mass flux which, with water density ≈ 1000 kg/m³
(1 kg/m² ≡ 1 mm of water), reads directly as **mm/hr**. Terrain, grid spacing (``dx``/``dy``) and the
water-vapour scale height ``H_w`` are in **metres**; wind speed in **m/s**; ``direction`` in degrees
(meteorological: 0 = wind from the north, 270 = from the west, i.e. a westerly blowing eastward). The
conversion to the biome map's cm/yr and the placement on the global lat×lon grid are 5A.2.

Cited implementation
--------------------
Formulae and pinned constants follow Smith & Barstad (2004), "A Linear Theory of Orographic
Precipitation", *J. Atmos. Sci.* **61**, 1377–1391, cross-checked against the PISM/QGIS reference
implementation (``pism/LinearTheoryOrographicPrecipitation``) — including the triangle-ridge exact
solution used as the tight anchor. See [[smith-barstad-orographic-source]].
"""
from __future__ import annotations

import numpy as np

# --------------------------------------------------------------------------- #
# Pinned Smith & Barstad (2004) constants ([[smith-barstad-orographic-source]]).
# Cited at build (cross-checked vs the PISM/QGIS reference), NOT from memory. These set the *absolute*
# mm/hr magnitude (loose / calibration-dependent); the tight anchors are structural and hold for any
# positive C_w. An idealized set that exposes the windward-wet / lee-dry rain shadow cleanly.
# --------------------------------------------------------------------------- #
TAU_C_S = 1000.0              # s      — cloud-water → hydrometeor conversion time
TAU_F_S = 1000.0              # s      — hydrometeor fallout time
NM_PER_S = 0.005              # s⁻¹    — moist stability (moist Brunt–Väisälä) frequency
HW_M = 2500.0                 # m      — water-vapour scale height (moist-layer depth)
RHO_SREF_KG_M3 = 7.4e-3       # kg/m³  — reference saturation water-vapour density
MOIST_LAPSE_K_PER_KM = -6.5   # K/km   — moist adiabatic lapse rate Θ_m (negative: T falls with height)
ENV_LAPSE_K_PER_KM = -5.8     # K/km   — environmental lapse rate γ (negative)
# Uplift sensitivity C_w = ρ_Sref·Θ_m/γ — condensation produced per unit forced ascent. The two lapse
# rates are both negative, so their ratio is positive → C_w > 0 (a positive rain response to uplift).
CW_KG_M3 = RHO_SREF_KG_M3 * MOIST_LAPSE_K_PER_KM / ENV_LAPSE_K_PER_KM   # ≈ 8.29e-3 kg/m³

U_REF_M_S = 15.0              # m/s    — reference background wind speed
DIRECTION_WESTERLY_DEG = 270.0  # °    — a westerly (wind FROM the west, blowing eastward, +x)
CORIOLIS_LAT_DEG = 45.0       # °      — default latitude for the Coriolis parameter f (wide-mountain term)

OMEGA_EARTH = 7.2921e-5       # s⁻¹    — Earth's rotation rate (for f = 2Ω sin φ)
RHO_WATER = 1000.0            # kg/m³  — so 1 kg/m² of condensate ≡ 1 mm of water
SECONDS_PER_HOUR = 3600.0


def coriolis_f(latitude_deg: float) -> float:
    """The Coriolis parameter ``f = 2Ω sin φ`` (s⁻¹). Zero at the equator; enters the wide-mountain term."""
    return 2.0 * OMEGA_EARTH * np.sin(np.radians(latitude_deg))


def wind_components(speed: float, direction_deg: float) -> tuple[float, float]:
    """``(u, v)`` wind components (m/s) from speed + meteorological ``direction`` (0 = from N, 270 = from W).

    Meteorological convention: ``direction`` is where the wind blows *from*. A **westerly**
    (``direction = 270``) therefore blows toward the **east**, giving ``u = +speed``, ``v = 0``.
    """
    u = -np.sin(np.radians(direction_deg)) * speed
    v = -np.cos(np.radians(direction_deg)) * speed
    return float(u), float(v)


def orographic_precip(orography: np.ndarray, dx: float, dy: float, *,
                      speed: float = U_REF_M_S, direction_deg: float = DIRECTION_WESTERLY_DEG,
                      tau_c: float = TAU_C_S, tau_f: float = TAU_F_S,
                      Nm: float = NM_PER_S, Hw: float = HW_M, Cw: float = CW_KG_M3,
                      latitude_deg: float = CORIOLIS_LAT_DEG,
                      background_mm_hr: float = 0.0, truncate: bool = True) -> np.ndarray:
    """Smith & Barstad (2004) linear orographic precipitation ``P(x, y)`` in **mm/hr**.

    Applies the wavenumber-space transfer function to a 2-D terrain ``orography`` (metres) on a regular
    grid of spacing ``dx``/``dy`` (metres), under a **prescribed uniform** background wind
    (``speed``/``direction_deg``). ``orography`` is indexed ``[y, x]`` (rows = latitude/northward,
    columns = longitude/eastward), matching the map grid.

    The transfer function (Smith & Barstad 2004, eq. 49)::

        P̂(k, l) = C_w · i σ · ĥ(k, l)
                  ───────────────────────────────────────────────
                  (1 − i m H_w)(1 + i σ τ_c)(1 + i σ τ_f)

    with the horizontal wind ``(U, V)``, intrinsic frequency ``σ = U k + V l``, and the vertical
    wavenumber ``m`` from the linear mountain-wave dispersion::

        m² = (N_m² − σ²)(k² + l²) / (σ² − f²)

    Three numerical subtleties (each an anchor's failure mode):

    * **Branch of** ``m`` **(the rain-shadow sign).** For propagating modes (``m² ≥ 0``) the root is
      chosen with the sign of ``σ`` so the wave tilts **upwind** and the drying sits in the **lee**;
      get this wrong and the shadow flips to the windward side. Evanescent modes (``m² < 0``) take the
      decaying (bounded) branch automatically.
    * **The** ``σ = 0`` **locus** (modes ⊥ the wind, including the mean): the numerator ``i σ ĥ`` is
      zero there, so ``P̂ = 0`` — set explicitly to avoid a 0/∞.
    * **FFT periodicity:** the terrain is zero-padded by a full domain width so lee-side drying cannot
      wrap around into the upwind edge.

    ``background_mm_hr`` is a spatially-uniform ambient rate added after the transfer (the ``P_∞`` of
    the linear theory, isolated to 0 for the pure orographic anomaly). With ``truncate`` the result is
    clamped to ``≥ 0`` (the model can produce small negative values in the lee that are unphysical rain).
    """
    orography = np.asarray(orography, dtype=float)
    eps = 1e-18

    # Zero-pad by a full domain width: keeps the lee shadow from wrapping into the windward edge.
    pad = max(orography.shape)
    h = np.pad(orography, pad, "constant")
    nrows, ncols = h.shape

    h_hat = np.fft.fft2(h)

    # Angular wavenumbers (rad/m) for the padded grid.
    kx = np.fft.fftfreq(ncols, dx / (2.0 * np.pi))
    ky = np.fft.fftfreq(nrows, dy / (2.0 * np.pi))
    kx, ky = np.meshgrid(kx, ky)

    u0, v0 = wind_components(speed, direction_deg)
    f = coriolis_f(latitude_deg)

    sigma = u0 * kx + v0 * ky                      # intrinsic frequency σ = Uk + Vl

    # m² = (N_m² − σ²)(k²+l²)/(σ² − f²); regularize the σ² = f² resonance so it never divides by zero.
    denom = sigma**2 - f**2
    denom[np.logical_and(np.fabs(denom) < eps, denom >= 0)] = eps
    denom[np.logical_and(np.fabs(denom) < eps, denom < 0)] = -eps
    m_squared = (Nm**2 - sigma**2) * (kx**2 + ky**2) / denom

    m = np.sqrt(np.asarray(m_squared, dtype=np.cdouble))
    # Propagating modes: pick the branch tilting upwind (sign of σ). Evanescent modes keep the
    # decaying branch numpy's complex sqrt already returns.
    propagating = np.logical_and(m_squared >= 0, sigma != 0)
    m[propagating] *= np.sign(sigma[propagating])

    P_hat = h_hat * (Cw * 1j * sigma
                     / ((1 - 1j * m * Hw) * (1 + 1j * sigma * tau_c) * (1 + 1j * sigma * tau_f)))

    P = np.real(np.fft.ifft2(P_hat))
    P = P[pad:-pad, pad:-pad] if pad > 0 else P    # strip the padding

    # The transfer function returns a condensation mass flux in kg/(m²·s). With water density
    # ρ_water = 1000 kg/m³, 1 kg/m² of water ≡ 1 mm depth, so kg/(m²·s) = mm/s → ×3600 = mm/hr.
    P *= SECONDS_PER_HOUR
    P += background_mm_hr                           # the ambient P_∞ (0 for the pure orographic anomaly)
    if truncate:
        P[P < 0.0] = 0.0                            # small negative lee values are unphysical rain
    return P


# --------------------------------------------------------------------------- #
# Idealized terrain — the spike geometries (a compact ridge for the FFT, a triangle for the exact anchor)
# --------------------------------------------------------------------------- #
def make_grid(half_width_m: float = 200e3, dx: float = 2e3) -> tuple[np.ndarray, float]:
    """A symmetric 1-D coordinate axis ``x ∈ [−half_width, +half_width]`` (m) at spacing ``dx`` (m)."""
    n = int(round(2 * half_width_m / dx)) + 1
    return np.linspace(-half_width_m, half_width_m, n), dx


def gaussian_ridge(x: np.ndarray, amplitude_m: float = 1000.0, sigma_m: float = 15e3) -> np.ndarray:
    """A compact Gaussian mountain ridge ``h(x) = A·exp(−x²/2σ²)`` (m) — smooth, decays inside the domain.

    Smoothness (no kink) makes it the clean geometry for the FFT model, and its compact support keeps
    the terrain ≈ 0 at the domain edges so the zero-padded FFT sees no wrap-around.
    """
    x = np.asarray(x, dtype=float)
    return amplitude_m * np.exp(-(x**2) / (2.0 * sigma_m**2))


def triangle_ridge(x: np.ndarray, amplitude_m: float = 500.0, half_width_m: float = 50e3) -> np.ndarray:
    """A symmetric triangular ridge ``h(x) = A·max(0, 1 − |x|/d)`` (m) — the shape with a *closed-form* rain.

    Used only for the exact anchor (:func:`triangle_ridge_exact`); its kink at the crest is what makes
    the reduced-limit precipitation integrable in closed form.
    """
    x = np.asarray(x, dtype=float)
    return np.maximum(amplitude_m * (1.0 - np.fabs(x) / half_width_m), 0.0)


def triangle_ridge_exact(x: np.ndarray, u: float, Cw: float, tau: float,
                         amplitude_m: float = 500.0, half_width_m: float = 50e3) -> np.ndarray:
    """The **exact** orographic precip (mm/hr) for a triangle ridge in the reduced limit ``H_w = τ_c = 0``.

    With no water-vapour scale height and no conversion delay, the Smith & Barstad transfer function
    collapses to advection–fallout of condensate along the wind, and the triangle ridge (wind blowing
    toward +x) has the closed-form solution (Smith & Barstad 2004; PISM reference):

    * **windward slope** ``−d ≤ x < 0``:  ``P = C·(1 − e^{−(x+d)/Uτ})`` — rain builds up the slope;
    * **lee** ``0 ≤ x ≤ x_c``:  ``P = C·(e^{−x/Uτ}·(2 − e^{−d/Uτ}) − 1)`` — the **rain shadow**, decaying
      to zero at the cutoff ``x_c = Uτ·ln(2 − e^{−d/Uτ})``;
    * elsewhere ``P = 0``,

    where ``C = C_w·u·A/d`` (the upslope rain rate) and ``Uτ = u·τ`` (the fallout advection length). This
    is the tight analytic anchor: the FFT :func:`orographic_precip` converges to it as the grid refines.
    """
    d, A = half_width_m, amplitude_m
    C = Cw * u * A / d
    Ut = u * tau
    xc = Ut * np.log(2.0 - np.exp(-d / Ut))

    def _p(xi: float) -> float:
        if -d <= xi < 0.0:
            return C * (1.0 - np.exp(-(xi + d) / Ut))
        if 0.0 <= xi <= xc:
            return C * (np.exp(-xi / Ut) * (2.0 - np.exp(-d / Ut)) - 1.0)
        return 0.0

    x = np.asarray(x, dtype=float)
    flat = np.array([_p(float(xi)) for xi in np.atleast_1d(x).ravel()])
    return (SECONDS_PER_HOUR * flat).reshape(x.shape) if x.ndim else float(SECONDS_PER_HOUR * flat[0])
