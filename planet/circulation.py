"""Planet's instantiation of the shallow-water engine — a midlatitude β-plane (Planet Phase 3).

This is the **consumer** of the shared :mod:`engines.fluid` (load
``engines/fluid/CONTRACT.md``, never the engine internals): it pins the *planetary*
constants the engine leaves to its consumer (Earth's rotation → ``f₀``, ``β``; an
equivalent depth ``H`` → gravity-wave speed and the deformation radius ``L_R``) and drives
the engine through the two banked Phase-3 artifacts — a **geostrophic adjustment** (an
unbalanced height anomaly radiating gravity waves and settling to a balanced vortex over
``L_R``) and a **westward Rossby wave** — with the conservation diagnostics (mass / energy /
potential enstrophy) tracked alongside.

This mirrors how Steel's ``jominy``/``carburize`` and Chip's ``diffusion_dopant`` consume the
diffusion spine: the *engine* validates the generic solver (geostrophic balance, wave speeds,
PV at finite amplitude — ``engines/fluid/tests/``); *here* the planetary numbers are pinned
and validated (the realistic ``L_R``, the adjusted-jet scale, the westward propagation).

Phase 3 is **one-layer dry dynamics in isolation** — there is no coupling to the EBM yet.
Phase 4 (``coupler.py``) forces this flow with the EBM's meridional temperature gradient so a
geostrophically-balanced jet emerges; that is where the interactive map registers its
``vector_overlay`` circulation layer (no renderer edit). The dry single layer carries no
thermodynamic variable, which is exactly *why* Phase-4 coupling is one-way (plan §3–4).

Units — SI ([[shallow-water-source]])
-------------------------------------
``f₀``, ``β`` from Earth's rotation (1/s, 1/(m·s)); lengths in m, speeds in m/s, time in s;
``H`` an equivalent depth in m. Latitudes in degrees on input.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from engines.fluid import ShallowWater, SWState, uniform_grid

# --------------------------------------------------------------------------- #
# Pinned planetary constants ([[shallow-water-source]] — Gill 1982 / Vallis 2017).
# Earth's rotation sets the β-plane; the equivalent depth is a calibration (loose) chosen so
# the midlatitude deformation radius matches the cited extratropical value (~1000 km).
# --------------------------------------------------------------------------- #
OMEGA_EARTH = 7.292e-5     # rad/s  — Earth's angular velocity
R_EARTH = 6.371e6          # m      — Earth's mean radius
G_EARTH = 9.81             # m/s²   — gravitational acceleration
PHI_REF_DEG = 45.0         # °      — reference latitude for the β-plane (mid-latitudes)
H_EQUIV = 1000.0           # m      — equivalent depth; L_R(45°) ≈ 960 km (calibrated, loose)


def coriolis_f0(phi_deg: float = PHI_REF_DEG) -> float:
    """Coriolis parameter ``f₀ = 2Ω sin φ`` (1/s) at latitude ``phi_deg``."""
    return 2.0 * OMEGA_EARTH * np.sin(np.radians(phi_deg))


def coriolis_beta(phi_deg: float = PHI_REF_DEG) -> float:
    """Meridional Coriolis gradient ``β = 2Ω cos φ / a`` (1/(m·s)) at latitude ``phi_deg``."""
    return 2.0 * OMEGA_EARTH * np.cos(np.radians(phi_deg)) / R_EARTH


def midlatitude_beta_plane(nx: int = 96, ny: int = 96, n_LR: float = 6.0,
                           phi_ref_deg: float = PHI_REF_DEG, H: float = H_EQUIV,
                           beta: float | None = None):
    """A planetary midlatitude β-plane: ``(grid, sw)`` sized to ``n_LR`` deformation radii.

    The domain spans ``n_LR × L_R`` in each direction (so the adjustment / Rossby features
    are well inside it), with ``f₀``/``β`` at ``phi_ref_deg`` and equivalent depth ``H``.
    ``beta=0.0`` gives an f-plane (for the adjustment demo, where β is not needed).
    """
    f0 = coriolis_f0(phi_ref_deg)
    b = coriolis_beta(phi_ref_deg) if beta is None else float(beta)
    L_R = np.sqrt(G_EARTH * H) / f0
    L = n_LR * L_R
    grid = uniform_grid(L, L, nx, ny)
    sw = ShallowWater(grid, G_EARTH, H, f0=f0, beta=b)
    return grid, sw


# --------------------------------------------------------------------------- #
# Banked artifact 1 — geostrophic adjustment (height anomaly → balanced jet over L_R)
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class AdjustmentRun:
    """A geostrophic-adjustment record (plain arrays — the loose-coupling currency).

    ``x``/``y`` cell-center coordinates (m); ``eta_init`` the initial unbalanced height anomaly,
    ``eta_balanced`` the late time-averaged balanced remnant, ``eta_helmholtz`` the analytic
    adjusted state ``(1−L_R²∇²)⁻¹η_init`` (m). ``times`` and ``mass``/``energy``/``enstrophy``
    are the conservation diagnostics over the run (relative to the initial value). ``L_R`` the
    deformation radius (m).
    """

    x: np.ndarray
    y: np.ndarray
    eta_init: np.ndarray
    eta_balanced: np.ndarray
    eta_helmholtz: np.ndarray
    times: np.ndarray
    mass: np.ndarray
    energy: np.ndarray
    enstrophy: np.ndarray
    L_R: float


def geostrophic_adjustment(nx: int = 96, ny: int = 96, n_periods: float = 25.0,
                           sigma_frac: float = 1.0 / 3.0, amp: float = 1.0) -> AdjustmentRun:
    """Run the geostrophic-adjustment demo: an f-plane height bump settling to a balanced vortex.

    Initializes a Gaussian height anomaly (width ``sigma_frac·L_R``, amplitude ``amp`` m) at rest,
    integrates ``n_periods`` inertial periods, and time-averages the late third to extract the
    balanced remnant (filtering the radiated gravity waves). Returns an :class:`AdjustmentRun`.
    """
    grid, sw = midlatitude_beta_plane(nx, ny, beta=0.0)         # adjustment uses the f-plane
    L_R = sw.rossby_radius
    X, Y = grid.center_mesh()
    sigma = sigma_frac * L_R
    eta_init = amp * np.exp(-(((X - grid.Lx / 2) ** 2 + (Y - grid.Ly / 2) ** 2)) / (2 * sigma ** 2))
    s = SWState(h=sw.H + eta_init, u=np.zeros((ny, nx)), v=np.zeros((ny, nx)))

    f0 = sw.f0
    dt = sw.max_dt(s) * 0.6
    n = int(n_periods * 2 * np.pi / f0 / dt)
    nwin = max(1, n // 3)
    m0, e0, z0 = sw.mass(s), sw.energy(s), sw.potential_enstrophy(s)
    times, mass, energy, enstrophy = [], [], [], []
    eta_acc = np.zeros((ny, nx))
    t = 0.0
    for i in range(n):
        s = sw.step(s, dt); t += dt
        times.append(t)
        mass.append(sw.mass(s) / m0 - 1.0)
        energy.append(sw.energy(s) / e0 - 1.0)
        enstrophy.append(sw.potential_enstrophy(s) / z0 - 1.0)
        if i >= n - nwin:
            eta_acc += (s.h - sw.H)
    eta_balanced = eta_acc / nwin

    # analytic adjusted state: (1 − L_R²∇²) η_adj = η_init  (spectral solve, periodic grid)
    kx = 2 * np.pi * np.fft.fftfreq(nx, d=grid.dx)
    ky = 2 * np.pi * np.fft.fftfreq(ny, d=grid.dy)
    KX, KY = np.meshgrid(kx, ky)
    eta_helm = np.real(np.fft.ifft2(np.fft.fft2(eta_init) / (1.0 + L_R ** 2 * (KX ** 2 + KY ** 2))))

    return AdjustmentRun(
        x=grid.x_centers(), y=grid.y_centers(),
        eta_init=eta_init, eta_balanced=eta_balanced, eta_helmholtz=eta_helm,
        times=np.array(times), mass=np.array(mass), energy=np.array(energy),
        enstrophy=np.array(enstrophy), L_R=L_R,
    )


# --------------------------------------------------------------------------- #
# Banked artifact 2 — a westward-propagating Rossby wave
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class RossbyRun:
    """A Rossby-wave record. ``x``/``y`` coordinates (m); ``snapshots`` a list of height-anomaly
    fields (m) at ``snapshot_times`` (s); ``c_measured``/``c_analytic`` the zonal phase speed
    (m/s, negative = westward); ``mk``/``mk_l`` the zonal/meridional mode numbers."""

    x: np.ndarray
    y: np.ndarray
    snapshots: list
    snapshot_times: np.ndarray
    c_measured: float
    c_analytic: float
    mk: int
    ml: int


def rossby_wave(nx: int = 96, ny: int = 96, mk: int = 1, ml: int = 1,
                n_snapshots: int = 4, frac_period: float = 0.5) -> RossbyRun:
    """Run the Rossby-wave demo: a balanced single mode propagating **westward** on the β-plane.

    Initializes a geostrophically-balanced ``cos(kx+ly)`` mode, integrates ``frac_period`` of its
    analytic Rossby period taking ``n_snapshots`` evenly-spaced frames, and measures the zonal
    phase speed from the tracked mode phase. Returns a :class:`RossbyRun` (``c < 0`` = westward).
    """
    grid, sw = midlatitude_beta_plane(nx, ny)
    f0 = sw.f0
    k = 2 * np.pi * mk / grid.Lx
    l = 2 * np.pi * ml / grid.Ly
    Xc, Yc = grid.center_mesh()
    xU = grid.x_corners()[None, :] * np.ones((ny, 1)); yU = grid.y_centers()[:, None] * np.ones((1, nx))
    xV = grid.x_centers()[None, :] * np.ones((ny, 1)); yV = grid.y_corners()[:, None] * np.ones((1, nx))
    eta = (f0 / sw.g) * np.cos(k * Xc + l * Yc)              # geostrophic streamfunction Ψ=1
    u = l * np.sin(k * xU + l * yU)
    v = -k * np.sin(k * xV + l * yV)
    s = SWState(h=sw.H + eta, u=u, v=v)

    om_an = -sw.beta * k / (k ** 2 + l ** 2 + 1.0 / sw.rossby_radius ** 2)
    dt = sw.max_dt(s) * 0.8
    T_end = frac_period * 2 * np.pi / abs(om_an)
    n = int(T_end / dt)
    snap_every = max(1, n // (n_snapshots - 1))
    snapshots, snap_times = [(s.h - sw.H).copy()], [0.0]
    phase, times = [], []
    t = 0.0
    for i in range(1, n + 1):
        s = sw.step(s, dt); t += dt
        de = s.h - sw.H
        a = np.mean(de * np.cos(k * Xc + l * Yc)); b = np.mean(de * np.sin(k * Xc + l * Yc))
        phase.append(np.arctan2(b, a)); times.append(t)
        if i % snap_every == 0 and len(snapshots) < n_snapshots:
            snapshots.append(de.copy()); snap_times.append(t)
    om_meas = np.polyfit(np.array(times), np.unwrap(phase), 1)[0]
    return RossbyRun(
        x=grid.x_centers(), y=grid.y_centers(),
        snapshots=snapshots, snapshot_times=np.array(snap_times),
        c_measured=om_meas / k, c_analytic=om_an / k, mk=mk, ml=ml,
    )
