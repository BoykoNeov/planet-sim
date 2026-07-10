"""The full-sphere MSE-diffusing EBM — and why moisture does NOT tighten the ITCZ sensitivity.

This module answers a backlog question with a **clean negative** (plan §12.2, the rung-2.x refinement
"tighten the ITCZ-migration sensitivity"). Rung 2.x's dry sphere EBM (:mod:`planet.sphere_ebm`) gives an
ITCZ-migration sensitivity ``≈ −6.3 deg/PW`` — a factor ~1.5–2 above the observed ``~3`` (Donohoe et al.
2013). The backlog line proposed *re-deriving ``D``* to close the gap; the dry module now shows why that
cannot work (the sensitivity is a **radiation** quantity, ``δ/AHT = −1/(2π a² NEI(0))``, that ``D`` cancels
out of). The natural next hypothesis is that the gap is **moisture**: the observed ``AHT`` is *moist* static
energy transport, which carries ~``(1 + β)`` × the sensible energy per degree at the warm equator, so a
given EFE shift should correspond to a larger ``PW`` and hence a smaller ``deg/PW``. This module builds that
model — the **rung-2.5 MSE-diffusing physics** (:mod:`planet.moist_ebm`) on the **rung-2.x full sphere**
(:mod:`planet.sphere_ebm`) — and measures the sensitivity. **The hypothesis fails, structurally.**

The result — moisture moves it ~10 %, not ~50 %, and the NEI identity says why
-----------------------------------------------------------------------------
Recalibrated so the moist present-day climate reproduces the dry equator-to-pole **contrast** (the rung-2.5
double-count discipline, :func:`recalibrate_sensible_D_sphere` → ``D_s ≈ 0.28``), the moist sensitivity is
``≈ −5.7 deg/PW`` — the dry ``−6.3`` moved by only ~10 %, and it *saturates* there across RH. The reason is
**not** a coincidence: the moisture-amplified diffusivity ``D_eff(0) = D_s(1 + β(T_eq))`` rises (~0.55 →
~0.93) but the moisture-flattened equatorial curvature ``T̄ₓₓ(0)`` shrinks in near-lockstep (~−64 → ~−43),
so their **product ``D_eff(0)·T̄ₓₓ(0)`` barely moves** — and that product is exactly the thing pinned by the
equatorial energy balance:

    D_eff(0) · T̄ₓₓ(0) = −NEI(0)         (the same identity the dry module banks)

Both the enhanced diffusivity **and** the flatter curvature are two faces of the *same* equatorial moisture
amplification, so they cancel in the product. The MSE upgrade changes the sensitivity **only** through the
~1.7 K by which the moist equator runs cooler than the dry one (a small ``−B·ΔT_eq`` shift of ``NEI(0)``) —
it does **not** change the fact that the sensitivity is set by ``NEI(0)``, a **radiation** quantity that
neither ``D`` nor ``D_eff`` reaches. So the moist model corroborates the dry module's headline rather than
overturning it: **the factor-2 gap to observed is not a moisture-transport effect** — it needs a stronger
equatorial radiative surplus (rung 4) or the gross-moist-stability dynamics the diffusive closure omits
(the deferred rung-3+ wall, :mod:`planet.moist` Hadley note). This is a negative result, banked at the
altitude the identity supports.

A SIBLING — sphere_ebm.py and moist_ebm.py are both untouched
-------------------------------------------------------------
Composition, not surgery (the every-rung-is-a-sibling discipline, ADR 0005): this reuses
:class:`planet.sphere_ebm.SphereEBM` for the pole-to-pole geometry, energetic-ITCZ diagnostics and the
``NEI`` reading, and :func:`planet.moist_ebm.effective_diffusivity` for the moisture-amplified ``D_eff(T)``
placed **inside** the conservative divergence. The steady solve is the dt-free **Picard** iteration of
rung 2.5 (:func:`planet.moist_ebm.moist_steady_direct`) lifted to ``x ∈ [−1, 1]`` — freeze ``D_eff(T_k)``,
solve the linear ``(L_T[D_eff] − B·I)T = A − S(1−α) − Q`` in one tridiagonal banded solve, repeat. No
time-stepping ⟹ **no operator-splitting bias** (the rung-2.5 gotcha). ``RH = 0`` and ``D_s = D_TRANSPORT``
reduce it **bit-for-bit** to :meth:`planet.sphere_ebm.SphereEBM.steady_linear` (the plumbing reduction).

Validation triad (plan §3)
--------------------------
* **Tight.** Bit-for-bit reduction to the dry :meth:`SphereEBM.steady_linear` at ``RH = 0`` (``β ≡ 0`` ⟹
  the dry operator); EFE ``= 0`` exactly for a symmetric climate; the **NEI identity**
  ``D_eff(0)·T̄ₓₓ(0) = −NEI(0)`` and hence the ``NEI``-form sensitivity matching the **measured** Q-sweep
  migration.
* **Real-but-loose (the unlock — a NEGATIVE).** The moist sensitivity ``≈ −5.7 deg/PW`` — the dry ``−6.3``
  moved only ~10 %, saturating across RH; the factor-2 gap to observed is **not** closed by MSE transport.
  The ~10 % is the ``−B·ΔT_eq`` equatorial-cooling shift of ``NEI(0)``, not a transport effect.
* **Plumbing.** The recalibration ``D_s`` matches the dry contrast; an imposed odd ``Q(x) = q·x`` shifts the
  EFE toward the warmed hemisphere (direction by-construction).

Units — SI (W m⁻², °C, ``x = sin φ`` on [−1, 1]); transport in PW. Sources: Bischoff & Schneider 2014
(the ``δ ≈ −AHT/NEI`` identity), Donohoe et al. 2013 (~3 deg/PW), Hwang & Frierson 2010 / Flannery 1984
(the fixed-RH MSE diffusion). See [[planet-rung2x-itcz]], [[moist-ebm-source]], [[seasonal-ebm-source]].
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np
from scipy.linalg import solve_banded
from scipy.optimize import brentq

from . import ebm, moist
from .ebm import A_OLR, B_OLR, D_TRANSPORT, T_FREEZE, WATER_DEPTH
from .moist_ebm import effective_diffusivity, moisture_amplification
from .sphere_ebm import (
    AREA_FACTOR, PW, SphereClimate, SphereEBM, efe_from_transport,
    itcz_sensitivity_from_nei,
)

RH_DEFAULT = moist.RH_DEFAULT


def constant_albedo_absorbed(a0: float = ebm.ALBEDO_A0) -> Callable[[np.ndarray, np.ndarray], np.ndarray]:
    """The ice-OFF absorbed shortwave ``S(x)·(1 − a₀)`` — the linear, splitting-free reference forcing.

    Identical to the dry module's ``_absorbed_const`` (constant albedo, so the direct Picard solve applies
    and the ``RH = 0`` reduction is bit-for-bit against :meth:`SphereEBM.steady_linear`).
    """
    def absorbed(x, T):
        return ebm.insolation(x) * (1.0 - float(a0))
    return absorbed


class SphereMoistEBM:
    """Full-sphere MSE-diffusing moist EBM; see the module docstring. Composes :class:`SphereEBM`.

    ``D_s`` is the dry **sensible** coefficient (the moist ``D_eff = D_s·(1 + β(T))`` amplifies it, inside
    the divergence); ``RH`` the fixed relative humidity (``RH = 0`` ⟹ the dry sphere, bit-for-bit). All
    geometry, the equatorial-``NEI`` reading and the parameters live on the wrapped :class:`SphereEBM`.
    """

    def __init__(self, D_s: float = 0.28, RH: float = RH_DEFAULT, A: float = A_OLR, B: float = B_OLR,
                 T_freeze: float = T_FREEZE, water_depth: float = WATER_DEPTH, n_cells: int = 360):
        self.base = SphereEBM(A=A, B=B, D=D_s, T_freeze=T_freeze, water_depth=water_depth, n_cells=n_cells)
        self.D_s = float(D_s)
        self.RH = float(RH)

    # -- geometry / parameter passthroughs --------------------------------- #
    @property
    def x(self) -> np.ndarray:
        return self.base.x

    @property
    def A(self) -> float:
        return self.base.A

    @property
    def B(self) -> float:
        return self.base.B

    # -- the moist energy transport (D_eff(T) inside) ---------------------- #
    def moist_transport(self, T: np.ndarray) -> np.ndarray:
        """Northward MOIST energy transport ``H(x) = −2π a² D_eff(T)(1−x²) ∂T/∂x`` (PW).

        The dry :meth:`SphereEBM.atmospheric_transport` with the moisture-amplified ``D_eff(T)`` in place of
        the constant ``D`` — this is the moist static energy flux, whose zero nearest the equator is the EFE.
        """
        T = np.asarray(T, dtype=float)
        Deff = effective_diffusivity(T, self.D_s, self.RH)
        return -AREA_FACTOR * Deff * (1.0 - self.x ** 2) * np.gradient(T, self.x) / PW

    def energy_flux_equator(self, T: np.ndarray) -> tuple[float, float]:
        """``(φ_EFE, AHT_eq)`` from the **moist** transport (the identical rule, :func:`efe_from_transport`)."""
        return efe_from_transport(self.x, self.moist_transport(T))

    def net_radiative_input_equator(self, absorbed_fn, T: np.ndarray) -> float:
        """``NEI(0)`` at the equator — reuses :meth:`SphereEBM.net_radiative_input_equator` (radiation only)."""
        return self.base.net_radiative_input_equator(absorbed_fn, T)

    def equatorial_curvature(self, T: np.ndarray, band: float = 0.15) -> float:
        """``T̄ₓₓ(0)`` — reuses the dry module's degree-2 equatorial fit."""
        return self.base.equatorial_curvature(T, band)

    def _result(self, T: np.ndarray, converged: bool, it: int) -> SphereClimate:
        phi_efe, aht_eq = self.energy_flux_equator(T)
        return SphereClimate(x=self.x, T=T, global_mean_T=self.base.global_mean(T),
                             phi_efe=phi_efe, aht_eq=aht_eq, converged=converged, iterations=it)

    # -- the dt-free Picard steady solve (D_eff frozen each iterate) ------- #
    def steady(self, absorbed_fn, Q: Optional[np.ndarray] = None,
               tol: float = 1e-11, max_iter: int = 800, T_init=None) -> SphereClimate:
        """The dt-free moist steady climate by Picard on the frozen-``D_eff`` linear solve (constant albedo).

        Each iterate freezes ``D_eff(T_k) = D_s(1 + β(T_k))`` inside the conservative sphere operator
        ``d/dx[D_eff(1−x²) d/dx]`` (harmonic-mean faces, the engine's idiom), subtracts ``B``, and solves
        ``(L_T − B·I)T = A − S(1−α) − Q`` in one banded solve — repeating to ``max|ΔT| < tol``. **No
        time-stepping ⟹ no operator-splitting bias** (the rung-2.5 headline path lifted to the sphere).
        ``Q`` is an optional odd cross-equatorial source (W m⁻²; ``q·x`` drives the EFE off the equator).
        **Raises** on a state-dependent (ice) ``absorbed_fn`` — that must go through a relaxation, out of
        scope for this diagnostic. ``RH = 0`` ⟹ constant ``D_eff`` ⟹ the dry :meth:`SphereEBM.steady_linear`
        bit-for-bit.
        """
        x = self.x
        a_cold = absorbed_fn(x, np.full(x.shape, -100.0))
        a_warm = absorbed_fn(x, np.full(x.shape, 100.0))
        if not np.allclose(a_cold, a_warm):
            raise ValueError("SphereMoistEBM.steady requires a state-independent (constant-albedo) absorbed "
                             "field; the ice feedback is out of scope for this diagnostic")
        absorbed = np.asarray(a_cold, dtype=float)
        Qarr = np.zeros_like(x) if Q is None else np.asarray(Q, dtype=float)
        rhs = self.A - absorbed - Qarr
        dx = self.base.grid.widths
        dxc = np.diff(self.base.grid.centers)
        n = x.size
        T = (15.0 - 40.0 * x ** 2) if T_init is None else np.array(np.broadcast_to(T_init, x.shape), float)
        converged, it = False, 0
        for it in range(1, max_iter + 1):
            Dcoef = effective_diffusivity(T, self.D_s, self.RH) * (1.0 - x ** 2)   # D_eff(T) inside divergence
            Dface = 2.0 * Dcoef[:-1] * Dcoef[1:] / (Dcoef[:-1] + Dcoef[1:])        # harmonic-mean faces
            Tt = Dface / dxc
            sub = np.zeros(n); diag = np.zeros(n); sup = np.zeros(n)
            sup[:-1] += Tt / dx[:-1]; diag[:-1] += -Tt / dx[:-1]
            sub[1:] += Tt / dx[1:];   diag[1:] += -Tt / dx[1:]
            diag = diag - self.B
            ab = np.zeros((3, n)); ab[0, 1:] = sup[:-1]; ab[1, :] = diag; ab[2, :-1] = sub[1:]
            T_new = solve_banded((1, 1), ab, rhs)
            if np.max(np.abs(T_new - T)) < tol:
                T = T_new
                converged = True
                break
            T = T_new
        return self._result(T, converged, it)

    # -- the ITCZ-migration sensitivity (measured + the NEI closed form) --- #
    def itcz_sensitivity(self, absorbed_fn, q_amps=(0.0, 1.0, 2.0, 4.0)) -> tuple[float, float]:
        """``(slope_measured, slope_nei)`` in deg/PW — the moist EFE migration rate and its ``NEI`` closed form.

        Sweeps an odd Q-flux ``Q(x) = q·x`` over ``q_amps``, fits ``φ_EFE`` against the **moist** ``AHT_eq``,
        and returns that slope beside :func:`itcz_sensitivity_from_nei` evaluated at the symmetric base
        state's ``NEI(0)``. They agree (the tight identity) — and both land ``≈ −5.7`` (the dry ``−6.3``
        barely moved; the negative result).
        """
        phis, ahts = [], []
        nei0 = None
        for q in q_amps:
            c = self.steady(absorbed_fn, Q=(q * self.x))
            phis.append(c.phi_efe); ahts.append(c.aht_eq)
            if q == 0.0:
                nei0 = self.net_radiative_input_equator(absorbed_fn, c.T)
        slope = float(np.polyfit(ahts, phis, 1)[0])
        return slope, itcz_sensitivity_from_nei(nei0)


def equator_pole_contrast(x: np.ndarray, T: np.ndarray) -> float:
    """The full-sphere equator-to-pole contrast ``T(0°) − T(90°)`` (°C) — equator by interpolation to x=0."""
    x = np.asarray(x, dtype=float); T = np.asarray(T, dtype=float)
    return float(np.interp(0.0, x, T) - T[-1])


def recalibrate_sensible_D_sphere(RH: float = RH_DEFAULT, a0: float = ebm.ALBEDO_A0, n_cells: int = 360,
                                  A: float = A_OLR, B: float = B_OLR, D_lo: float = 0.05,
                                  D_hi: float = 0.9) -> float:
    """The sensible ``D_s`` whose moist sphere climate matches the DRY sphere equator-pole contrast (W m⁻² K⁻¹).

    The rung-2.5 double-count discipline (:func:`planet.moist_ebm.recalibrate_sensible_D`) on the full
    sphere: rung-0's ``D = 0.555`` already lumps in latent transport, so diffusing MSE explicitly on top
    would double-count it — re-derive a smaller sensible ``D_s`` that reproduces the dry contrast. For Earth
    defaults + RH 0.8 this gives ``D_s ≈ 0.28`` (the same value the hemisphere recalibration lands on).
    """
    absorbed = constant_albedo_absorbed(a0)
    dry = SphereEBM(A=A, B=B, D=D_TRANSPORT, n_cells=n_cells).steady_linear(absorbed)
    dry_contrast = equator_pole_contrast(dry.x, dry.T)

    def gap(D_s: float) -> float:
        c = SphereMoistEBM(D_s=D_s, RH=RH, A=A, B=B, n_cells=n_cells).steady(absorbed)
        return equator_pole_contrast(c.x, c.T) - dry_contrast

    return float(brentq(gap, D_lo, D_hi, xtol=1e-4))


@dataclass(frozen=True)
class MoistSensitivity:
    """The moist-vs-dry ITCZ-sensitivity comparison and the NEI decomposition of the (non-)tightening.

    ``D_s``/``RH`` the recalibrated moist closure; ``slope_dry``/``slope_moist`` the ITCZ sensitivities
    (deg/PW); ``nei_dry``/``nei_moist`` the equatorial net radiative inputs (W m⁻², the sensitivity's true
    denominator); ``Teq_dry``/``Teq_moist`` the equatorial temperatures (the ~1.7 K cooling that is the
    *entire* moist effect); ``Deff0``/``Txx_moist`` the moist equatorial diffusivity and curvature (whose
    product is pinned to ``−nei_moist``). ``observed_deg_per_pw`` the Donohoe target. Plain scalars.
    """

    D_s: float
    RH: float
    slope_dry: float
    slope_moist: float
    nei_dry: float
    nei_moist: float
    Teq_dry: float
    Teq_moist: float
    Deff0: float
    Txx_moist: float
    observed_deg_per_pw: float


def moist_vs_dry_sensitivity(RH: float = RH_DEFAULT, n_cells: int = 360,
                             observed: float = -3.0) -> MoistSensitivity:
    """Compute the dry-vs-moist ITCZ sensitivity and the ``NEI`` decomposition — the module headline.

    Recalibrates ``D_s`` (matches the dry contrast), builds the symmetric dry and moist base states, and
    measures both ITCZ sensitivities. Returns a :class:`MoistSensitivity` exposing that the moist number
    ``≈ −5.7`` barely moved the dry ``≈ −6.3`` and that the whole difference is the ``−B·ΔT_eq`` shift of
    ``NEI(0)`` (a radiation effect), not the moisture-amplified transport.
    """
    absorbed = constant_albedo_absorbed()
    D_s = recalibrate_sensible_D_sphere(RH=RH, n_cells=n_cells)

    dry = SphereEBM(D=D_TRANSPORT, n_cells=n_cells)
    cdry = dry.steady_linear(absorbed)
    slope_dry, _ = dry.itcz_sensitivity(absorbed, linear=True)
    nei_dry = dry.net_radiative_input_equator(absorbed, cdry.T)
    Teq_dry = float(np.interp(0.0, cdry.x, cdry.T))

    m = SphereMoistEBM(D_s=D_s, RH=RH, n_cells=n_cells)
    cm = m.steady(absorbed)
    slope_moist, _ = m.itcz_sensitivity(absorbed)
    nei_moist = m.net_radiative_input_equator(absorbed, cm.T)
    Teq_moist = float(np.interp(0.0, cm.x, cm.T))
    Deff0 = float(effective_diffusivity(Teq_moist, D_s, RH))
    Txx_moist = m.equatorial_curvature(cm.T)

    return MoistSensitivity(
        D_s=float(D_s), RH=float(RH), slope_dry=float(slope_dry), slope_moist=float(slope_moist),
        nei_dry=float(nei_dry), nei_moist=float(nei_moist), Teq_dry=Teq_dry, Teq_moist=Teq_moist,
        Deff0=Deff0, Txx_moist=float(Txx_moist), observed_deg_per_pw=float(observed),
    )
