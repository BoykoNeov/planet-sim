"""Two-layer quasi-geostrophic turbulence — the saturated baroclinic eddy-flux engine (rung 3 Phase B).

This is the **rung-3 Phase-B engine**: the standard tool for the experiment rung 1 named but could
not run — does a *nonlinear, saturated* baroclinic eddy field produce a **down-gradient, irreversible,
realistically-scaled** meridional thickness (heat) flux, so that the rung-1 reduction-to-diffusive-EBM
is finally **non-vacuous**? Rung 1's *barotropic* flux was ~1000× too weak and ~90 % reversible
(``irr_frac ~ 0.1``); the prediction is a strong, persistent baroclinic flux with ``irr ~ O(1)``.

Why QG, not the free-surface shallow-water engine (the spike's load-bearing finding)
------------------------------------------------------------------------------------
Phase A built the linear baroclinic growth rate in a free-surface two-layer SW engine
(:class:`engines.fluid.LayeredShallowWater`) and validated it to ~4 %. But the *saturated* runs
that Phase B needs **outcrop** in that engine: at finite amplitude the interface displacement
reaches the layer depth ``H`` ⟹ the thickness ``h → 0`` ⟹ ``PV = (f+ζ)/h`` detonates. The control
number is the Froude ratio ``Fr = U_s/√(g'H) ≈ η_sat/H`` and the first-saturation overshoot drives
``η/H ≈ 12·Fr``; avoiding it costs ``Fr ≲ 0.04`` (e-fold ≳ 370 h) — the full QG-regime cost in the
tool worst-suited to it. This is the empirical reason **two-layer turbulence is done in QG**
(Held & Larichev 1996 *is* two-layer QG): QG is **linearized in the thickness**, so the interface
displacement has **no layer-depth floor** (saturation is well-posed, no outcropping), and the
**rigid lid filters the fast external gravity wave** (no external-mode CFL — the timestep is set by
the slow advective speed). Both stiffnesses the free-surface model paid are gone.

The model (Phillips 1954 / Held & Larichev 1996)
------------------------------------------------
Two layers of QG potential vorticity on a doubly-periodic β-plane, the **perturbation** PV
(anomaly about a fixed mean zonal shear) advected by its own flow:

    ∂q_k/∂t + U_k ∂q_k/∂x + (∂q̄_k/∂y)·∂ψ_k/∂x + J(ψ_k, q_k) = D_k          (k = 1 upper, 2 lower)
    q_k = ∇²ψ_k + (−1)^k F_k (ψ_1 − ψ_2)            (the PV anomaly ↔ streamfunction relation)
    F_k = f₀² / (g' H_k)                            (layer coupling / inverse deformation length²)
    ∂q̄_1/∂y = β + F_1 (U_1 − U_2),  ∂q̄_2/∂y = β − F_2 (U_1 − U_2)   (the mean PV gradients)

The fields ``q_k, ψ_k`` are the **eddies** (perturbations); the mean zonal flows ``U_k`` and the
mean PV gradients ``∂q̄_k/∂y`` enter as **constant background coefficients** — exactly the
"fixed mean shear + prognostic eddies on a doubly-periodic plane" geostrophic-turbulence setup
(the same background design as the SW engine's optional mean shear). ``J(ψ, q) = ψ_x q_y − ψ_y q_x``
is the eddy–eddy advection (the nonlinearity); dropping it leaves the **linear** instability
operator (the tight anchor, :meth:`TwoLayerQG.growth_rate`). ``D_k`` is the dissipation
(hyperviscosity + bottom Ekman drag, both default-off — :meth:`TwoLayerQG.step`).

The mean-flow shear ``U_1 − U_2`` is the **background available-potential-energy reservoir** that
feeds the instability; **β returns** (the SW stability solver was f-plane) and sets both the
**Rhines arrest** of the inverse cascade and a **finite critical shear** ``U_crit = β/F`` — the
Charney–Stern condition that the lower-layer mean PV gradient ``β − F·U_s`` reverse sign.

The spectral PV inversion (and the K=0 trap)
--------------------------------------------
ψ is recovered from q by a **2×2 spectral inversion**: at each wavenumber ``K² = k² + l²`` the
coupled relation ``q̂ = A(K²)·ψ̂`` with

    A = [ −(K²+F_1)      F_1      ]            det A = K²·(K² + F_1 + F_2)
        [   F_2       −(K²+F_2)   ]

is a 2×2 solve. The determinant is well-conditioned for every ``K > 0`` but **zero at ``K = 0``**:
the domain-mean streamfunction is undetermined from the PV (a gauge), so the ``K = 0`` mode is set
to ``ψ̂ = 0`` (the standard convention — the domain-mean flow carries no available energy here).

Validation (the rung discipline)
--------------------------------
* *tight* — the **linear stability** is the analytic two-layer QG (Phillips) dispersion: zero shear
  ⟹ neutral; a short-wave cutoff at ``K² = 2F``; a finite critical shear ``U_crit = β/F``; the
  most-unstable wavelength a few × the deformation radius. The rooted 2×2 dispersion matches the
  independent Phillips closed form to ~1e-9 (equal layers), and the SW 6×6 solver
  (:class:`engines.fluid.TwoLayerStability`) in the rigid-lid limit (``g → ∞``) agrees to <0.5 %
  — the cross-model bridge between the Phase-A (SW) and Phase-B (QG) engines.
* *real-but-loose (the unlock)* — the **emergent saturated eddy thickness diffusivity** feeding the
  reduction-to-EBM: direction (down-gradient) and irreversibility are the banked claims; the
  magnitude is dimensionless / config-tuned (idealized ``κ ~ v'·L_d`` is intrinsically 15–60× below
  Earth's observed ``κ₀ ≈ 2.2×10⁶ m²/s``, so the discriminators are validated **dimensionless**).
* *plumbing* — ``q ↔ ψ`` round-trips to machine precision; zero shear ⟹ no eddies (decay).

Scope edges (named)
-------------------
* **A new model outside ``engines/fluid``** (pseudospectral, not the C-grid) → there is **no
  bit-for-bit single-layer reduction**; Phase A (SW) and Phase B (QG) validate *different* models,
  bridged only by the shared two-layer linear instability (the <0.5 % rigid-lid cross-check above).
* **Homogeneous box → a domain-bulk ``κ``**, not a latitude-resolved ``κ(y)``; a meridional channel
  (which would test the *operator shape*) is the named BC extension, not this increment.
* **Idealized ``(f₀, β, g', H, U_s)``** chosen for a resolvable deformation radius and a dealiased,
  affordable run — honest at rung 3 (validates the *mechanism + scaling*, not Earth jet speeds), the
  same config-tuned honesty banked at rungs 1–2.

Units — SI ([[shallow-water-source]], [[ebm-radiation-source]])
---------------------------------------------------------------
``ψ`` in m²/s; ``q`` in 1/s; ``f₀`` in 1/s; ``β`` in 1/(m·s); ``g'`` in m/s²; ``H_k`` in m;
``U_k`` in m/s; lengths in m; time in s. ``F_k = f₀²/(g'H_k)`` in 1/m².

Sources (extending [[shallow-water-source]]): **Held & Larichev 1996** (two-layer QG turbulence and
the eddy-diffusivity scaling); Phillips 1954 / Eady 1949 (two-layer baroclinic instability);
Vallis 2017 *AOFD* (QG formulation, spectral inversion, the 2/3 dealias rule).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np


# --------------------------------------------------------------------------- #
# State — the stable data boundary: the two-layer PV anomaly field
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class QGState:
    """Two-layer QG prognostic state: the PV anomaly ``q`` stacked over a leading layer axis.

    ``q`` has shape ``(2, ny, nx)`` (layer 1 = upper, layer 2 = lower) — the eddy PV in physical
    space at cell centres, the single stable array boundary. The streamfunction ``ψ`` is *not*
    stored: it is recovered from ``q`` by the spectral inversion (:meth:`TwoLayerQG.invert`)
    whenever a velocity or flux is needed, so ``q`` alone is the prognostic carrier.
    """

    q: np.ndarray

    def copy(self) -> "QGState":
        return QGState(q=np.array(self.q, dtype=float))


# --------------------------------------------------------------------------- #
# The two-layer QG model
# --------------------------------------------------------------------------- #
class TwoLayerQG:
    """Doubly-periodic two-layer quasi-geostrophic model; see the module docstring.

    Parameters
    ----------
    nx, ny : int
        Grid points in x (zonal) and y (meridional); the domain is ``Lx × Ly`` doubly periodic.
    Lx, Ly : float
        Domain extent (m).
    f0 : float
        Coriolis parameter (1/s) — sets the layer coupling ``F_k = f₀²/(g'H_k)``.
    beta : float
        Meridional PV gradient ``β = df/dy`` (1/(m·s)); re-enters at rung 3 (the SW solver was
        f-plane). ``0.0`` is an f-plane (Eady-like: no critical shear).
    gp : float
        Reduced gravity ``g'`` (m/s²) at the internal interface — sets the deformation radius.
    H1, H2 : float
        Rest thicknesses of the upper / lower layer (m).
    U1, U2 : float
        Mean zonal flow in each layer (m/s); the shear ``U_1 − U_2`` is the APE reservoir.
        Use :meth:`symmetric` for the pure-baroclinic ``U_1 = +U_s/2, U_2 = −U_s/2``.
    nu4 : float
        Biharmonic hyperviscosity ``ν₄`` (m⁴/s) on the PV (``−ν₄∇⁴q``); ``0.0`` is inviscid.
        The grid-scale enstrophy sink (turbulent runs cascade to the grid); has its own CFL.
    r_drag : float
        Linear bottom Ekman drag ``r`` (1/s) on the **lower** layer (``−r∇²ψ_2`` added to
        ``∂q_2/∂t``); ``0.0`` is frictionless. Arrests the inverse energy cascade (Held–Larichev).
    """

    def __init__(
        self,
        nx: int,
        ny: int,
        Lx: float,
        Ly: float,
        f0: float,
        gp: float,
        H1: float,
        H2: float,
        beta: float = 0.0,
        U1: float = 0.0,
        U2: float = 0.0,
        nu4: float = 0.0,
        r_drag: float = 0.0,
    ) -> None:
        if nx < 2 or ny < 2:
            raise ValueError("need at least 2 grid points in each direction")
        if gp <= 0.0 or H1 <= 0.0 or H2 <= 0.0:
            raise ValueError("g', H1, H2 must be positive")
        self.nx, self.ny = int(nx), int(ny)
        self.Lx, self.Ly = float(Lx), float(Ly)
        self.dx, self.dy = self.Lx / self.nx, self.Ly / self.ny
        self.f0, self.beta, self.gp = float(f0), float(beta), float(gp)
        self.H1, self.H2 = float(H1), float(H2)
        self.U1, self.U2 = float(U1), float(U2)
        self.nu4, self.r_drag = float(nu4), float(r_drag)

        # Layer coupling F_k = f0^2/(g'H_k) (1/m^2).
        self.F1 = self.f0 ** 2 / (self.gp * self.H1)
        self.F2 = self.f0 ** 2 / (self.gp * self.H2)

        # Spectral grid: wavenumbers (rad/m). Field layout f[j, i] = f(y_j, x_i); fft2 over (y, x).
        kx = 2.0 * np.pi * np.fft.fftfreq(self.nx, d=self.dx)
        ky = 2.0 * np.pi * np.fft.fftfreq(self.ny, d=self.dy)
        self.KX, self.KY = np.meshgrid(kx, ky)               # both (ny, nx)
        self.K2 = self.KX ** 2 + self.KY ** 2
        self.K4 = self.K2 ** 2

        # Inversion determinant det A = K^2 (K^2 + F1 + F2); zero only at K=0 (handled).
        self._det = self.K2 * (self.K2 + self.F1 + self.F2)
        self._det_safe = np.where(self._det == 0.0, 1.0, self._det)
        self._k0 = (self.K2 == 0.0)                          # the undetermined domain-mean mode

        # 2/3 dealias mask (Orszag): keep |k| <= 2/3 Nyquist in each direction.
        kx_cut = (2.0 / 3.0) * np.pi / self.dx
        ky_cut = (2.0 / 3.0) * np.pi / self.dy
        self._dealias = (np.abs(self.KX) <= kx_cut) & (np.abs(self.KY) <= ky_cut)

    # -- convenience constructors ------------------------------------------- #
    @classmethod
    def symmetric(cls, nx, ny, Lx, Ly, f0, gp, H1, H2, Us, **kw) -> "TwoLayerQG":
        """Build with a **symmetric** shear ``U_1 = +U_s/2, U_2 = −U_s/2`` (zero barotropic mean)."""
        return cls(nx, ny, Lx, Ly, f0, gp, H1, H2, U1=0.5 * Us, U2=-0.5 * Us, **kw)

    # -- derived physical scales -------------------------------------------- #
    @property
    def Us(self) -> float:
        """The vertical shear ``U_1 − U_2`` (m/s) — the available-potential-energy reservoir."""
        return self.U1 - self.U2

    @property
    def Ld(self) -> float:
        """Two-layer deformation radius ``L_d = √(g'H_e)/|f₀|`` (m), ``H_e = H₁H₂/(H₁+H₂)``."""
        He = self.H1 * self.H2 / (self.H1 + self.H2)
        return float(np.sqrt(self.gp * He) / abs(self.f0))

    @property
    def mean_pv_gradients(self) -> tuple[float, float]:
        """The mean PV gradients ``(∂q̄_1/∂y, ∂q̄_2/∂y) = (β + F_1·U_s, β − F_2·U_s)`` (1/(m·s)).

        Derived (module docstring): ``ψ̄_k = −U_k y`` ⟹ ``q̄_k = [β − (−1)^k F_k U_s]·y``. The
        **lower-layer reversal** ``β − F_2 U_s < 0`` (with ``β + F_1 U_s > 0`` always) is the
        Charney–Stern necessary condition for instability ⟹ critical shear ``U_crit = β/F`` — the
        sign convention the linear anchor pins (a flip here is the #1 day-one bug).
        """
        Us = self.Us
        return self.beta + self.F1 * Us, self.beta - self.F2 * Us

    @property
    def critical_shear(self) -> float:
        """Charney–Stern critical shear ``U_crit = β/F_2`` (m/s); ``0`` on an f-plane (β=0)."""
        return self.beta / self.F2

    # -- the spectral PV inversion (q ↔ ψ) ---------------------------------- #
    def invert(self, q: np.ndarray) -> np.ndarray:
        """Recover the streamfunction ``ψ`` (shape ``(2, ny, nx)``) from the PV anomaly ``q``.

        The 2×2 spectral solve ``ψ̂ = A⁻¹ q̂`` (module docstring); the undetermined domain-mean
        (``K=0``) mode is set to zero. Returns a real field.
        """
        qh = np.fft.fft2(q, axes=(-2, -1))
        K2, F1, F2 = self.K2, self.F1, self.F2
        # A^{-1} = (1/det)·[[ -(K2+F2),  -F1     ],
        #                   [ -F2,       -(K2+F1)]]
        psih1 = (-(K2 + F2) * qh[0] - F1 * qh[1]) / self._det_safe
        psih2 = (-F2 * qh[0] - (K2 + F1) * qh[1]) / self._det_safe
        psih1[self._k0] = 0.0
        psih2[self._k0] = 0.0
        psih = np.stack([psih1, psih2])
        return np.real(np.fft.ifft2(psih, axes=(-2, -1)))

    def pv_from_psi(self, psi: np.ndarray) -> np.ndarray:
        """Forward map ``q = A·ψ`` (the PV anomaly from a streamfunction) — the inverse of :meth:`invert`.

        ``q_1 = ∇²ψ_1 + F_1(ψ_2 − ψ_1)``, ``q_2 = ∇²ψ_2 + F_2(ψ_1 − ψ_2)`` (the anomaly part; the
        ``βy`` mean is carried as a background coefficient, not here). Used by the ``q ↔ ψ``
        round-trip plumbing test and for seeding a state from a streamfunction.
        """
        psih = np.fft.fft2(psi, axes=(-2, -1))
        K2, F1, F2 = self.K2, self.F1, self.F2
        qh1 = -(K2 + F1) * psih[0] + F1 * psih[1]
        qh2 = F2 * psih[0] - (K2 + F2) * psih[1]
        qh = np.stack([qh1, qh2])
        return np.real(np.fft.ifft2(qh, axes=(-2, -1)))

    def velocities(self, psi: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """``(u, v)`` from ``ψ``: geostrophic ``u = −∂ψ/∂y``, ``v = ∂ψ/∂x`` (shape ``(2, ny, nx)``)."""
        psih = np.fft.fft2(psi, axes=(-2, -1))
        u = np.real(np.fft.ifft2(-1j * self.KY * psih, axes=(-2, -1)))
        v = np.real(np.fft.ifft2(1j * self.KX * psih, axes=(-2, -1)))
        return u, v

    # -- the LINEAR instability operator (the tight analytic anchor) -------- #
    def linear_frequencies(self, k: float, l: float = 0.0) -> np.ndarray:
        """The two complex QG wave frequencies ``ω`` at zonal/meridional wavenumber ``(k, l)``.

        Linearizing (drop the eddy–eddy ``J``) gives ``ω q̂ = k(U·A + Q)ψ̂`` with ``q̂ = A ψ̂`` ⟹
        the 2×2 eigenproblem ``ω ψ̂ = k·A⁻¹(U A + Q) ψ̂``, ``A`` the inversion matrix, ``U =
        diag(U_1,U_2)``, ``Q = diag(∂q̄_1/∂y, ∂q̄_2/∂y)``. Built directly from the equations.
        """
        K2 = k * k + l * l
        if K2 == 0.0:
            return np.zeros(2, dtype=complex)
        F1, F2 = self.F1, self.F2
        A = np.array([[-(K2 + F1), F1], [F2, -(K2 + F2)]], dtype=complex)
        U = np.diag([self.U1, self.U2]).astype(complex)
        Q1, Q2 = self.mean_pv_gradients
        Q = np.diag([Q1, Q2]).astype(complex)
        M = k * np.linalg.solve(A, U @ A + Q)
        return np.linalg.eigvals(M)

    def growth_rate(self, k: float, l: float = 0.0) -> float:
        """Linear growth rate ``σ = max Im(ω)`` (1/s) at ``(k, l)`` for the configured shear."""
        return float(np.max(self.linear_frequencies(k, l).imag))

    def most_unstable(self, n: int = 400, k_lo: float = 1e-3,
                      k_hi: float = 1.5) -> tuple[float, float]:
        """Scan zonal wavenumber ``k`` (``l = 0``) and return ``(k*, σ_max)``.

        ``k`` is scanned over ``[k_lo, k_hi]·√(2F)`` where ``√(2F)`` (``F = f₀²/(g'H₁)``) is the
        β=0 short-wave cutoff ``K² = 2F``; the default bracket spans the unstable band.
        """
        k_cut = np.sqrt(2.0 * self.F1)
        ks = np.linspace(k_lo, k_hi, n) * k_cut
        sig = np.array([self.growth_rate(k, 0.0) for k in ks])
        i = int(np.argmax(sig))
        return float(ks[i]), float(sig[i])

    # -- the nonlinear right-hand side (pseudospectral) --------------------- #
    def _invert_spectral(self, qh: np.ndarray) -> np.ndarray:
        """Spectral inversion ``ψ̂ = A⁻¹ q̂`` (the K=0 mode zeroed) — the inner kernel of :meth:`invert`."""
        K2, F1, F2 = self.K2, self.F1, self.F2
        psih1 = (-(K2 + F2) * qh[0] - F1 * qh[1]) / self._det_safe
        psih2 = (-F2 * qh[0] - (K2 + F1) * qh[1]) / self._det_safe
        psih1[self._k0] = 0.0
        psih2[self._k0] = 0.0
        return np.stack([psih1, psih2])

    def _rhs(self, q: np.ndarray) -> np.ndarray:
        """Tendency ``∂q/∂t`` (physical, shape ``(2, ny, nx)``).

        Assembled in spectral space: the mean-flow Doppler ``−U_k ∂_x q``, the baroclinic/β source
        ``−(∂q̄_k/∂y)∂_x ψ`` (the linear instability terms), the **eddy–eddy** Jacobian
        ``−J(ψ, q) = −(u q_x + v q_y)`` computed pseudospectrally and **2/3-dealiased**, the
        hyperviscosity ``−ν₄∇⁴q → −ν₄K⁴q̂``, and the lower-layer bottom drag ``−r∇²ψ_2 → +rK²ψ̂_2``.
        Dropping the Jacobian leaves exactly the linear operator of :meth:`linear_frequencies`.
        """
        KX, KY, K2, K4 = self.KX, self.KY, self.K2, self.K4
        qh = np.fft.fft2(q, axes=(-2, -1))
        psih = self._invert_spectral(qh)
        # velocities and PV gradients in physical space for the Jacobian J = u·∇q
        u = np.real(np.fft.ifft2(-1j * KY * psih, axes=(-2, -1)))
        v = np.real(np.fft.ifft2(1j * KX * psih, axes=(-2, -1)))
        qx = np.real(np.fft.ifft2(1j * KX * qh, axes=(-2, -1)))
        qy = np.real(np.fft.ifft2(1j * KY * qh, axes=(-2, -1)))
        adv_h = np.fft.fft2(u * qx + v * qy, axes=(-2, -1)) * self._dealias
        # background coefficients (broadcast over the leading layer axis)
        Q1, Q2 = self.mean_pv_gradients
        Q = np.array([Q1, Q2])[:, None, None]
        U = np.array([self.U1, self.U2])[:, None, None]
        dqh = -(U * (1j * KX) * qh) - (Q * (1j * KX) * psih) - adv_h - self.nu4 * K4 * qh
        if self.r_drag > 0.0:
            dqh[1] = dqh[1] + self.r_drag * K2 * psih[1]      # −r∇²ψ_2 = +rK²ψ̂_2 on the lower layer
        return np.real(np.fft.ifft2(dqh, axes=(-2, -1)))

    # -- time stepping (SSP-RK3) -------------------------------------------- #
    def max_dt(self, state: QGState, safety: float = 0.3) -> float:
        """CFL-stable step (s): ``safety · min(min(dx,dy)/(π·speed), 0.035·min(dx,dy)⁴/ν₄)``.

        ``speed = |u|_max + max|U_k|`` is the fastest advective signal; the ``π`` factor is the
        spectral max-wavenumber (``k_max = π/dx``). With ``ν₄ > 0`` the biharmonic term adds its own
        ``dx⁴/ν₄`` limit (the same fraction the de-risking spike used). QG carries **no fast external
        gravity wave** (the rigid lid filtered it), so this advective step is far larger than the
        free-surface engine's — the cost win that made the saturated runs affordable.
        """
        u, v = self.velocities(self.invert(state.q))
        speed = float(np.max(np.sqrt(u ** 2 + v ** 2))) + max(abs(self.U1), abs(self.U2))
        lim = [min(self.dx, self.dy) / (np.pi * max(speed, 1e-30))]
        if self.nu4 > 0.0:
            lim.append(0.035 * min(self.dx, self.dy) ** 4 / self.nu4)
        return float(safety * min(lim))

    def step(self, state: QGState, dt: float) -> QGState:
        """Advance ``state`` by one SSP-RK3 step ``dt`` (s); returns a new state (no mutation).

        Raises if ``dt`` exceeds the CFL limit (the explicit-solver analogue of the diffusion
        engine's stability guarantee — conditional, so enforced). The convex-combination form is the
        Shu–Osher SSP-RK3, mirroring :class:`engines.fluid.ShallowWater`.
        """
        if dt <= 0.0:
            raise ValueError("dt must be positive")
        cfl = self.max_dt(state, safety=1.0)
        if dt > cfl:
            raise ValueError(f"dt={dt:g} exceeds the CFL stability limit {cfl:g}; reduce dt (use max_dt())")
        q = state.q
        k1 = self._rhs(q)
        q1 = q + dt * k1
        k2 = self._rhs(q1)
        q2 = 0.75 * q + 0.25 * (q1 + dt * k2)
        k3 = self._rhs(q2)
        qn = (1.0 / 3.0) * q + (2.0 / 3.0) * (q2 + dt * k3)
        return QGState(q=qn)

    def solve(self, state: QGState, t_end: float, dt: Optional[float] = None,
              safety: float = 0.3, recompute_dt_every: int = 50) -> QGState:
        """March from ``0`` to ``t_end`` (s). With ``dt=None`` adapts the CFL step as the eddies spin
        up (recomputed every ``recompute_dt_every`` steps, since a growing/saturating flow accelerates).
        """
        if t_end < 0.0:
            raise ValueError("t_end must be non-negative")
        s = state
        t = 0.0
        step_dt = dt if dt is not None else self.max_dt(s, safety=safety)
        n = 0
        while t < t_end - 1e-12 * max(1.0, t_end):
            if dt is None and n % recompute_dt_every == 0:
                step_dt = self.max_dt(s, safety=safety)
            h = min(step_dt, t_end - t)
            s = self.step(s, h)
            t += h
            n += 1
        return s

    # -- seeding ------------------------------------------------------------ #
    def random_state(self, amplitude: float = 1e-3, seed: int = 0,
                     kmax_factor: float = 2.5) -> QGState:
        """A small band-limited random **baroclinic** seed (PV anomaly) to trigger the instability.

        Builds a random baroclinic streamfunction ``τ = ψ_1 − ψ_2`` band-limited to ``|K| <
        kmax_factor·√(2F)`` (``√(2F)`` the deformation wavenumber, the cutoff scale — robust at any
        shear, including the zero-shear decay test), scaled so the initial eddy meridional velocity
        has rms ``amplitude`` (m/s), then maps it to ``q`` via :meth:`pv_from_psi`.
        """
        rng = np.random.default_rng(seed)
        k_ref = np.sqrt(2.0 * self.F1)
        fld = rng.standard_normal((self.ny, self.nx))
        fh = np.fft.fft2(fld)
        fh[np.sqrt(self.K2) > kmax_factor * k_ref] = 0.0
        fh[self._k0] = 0.0
        tau = np.real(np.fft.ifft2(fh))
        psi = np.stack([0.5 * tau, -0.5 * tau])
        _, v = self.velocities(psi)
        vrms = float(np.sqrt(np.mean(v ** 2)))
        psi = psi * (amplitude / (vrms + 1e-300))
        return QGState(q=self.pv_from_psi(psi))

    # -- diagnostics -------------------------------------------------------- #
    def eddy_kinetic_energy(self, state: QGState) -> float:
        """Eddy kinetic energy ``½⟨u² + v²⟩`` (m²/s²; mean over both layers and the domain).

        The fields are the eddies (the mean flow is a background coefficient), so the prognostic
        ``(u, v)`` are already perturbations — the spin-up/saturation/plateau diagnostic.
        """
        u, v = self.velocities(self.invert(state.q))
        return float(0.5 * np.mean(u ** 2 + v ** 2))

    def v_rms(self, state: QGState) -> float:
        """RMS eddy meridional velocity ``√⟨v²⟩`` over both layers (m/s) — the eddy velocity scale."""
        _, v = self.velocities(self.invert(state.q))
        return float(np.sqrt(np.mean(v ** 2)))

    def bulk_eddy_flux(self, state: QGState) -> tuple[float, float]:
        """Per-layer domain-bulk meridional eddy flux of the baroclinic streamfunction ``⟨v_k'·τ'⟩``.

        ``τ = ψ_1 − ψ_2`` is proportional to the interface displacement ``η = (f₀/g')τ`` (the
        dynamical temperature). The down-gradient thickness flux ``⟨v_k'·h_k'⟩`` (``h_1' ∝ −τ'``,
        ``h_2' ∝ +τ'``) reduces — once the common ``f₀/g'`` constant is dropped (it cancels in
        :meth:`kappa_eff`, advisor) — to ``⟨v_k'·τ'⟩`` for **both** layers (the layer-1 sign and the
        layer-1 mean-gradient sign both flip, so the ratio is the same). Primes = domain-mean removed.
        ``> 0`` ⟺ **down-gradient** in both layers. Returned as a domain-bulk scalar per layer (a
        meridional channel — the named BC extension — would be needed for a latitude-resolved κ(y));
        its **time series** is the irreversibility diagnostic (temporal persistence, not spatial).
        """
        psi = self.invert(state.q)
        _, v = self.velocities(psi)
        vp = v - v.mean(axis=(-2, -1), keepdims=True)
        tau = psi[0] - psi[1]
        tau = tau - tau.mean()
        return float(np.mean(vp[0] * tau)), float(np.mean(vp[1] * tau))

    def kappa_eff(self, state: QGState) -> tuple[float, float]:
        """Per-layer domain-bulk down-gradient eddy thickness diffusivity ``(κ_1, κ_2)`` (m²/s).

        ``κ_k = −⟨v_k'h_k'⟩/(∂h̄_k/∂y)``. With symmetric shear ``∂h̄_1/∂y ∝ +U_s``,
        ``∂h̄_2/∂y ∝ −U_s`` (the proportionality constant cancels against :meth:`bulk_eddy_flux`),
        so ``κ_k = ⟨v_k'·τ'⟩/U_s`` for **both** layers — and both must come out **positive and
        comparable** if the flux is genuinely down-gradient (the cross-layer consistency check). The
        physical (dimensional) ``κ`` carries the dropped ``f₀/g'`` constant; the **dimensionless**
        ``κ_eff/(v'_rms·L_d)`` (the pre-registered discriminator) does not, so it is reported on the
        constant-free form (idealized ``κ`` is intrinsically far below Earth's — see the module docstring).
        """
        if self.Us == 0.0:
            raise ValueError("no mean shear (U_s = 0) — the thickness diffusivity is undefined")
        f1, f2 = self.bulk_eddy_flux(state)
        return f1 / self.Us, f2 / self.Us
