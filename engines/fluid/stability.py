"""Two-layer shallow-water **linear stability** — the analytic baroclinic growth-rate anchor.

This is the *tight* anchor for the N-layer engine extension (rung 3 of the GCM climb):
the growth rate :class:`LayeredShallowWater` must reproduce on an unstable thermal-wind
basic state. Following the rung discipline (match a first-principles analytic relation, not
a recalled formula), the rate comes from **linearizing the two-layer SW equations and
rooting the 6×6 dispersion matrix** — *not* from a recalled QG-Phillips quartic. The matrix
is **self-validating** (and the test suite asserts it):

* at **zero shear** every mode is **neutral to machine precision** (``max|Im ω| ~ 1e-19``), and
* the gravity branches sit exactly on the **two-layer Poincaré dispersions** — external
  ``ω² = f₀² + g·H_tot·k²`` and internal ``ω² = f₀² + g'·H_e·k²`` (``H_e = H₁H₂/(H₁+H₂)``).

Above critical shear it gives ``σ(k) = max Im ω`` with a **short-wave cutoff** and a
most-unstable wavelength a few × the internal deformation radius. The baroclinic terms are
externally anchored (the zero-shear check leaves them untested): the max-growth coefficient
``σ_max ≈ 0.30·U_s/L_d`` matches the **Eady** model (a wholly independent derivation, 0.310)
to ~2 %.

**Scope (the honest edges, banked in the spike).** The operator is **f-plane** — ``β`` is
*not* in the perturbation matrix, so the model is **Eady-like: no critical shear, unstable
for all shear** (a *finite* critical shear ``U_s,crit = β·g'H/f₀²`` needs a β-capable
PV-gradient treatment = the named within-rung extension). Symmetric shear ``U₁ = +U_s/2``,
``U₂ = −U_s/2`` (zero barotropic mean → pure baroclinic). All SI.

Sources (extending ``[[shallow-water-source]]``): two-layer SW formulation — Vallis 2017
*AOFD* / Cushman-Roisin & Beckers; baroclinic instability — Phillips 1954 / Eady 1949 /
Charney 1947.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class TwoLayerStability:
    """Linear-stability operator for two stacked free-surface SW layers (densities ρ₁<ρ₂).

    Parameters
    ----------
    f0 : float
        Coriolis parameter (1/s); the operator is f-plane (β absent — see the module docstring).
    g : float
        Gravity (m/s²) — sets the external/barotropic mode.
    gp : float
        Reduced gravity ``g' = g·Δρ/ρ`` (m/s²) at the internal interface — sets the baroclinic mode.
    H1, H2 : float
        Rest thicknesses of the upper / lower layer (m).
    """

    f0: float
    g: float
    gp: float
    H1: float
    H2: float

    # -- derived scales ----------------------------------------------------- #
    @property
    def He(self) -> float:
        """Equivalent depth of the internal mode ``H₁H₂/(H₁+H₂)`` (m)."""
        return self.H1 * self.H2 / (self.H1 + self.H2)

    @property
    def Htot(self) -> float:
        return self.H1 + self.H2

    @property
    def Ld_int(self) -> float:
        """Internal (baroclinic) deformation radius ``√(g'H_e)/|f₀|`` (m)."""
        return np.sqrt(self.gp * self.He) / abs(self.f0)

    def poincare_external(self, k: float) -> float:
        """External (barotropic) Poincaré frequency ``√(f₀² + g·H_tot·k²)`` (1/s)."""
        return float(np.sqrt(self.f0 ** 2 + self.g * self.Htot * k ** 2))

    def poincare_internal(self, k: float) -> float:
        """Internal (baroclinic) Poincaré frequency ``√(f₀² + g'·H_e·k²)`` (1/s)."""
        return float(np.sqrt(self.f0 ** 2 + self.gp * self.He * k ** 2))

    # -- thermal-wind basic state ------------------------------------------- #
    def basic_state_gradients(self, U1: float, U2: float) -> tuple[float, float]:
        """Basic-state thickness gradients ``G_k = d(H̄_k)/dy`` from thermal wind (1/–).

        The mean zonal flows ``U_k`` are in geostrophic balance with sloping interfaces:
        the top surface slope ``dη₀/dy = −f₀U₁/g`` and the internal interface slope
        ``dη₁/dy = f₀(U₁−U₂)/g'``. With ``h₁ = H₁ + η₀ − η₁`` and ``h₂ = H₂ + η₁``,

            ``G₁ = dη₀/dy − dη₁/dy``,   ``G₂ = dη₁/dy``.

        These are the coefficients of the baroclinic term ``−G_k·v'`` in the linearized
        continuity (advection of the mean thickness by the perturbation meridional flow) —
        the engine injects exactly the same ``G_k`` as its optional background.
        """
        deta0 = -self.f0 * U1 / self.g
        deta1 = self.f0 * (U1 - U2) / self.gp
        return deta0 - deta1, deta1

    # -- the dispersion matrix ---------------------------------------------- #
    def dispersion_matrix(self, k: float, l: float, U1: float, U2: float) -> np.ndarray:
        """The 6×6 complex matrix ``M(k,l)`` whose eigenvalues are the wave frequencies ``ω``.

        Perturbation amplitudes ``X = (u₁,v₁,h₁, u₂,v₂,h₂) ~ exp[i(kx+ly−ωt)]`` of the
        linearized two-layer SW equations about a uniform-flow / sloped-interface basic state:

            (−iω + ikU_k) u_k − f v_k                       = −ik P_k
            (−iω + ikU_k) v_k + f u_k                       = −il P_k
            (−iω + ikU_k) h_k + G_k v_k + H_k(ik u_k+il v_k) = 0
            P₁ = g(h₁+h₂),   P₂ = g(h₁+h₂) + g' h₂

        rearranged via ``−iω·φ = R  ⇒  ω·φ = i·R`` into ``ω X = M X``. ``σ = max Im ω``.
        Built directly from the equations (audit each row against the system above) — *not*
        a recalled quartic.
        """
        g, gp, f = self.g, self.gp, self.f0
        H1, H2 = self.H1, self.H2
        G1, G2 = self.basic_state_gradients(U1, U2)
        M = np.zeros((6, 6), dtype=complex)
        u1, v1, h1, u2, v2, h2 = 0, 1, 2, 3, 4, 5
        # --- layer 1 ---  P₁ = g·h1 + g·h2
        M[u1, u1] += k * U1
        M[u1, v1] += 1j * f
        M[u1, h1] += k * g
        M[u1, h2] += k * g
        M[v1, v1] += k * U1
        M[v1, u1] += -1j * f
        M[v1, h1] += l * g
        M[v1, h2] += l * g
        M[h1, h1] += k * U1
        M[h1, u1] += H1 * k
        M[h1, v1] += H1 * l - 1j * G1
        # --- layer 2 ---  P₂ = g·h1 + (g+g')·h2
        M[u2, u2] += k * U2
        M[u2, v2] += 1j * f
        M[u2, h1] += k * g
        M[u2, h2] += k * (g + gp)
        M[v2, v2] += k * U2
        M[v2, u2] += -1j * f
        M[v2, h1] += l * g
        M[v2, h2] += l * (g + gp)
        M[h2, h2] += k * U2
        M[h2, u2] += H2 * k
        M[h2, v2] += H2 * l - 1j * G2
        return M

    def frequencies(self, k: float, l: float, U1: float = 0.0, U2: float = 0.0) -> np.ndarray:
        """The six complex wave frequencies ``ω`` at ``(k, l)`` for layer flows ``U1, U2``."""
        return np.linalg.eigvals(self.dispersion_matrix(k, l, U1, U2))

    def growth_rate(self, k: float, l: float, Us: float) -> tuple[float, complex]:
        """``(σ, ω)`` — max growth ``σ = max Im ω`` and its frequency, for **symmetric** shear.

        Symmetric ``U₁ = +U_s/2``, ``U₂ = −U_s/2`` (zero barotropic mean → pure baroclinic).
        """
        w = self.frequencies(k, l, 0.5 * Us, -0.5 * Us)
        i = int(np.argmax(w.imag))
        return float(w[i].imag), complex(w[i])

    def most_unstable(self, Us: float, l: float = 0.0, n: int = 200,
                      k_lo: float = 0.02, k_hi: float = 1.5) -> tuple[float, float]:
        """Scan zonal wavenumber ``k`` and return ``(k*, σ_max)`` for symmetric shear ``Us``.

        ``k`` is scanned over ``[k_lo, k_hi]·k_cut`` where the short-wave cutoff scale is set by
        the internal deformation wavenumber ``k_cut = √2·f₀/√(g'H₁)`` (the β=0 Phillips cutoff
        ``K² = 2F``, ``F = f₀²/(g'H₁)``). Defaults bracket the most-unstable mode for the
        idealized rung-3 parameters.
        """
        k_cut = np.sqrt(2) * abs(self.f0) / np.sqrt(self.gp * self.H1)
        ks = np.linspace(k_lo, k_hi, n) * k_cut
        sig = np.array([self.growth_rate(k, l, Us)[0] for k in ks])
        i = int(np.argmax(sig))
        return float(ks[i]), float(sig[i])
