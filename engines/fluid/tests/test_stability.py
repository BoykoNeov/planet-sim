"""Two-layer SW linear-stability solver — the analytic baroclinic growth-rate anchor.

:mod:`engines.fluid.stability` linearizes the two-layer shallow-water equations and roots the
6×6 dispersion matrix for the growing mode (built from the equations, **not** a recalled
quartic). These tests are the solver's *self-validation* — the same first-principles checks the
de-risking spike used, now pinned — so that the engine's measured growth rate
(``test_layered.py``) can be trusted against it:

* **Zero shear → neutral to machine precision.** No basic-state shear ⇒ no available potential
  energy ⇒ every mode must be exactly neutral (``Im ω = 0`` to round-off). The decisive check
  that the linearization carries no spurious growth.
* **Recovers both two-layer Poincaré dispersions.** At zero shear the gravity branches sit on the
  external ``ω² = f₀² + g·H_tot·k²`` and internal ``ω² = f₀² + g'·H_e·k²`` relations — the
  rotation/coupling check (a wrong Coriolis or Montgomery stack would miss them).
* **Short-wave cutoff.** Baroclinic instability has a high-``k`` cutoff (``σ → 0`` past it).
* **Eady coefficient — the external anchor on the baroclinic terms.** The zero-shear checks leave
  the baroclinic ``G_k`` terms untested, so they are anchored against the **Eady** model (a wholly
  independent derivation): ``σ_max ≈ 0.31·U_s/L_d`` with ``L_d = √(g'H)/f₀``.
* **f-plane ⇒ no critical shear.** ``β`` is absent from the operator, so the model is Eady-like:
  growth ∝ shear with **no threshold** (a finite critical shear needs a β-capable PV-gradient
  treatment = the named within-rung extension).
"""
import numpy as np
import pytest

from engines.fluid import TwoLayerStability


# Idealized rung-3 parameters (the spike's: modest √(gH_tot), resolvable internal L_d).
PARAMS = dict(f0=1.0e-4, g=10.0, gp=0.2, H1=500.0, H2=500.0)


@pytest.fixture
def st():
    return TwoLayerStability(**PARAMS)


# --------------------------------------------------------------------------- #
# Zero shear → neutral to machine precision (the decisive first-principles check)
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("lam_km", [100.0, 300.0, 1000.0])
def test_zero_shear_is_neutral_to_machine_precision(st, lam_km):
    k = 2 * np.pi / (lam_km * 1e3)
    w = st.frequencies(k, 0.0, U1=0.0, U2=0.0)
    assert np.max(np.abs(w.imag)) < 1e-15


# --------------------------------------------------------------------------- #
# Recovers both two-layer Poincaré dispersions (rotation + Montgomery coupling)
# --------------------------------------------------------------------------- #
def test_recovers_two_layer_poincare_dispersions(st):
    """The gravity branches sit on the two-layer Poincaré dispersions. The match is to ~O(Δρ/ρ),
    **not** exact: the Poincaré formulas are the *decoupled / rigid-lid* limits, and the exact
    free-surface eigenvalues differ by the small free-surface coupling (the same extra degree of
    freedom that makes the free-surface SW model ~4 % more unstable than rigid-lid Phillips). The
    decisive *exact* check is the zero-shear neutrality above (Im ω → machine zero)."""
    k = 2 * np.pi / 200e3
    freqs = np.sort(np.abs(st.frequencies(k, 0.0).real))   # zero shear
    # six roots: two ~zero (geostrophic), two internal (±), two external (±)
    assert freqs[0] == pytest.approx(0.0, abs=1e-9)
    assert freqs[2] == pytest.approx(st.poincare_internal(k), rel=1e-2)
    assert freqs[4] == pytest.approx(st.poincare_external(k), rel=1e-2)


# --------------------------------------------------------------------------- #
# Short-wave cutoff — σ → 0 past the cutoff wavenumber
# --------------------------------------------------------------------------- #
def test_short_wave_cutoff(st):
    Us = 4.0
    kstar, smax = st.most_unstable(Us)
    assert smax > 0.0
    # β=0 cutoff K² = 2F, F = f₀²/(g'H₁); past it the mode is neutral
    k_cut = np.sqrt(2) * st.f0 / np.sqrt(st.gp * st.H1)
    assert st.growth_rate(1.5 * k_cut, 0.0, Us)[0] < 1e-12 * smax + 1e-18


# --------------------------------------------------------------------------- #
# Eady coefficient — the independent anchor on the baroclinic G_k terms
# --------------------------------------------------------------------------- #
def test_eady_growth_coefficient(st):
    """σ_max ≈ 0.31·U_s/L_d (L_d = √(g'H₁)/f₀) — within a few % of the Eady model, which is a
    *wholly independent* derivation (continuous stratification). This anchors the baroclinic
    terms that the zero-shear Poincaré checks leave untested."""
    Us = 4.0
    _, smax = st.most_unstable(Us)
    Ld = np.sqrt(st.gp * st.H1) / st.f0
    coeff = smax / (Us / Ld)
    assert coeff == pytest.approx(0.31, abs=0.02)


# --------------------------------------------------------------------------- #
# f-plane ⇒ Eady-like: growth ∝ shear with no critical-shear threshold
# --------------------------------------------------------------------------- #
def test_fplane_has_no_critical_shear(st):
    """The operator is f-plane (β absent), so it is unstable for *all* shear and the most-unstable
    growth scales ~linearly with shear: σ/U_s is roughly constant (no threshold). A finite critical
    shear is the named β-capable extension."""
    ratios = [st.most_unstable(U)[1] / U for U in (0.1, 1.0, 4.0)]
    # all positive (unstable at every shear) and within ~10 % of each other (≈ linear)
    assert all(r > 0 for r in ratios)
    assert max(ratios) / min(ratios) < 1.10


# --------------------------------------------------------------------------- #
# Most-unstable wavelength is a few × the deformation radius (sanity)
# --------------------------------------------------------------------------- #
def test_most_unstable_wavelength_scale(st):
    kstar, _ = st.most_unstable(4.0)
    lam = 2 * np.pi / kstar
    Ld = np.sqrt(st.gp * st.H1) / st.f0
    assert 4.0 < lam / Ld < 9.0          # ~7× at these params (≈ 700 km)


# --------------------------------------------------------------------------- #
# Thermal-wind gradients are the coefficients the engine injects (cross-link)
# --------------------------------------------------------------------------- #
def test_basic_state_gradients_symmetric_shear(st):
    """Symmetric shear U₁=+Us/2, U₂=−Us/2. The lower-layer gradient is the internal-interface
    slope; the upper-layer gradient adds the free-surface slope. These G_k are exactly what
    LayeredShallowWater.thermal_wind injects (asserted in test_layered.py)."""
    Us = 4.0
    G1, G2 = st.basic_state_gradients(0.5 * Us, -0.5 * Us)
    # interface slope dη₁/dy = f₀·Us/g'; G₂ equals it
    assert G2 == pytest.approx(st.f0 * Us / st.gp, rel=1e-12)
    # surface slope dη₀/dy = −f₀·(Us/2)/g; G₁ = dη₀/dy − dη₁/dy
    assert G1 == pytest.approx(-st.f0 * (0.5 * Us) / st.g - G2, rel=1e-12)
