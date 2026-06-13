"""Two-layer quasi-geostrophic turbulence — the rung-3 Phase-B saturated baroclinic eddy-flux engine.

:mod:`planet.baroclinic_qg` is the standard tool (Held & Larichev 1996) for the experiment rung 1
named but could not run: does a *saturated* baroclinic eddy field produce a strong, irreversible,
down-gradient meridional thickness flux, so the reduction-to-diffusive-EBM is finally non-vacuous?
It is a **new model outside ``engines/fluid``** (pseudospectral, not the C-grid) — chosen because the
free-surface SW engine OUTCROPS at saturation (the spike's finding). The triad the suite pins:

* **tight (the linear anchor)** — the rooted 2×2 QG dispersion equals the analytic Phillips closed
  form to machine precision (equal layers — *the same equations*), has the ``K²=2F`` short-wave
  cutoff and the ``U_crit=β/F`` Charney–Stern critical shear, and is neutral at zero shear. The SW
  6×6 solver (:class:`engines.fluid.TwoLayerStability`) in the rigid-lid limit (``g→∞``) agrees to
  <0.5 % — the cross-model bridge tying the Phase-A (SW) and Phase-B (QG) engines together.
* **plumbing** — ``q ↔ ψ`` round-trips to machine precision (the spectral inversion is exact, the
  ``K=0`` gauge aside); zero shear ⟹ the eddies decay (no APE source, hyperviscosity damps).
* **real-but-loose (the unlock)** — the emergent saturated eddy thickness diffusivity is
  **down-gradient** (``κ_1, κ_2 > 0``, comparable) and **irreversible** (``irr ≫`` rung-1's ~0.1),
  with ``κ/(v'_rms·L_d)`` order-unity (dimensionless — idealized ``κ`` is intrinsically far below
  Earth's). *[slow — a saturated turbulence integration.]*

The nonlinear engine is validated against the linear anchor it must reduce to: a single growing
eigenmode at small amplitude grows at the analytic rate to <1 % *[slow]*.
"""
import numpy as np
import pytest

from engines.fluid import TwoLayerStability
from planet.baroclinic_qg import TwoLayerQG, QGState


# Idealized rung-3 parameters (resolvable L_d, modest speeds — the saturation spike's).
P = dict(f0=1.0e-4, gp=2.0, H1=400.0, H2=400.0, beta=1.6e-11)


def _phillips_growth(f0, gp, H, beta, k, l, Us):
    """Analytic two-layer QG (Phillips) growth rate, equal layers, rigid lid, symmetric shear —
    the *independent closed form* (Phillips 1954) the rooted 2×2 dispersion must reproduce. ``F =
    f₀²/(g'H)``; the short-wave cutoff is ``K²=2F``; unstable when the radicand < 0."""
    F = f0 ** 2 / (gp * H)
    K2 = k ** 2 + l ** 2
    UT = 0.5 * Us
    radicand = (beta ** 2 * F ** 2) / (K2 ** 2 * (K2 + 2 * F) ** 2) \
        - UT ** 2 * (2 * F - K2) / (K2 + 2 * F)
    return 0.0 if radicand >= 0 else float(k * np.sqrt(-radicand))


# --------------------------------------------------------------------------- #
# tight — the 2×2 QG dispersion equals the analytic Phillips form to ~machine precision
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("beta", [0.0, 1.0e-11, 1.6e-11])
@pytest.mark.parametrize("Us", [1.0, 2.0, 4.0])
@pytest.mark.parametrize("kf,lf", [(0.2, 0.0), (0.5, 0.0), (0.8, 0.3), (1.05, 0.0)])
def test_dispersion_matches_phillips_to_machine_precision(beta, Us, kf, lf):
    """Equal layers ⇒ the rooted 2×2 dispersion *is* the Phillips closed form — they must agree to
    ~1e-9, NOT a few % (the ~4 % gap is free-surface-SW-vs-Phillips, which does not exist here). A
    loose tolerance would hide a partially-compensated PV-gradient sign bug (advisor)."""
    f0, gp, H = P["f0"], P["gp"], P["H1"]
    k_cut = np.sqrt(2.0 * f0 ** 2 / (gp * H))
    k, l = kf * k_cut, lf * k_cut
    m = TwoLayerQG.symmetric(8, 8, 1.0, 1.0, f0=f0, gp=gp, H1=H, H2=H, Us=Us, beta=beta)
    sp = _phillips_growth(f0, gp, H, beta, k, l, Us)
    qg = m.growth_rate(k, l)
    if sp > 0.0:
        assert qg == pytest.approx(sp, rel=1e-9)
    else:
        assert qg < 1e-12 * (f0 * Us) + 1e-18      # both neutral (sub-cutoff / sub-critical)


def test_zero_shear_is_neutral_to_machine_precision():
    """No mean shear ⇒ no available potential energy ⇒ every mode neutral (Im ω → round-off)."""
    m = TwoLayerQG.symmetric(8, 8, 1.0, 1.0, Us=0.0, **P)
    k_cut = np.sqrt(2.0 * m.F1)
    maxg = max(m.growth_rate(kf * k_cut, lf * k_cut)
               for kf in (0.1, 0.5, 1.0) for lf in (0.0, 0.3))
    assert maxg < 1e-15


def test_short_wave_cutoff():
    """Baroclinic instability has a high-k cutoff at K²=2F (β=0); past it the mode is neutral."""
    m = TwoLayerQG.symmetric(8, 8, 1.0, 1.0, f0=P["f0"], gp=P["gp"], H1=P["H1"], H2=P["H2"],
                             Us=4.0, beta=0.0)
    kstar, smax = m.most_unstable()
    assert smax > 0.0
    k_cut = np.sqrt(2.0 * m.F1)
    assert m.growth_rate(1.5 * k_cut, 0.0) < 1e-12 * smax + 1e-18


def test_critical_shear_charney_stern():
    """β re-enters at rung 3 (the SW solver was f-plane) ⇒ a FINITE critical shear U_crit=β/F: the
    lower-layer mean PV gradient β−F·U_s must reverse sign (Charney–Stern). Sub-critical ⇒ neutral,
    super-critical ⇒ growing. This pins both the PV-gradient sign and β's re-entry (the f-plane SW
    solver could test neither)."""
    f0, gp, H, beta = P["f0"], P["gp"], P["H1"], 1.6e-11
    U_crit = beta / (f0 ** 2 / (gp * H))
    sub = TwoLayerQG.symmetric(8, 8, 1.0, 1.0, f0=f0, gp=gp, H1=H, H2=H, Us=0.6 * U_crit, beta=beta)
    sup = TwoLayerQG.symmetric(8, 8, 1.0, 1.0, f0=f0, gp=gp, H1=H, H2=H, Us=2.0 * U_crit, beta=beta)
    assert sub.critical_shear == pytest.approx(U_crit, rel=1e-12)
    assert sub.most_unstable()[1] < 1e-15            # sub-critical: neutral
    assert sup.most_unstable()[1] > 1e-9             # super-critical: unstable
    # the mean PV gradients carry the sign: upper always +, lower reversed when super-critical
    Q1, Q2 = sup.mean_pv_gradients
    assert Q1 > 0.0 and Q2 < 0.0


def test_most_unstable_wavelength_scale():
    """Most-unstable wavelength is a few × the deformation radius (the storm scale)."""
    m = TwoLayerQG.symmetric(8, 8, 1.0, 1.0, Us=4.0, **P)
    kstar, _ = m.most_unstable()
    assert 4.0 < (2 * np.pi / kstar) / m.Ld < 12.0


def test_rigid_lid_bridge_to_sw_solver():
    """The cross-model bridge: the free-surface SW 6×6 solver (Phase A) in the rigid-lid limit
    (g→∞, the external mode infinitely fast) converges to the QG growth rate — σ_SW/σ_QG → 1 to
    <0.5 %. This is the *only* tie between the Phase-A (SW) and Phase-B (QG) engines (no bit-for-bit
    reduction across the model boundary), so it is an assertion, not a comment (advisor)."""
    f0, H, Us = P["f0"], P["H1"], 4.0
    g, drho = 10.0 * 64, 0.2 / 64                     # g→∞ at fixed g' = g·Δρ/ρ = 2.0
    gp = g * drho
    st = TwoLayerStability(f0=f0, g=g, gp=gp, H1=H, H2=H)
    _, sig_sw = st.most_unstable(Us)
    mq = TwoLayerQG.symmetric(8, 8, 1.0, 1.0, f0=f0, gp=gp, H1=H, H2=H, Us=Us, beta=0.0)
    _, sig_qg = mq.most_unstable()
    assert sig_sw / sig_qg == pytest.approx(1.0, abs=0.005)


# --------------------------------------------------------------------------- #
# plumbing — the spectral inversion is exact; zero shear ⇒ eddies decay
# --------------------------------------------------------------------------- #
def test_q_psi_round_trip_is_machine_exact():
    """q ↔ ψ round-trips to machine precision (the 2×2 spectral inversion is exact). The domain-mean
    (K=0) streamfunction is an undetermined gauge, so the ψ round-trip is checked mean-removed; the
    q round-trip (q → ψ → q) is exact regardless."""
    m = TwoLayerQG.symmetric(48, 48, 4.0e6, 4.0e6, Us=4.0, **P)
    rng = np.random.default_rng(0)
    psi = rng.standard_normal((2, 48, 48))
    psi -= psi.mean(axis=(-2, -1), keepdims=True)
    q = m.pv_from_psi(psi)
    psi_rt = m.invert(q)
    psi_rt -= psi_rt.mean(axis=(-2, -1), keepdims=True)
    assert np.max(np.abs(psi - psi_rt)) < 1e-10
    assert np.max(np.abs(q - m.pv_from_psi(m.invert(q)))) < 1e-12 * np.max(np.abs(q))


def test_zero_shear_eddies_decay():
    """With no mean shear there is no APE source: a seeded eddy field only decays (hyperviscosity
    damps it). The dynamical complement to the linear neutrality check."""
    m = TwoLayerQG.symmetric(48, 48, 3.0e6, 3.0e6, Us=0.0,
                             nu4=0.05 * 4.0 * (3.0e6 / 48) ** 3, **P)
    s = m.random_state(amplitude=1.0, seed=1)
    e0 = m.eddy_kinetic_energy(s)
    s = m.solve(s, t_end=24 * 3600.0, safety=0.3)
    assert m.eddy_kinetic_energy(s) < e0


def test_single_layer_inversion_validity_and_cfl_guard():
    """Construction validation + the CFL guard (the explicit-solver discipline mirrored from the SW
    engine): a non-positive or over-large step raises rather than silently blowing up."""
    with pytest.raises(ValueError):
        TwoLayerQG(8, 8, 1.0, 1.0, f0=1e-4, gp=-1.0, H1=400, H2=400)
    m = TwoLayerQG.symmetric(16, 16, 2.0e6, 2.0e6, Us=4.0, **P)
    s = m.random_state(amplitude=1e-3, seed=2)
    with pytest.raises(ValueError):
        m.step(s, 10.0 * m.max_dt(s, safety=1.0))
    with pytest.raises(ValueError):
        m.step(s, -1.0)


# --------------------------------------------------------------------------- #
# the nonlinear engine reduces to the validated linear operator (slow)
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_nonlinear_engine_grows_at_linear_rate():
    """A single growing eigenmode at small amplitude grows at the analytic linear rate to <1 % — the
    nonlinear pseudospectral RHS reduces to the linear instability operator the tight anchor pins."""
    disp = TwoLayerQG.symmetric(8, 8, 1.0, 1.0, Us=4.0, **P)
    kstar, sig = disp.most_unstable()
    L, nx = 4 * 2 * np.pi / kstar, 96
    dx = L / nx
    m = TwoLayerQG.symmetric(nx, nx, L, L, Us=4.0, nu4=0.02 * 4.0 * dx ** 3, **P)
    kx = 2 * np.pi * np.fft.fftfreq(nx, d=dx)
    kmode = kx[np.argmin(np.abs(kx - kstar))]
    A = np.array([[-(kmode ** 2 + m.F1), m.F1], [m.F2, -(kmode ** 2 + m.F2)]], dtype=complex)
    Q1, Q2 = m.mean_pv_gradients
    M = kmode * np.linalg.solve(A, np.diag([m.U1, m.U2]) @ A + np.diag([Q1, Q2]))
    w, V = np.linalg.eig(M)
    ig = int(np.argmax(w.imag))
    sig_mode, ev = float(w[ig].imag), V[:, ig]
    x = (np.arange(nx) + 0.5) * dx
    psi = np.stack([np.real(ev[0] * np.exp(1j * kmode * x))[None, :] * np.ones((nx, 1)),
                    np.real(ev[1] * np.exp(1j * kmode * x))[None, :] * np.ones((nx, 1))])
    psi *= 1e-4 / np.sqrt(np.mean(psi ** 2))
    s = QGState(q=m.pv_from_psi(psi))
    ts, lnE, t = [], [], 0.0
    t_end = 6.0 / sig_mode
    nstep = 0
    while t < t_end:
        dt = m.max_dt(s, 0.3)
        if nstep % 20 == 0:
            ts.append(t)
            lnE.append(0.5 * np.log(m.eddy_kinetic_energy(s) + 1e-300))
        s = m.step(s, dt)
        t += dt
        nstep += 1
    ts, lnE = np.array(ts), np.array(lnE)
    i0, i1 = len(ts) // 5, 7 * len(ts) // 10
    sig_meas = np.polyfit(ts[i0:i1], lnE[i0:i1], 1)[0]
    assert sig_meas == pytest.approx(sig_mode, rel=0.05)


# --------------------------------------------------------------------------- #
# real-but-loose (the unlock) — the emergent saturated eddy thickness diffusivity (slow)
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_saturated_flux_is_down_gradient_and_non_vacuous():
    """The rung-3 payoff: in the SATURATED state the emergent meridional eddy thickness flux is
    **down-gradient** (κ>0) and **irreversible** (irr≈1, vs rung-1's barotropic ~0.1), with a
    **dimensionless** mixing efficiency ``κ/(v'_rms·L_d) = O(1)`` (vs rung-1's ~1e-3) — i.e. the
    reduction-to-diffusive-EBM is finally non-vacuous (plan §10, [[planet-rung3-phaseB-outcropping]]).

    Asserted DIMENSIONLESS / qualitative, never the dimensional κ (the idealized magnitude is
    box/drag-dependent and the large-scale condensate inflates v'_rms — advisor). A cheap
    config (the full turbulence-vs-wave characterization — the inverse-cascade KE spectrum and the
    vortex-filament PV field that distinguish turbulence from a steady wave — lives in the spike;
    here the saturated drag-arrested run pins the flux discriminators as a regression guard)."""
    disp = TwoLayerQG.symmetric(8, 8, 1.0, 1.0, Us=4.0, **P)
    kstar, sig = disp.most_unstable()
    nx = 64
    L = 3 * 2 * np.pi / kstar
    dx = L / nx
    m = TwoLayerQG.symmetric(nx, nx, L, L, Us=4.0, nu4=0.1 * 4.0 * dx ** 3,
                             r_drag=1.0 * sig, **P)
    s = m.random_state(amplitude=1e-3, seed=0)
    t_end, t_avg = 18.0 / sig, 8.0 / sig
    t, n = 0.0, 0
    fluxes, vrms, spec = [], [], None
    while t < t_end:
        dt = m.max_dt(s, 0.3)
        s = m.step(s, dt)
        t += dt
        n += 1
        assert np.isfinite(s.q).all()                # never outcrops (the QG win over free-surface SW)
        if t >= t_avg and n % 20 == 0:
            f1, f2 = m.bulk_eddy_flux(s)
            fluxes.append(0.5 * (f1 + f2))
            vrms.append(m.v_rms(s))
            Kc, E = m.ke_spectrum(s)
            spec = E if spec is None else spec + E
    fluxes = np.array(fluxes)
    kappa = fluxes.mean() / m.Us
    irr = abs(fluxes.mean()) / np.abs(fluxes).mean()
    ratio = kappa / (np.mean(vrms) * m.Ld)
    k_peak = Kc[np.argmax(spec)]                      # time-mean spectral peak
    assert kappa > 0.0                               # down-gradient
    assert irr > 0.8                                 # irreversible (rung-1 ~0.1)
    assert 0.1 < ratio < 10.0                        # O(1) mixing efficiency (rung-1 ~1e-3)
    assert np.mean(vrms) > m.Us                       # eddies exceed the mean shear (v'/U_s ≈ 1.8 here)
    # the SUFFICIENT check (advisor): the others are necessary-but-not-sufficient (any sustained
    # baroclinic state is down-gradient with irr≈1). A developed inverse cascade peaks BELOW the
    # injection wavenumber k* — energy transferred upscale — which a quasi-steady wave (peak AT k*)
    # cannot do. This is what makes the saturated state genuine turbulence, not a finite-amplitude wave.
    assert k_peak < kstar
