"""Planet rung-3 Phase-B integration: the banked QG-turbulence evidence figure end to end (§10).

Runs the demo that banks the PV-fields + KE-spectrum figure and renders it. ``slow``-marked (it spins
the QG engine up to saturation); the *physics* — the inverse-cascade spectral peak ``k_peak < k*`` and
the down-gradient/irreversible flux discriminators — is asserted in ``test_baroclinic_qg.py``. This
render is an *execution* smoke-test (ADR 0002), ``importorskip``-gated on ``[viz]``, not a physics one.
"""
import pytest


@pytest.mark.slow
def test_demo_baroclinic_qg_renders(tmp_path):
    # ADR 0002: an execution smoke-test (test_baroclinic_qg validates the turbulence + flux numbers).
    pytest.importorskip("matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from planet import demo_baroclinic_qg as demo

    # a small, short run — just enough to exercise compute() + the figure path (not the full cascade)
    r = demo.compute(nx=32, n_lam=2, r_fac=1.0, n_efold_total=6.0, n_efold_avg=2.0)
    assert r.q.shape[0] == 2                               # two layers
    assert r.E_norm.max() == pytest.approx(1.0)            # spectrum normalised to its peak
    assert r.K_over_kstar.shape == r.E_norm.shape

    import planet.demo_baroclinic_qg as d
    d.DOCS_FIGURE = tmp_path / "qg-docs.png"
    d.OUTPUT_FIGURE = tmp_path / "qg-out.png"
    saved = demo.save_figure(r)
    assert saved.exists() and saved.stat().st_size > 0
    plt.close("all")
