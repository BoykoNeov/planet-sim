"""Frozen climlab / North reference table — the Phase-1 benchmark leg without the ``[climate]`` extra.

The validation triad's **benchmark** leg compares the EBM against climlab's `EBM` (the present-day
ice line, the Snowball threshold, the hysteresis). climlab is consumed as a **reference tool** (the
pycalphad pattern, ARCHITECTURE.md §7): it is *opt-in* behind the ``[climate]`` extra and **never
copied**, so to keep the committed triad green without it this module ships a **frozen table** of the
reference facts (cited [[ebm-radiation-source]]), and the *live* climlab cross-check
(:func:`climlab_present_day`) is a ``slow`` / ``importorskip`` test that runs only where climlab is
installed.

The numbers are **loose bands**, not point values: the exact ice line / threshold depend on the
calibrated radiation/albedo constants (the non-circularity split — the structure is asserted tight,
the calibrated thresholds only in bands), the way Steel's 1045 knee and Chip's contrast curve were.
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ClimateReference:
    """Frozen reference facts for the Phase-1 benchmark (climlab `EBM` defaults / North 1975).

    Cited [[ebm-radiation-source]]; used for comparison, never redistributed. All ranges are loose
    bands (calibration-dependent — the benchmark leg asserts structure tightly, thresholds loosely).
    """

    # Present-day climate (the finite-cap branch — Earth's branch). The ~70° edge is the OBSERVED
    # perennial snow/ice line (≈72°, North 1975 / Sellers), NOT this model's output (which gives ~73°,
    # so the reference is a genuine independent target, not the model relabeled — the non-circularity).
    present_ice_line_deg: float = 70.0                 # observed subpolar ice/snow edge (North 1975)
    present_ice_line_band: tuple[float, float] = (60.0, 80.0)
    present_global_mean_C: float = 14.5                # global-mean surface temperature ~14–15 °C
    present_global_mean_band: tuple[float, float] = (12.0, 18.0)

    # Snowball bifurcation (the dimming that triggers runaway glaciation = the large-ice-cap
    # instability). Budyko/North put it at a few-percent dimming; Voigt & Marotzke 2010 (a modern
    # coupled GCM) find 6–9 % — this model's ~8 % sits inside that range.
    snowball_dimming_pct_band: tuple[float, float] = (2.0, 12.0)
    snowball_global_mean_max_C: float = -20.0         # a Snowball planet is deeply frozen (T̄ well below 0)

    # The hysteresis is *wide*: the white planet re-melts only at a much brighter sun than it froze.
    hysteresis_positive: bool = True                  # melt_S0 > freeze_S0 (the loop has positive width)

    # The climlab EBM defaults these benchmarks come from (pinned; == planet.ebm constants).
    climlab_S0: float = 1365.2
    climlab_A: float = 210.0
    climlab_B: float = 2.0
    climlab_D: float = 0.555
    climlab_Tf: float = -10.0
    climlab_a0: float = 0.30
    climlab_a2: float = 0.078
    climlab_ai: float = 0.62
    climlab_s2: float = -0.48


REFERENCE = ClimateReference()


def climlab_present_day(num_lat: int = 180, years: float = 5.0):
    """Live climlab cross-check: the present-day global-mean T and ice-line latitude from climlab's `EBM`.

    Builds climlab's annual-mean ``EBM`` at its defaults and **seeds an Earth-like (capped) initial
    condition** — warm equator, frozen pole — so it settles on the *same* finite-cap branch this
    model's :func:`~planet.albedo.present_day_climate` targets (the system is bistable at
    present S₀, so the cross-check must compare like branch with like). Integrates to equilibrium and
    reads the global-mean surface temperature and the ice-line latitude (where ``Ts`` crosses ``Tf``).
    Used **only** by the ``slow`` / ``importorskip`` benchmark test (climlab is the opt-in
    ``[climate]`` extra) — returns ``(global_mean_C, ice_line_deg)``. Kept defensive: any climlab
    API/availability issue is the test's to skip on, not this module's to depend on.
    """
    import numpy as np
    import climlab                                     # the [climate] extra; the test importorskips

    model = climlab.EBM(num_lat=num_lat)
    try:                                               # seed the capped branch (best-effort across APIs)
        lat = np.asarray(model.lat).squeeze()
        model.Ts[:] = (30.0 - 60.0 * np.abs(np.sin(np.deg2rad(lat)))).reshape(model.Ts.shape)
    except Exception:                                  # fall back to climlab's own default IC
        pass
    model.integrate_years(years)
    Ts = np.asarray(model.Ts).squeeze()
    lat = np.asarray(model.lat).squeeze()
    nh = lat >= 0.0                                    # northern hemisphere (symmetric annual mean)
    lat_nh, Ts_nh = lat[nh], Ts[nh]
    global_mean = float(np.average(Ts, weights=np.cos(np.deg2rad(lat))))
    Tf = REFERENCE.climlab_Tf
    if np.all(Ts_nh > Tf):
        ice_line = 90.0
    elif np.all(Ts_nh <= Tf):
        ice_line = 0.0
    else:
        order = np.argsort(Ts_nh)                      # interpolate the Tf crossing
        ice_line = float(np.interp(Tf, Ts_nh[order], lat_nh[order]))
    return global_mean, ice_line
