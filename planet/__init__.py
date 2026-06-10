"""Earth-system / Planet simulator — *planetary knobs in, climate & habitability out*.

Project #3 of the BigSim program and its **capstone** (plan: ``docs/plans/planet-earth-system.md``).
It reuses the shared diffusion/heat spine (:mod:`engines.diffusion`), unchanged, a **third** time — as a
sphere's latitudinal heat transport — and, in later phases, builds the program's one remaining
shared engine (``engines/fluid``, the shallow-water solver).

Phase 1 public API — the latitudinal EBM & the Snowball bifurcation::

    from planet.ebm import (
        EnergyBalanceModel, ClimateState, insolation, legendre_P2,
        equilibrium_temperature_0d, two_mode_solution, ice_line_latitude,
    )
    from planet.albedo import (
        EBMParams, planetary_albedo, absorbed_shortwave,
        present_day_climate, snowball_hysteresis, HysteresisLoop,
    )

> **UNIT SYSTEM — SI / climlab-conventional** (W m⁻², °C, ``x = sin φ`` dimensionless), unlike
> Chip's per-module native units: the EBM constants (climlab/North ``A, B, D, α, Tf``) are
> tabulated in W m⁻²/°C, and the engine is fed the transport in those units directly.
> Latitudes are reported in degrees. See :mod:`~planet.ebm` for the full banner.
"""
