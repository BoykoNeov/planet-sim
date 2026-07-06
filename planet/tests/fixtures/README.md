# Test fixtures

## `oscar_subsample.npz` (~14 KB)

A 5° (stride-20) subsample of one real OSCAR surface-current granule, used by
`test_ocean_currents.py` so the O2 producer + the R1 round-trip run against **real** data without
committing the raw 33 MB netCDF (plan §9.6 data discipline).

* Source granule: `oscar_currents_final_20200601.nc` — OSCAR L4 v2.0, daily mean 2020-06-01,
  0.25°, NASA PO.DAAC, DOI `10.5067/OSCAR-25F20`. The product is freely distributable
  ("OSCAR products are supported by NASA and may be freely distributed" — granule
  `acknowledgment` attribute); this subsample keeps that attribution.
* Conventions are kept **raw** on purpose (0–360° longitude, NaN land, cell-centred ±89.75° lat)
  so the tests exercise the producer's rewrap / mask / pole handling, not a pre-cleaned copy.
* Regenerate (needs the granule + the `[ocean]` extra):

  ```python
  import numpy as np
  from planet.ocean_currents import load_oscar
  s = load_oscar("oscar_currents_final_20200601.nc", stride=20)
  np.savez_compressed("planet/tests/fixtures/oscar_subsample.npz",
                      lat=s.lat, lon=s.lon, u=s.u.astype(np.float32), v=s.v.astype(np.float32),
                      product=np.array(s.product), doi=np.array(s.doi), credit=np.array(s.credit),
                      date=np.array(s.date), depth_note=np.array(s.depth_note))
  ```
