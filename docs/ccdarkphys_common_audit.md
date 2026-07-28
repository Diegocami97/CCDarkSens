# ccdarkphys.common audit vs QEdark / pydme

Checked `python/ccdarkphys/common/` (constants, halo, io) against QEdark_f2.ipynb, QEdark_constants.py, and pydme mqedark.

---

## constants.py

| Constant | ccdarkphys | pydme QEdark_constants | QEdark notebook |
|----------|------------|------------------------|-----------------|
| **alpha** | 1/137.035999084 | 1/137.03599908 | 1/137 |
| **me_eV** | 510998.95 | 5.1099894e5 | 0.511e6 |
| **amu_kg** | 1.66053906660e-27 | 1.660538782e-27 | 1.660538782e-27 |
| **rho_X_eVcm3** | 0.3e9 | (in class, 0.3e9) | 0.4e9 (Cell 8) |
| **ccms** | 2.99792458e10 | c_light*1e2 (same) | — |
| **sec_per_year** | 365.25*86400 | sec2year (same) | — |

- Values match pydme; alpha/me_eV are slightly more precise than the notebook.
- **rho_X**: We use 0.3e9 (pydme default). Notebook uses 0.4e9; change here only if you want to match the notebook.

---

## halo.py

- **eta_shm_analytic**: Closed-form SHM (truncated Maxwellian + vE). Returns η in **(cm/s)^-1**. Same formula as standard DM-e literature; hard cutoff at vmin > vesc + vE; tiny negatives set to 0.
- **eta_shm_numeric**: Same physics as pydme `etaSHM` and `etaSHM_parallel`:
  - **KK** = v0³ [ -2 exp(-vesc²/v0²) π vesc/v0 + π^1.5 erf(vesc/v0) ] (same as pydme).
  - Regions (B4/B5 from 1509.01598): vmin ≤ vesc−vE, vesc−vE < vmin ≤ vesc+vE, vmin > vesc+vE.
  - Integrand (2π/KK)*vx*exp(-(vx²+vE²+2 vx vE cosq)/v0²), same as pydme.
- **kms_to_cms**: x * 1e5 (1 km/s = 10⁵ cm/s). Correct.

Docstring updated: module now states that η is in **(cm/s)^-1**, not dimensionless.

---

## io.py

- **data_path**: Resolves package-relative paths (e.g. for Si_f2.txt). No external dependency.
- **write_csv**: Writes E, R with header; **meta** must include material, mediator, table_path, table_sha1, v0_cm_s, vE_cm_s, vesc_cm_s, mchi_eV, sigma_e_cm2. Header states **events / kg / year / eV**. Matches what `entry.py` returns and what the C++ pipeline expects.

---

## Summary

- **constants**: Aligned with pydme; rho_X = 0.3e9 (notebook uses 0.4e9).
- **halo**: SHM η(vmin) matches pydme (KK, regions, integrand); units (cm/s)^-1 documented.
- **io**: CSV format and units correct for the pipeline.

No code changes required for consistency with pydme; only the halo module docstring was corrected (dimensionless → (cm/s)^-1).
