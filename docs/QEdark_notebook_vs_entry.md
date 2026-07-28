# QEdark notebook vs ccdarkphys.qedark.entry

This document compares the logic and units in **QEdark_f2.ipynb** (QEdark-python) with **python/ccdarkphys/qedark/entry.py** so that implementation and units match.

**Reference notebook:** `QEdark_f2.ipynb` (e.g. `/Users/.../QEdark-python/QEdark_f2.ipynb`).

---

## 1. dRdE return value and units

- **Notebook** (Cell 14): `dRdE(material, mX, Ee, FDMn, halo, params)` docstring says *"returns dR/dE [events/kg/year]"* and the code comment says `# [(kg-year)^-1]`. In **dRdnearray** they do `array_[ne] = np.sum(tmpdRdE)` over sub-bins **without** multiplying by dE. For that sum to be the integrated rate in the bin in **events/(kg·year)**, each `dRdE(...)` must return the **integrated rate in a bin of width dE** (0.1 eV), i.e. **events/(kg·year)**. So the notebook returns at each Ee the **rate in that 0.1 eV bin**, not the differential per eV. True dR/dE = (return value) / dE → **events/(kg·year·eV)**.

- **entry.py**: We compute the same prefactor×sum as the notebook, then **divide by dE** to get the differential: `dRdE_kg_year_eV = dRdE_kg_year / dE`. So the CSV is in **events/(kg·year·eV)** as required by the pipeline.

- **Pipeline**: CSV and C++ `RateTable` expect **events/(kg·year·eV)**. The division by dE is correct.

**How to verify:** Run the notebook’s `dRdnearray("Si", mX=0.5e6, Ebin=3.8, nFDM, "shm", vparams)` with sigma_e=1 and your halo; note the first non-zero bin value (integrated rate in [1.2, 5.0] eV in events/(kg·year)). Then run `python utils/verify_dRdE_units.py` and compare “Integrated rate (sum of dRdE_kg_year_eV * dE)” for the same bin. They should match (up to halo/constants). If they match when we use `dRdE_kg_year_eV = dRdE_kg_year / dE`, the convention is correct.

---

## 2. Formula alignment (entry.py vs notebook)

| Item | Notebook (Cell 14) | entry.py |
|------|--------------------|----------|
| **dQ** | `0.02*alpha*me_eV` | `0.02 * QEC.alpha * QEC.me_eV` |
| **dE** | `0.1` | `binsize_eV` (default 0.1) |
| **nq, nE** | 900, 500 | 900, 500 |
| **fcrys load** | `transpose(resize(loadtxt(..., skiprows=1), (nE,nq)))` | `_load_si_table_as_notebook(nE, nq)` same logic |
| **materials** | `[2*28.0855*amu2kg, 2.0, 1.2, 3.8, wk/4*fcrys['Si']]` | Same (Eprefactor=2, Egap=band_gap_eV, epsilon=eh_pair_eV, wk/4*fcrys) |
| **FDM(q,n)** | `(alpha*me_eV/q)**n`, n=0,1,2 | Same (n=0 → 1; n=2 → (α m_e/q)²) |
| **mu_Xe(mX)** | `mX*me_eV/(mX+me_eV)` | Same |
| **Prefactor** | `ccms**2*sec2year*rho_X/mX*1/Mcell*alpha*me_eV**2/mu_Xe(mX)**2` | Same with QEC (ccms, sec_per_year, rho_X_eVcm3, amu_kg) |
| **Ei** | `int(floor(Ee*10))` | `int(floor(Ee*10))` |
| **vmin** | `(q/(2*mX)+Ee/q)*ccms` | Same (with qsafe for q→0) |
| **Kinematic cut** | `if vmin > (vesc+vE)*1.1` | `if vmin > (vesc_cm_s + vE_cm_s) * 1.1` |
| **η(vmin)** | `etaSHM(vmin, params)` from DM_halo_dist | `HALO.eta_shm_numeric(...)` (same (cm/s)^-1 units) |
| **Summand** | `Eprefactor*(1/q)*eta*FDM(q,n)**2*materials[...][qi-1,Ei-1]` | Same |
| **Return** | `prefactor*np.sum(array_)` → events/(kg·year·eV) | Same; then × sigma_e, no ÷ dE |

---

## 3. Constants

| Constant | QEdark_constants.py / notebook | ccdarkphys.common.constants |
|----------|-------------------------------|-----------------------------|
| **sec2year / sec_per_year** | 60*60*24*365.25 | 365.25*86400 (same) |
| **ccms** | c_light*1e2 | 2.99792458e10 (same) |
| **alpha** | 1/137 | 1/137.035999084 |
| **me_eV** | 0.511e6 | 510998.95 |
| **amu2kg** | 1.660538782e-27 | 1.66053906660e-27 |
| **rho_X** | Notebook Cell 8: **0.4e9** eV/cm³ | **0.3e9** eV/cm³ |

So **rho_X** differs: notebook uses **0.4e9**, we use **0.3e9**. Rates scale linearly with rho_X, so our rates are **0.3/0.4 = 0.75** of the notebook if all else is equal. To match the notebook exactly, set `rho_X_eVcm3 = 0.4e9` in `python/ccdarkphys/common/constants.py` (or make it a parameter).

---

## 4. Halo η(vmin)

- **Notebook**: `etaSHM(vmin, params)` in `DM_halo_dist.py`; params = [v0, vE, vesc] in **cm/s**; returns a scalar with dimension **(cm/s)^-1** (integral of f/v over the allowed region).
- **entry.py**: `HALO.eta_shm_numeric(vmin, v0_cm_s, vE_cm_s, vesc_cm_s)` with velocities in **cm/s**; return is **(cm/s)^-1**. The analytic form in `eta_shm_analytic` follows the usual truncated Maxwellian; the numeric version is used for consistency and to avoid edge-case cancellations.

---

## 5. Summary

- **Implementation**: entry.py matches the notebook’s dRdE formula, indexing, kinematic cut (×1.1), and interpretation of the return value as **events/(kg·year·eV)**. No division by dE is applied when building the output array.
- **Units**: Pipeline output is **events/(kg·year·eV)** end-to-end.
- **Difference**: **rho_X** is 0.3e9 here vs 0.4e9 in the notebook; change `rho_X_eVcm3` to 0.4e9 if you want exact numerical agreement with the notebook.
