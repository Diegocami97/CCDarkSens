# pydme mqedark rate implementation vs CCDarkSens

Checked `collab_frameworks/pydme/pydme/mqedark/` to align with pydme for limit comparison.

---

## 1. ρ_X (DM density)

| Location | Value |
|----------|--------|
| **qedark4dm.py** `QEDark4DM.__init__` | `rho_X=0.3e9` (default), in eV/cm³ |
| **qcdark4dm.py** `QCDark4DM.__init__` | `rho_X=0.3e9` (default) |
| **qedark.py** | Uses `rho_X` in prefactor (commented at top: `#rho_X = 0.3e9`); low-level API gets it from caller/material |

So **pydme uses 0.3e9** by default — **same as our** `ccdarkphys.common.constants.rho_X_eVcm3 = 0.3e9`. No change needed for ρ_X to match pydme.

---

## 2. Prefactor and formula (QEDark4DM)

- **Prefactor** (qedark4dm.py ~L208):  
  `(ccms**2*sec2year) * (rhoX_over_mX) * (1/Mcell)*alpha*self.xsec_e*me_eV**2/mu_Xe**2`  
  with `rhoX_over_mX = 1 if self.SRDMFiles else self.rho_X/self.mX`.
- **Kernel**: `Eprefactor * (1/q) * eta * FDM(q,n)**2 * fcrys_atEe`, then `_dRdE = prefactor*np.sum(int_vect, axis=0)`.
- **Return**: Comment says "events/kg/year"; dimensionally this is **dR/dE in events/(kg·year·eV)** at each Ee (same as our interpretation).

So the **rate formula and units** match what we use in `entry.py` (prefactor × sum over q, output as differential rate per eV).

---

## 3. Default halo (velocity)

| Parameter | pydme (qedark4dm default) | Our config / entry |
|-----------|---------------------------|---------------------|
| v0        | **238e5** cm/s (238 km/s) | 220 km/s            |
| vE        | **263e5** cm/s (263 km/s) | 232 km/s            |
| vesc      | **544e5** cm/s (544 km/s)| 544 km/s            |

pydme’s **default** v0 and vE are **higher** than the usual SHM (220, 232). If the reference limit was produced with pydme defaults, our rates use a different halo (e.g. 220, 232) and can differ. To match pydme’s limit curve, use the **same** v0, vE, vesc as in the pydme run (e.g. 238, 263, 544 if they used defaults).

---

## 4. Kinematic cut

| Code      | Cut |
|-----------|-----|
| **qedark.py** (dRdE_nonvect) L123 | `if vmin > (vesc+vE):` — **no** 1.1 factor |
| **Our entry.py** (after fix)      | `if vmin > (vesc_cm_s + vE_cm_s) * 1.1` |

So **pydme** uses a **stricter** cut (smaller allowed region in q). We use a 1.1 safety margin like the QEdark notebook. That makes our rate **slightly higher** than pydme’s for the same (ρ_X, halo), not lower — so it does **not** explain a **weaker** limit on our side.

---

## 5. Export units (get_rates_SHM)

- `rates = self.dRdE(Ee, ...)/1000/365.25` → converts to **events/(g·day)** per energy (or per 0.1 eV bin, depending on Ee spacing).
- They store `"dRdE" : rates/self.xsec_e` (rate per unit σ_e).
- Upper limit script uses **exposure in g·day** (e.g. 1250 g·day) and multiplies.

We use **events/(kg·year·eV)** in CSVs and exposure in **kg·year**; the C++ pipeline applies exposure once. So our **convention** differs (kg·year vs g·day) but the **physics** (rate × exposure → counts) is the same.

---

## 6. Summary for “weaker limit” vs pydme

- **ρ_X**: pydme = 0.3e9, we = 0.3e9 → **aligned**.
- **Halo**: pydme defaults **238, 263, 544** km/s; we often use **220, 232, 544**. Using **238, 263, 544** in our rate generation (and possibly in the scan if pydme used that for the reference) can help match the reference.
- **Kinematic cut**: We use 1.1, pydme does not → we are slightly more permissive; not the cause of a weaker limit.
- **Other**: Exposure (livetime × mass), Bp/Br template, and CL/test statistic must also match the reference.

**Practical step:** Regenerate rate CSVs and run the limit with **v0=238, vE=263, vesc=544** (km/s) in the grid config halo to match pydme defaults, and ensure exposure and Bp/Br match the pydme analysis.
