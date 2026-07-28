# pydme: background, data, and rates scaling

Summary of how pydme scales background, data, and signal rates in the likelihood and upper-limit analysis (from `collab_frameworks/pydme`).

---

## 1. Background (Bp, Br)

**File:** `pydme/detector/background_models/background_pattern.py`

- Template values are the same as in CCDarkSens (e.g. Bp = [141.4, 0.111, 0.042, ...], Br = [0.039, 0.039, ...]).
- **Scaling:** Bp and Br are **divided by `len(gamma)`** (number of bins) before forming the per-bin expectation:
  ```python
  Bp = df['Bp'].values / len(gamma)
  Br = df['Br'].values / len(gamma)
  B = Bp[row] + theta[0]*Br[row]   # per-bin expectation
  B_array = np.full_like(gamma, B)
  ```
- So **per-bin** expected count = (Bp_total + θ·Br_total) / n_bins; the **total** per pattern is still Bp_total + θ·Br_total. No extra global scaling of the background total.
- **Constraint:** The prior term is applied **per bin**: `L_array = -theta[0]*Br[row] + 98*log(theta[0]*Br[row])` where `Br[row]` is the **per-bin** Br (i.e. Br/len(gamma)). So the effective constraint strength depends on the number of bins (n_bins terms are summed in the NLL). In CCDarkSens we have one bin per pattern and a single constraint term per pattern; if pydme uses many gamma bins, the comparison is still consistent in total B = Bp + θ·Br, but the constraint shape can differ slightly.

---

## 2. Data

- **No additional scaling.** Data (e.g. `Count_Candidate_11`, etc.) are passed in and used as-is in the Poisson NLL.
- In the SRDM script, `D` is an array (one entry per time/gamma bin); `gamma`, `texp`, `Npix` have the same length. So the likelihood is per (pattern, bin).
- `self.exposure` is set from data (e.g. `exposure = np.sum(Npix * texp * Mpix)` or hardcoded 1250) and is used for **printing only** (“Exposure used”); it is **not** used inside the NLL. Exposure enters the **model** only via `t_exp`, `mass_pix`, `N_pix` in the signal rate → counts conversion.

---

## 3. Signal rates

**Files:** `pydme/dme_nLL_exclusion.py` (Model_pattern), rate tables / `get_ne_from_energy_rates`

- Signal rates from the input files are in **events/gram/day** (see docstring in `NLLAnalysis`).
- In `Model_pattern` the conversion to counts is:
  ```python
  s.append(S[idx] * t_exp * mass_pix * N_pix)
  ```
  So: (events/g/day) × (days) × (g) × (pixels) → counts. **No extra scaling** beyond this exposure factor.
- In the rate computation (`qedark4dm.py` / `qcdark4dm.py`): `rates = self.dRdE(...)/1000/365.25` → converts to **events/kg/day**. Standard unit conversion, no further scaling.

---

## 4. Summary table

| Item        | pydme scaling | Notes |
|------------|----------------|--------|
| **Bp, Br** | ÷ len(gamma) per bin | Total per pattern = Bp + θ·Br unchanged. |
| **Data**   | None           | Used as-is in Poisson. |
| **Signal** | × t_exp × mass_pix × N_pix | Exposure only; rates in events/g/day. |
| **Rates**  | /1000/365.25   | dRdE → events/kg/day. |
| **Exposure** | Not in NLL   | Only in model via t_exp, mass_pix, N_pix; also printed. |

For **CCDarkSens** using one bin per pattern and the same Bp, Br **totals**, the background model is consistent with pydme’s **total** (Bp + θ·Br) per pattern. The only structural difference is pydme’s per-bin splitting and per-bin constraint when there are many gamma bins.
