# Model_pattern: pydme vs CCDarkSens

Comparison of the pattern-space likelihood model in pydme (`Model_pattern` + `Background_pattern` + NLL) and CCDarkSens (profile likelihood with pattern bins).

---

## 1. Model structure (same idea)

| Aspect | pydme | CCDarkSens |
|--------|--------|------------|
| **Observable** | Counts per pattern (and per γ bin) | Counts per pattern (6 bins) |
| **Expected rate** | R = S + B per (pattern, bin) | μ_i = S_pat[i] + B_pat[i] per pattern i |
| **Background** | B = Bp + θ·Br (same totals as here) | B_i = Bp[i] + θ·Br[i] |
| **Signal** | S from rate tables × t_exp × mass_pix × N_pix | S_pat from pipeline × exposure (kg·year) |
| **Likelihood** | Poisson: ∑ (μ − D ln μ) | Same: ∑_i (μ_i − n_i ln μ_i) |
| **Constraint** | −θ·Br + 98·ln(θ·Br) per term | Same form per pattern |

So at the level of “one expected count and one observed count per pattern,” the two implementations are the same: Poisson with μ = S + Bp + θ·Br and the same constraint on θ·Br.

---

## 2. Binning (main difference)

**pydme**

- One “dataset” per pattern (6 datasets: 11, 21, 111, 31, 22, 211).
- Each dataset has **several bins** (one per time/γ bin): `gamma` (or equivalent) has length `n_bins`.
- Per-bin expectation for pattern `patt`:
  - Background: `(Bp[patt] + θ·Br[patt]) / len(gamma)` (so total B = Bp + θ·Br).
  - Signal: `S[patt] * t_exp * mass_pix * N_pix`; if these are arrays, signal is distributed over bins; total S still = S_pat × exposure.
- NLL: one Poisson term per (pattern, γ bin):  
  `∑_pattern ∑_bin ( μ_bin − D_bin·ln(μ_bin) + constraint_bin )`.
- Constraint is applied **per γ bin**: each bin contributes  
  `−θ·Br_bin + 98·ln(θ·Br_bin)` with `Br_bin = Br_total / len(gamma)`.  
  So for one pattern the total constraint is  
  `−θ·Br_total + 98·len(gamma)·ln(θ·Br_total / len(gamma))`.

**CCDarkSens**

- **One bin per pattern** (6 bins total). No γ or time splitting.
- Data and model are vectors of length 6 (one count per pattern).
- μ_i = S_pat[i] + Bp[i] + θ·Br[i]; one Poisson term per pattern.
- Constraint: one term per pattern,  
  `−θ·Br[i] + 98·ln(θ·Br[i])`,  
  i.e. **no** division by `len(gamma)`.

So:

- **Same:** Total expected counts per pattern (S + Bp + θ·Br), same Poisson formula, same Bp/Br numbers and same constraint **form** (−θ·Br + 98·ln(θ·Br)).
- **Different:** pydme uses (pattern × γ) fine bins; CCDarkSens uses pattern-only bins. So pydme’s constraint is “per γ bin” and becomes  
  `98·len(gamma)·ln(θ·Br/len(gamma))` per pattern instead of `98·ln(θ·Br)`.

---

## 3. Where things live in code

**pydme**

- `dme_nLL_exclusion.Model_pattern`: builds R = s_ + b_ per pattern; s from `Signal(…, DoPatterns=True)` × exposure; b from `Background(…, pattern=patt, Npix=N_pix)`.
- `background_pattern.Background_pattern`: B = (Bp + θ·Br) / len(gamma) per bin; returns B_array, L_array (constraint per bin).
- `nLL_Poisson`: loops over datasets (patterns) and bins; `Nim[f]` = model counts per bin; `self.data[i]` = data per bin.

**CCDarkSens**

- Scan app: builds `B_pat` (from config Bp/Br when `background_source == "bp_br_template"`), `S_pat` via `FoldNeToPatternRates(S_true, …, pattern_eff_map)`.
- `ProfileLikelihood::NLL(S, param)`: μ_i = S[i] + Bp[i] + param*Br[i]; one term per pattern; constraint uses `constrain_prior_strength_` (98) and param*Br[i].

---

## 4. Summary

- **Same:** Poisson likelihood, B = Bp + θ·Br totals per pattern, same Bp/Br values, same type of constraint (−θ·Br + 98·ln(θ·Br)), signal × exposure → counts.
- **Different:** pydme has (pattern × γ) bins and divides Bp/Br by `len(gamma)` per bin; constraint is summed over γ so it becomes `98·len(gamma)·ln(θ·Br/len(gamma))` per pattern. CCDarkSens has one bin per pattern and uses 98·ln(θ·Br) per pattern with no γ splitting.

So **Model_pattern in pydme is very similar to CCDarkSens**: same physical model and same likelihood type. The only structural difference is **binning** (pattern×γ vs pattern-only) and the resulting **constraint normalization** when `len(gamma) > 1`. For Asimov or single-bin data, results should be close; small differences can come from the constraint when pydme uses many γ bins.

---

## 5. What actually sets len(gamma): data + do_single_bin_likelihood

**len(gamma) is not fixed by “analysis type”.** It comes from the **data**:

- **SRDM** (`analysis/SRDM/upper_limit/lbc_dmanalysis_upperlimits.py`):  
  `gamma = conversion_helpers.get_isoreflection_angles_vectorized(utc_dates) / 180`  
  with `utc_dates` from the CSV (e.g. `Final_Combined_Image_Data.csv`). So **len(gamma) = number of rows** in that file (one per image/time bin).

- **Daily modulation** (`analysis/DailyModulation/.../lbc_dmanalysis_upperlimits.py`):  
  `gamma = data[col_naming['gamma']].values / np.pi` (or from `get_isoreflection_angles_vectorized`). Again **len(gamma) = number of rows** in the data.

So **len(gamma) = 1** only if the input table has a single row (one aggregate bin). Otherwise it’s the number of time/γ bins in the dataset.

The flag that changes how the likelihood uses those bins is **`do_single_bin_likelihood`** in `dme_nLL_exclusion.NLLAnalysis` (default `False`):

| Script              | `do_single_bin_likelihood` | Effect |
|---------------------|----------------------------|--------|
| **SRDM** upper limit| **True** (set at line 166)  | Poisson uses **totals** per pattern: `_model = np.sum(Nim[f])`, `_data = np.sum(self.data[i])`. One Poisson term per pattern, like CCDarkSens. |
| **Daily modulation**| **False** (line 224)        | Poisson uses **per-bin** model and data; many terms per pattern. |

Important: when `do_single_bin_likelihood = True`, only the **Poisson** part is collapsed to one bin per pattern. The **background and constraint** are still built from the full **gamma** array:

- `Model_loader` → `Background_pattern(theta, pattern, gamma, …)` is called with **full** `self.gamma[i]`.
- So **Bp/Br are still divided by len(gamma)** and the **constraint** is still **per bin**, then summed: `_constrain = np.sum(constrain[f])` ⇒ total constraint = **−θ·Br + 98·len(gamma)·ln(θ·Br/len(gamma))** per pattern.

So for **SRDM** (single-bin likelihood):

- **Poisson:** matches CCDarkSens (one term per pattern, totals).
- **Constraint:** still depends on **len(gamma)** from the data. Only if the SRDM CSV has **one row** do you get len(gamma)=1 and the same constraint as CCDarkSens. If the SRDM file has many rows, the constraint in pydme is still “multi-bin” and can disagree with CCDarkSens.

**Why constrain_n_bins=4450 often does not move the limit:** With n_bins=4450 the constraint term is huge and dominates the Poisson term, so the minimizer pushes theta to the lower bound (0.5). Then B(theta)=Bp+0.5*Br; with Bp >> Br we have B(0.5) approx B(1), so the limit curve barely changes. Use constrain_n_bins=1 for the usual single-bin constraint.
