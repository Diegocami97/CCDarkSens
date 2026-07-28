# pydme vs CCDarkSens: Comparison and How to Resolve Limit Discrepancy

This document compares the two frameworks step-by-step and lists concrete checks to align the limit curves.

---

## 1. Exposure

| | pydme | CCDarkSens |
|---|--------|-------------|
| **Units** | g·day | kg·year (internally); config targets kg·day |
| **Definition** | `exposure = np.sum(Npix * texp * Mpix)` or override `1250` | `exposure_kg_year = livetime_days * duty_cycle * mass_kg / 365.25` |
| **Where** | `lbc_dmanalysis_upperlimits.py:64,67`; then `dmeLL.exposure = exposure` | `ExperimentSetup.cc:33-34`; `summary.exposure_kg_year` |
| **Signal use** | `S[idx]*t_exp*mass_pix*N_pix` per row (events/g/day × g·day) | `FoldToNe(dRdE, exposure_kg_year, ...)` then pattern fold |

**Alignment:** Use 1.3 kg·day = 1300 g·day. CCDarkSens: set `livetime_days` so `livetime_days * mass_kg = 1.3` (e.g. 85.356 × 0.01523). pydme: set `exposure = 1300` (g·day) if using the override.

---

## 2. Background (B_pat)

| | pydme | CCDarkSens |
|---|--------|-------------|
| **Source** | `Background_pattern(theta, pattern, N_pix)` | `background_Bp`, `background_Br` when `background_source: "bp_br_template"` |
| **Formula** | Bp, Br from dict; **Bp/len(gamma), Br/len(gamma)** then `B = Bp[row]+θ*Br[row]` per pattern. (Third arg is actually N_pix; len = number of rows.) | `B_pat[i] = Bp[i] + Br[i]` (θ=1 nominal); θ profiled in NLL. |
| **Total** | Sum over rows gives Bp+θ·Br per pattern → same total as CCDarkSens when one row per pattern. | B_pat = Bp + θ·Br (same numbers) ✓ |
| **Flat / DC** | Not added in `Background_pattern`; template is pre-derived. | With `bp_br_template`, flat and DC are **not** applied; only Bp, Br from config. |
| **Migration / background efficiency** | `Background_pattern_()` can use `Background_efficiencies.csv`; SRDM pattern analysis uses `Background_pattern()` (no CSV). | `background_efficiency_csv` only used when `background_source != "bp_br_template"`. With template: **not used**. |

**Alignment:** Bp, Br and pattern order already match. No extra background efficiency or flat is applied in either when using the template.

---

## 3. Signal (S_pat)

| | pydme | CCDarkSens |
|---|--------|-------------|
| **Source** | Pre-computed rates in CSVs (events/g/day per pattern) | QEDark dR/dE tables → ionization → pattern efficiency |
| **Path** | `event_interpol_mX` (from `load_data_file`) → `Signal(..., DoPatterns=True)` → `S[idx]*t_exp*mass_pix*N_pix` | `DMElectronModel::MakeSpectrum_E()` → `ChargeIonization::FoldToNe(dRdE, exposure_kg_year)` → `FoldNeToPatternRates(..., pattern_eff_map)` |
| **Halo** | **SRDMModulation** in LBC analysis: signal depends on (γ, log10(σ_e)); interpolated from files. | **SHM** (standard halo): dR/dE from QEDark, no γ dependence. |
| **Efficiency** | Baked into the rate files (e.g. from `Efficiencies_patterns_Nsims100000_...` when generating rates). | Applied at run time via `efficiency_csv` / `pattern_eff_map`. |

**Critical difference:**  
- **Halo:** pydme SRDM uses **SRDMModulation** (velocity distribution and possibly modulation); CCDarkSens uses **SHM** with QEDark. Different halo → different S(m_χ, σ_e) → different limits.  
- **Signal pipeline:** pydme = pre-folded rates (events/g/day); CCDarkSens = dR/dE + ionization + pattern ε. To match, the same mediator, material, ionization, and efficiency must be used when generating pydme’s rate files as in CCDarkSens.

---

## 4. Efficiency

| | pydme | CCDarkSens |
|---|--------|-------------|
| **Signal** | Encoded in the pattern rate CSVs (e.g. from 100k or 1M efficiency table when building rates). | `efficiency_csv` → `pattern_eff_map`; optional `efficiency_csv_reference` overlay. |
| **Background** | No separate table; B = Bp + θ·Br. | With `bp_br_template`: no background efficiency applied. |

**Alignment:** Use the **same** efficiency file as the reference (e.g. `Efficiencies_patterns_Nsims1000000_DCTrue_alpha1.csv`). If the reference was produced with 100k efficiencies, use that; if 1M, use 1M. Mismatch (e.g. 22@4, 31@5) shifts limits.

---

## 5. Likelihood and upper limit

| | pydme | CCDarkSens |
|---|--------|-------------|
| **Per-bin** | Poisson: `μ - D*ln(μ)` | Same |
| **Bins** | 6 pattern bins; `do_single_bin_likelihood=True` → sum data and model per pattern, then 6 terms. | `pattern_roi` (6 patterns); one Poisson term per pattern. ✓ |
| **Constraint** | `-θ·Br + 98·ln(θ·Br)` per pattern | Same (`constrain_prior_strength: 98`) ✓ |
| **θ bounds** | `[[0.5, 10]]` | `theta_lo: 0.5`, `theta_hi: 10` ✓ |
| **q_μ, q0, CL** | `2*(NLL(σ)-NLL_min)`; target_q = (norm.ppf(CL))² | Same |
| **UL algorithm** | Bracketing + bisection in log10(σ_e) | Same when `profile_minimizer: "pydme"` or Minuit 2D |

**Alignment:** Likelihood and UL procedure are already aligned when using the same CL and θ range.

---

## 6. Checklist to resolve the discrepancy

Use this order.

1. **Exposure**
   - CCDarkSens: `livetime_days * mass_kg = 1.3` (e.g. 85.356 × 0.01523).
   - If comparing to pydme with 1250 g·day, that’s 1.25 kg·day; use 1.3 kg·day for the 2025 reference if that’s what they quote.

2. **Efficiency**
   - Set `efficiency_csv` (and `efficiency_csv_reference` if used) to the **exact** file used for the reference curve (1M or 100k, same DC/alpha).
   - Compare a few (pattern, n_e) bins between that file and pydme’s rate-generation efficiency if possible.

3. **Halo and signal source**
   - If the reference is **SHM**: CCDarkSens (SHM + QEDark) is the right setup; ensure rate tables (mediator, material, path) match the reference.
   - If the reference is **SRDM**: pydme’s SRDM pipeline (SRDMModulation + pre-folded rates) is the right setup; CCDarkSens with SHM will **not** match; you’d need SRDM signal in CCDarkSens or compare to pydme run with SHM.

4. **Background**
   - With `bp_br_template`, Bp/Br and θ treatment already match; no flat or background efficiency is applied in either. Ensure Bp/Br were derived for the **same** exposure (e.g. 1.3 kg·day) and same patterns.

5. **Rate tables**
   - Same material (Si), mediator (e.g. heavy), and dR/dE set (e.g. after_eta_fix). pydme signal files are pre-folded; CCDarkSens folds at runtime — only the combined (dR/dE × ionization × efficiency) should match.

6. **θ bounds**
   - Keep `[0.5, 10]` and ensure the best-fit θ is not hitting the boundary (e.g. increase `theta_hi` if needed for patterns with larger Br).

---

## 7. Quick reference: where things live

| Component | pydme | CCDarkSens |
|-----------|--------|------------|
| Exposure | `lbc_dmanalysis_upperlimits.py:64,67,157` | `ExperimentSetup.cc:33-34` |
| B_pat | `pydme/detector/background_models/background_pattern.py` | `apps/ccdarksens_scan_dmelectron_pattern.cc:856-865` |
| S_pat | `dme_nLL_exclusion.py` (`Signal`, `Model_pattern`); `data/__init__.py` (load) | `FoldToNe` + `FoldNeToPatternRates` in scan app |
| NLL / UL | `dme_nLL_exclusion.py` (`nLL_Poisson`, `compute_upper_limit_fast`) | `ProfileLikelihood.cc`, scan app pydme/Minuit branch |

---

**Bottom line:** The main knobs that change the limit curve are **exposure**, **efficiency table**, and **halo/signal source** (SHM vs SRDM). Match those to the reference first; then Bp/Br and θ treatment are already consistent between pydme and CCDarkSens when using the template.
