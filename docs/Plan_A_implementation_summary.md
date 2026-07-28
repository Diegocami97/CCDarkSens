# Plan A implementation summary

Summary of code changes made to add **pattern-space observable** (Part A) to CCDarkSens, so the likelihood can be computed in pattern-ID bins (e.g. 11, 21, 111, 31, 22, 211) instead of n_e bins.

---

## 1. Experiment config (pattern_roi, observable_bins)

**Files:** `include/ccdarksens/experiment/ExperimentSetup.hh`, `src/experiment/ExperimentSetup.cc`

- **ExperimentConfig**
  - `pattern_roi`: `std::vector<int>` — list of pattern IDs used for the likelihood when observable is pattern (e.g. `[11, 21, 111, 31, 22, 211]`).
  - `observable_bins`: `std::string` — `"n_e"` (default) or `"pattern"`.
- **ExperimentSummary**
  - Same fields copied in `prepare_summary()` for use in apps.

---

## 2. Config parsing

**File:** `src/io/ConfigManager.cc`

- In `parse_experiment_()`:
  - If key `"pattern_roi"` present → `exp_cfg_.pattern_roi = j.at("pattern_roi").get<std::vector<int>>()`.
  - If key `"observable_bins"` present → `exp_cfg_.observable_bins = j.at("observable_bins").get<std::string>()`.

---

## 3. FoldNeToPatternRates (pattern rates from S(n_e))

**New files:** `include/ccdarksens/response/PatternRates.hh`, `src/response/PatternRates.cc`  
**Build:** `src/response/PatternRates.cc` added to `ccdarksens_core` in `CMakeLists.txt`.

- **Function:** `FoldNeToPatternRates(h_ne, ne_min, ne_max, pattern_roi, pattern_eff_map)`
  - Inputs: TH1D of rates in n_e, n_e range, list of pattern IDs, map `(pattern_id, ne) -> efficiency` (from CSV).
  - For each pattern `p` in `pattern_roi`:  
    `rate_p = sum_{ne = ne_min}^{ne_max} h_ne(ne) * epsilon(p, ne)`  
  - Missing `(p, ne)` in the map are treated as 0.
  - Returns: `std::vector<double>` of rates per pattern, same order as `pattern_roi`.

---

## 4. Scan app branch (observable_bins)

**File:** `apps/ccdarksens_scan_dmelectron_pattern.cc`

- **Include:** `#include "ccdarksens/response/PatternRates.hh"`.
- **After background:** If `summary.observable_bins == "pattern"`:
  - Require non-empty `pattern_roi`. If `pattern_eff_map` is empty, fill it from `emc_ptr->GetPatternTable()` (see §8); then require non-empty map.
  - Precompute `B_pat = FoldNeToPatternRates(*B_tot, ne_min, ne_max, pattern_roi, pattern_eff_map)` once.
- **Per grid point:**
  - If **pattern**: `S_pat = FoldNeToPatternRates(*S_obs, ...)`, then `data = B_pat`, `model_null = B_pat`, `model_test[i] = S_pat[i] + B_pat[i]`; run PLR on these vectors.
  - If **n_e**: unchanged; build `data` / `model_null` / `model_test` from `roi_bins` and S_obs/B_tot in n_e; run PLR.

---

## 5. Config JSON

**File:** `configs/ccdarksens_scan_dmelectron_pattern.json`

- Under `experiment`: added optional `"pattern_roi": [11, 21, 111, 31, 22, 211]` and `"observable_bins": "n_e"`.
- To run in pattern space: set `"observable_bins": "pattern"` and keep `response.pattern_mc.efficiency_csv` and `experiment.pattern_roi`.

---

## 6. Notebook-style pattern efficiency logic (Plan A)

To **calculate pattern efficiencies in C++ using the same logic as the notebook** (e.g. `pattern_efficiency.py` / `efficiencies.py`), the following was added.

**PatternClassifier**
- Pattern statistics use the same formulas as the notebook by default:
  - **Pm(q,m)** = -log(norm.cdf(q, m, σ))
  - **Pmn(q0,q1,m,n)** = min over permutations (m,n)/(n,m) of -log(cdf(q0,m)·cdf(q1,n))
  - **Pmnl(q0,q1,q2,m,n,l)** = -log(cdf(q0,m)·cdf(q1,n)·cdf(q2,l))
- Same threshold comparison: pass when stat < thr_M / thr_MN / thr_MNL. CDF is implemented with `std::erf`; probabilities are clamped to avoid log(0).

Image generation remains the row-segment path (ChargeTransport + PixelSimulator).

**Alignment with `collab_frameworks/efficiencies.py` and `How_to_calculate_efficiencies_for_patterns.ipynb`:**
- The same **Pm, Pmn, Pmnl** formulas (normal CDF, −log, min over permutations for Pmn) and **thresholds** (thr_M, thr_MN, thr_MNL, configurable to match e.g. 3.5, 4, 5.5) are used in `PatternClassifier`, so the pattern-identification logic matches the notebook.
- Diffusion (σ_xy from A, b, z, E), readout noise, and optional dark current are applied in PatternMC’s pipeline (ChargeTransport + PixelSimulator) before classification; the notebook uses the same physics in `generate_image_E`.
- **Making the same plots:** The notebook plots efficiency vs n_e (and vs E), charge distributions per n_e, and sanity-check scatter (e.g. q1 vs q2 for pattern 11). To reproduce them from CCDarkSens: (1) Run the scan or example app to get the pattern table (or use an efficiency CSV produced by the C++ side); (2) Write efficiencies (pattern, ne, efficiency) to a CSV in the same format as `efficiencies.generate_file` if desired; (3) Use the notebook’s `Plot_efficiencies_from_file`, `Ploting_simulation_specific_E`, or `sanity_checks` on that data, or a small Python script that reads the ROOT output (e.g. `S_pat_validation`, `B_pat_validation`) or the efficiency CSV. No extra C++ plot code is required; the existing notebook and `efficiencies.py` can consume the exported tables.

---

## 7. ROOT const-correctness workaround

**Files:** `src/response/DetectorResponsePipeline.cc`, `src/response/PatternRates.cc`

- ROOT’s `TH1D::FindBin()` is not marked const. Calls on `const TH1D&` use `const_cast<TH1D&>(h).FindBin(...)` so the code compiles.

---

## Flow when observable_bins == "pattern"

1. Load efficiency CSV (or fill from PatternMC) → `pattern_eff_map` (pattern, ne) → P(pattern|n_e).
2. **Option A (no double application):** Pipeline and BackgroundBuilder do **not** apply total ε(n_e). Pipeline returns S_rec(n_e), background is B_raw(n_e).
3. B_pat = FoldNeToPatternRates(B_raw, …); once.
4. For each (m_χ, σ_e): S_rec = pipeline.Apply(dR/dE, …); S_pat = FoldNeToPatternRates(S_rec, …); PLR(data=B_pat, model_test=S_pat+B_pat, model_null=B_pat) → q; fill TH2.
5. Write TH2 and (with validation) pattern-rate histograms for the first grid point.

---

## Validation output (single-case examples)

When running in pattern mode, for the **first grid point** only the app now:

- Prints a short validation block: S_obs(n_e) for ne_min..ne_max, then for each pattern bin: pattern ID, S_pat, B_pat, model_test.
- Writes to the output ROOT file:
  - `S_pat_validation`: TH1D of signal rate per pattern bin (bin index = pattern index, x-axis labels = pattern IDs).
  - `B_pat_validation`: TH1D of background rate per pattern bin (same binning).

Use these to check that folding and pattern IDs match your expectations (e.g. pydme or a hand calculation).

### Viewing the histograms

Run the single-point example app with your config JSON (e.g. `configs/ccdarksens_scan_dmelectron_pattern.json`) so that `experiment.observable_bins` is `"pattern"`:

```bash
ccdarksens_example_one_point_pattern configs/ccdarksens_scan_dmelectron_pattern.json
```

The app writes into `run.outdir` (e.g. `outputs/scan_pattern/`):

- **ROOT file** `example_one_point_pattern.root`: contains `S_pat_validation`, `B_pat_validation` (rate per pattern), plus `S_obs_ne`, `B_tot_ne`.
- **PNG plot** `example_one_point_pattern_rates_per_pattern.png`: bar chart of signal and background rate per pattern (same data as the ROOT histograms), so you can inspect the rate per pattern without opening ROOT.

To inspect or replot from the ROOT file in ROOT:

```bash
root -l outputs/scan_pattern/example_one_point_pattern.root
```

Then e.g. `S_pat_validation->Draw("BAR")`, `B_pat_validation->Draw("BAR SAME")`, or `S_obs_ne->Draw()`, `B_tot_ne->Draw("SAME")`.

### Flow of the example one-point app (with new features)

**Usage:** `ccdarksens_example_one_point_pattern config.json [mchi_MeV] [sigma_e_cm2]`

**Config (from JSON):**

- **run:** outdir, verbosity, rng_seed  
- **experiment:** binning (ne_min, ne_max), pattern_roi, observable_bins  
- **response.pattern_mc:** n_events_per_ne, row_length (or 2×half_window_pix+1), sigma_readout_e, diffusion (A, b, α, β), efficiency_csv (optional), etc.  
- **response.pattern_image:** nrows_binned, ncols, row_binning (for 2D image / pattern efficiency; parsed and available as `cfg.response().pattern_image`).  
- **response.pattern_classifier:** thr_M, thr_MN, thr_MNL, Qmin_e, Qmax_e, etc.

**Flow (pattern mode, `observable_bins == "pattern"`):**

1. **Parse config** → ConfigManager + raw JSON (for keys not in ConfigManager, e.g. debug, efficiency_csv path).
2. **Setup** → ExperimentSetup (summary), detector, ChargeTransport, PatternMC (row_length from JSON or half_window), PatternClassifier.
3. **Efficiencies:**
   - If **efficiency_csv** in JSON: load (pattern_id, n_e, efficiency) → `pattern_eff_map`; set `efficiency_source_desc = "CSV: <path>"`.
   - Else: leave `pattern_eff_map` empty for now.
4. **Pipeline** → DetectorResponsePipeline (Pattern, skip pattern efficiency when pattern mode), PatternEfficiency from h_eps_ne (PrecomputeEpsilon or PrecomputeEpsilonWithPatternEff).
5. **Precompute ε(n_e)** → PatternMC builds table if needed; h_eps_ne set on pipeline (used only when not in pattern mode).
6. **Background** → BackgroundBuilder (no pattern efficiency in pattern mode) → B_dc_ne; flat dR/dE → pipe.Apply → B_flat_ne; B_tot = B_dc_ne + B_flat_ne.
7. **Pattern efficiencies (pattern mode):**
   - If `pattern_eff_map` still empty → fill from `PatternMC::GetPatternTable()`; set `efficiency_source_desc = "PatternMC (BuildPatternTable, N trials per n_e)"`.
   - Print: `Efficiencies computed from: <efficiency_source_desc>`.
   - Write **efficiency_per_pattern.csv** (pattern_id, n_e, efficiency) in run.outdir.
   - If **verbosity ≥ 2:** print table P(pattern|n_e) (rows = n_e, cols = pattern_roi).
   - **B_pat** = FoldNeToPatternRates(B_tot, pattern_roi, pattern_eff_map).
8. **Signal** → dRdE_sig → pipe.Apply → **S_obs** (in pattern mode this is S_rec(n_e), no ε applied in pipeline).
9. **Pattern space** → S_pat = FoldNeToPatternRates(S_obs, …); model_test = S_pat + B_pat.
10. **Print** S_obs(n_e), then per-pattern S_pat, B_pat, model_test.
11. **PLR** → q_ts = EvaluateRatio(B_pat, model_test, B_pat) (Asimov).
12. **ROOT output** → example_one_point_pattern.root: S_obs_ne, B_tot_ne, S_pat_validation, B_pat_validation (bin labels = pattern_roi).
13. **PNG** → example_one_point_pattern_rates_per_pattern.png (bar chart S_pat + B_pat).

**Outputs in run.outdir:**

| Output | Description |
|--------|-------------|
| example_one_point_pattern.root | S_obs_ne, B_tot_ne, S_pat_validation, B_pat_validation |
| example_one_point_pattern_rates_per_pattern.png | Rate per pattern bar chart |
| efficiency_per_pattern.csv | pattern_id, n_e, efficiency (P(pattern\|n_e)) |

**When `observable_bins == "n_e"`:** Same pipeline and background; no fold to pattern; q_ts from n_e ROI; no ROOT/PNG/CSV pattern outputs.

### Efficiency-per-pattern table and printouts (notebook-style)

When running the single-point example in pattern mode, the app also:

- **Prints where efficiencies are computed:** One line stating `Efficiencies computed from: CSV: <path>` or `PatternMC (BuildPatternTable, N trials per n_e)`.
- **Writes** `run.outdir/efficiency_per_pattern.csv` with columns `pattern_id,n_e,efficiency` (P(pattern|n_e)) for each pattern in `pattern_roi` and each n_e in the binning range.
- **Verbosity ≥ 2:** Prints a table of P(pattern|n_e) (rows = n_e, columns = pattern_roi), like the notebook efficiency output.
- **Image visualization:** The app will plot the full **2D binned image** (from `response.pattern_image` when the 2D generator is implemented); the 1D row printout (PrintExampleRows) is not used.

---

## 8. Pattern efficiency from PatternMC when no CSV

**File:** `apps/ccdarksens_scan_dmelectron_pattern.cc`

- When `observable_bins == "pattern"` and **no** `response.pattern_mc.efficiency_csv` is set, the app fills `pattern_eff_map` from PatternMC’s internal table (same binning as the pipeline).
- Right after requiring non-empty `pattern_roi`, if `pattern_eff_map.empty()`:
  - Read `emc_ptr->GetPatternTable()` (map: n_e → (PatternLabel → P)).
  - For each (n_e, label → P), encode label as integer code: `code = 0; for (int d : label.q) code = code*10 + d` (e.g. [1,1] → 11).
  - Set `pattern_eff_map[{code, ne}] = P`.
- If the map is still empty after this (e.g. table not built), the app throws.  
- This allows running in pattern space without an external CSV when you want efficiencies computed in C++ with the same detector and binning.

---

## 9. Single-point example app (one mass, one cross section)

**New app:** `apps/ccdarksens_example_one_point_pattern.cc`  
**Build:** `ccdarksens_example_one_point_pattern` in `CMakeLists.txt`.

- **Purpose:** Run the full pipeline for **one** (m_χ, σ_e) to validate S_obs → (optional) S_pat/B_pat and PLR without doing a full grid scan.
- **Usage:**  
  `ccdarksens_example_one_point_pattern config.json [mchi_MeV] [sigma_e_cm2]`  
  If `mchi_MeV` and `sigma_e_cm2` are omitted, the first grid point from `model.grid` in the config is used.
- **Behaviour:** Same config/setup as the pattern scan app (detector response, PatternMC, backgrounds). For the single point it:
  - Builds dR/dE → pipeline → S_obs(n_e).
  - If `observable_bins == "pattern"`: fills `pattern_eff_map` from PatternMC when no CSV; folds to S_pat, B_pat; prints S_obs, S_pat, B_pat, model_test; computes q; writes `example_one_point_pattern.root` with `S_pat_validation`, `B_pat_validation`, `S_obs_ne`, `B_tot_ne`; and writes `example_one_point_pattern_rates_per_pattern.png` (rate per pattern bar chart) in the same directory.
  - If `observable_bins == "n_e"`: prints S_obs and q in n_e ROI.
- **Config:** Use the same JSON as the pattern scan (e.g. `configs/ccdarksens_scan_dmelectron_pattern.json`). Set `experiment.observable_bins` to `"pattern"` to exercise the pattern path and (when no CSV) the fill-from-PatternMC logic.

---

## 10. Option A: pattern acceptance applied only in the fold

Efficiencies must not be applied twice in either mode. When `observable_bins == "n_e"`, the pipeline and BackgroundBuilder apply ε_total(n_e) once; we do not fold to pattern, so there is no second application. When `observable_bins == "pattern"`, we apply acceptance only in the fold.

- **DetectorResponsePipeline:** `SetSkipPatternEfficiency(bool)`. When true (set by apps when `observable_bins == "pattern"`), `Apply` does **not** call `pe_->Apply`; it returns S_rec(n_e) instead of S_obs(n_e) = S_rec × ε_total.
- **BackgroundBuilder:** When `observable_bins == "pattern"`, the apps do **not** call `bld.SetPatternEfficiency(pe)`, so `BuildBkgAsimov()` returns B_raw(n_e) (no ε_total applied).
- **Fold:** S_pat[p] = Σ_n_e S_rec(n_e)×P(p|n_e), B_pat[p] = Σ_n_e B_raw(n_e)×P(p|n_e). Total accepted rate = Σ_p S_pat[p] = Σ_n_e S_rec(n_e)×ε_total(n_e).

**PatternMC and detector effects:** The pattern table P(pattern|n_e) from PatternMC is built with **ChargeTransport** (diffusion), **PixelSimulator** (readout noise, pixel geometry/binning, optional dark current via `lambda_dc`), and **PatternClassifier**. So readout noise, binning (row segment), diffusion, and dark current (if enabled in the MC config) are already included in the table; the fold then applies these probabilities once.
