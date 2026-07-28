# Pattern app flow walkthrough

Walkthrough of **ccdarksens_scan_dmelectron_pattern** and **ccdarksens_example_one_point_pattern**: config → response → pattern efficiency → outputs. Includes a short gap checklist.

---

## 1. Config and experiment setup

- **ConfigManager** parses JSON: `run`, `detector`, `experiment`, `response`, `model`, `backgrounds`, `timing`.
- **Detector**: `rows`, `cols`, `pixel_size_um`, `thickness_mm`, etc. Image size for pattern 2D is taken from here when `pattern_image` uses detector (raw_rows/raw_cols).
- **Experiment**: `binning` (ne_min, ne_max), `roi_bins`, `pattern_roi`, `observable_bins` ("pattern" or n_e bins).
- **Response**: `mode` = "pattern", `pattern_mc`, `pattern_image`, `pattern_classifier`, `pcd`.

---

## 2. Detector response components (both apps)

1. **ChargeIonization**  
   - Table from `data/p100K_table.csv`.  
   - Used to fold dR/dE → S_true(n_e).

2. **ChargeTransport**  
   - From detector thickness + `pattern_mc`: A_um2, b_umInv, alpha, beta_per_keV.  
   - Used by PatternMC, PatternImageGenerator, PCD.

3. **PatternMC**  
   - Config: `ne_trials`, `row_length` (or 2×half_window_pix+1), pixel/readout from detector and `pattern_mc`.  
   - **Accepted labels**: from `experiment.pattern_roi` (pattern codes → `PatternLabel` with `q` vector).  
   - **Efficiency source** (mutually exclusive in practice):
     - **CSV**: if `response.pattern_mc.efficiency_csv` is set, app loads (pattern_id, n_e, efficiency) into `pattern_eff_map`.
     - **PatternMC table**: if no CSV (or as fallback when `observable_bins == "pattern"` and map empty), table comes from `PatternMC::BuildPatternTable`.
   - **2D image efficiency**: if `pattern_mc.use_2d_image_efficiency` is true and `pattern_image` has valid binning (and detector size or nrows_binned/ncols), app builds **PatternImageGenerator** (image size from detector when possible, binning from `pattern_image`), calls **SetPatternImageGenerator**. Then **BuildPatternTable** (called from PrecomputeEpsilon / PrecomputeEpsilonWithPatternEff) uses 2D: generate image → ScanImage2DWithIsolation → one pattern per event → P(pattern|n_e).

4. **PatternClassifier**  
   - From `response.pattern_classifier`: Qmin_e, Qmax_e, thresholds, enable_MN/MNL.  
   - Used by PatternMC (1D row or 2D scan) and by the example app for middle-row scan / simulate_cluster.

5. **PatternImageGenerator** (when used)  
   - **Image size**: from detector (`raw_rows` = det.rows, `raw_cols` = det.cols) when available; else from `pattern_image.nrows_binned` / `ncols`.  
   - **Binning**: `pattern_image.row_binning`, `col_binning`.  
   - Produces 2D binned image (diffusion, row/col binning, readout noise, optional DC).  
   - Used for: (a) 2D efficiency in PatternMC when `use_2d_image_efficiency`, (b) example app 2D plot and middle-row scan.

6. **PCD**  
   - PCDBasedResponse + PCDCalculator.  
   - Pipeline uses them to build P(q|n_e) and P(n_obs|n_true), then S_true → S_rec(n_e).

7. **DetectorResponsePipeline**  
   - **Pattern mode**:  
     - dR/dE → **S_true(n_e)** (ChargeIonization).  
     - If PCD: build kernel P(n_obs|n_true), fold S_true → **S_rec(n_e)**.  
     - If no PCD: optional Diffusion in n_e.  
     - Then **PatternEfficiency** ε(n_e) applied (unless `SetSkipPatternEfficiency(true)` for pattern bins) → **S_obs(n_e)**.  
   - Pipeline does **not** call PatternMC itself; it only uses the **precomputed** ε(n_e) histogram (PatternEfficiency).  
   - `SetPatternMC` is commented out in the scan app: no E-dependent ε(E,n_e) in the pipeline for this flow.

---

## 3. Pattern efficiency ε(n_e)

- **PrecomputeEpsilon(ne_min, ne_max, Ee_ref_eV)**  
  - Calls **BuildPatternTable** (1D row or 2D image path if generator set).  
  - Builds a single TH1D ε(n_e) at Ee_ref_eV (e.g. 50 eV).  
  - Used when no CSV is provided.

- **PrecomputeEpsilonWithPatternEff(..., pattern_eff_map)**  
  - If CSV was loaded, uses **pattern_eff_map** to fill ε(n_e) (and still calls BuildPatternTable for internal table if needed elsewhere).  
  - So when `efficiency_csv` is set, the **CSV** drives the histogram; when not set, **PatternMC** (1D or 2D) drives it.

- **PatternEfficiency**  
  - Wraps the TH1D; pipeline calls `pe_->Apply(*h_ne_rec)` to get S_obs.

- **Observable = pattern bins**  
  - When `observable_bins == "pattern"`, app folds S_obs(n_e) and B(n_e) into **pattern rates** via **FoldNeToPatternRates** using `pattern_roi` and `pattern_eff_map`.  
  - So the same efficiency map used to build ε(n_e) (or from CSV) must be used for folding to pattern space.

---

## 4. Backgrounds

- **BackgroundBuilder**: detector rows/cols, active_fraction, ne_min, ne_max, timing, dark current.  
- **B_dc_ne**: dark current in n_e (BuildBkgAsimov).  
- **B_flat_ne**: flat dR/dE through pipeline → n_e.  
- **B_tot** = B_dc_ne + B_flat_ne.  
- If pattern bins: **B_pat** = FoldNeToPatternRates(B_tot, pattern_roi, pattern_eff_map).

---

## 5. Scan app (ccdarksens_scan_dmelectron_pattern) flow

1. Parse config; create outdir.  
2. Experiment setup (binning, pattern_roi, roi_bins).  
3. Build ion, ChargeTransport, PatternMC config, classifier; **optional** efficiency CSV load.  
4. **Optional 2D efficiency**: if use_2d_image_efficiency and pattern_image + detector ok → PatternImageGenerator from detector size + binning → SetPatternImageGenerator.  
5. **ε(n_e)**:
   - If CSV loaded → PrecomputeEpsilonWithPatternEff(..., pattern_eff_map);  
   - else → PrecomputeEpsilon(...).  
6. **Hardcoded overrides** (see Gaps): bins 1–5 of h_eps_ne are overwritten (1.0, 0.38, 0.65, 0.79, 0.86).  
7. PCD + Diffusion; pipeline with PatternEfficiency(ε), SetSkipPatternEfficiency(use_pattern_bins).  
8. Backgrounds: B_dc_ne, B_flat_ne, B_tot; if pattern bins fill pattern_eff_map from PatternMC table if empty, then B_pat.  
9. Grid over (mchi, sigma): for each point, dR/dE → pipeline.Apply → S_obs; fold to S_pat if pattern bins; PLR with data=B_pat, model_null=B_pat, model_test=S_pat+B_pat; fill TH2D q.  
10. Write ROOT: B_tot_ne, q_mchi_sigma_pattern, validation histos.

---

## 6. Example app (ccdarksens_example_one_point_pattern) flow

1. Parse config; optional (mchi, sigma) from CLI or first grid point.  
2. Same detector response chain as above; optional CSV; optional 2D generator for PatternMC.  
3. ε(n_e) as in scan (with hardcoded overrides for bins 1–5).  
4. Pipeline, backgrounds, B_tot, B_flat_ne.  
5. If **observable_bins == "pattern"**:  
   - Fill pattern_eff_map from PatternMC table if no CSV.  
   - Optional: efficiency CSV out, verbose table, **E-dependent ε CSV** (debug.output_epsilon_E_csv), **2D image** (generate, plot PNG, middle-row scan, simulate_cluster, charge histos, sanity CSV).  
   - B_pat = FoldNeToPatternRates(B_tot, ...).  
6. Single signal point: dR/dE → pipeline → S_obs; fold to S_pat; PLR; optional ROOT/PNG outputs.

---

## 7. Potential gaps / things to double-check

| # | Item | Where | Note |
|---|------|--------|------|
| 1 | **Hardcoded ε(n_e) overrides** | Scan app ~L419–437; example ~L302–308 | Bins 1–5 of h_eps_ne are overwritten (1.0, 0.38, 0.65, 0.79, 0.86). If intentional (e.g. tuning), consider moving to config or a debug flag; otherwise remove so CSV or PatternMC result is used. |
| 2 | **E-dependent pattern efficiency** | Pipeline | `EnableEDependentPattern(false)` and `SetPatternMC` commented out. So ε is at a single Ee_ref_eV. If you want ε(E,n_e), pipeline would need to use PatternMC and E-dependent ε. |
| 3 | **2D image size vs full sensor** | pattern_image | With detector 1300×6300 and row_binning=100, col_binning=1, binned image is 13×6300. For a small “pattern window” (e.g. 3×50) you’d need different binning or a separate pattern-region size (not yet in config). |
| 4 | **Consistency: 2D generator vs CSV** | Config | If use_2d_image_efficiency is true but efficiency_csv is also set, the **CSV** is used for the ε(n_e) histogram; the 2D generator is only used for PatternMC’s internal table (e.g. if something else used it). So for “efficiency from 2D” you’d set use_2d_image_efficiency and omit efficiency_csv. |
| 5 | **pattern_eff_map and pattern_roi** | Both apps | When observable_bins is "pattern", pattern_eff_map must cover all (pattern_id, n_e) for pattern_roi × ne_min..ne_max. CSV and PatternMC table must use same pattern coding (e.g. 11, 21, 111) as pattern_roi. |
| 6 | **Skip pattern efficiency** | Pipeline | SetSkipPatternEfficiency(use_pattern_bins): when true (pattern bins), pipeline does not apply ε again because the fold to pattern space applies it. Ensure this matches your likelihood (single application of ε). |
| 7 | **Dark current in pipeline vs EfficiencyMC** | Config | lambda_per_exp computed but pix_cfg_pcd.lambda_dc = 0 and emc_cfg.pix_cfg.lambda_dc = 0 in scan app. DC is in BackgroundBuilder. If you want DC in pattern simulation (pileup), enable pileup_with_dc / pattern_image.lambda_dc and generator include_dark_current. |
| 8 | **Divisibility** | PatternImageGenerator | When using detector size, raw_rows % row_binning and raw_cols % col_binning must be 0; otherwise constructor throws. Config should ensure compatible detector.rows/cols and binning. |

---

## 8. Summary diagram (pattern path)

```
Config (JSON)
  → detector (rows, cols, ...)
  → response.pattern_mc (efficiency_csv?, use_2d_image_efficiency?)
  → response.pattern_image (row_binning, col_binning; optional nrows_binned, ncols)
  → experiment (pattern_roi, observable_bins)

App
  → ChargeIonization, ChargeTransport, PatternClassifier
  → PatternMC (accepted_labels from pattern_roi)
  → [Optional] PatternImageGenerator (size from detector, binning from pattern_image) → SetPatternImageGenerator
  → [Optional] Load efficiency_csv → pattern_eff_map
  → h_eps_ne = PrecomputeEpsilon(...) or PrecomputeEpsilonWithPatternEff(..., pattern_eff_map)
  → [Scan app] Overwrite h_eps_ne bins 1–5 (check if desired)
  → PatternEfficiency(ε) → pipe.SetPatternEfficiency(pe), SetSkipPatternEfficiency(use_pattern_bins)
  → PCD + Diffusion → pipeline.Apply(dR/dE) → S_obs(n_e)
  → If pattern bins: FoldNeToPatternRates(S_obs, pattern_roi, pattern_eff_map) → S_pat
  → Backgrounds: B_tot → B_pat if pattern bins
  → PLR(data=B_pat, model_test=S_pat+B_pat, model_null=B_pat) → q
  → Outputs: ROOT, optional CSV/PNG (example app).
```

If you want, we can next (a) remove or make configurable the hardcoded ε overrides, (b) add a small “pattern region” (rows/cols) separate from full detector, or (c) wire E-dependent ε in the pipeline.
