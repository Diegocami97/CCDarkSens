# ccdarksens_example_one_point_pattern — App flow

Single-point run: one (m_χ, σ_e) through the full pattern pipeline → S_obs(n_e) → pattern rates → PLR.

**Usage:** `ccdarksens_example_one_point_pattern config.json [mchi_MeV] [sigma_e_cm2]`  
If mchi/sigma omitted, the first grid point from `config.model.grid` is used.

---

## 1. Startup and config

- Parse **config.json** (ConfigManager + raw json `jroot`).
- **Run**: outdir, rng_seed, verbosity; create outdir.
- **ExperimentSetup**: binning (ne_min, ne_max), pattern_roi, roi_bins, **observable_bins** ("pattern" or n_e), exposure.
- **(mchi, sigma_e)**:
  - If CLI has 4 args → use `argv[2]`, `argv[3]`.
  - Else → expand `model.grid.mchi_MeV` and `model.grid.sigma_e_cm2`, take first point.

---

## 2. Detector response chain

- **ChargeIonization**: table `data/p100K_table.csv` (dR/dE → n_e).
- **ChargeTransport**: from detector thickness + response.pattern_mc (A_um2, b_umInv, alpha, beta, rng_seed).
- **PatternMCConfig**: ne_trials, row_length (from config or 2×half_window_pix+1), pixel/readout from detector and pattern_mc. **Accepted labels**: from experiment.pattern_roi (pattern codes → `PatternLabel.q`).
- **Efficiency CSV (optional)**:
  - If `response.pattern_mc.efficiency_csv` present → load (pattern_id, n_e, efficiency) into **pattern_eff_map**.
  - Else → pattern_eff_map stays empty; efficiency will come from PatternMC (1D or 2D).
- **PatternClassifier**: from response.pattern_classifier (Qmin_e, Qmax_e, thr_M, thr_MN, thr_MNL, enable_MN/MNL).
- **PatternMC** = PatternMC(pmc_config, ct, classifier); `include_dc_pileup` from config.
- **2D image efficiency (optional)**:
  - If `pattern_mc.use_2d_image_efficiency` and valid pattern_image (detector size or nrows_binned/ncols + binning):
  - Build **PatternImageConfig** (raw from detector or binned from JSON), then **PatternImageGenerator**, call **emc_ptr->SetPatternImageGenerator(img_gen_eff)**.
  - Then when ε(n_e) is built, PatternMC will use 2D: generate image → ScanImage2DWithIsolation → P(pattern|n_e).

---

## 3. Pattern efficiency ε(n_e)

- **Ee_ref_eV** = 50 eV.
- If **pattern_eff_map** not empty (CSV loaded) → **PrecomputeEpsilonWithPatternEff**(ne_min, ne_max, Ee_ref_eV, pattern_eff_map) → fills ε from CSV.
- Else → **PrecomputeEpsilon**(ne_min, ne_max, Ee_ref_eV) → calls **BuildPatternTable** (1D row or 2D image if generator set), then builds TH1D ε(n_e).
- **Hardcoded overrides**: bins for n_e=1..5 are overwritten (1.0, 0.38, 0.65, 0.79, 0.86).

---

## 4. Pipeline and PCD

- **PCD**: PixelSimulatorConfig (5×5 patch), PCDResponseConfig, PCDBasedResponse, PCDCalculator.
- **DetectorResponsePipeline**: ion, Diffusion, PCD response, PCD calculator, **AnalysisSpace::Pattern**.
- **PatternEfficiency** wraps h_eps_ne; **SetPatternEfficiency(pe)**, **SetSkipPatternEfficiency(use_pattern_bins)** (skip when observable = pattern bins, because fold applies ε).
- **EnableEDependentPattern(false)** — single E.

---

## 5. Backgrounds

- **BackgroundBuilder**: detector rows/cols, active_fraction, ne range, timing (livetime, duty cycle, exposure_time_s), dark current (lambda per year, norm_scale). If not pattern bins, set PatternEfficiency for bkg.
- **B_dc_ne** = BuildBkgAsimov() (dark current in n_e).
- **Flat bkg**: dummy DM model → dRdE_flat (constant rate); **B_flat_ne** = pipe.Apply(dRdE_flat, exposure, ne_min, ne_max, Ee_ref_eV).
- **B_tot** = B_dc_ne + B_flat_ne.

---

## 6. When observable_bins == "pattern"

- Require **pattern_roi** non-empty.
- If **pattern_eff_map** empty → fill from **PatternMC::GetPatternTable()** (pattern code, n_e → efficiency).
- Require pattern_eff_map non-empty.
- **efficiency_per_pattern.csv**: write (pattern_id, n_e, efficiency) for pattern_roi × ne.
- **Verbosity ≥ 2**: print table P(pattern|n_e).
- **Optional debug** (from `debug` in JSON):
  - **output_epsilon_E_csv**: PrecomputeEpsilonVsEnergy(E_grid), write epsilon_vs_E.csv.
  - **2D image block** (if pattern_image valid):
    - Build **PatternImageGenerator** (same logic as for efficiency: detector or nrows_binned/ncols).
    - Generate one 2D image (e.g. n_e=5), plot **example_2d_pattern_image.png**, scan middle row with classifier, print best pattern.
    - **simulate_cluster(b,c,d)** (if debug.simulate_cluster): 3×5 cluster + noise, classify, print.
    - **output_charge_histos**: charge distribution per n_e → ROOT.
    - **output_sanity_csv**: event, n_e, pattern_id, total_charge_e → sanity_check.csv.
- **B_pat** = FoldNeToPatternRates(B_tot, ne_min, ne_max, pattern_roi, pattern_eff_map).

---

## 7. Single signal point

- **DMElectronModel**: configure (mchi_MeV, sigma_e_cm2) from chosen point.
- **dRdE_sig** = MakeSpectrum_E().
- **S_obs** = pipe.Apply(dRdE_sig, exposure_kg_year, ne_min, ne_max, Ee_ref_eV)  
  → dR/dE → S_true(n_e) → [PCD kernel] → S_rec(n_e) → ε(n_e) → **S_obs(n_e)**.
- Print S_obs(n_e) and integral.

---

## 8. Test statistic and outputs (pattern bins)

- **S_pat** = FoldNeToPatternRates(S_obs, ne_min, ne_max, pattern_roi, pattern_eff_map).
- **model_test** = S_pat + B_pat; **model_null** = B_pat; **data** = B_pat (Asimov).
- **PLR** → **q_ts** = EvaluateRatio(data, model_test, model_null).
- **ROOT**: example_one_point_pattern.root — S_obs_ne, B_tot_ne, S_pat_validation, B_pat_validation, optional charge histos.
- **PNG**: example_one_point_pattern_rates_per_pattern.png — bar chart S_pat vs B_pat per pattern bin.

---

## 9. If observable_bins ≠ "pattern" (n_e ROI)

- data / model_null / model_test built from ROI n_e bins (S_obs and B_tot per bin).
- q_ts = PLR on those vectors; no pattern fold, no ROOT/PNG pattern outputs.

---

## Flow diagram (pattern path)

```
config.json
    ↓
(mchi, sigma_e) [CLI or first grid point]
    ↓
ChargeIonization, ChargeTransport, PatternMC(+ optional 2D generator), PatternClassifier
    ↓
pattern_eff_map [CSV or empty]
    ↓
h_eps_ne = PrecomputeEpsilon(...) or PrecomputeEpsilonWithPatternEff(..., pattern_eff_map)
    ↓ [overwrite bins 1–5]
PatternEfficiency(ε) → Pipeline (PCD, ε) → pipe.Apply(dRdE) → S_obs(n_e)
    ↓
Backgrounds: B_dc_ne, B_flat_ne → B_tot
    ↓
If pattern bins: pattern_eff_map from table if empty → B_pat = FoldNeToPatternRates(B_tot,...)
    ↓
S_obs = pipe.Apply(dRdE_sig,...) → S_pat = FoldNeToPatternRates(S_obs,...)
    ↓
PLR(B_pat, S_pat+B_pat, B_pat) → q
    ↓
ROOT + PNG (pattern validation and rates)
```
