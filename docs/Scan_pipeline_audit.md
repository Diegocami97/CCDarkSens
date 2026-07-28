# Scan pipeline audit: `ccdarksens_scan_dmelectron_pattern` + `configs/scan_dmelectron_pattern_pydme_minuit.json`

Step-by-step trace of where scaling/normalization is applied. Goal: confirm there is **no double or hidden scaling** of signal or exposure.

---

## Step 1: Config and experiment summary

- **Config**: `ConfigManager::parse()` reads JSON; `ExperimentSetup(setup(cfg.experiment_cfg(), det.mass_kg(), run.rng_seed)`.
- **Detector mass**: `det.mass_kg()` = config `detector.mass_kg` (0.01523 kg) if set; else `compute_mass_from_geometry_kg()` (rows×cols×pixel_size×thickness×active_fraction×density).
- **Exposure** (single place it is defined):
  ```text
  exposure_kg_year = livetime_days * duty_cycle * mass_kg / 365.25
  ```
  Implemented in `ExperimentSetup::prepare_summary()` (src/experiment/ExperimentSetup.cc).
- With livetime_days=85.356, duty_cycle=1, mass_kg=0.01523 → **exposure_kg_year ≈ 0.00356 kg·year** (~1.3 kg·day).  
- **No other factor** (e.g. active_fraction) is applied to exposure; active_fraction is already inside mass if mass is computed from geometry.

---

## Step 2: Signal rate table (dR/dE)

- **Model**: `DMElectronModel::Configure(mc)` resolves path from `rates_dir` + `filename_template` with `{mchi_MeV}`, `{sigma_e_cm2}`.
- **Load**: `RateTable::LoadCSV(path)` reads CSV; second column stored as **R_kg_year_eV_** (events/(kg·year·eV)); no conversion (RateTable.cc line 95 uses `y` as-is; the old `y * g_to_kg` is commented out).
- **Spectrum**: `MakeSpectrum_E()` = `table_->MakeTH1D(..., Emin_eV, Emax_eV, nbins)`; bin content = interpolated R_kg_year_eV. **No scaling**; each CSV is for one (mχ, σ_e), so rate is already for that cross section.

---

## Step 3: Rate → counts (n_e): FoldToNe

- **Call**: `ion->FoldToNe(*dRdE_sig, summary.exposure_kg_year, ne_min, ne_max)` (ChargeIonization.cc).
- **Formula** (per bin of dRdE):
  ```text
  rate = dRdE.GetBinContent(i)   // events/(kg·year·eV)
  counts = rate * exposure_kg_year * dE
  ```
  Then `counts` is distributed in n_e using P(n_e|E) from the ionization table.  
- **Exposure is applied exactly once** here for signal. No other multiplication by exposure or mass in the signal path.

---

## Step 4: n_e counts → pattern counts: FoldNeToPatternRates

- **Call**: `FoldNeToPatternRates(*S_true, ne_min, ne_max, pattern_roi, pattern_eff_map)` (PatternRates.cc).
- **Formula**: For each pattern in pattern_roi,  
  `S_pat[i] = sum_{n_e} S_true(n_e) * epsilon(pattern_i | n_e)`.  
  So **counts in → counts out**; efficiencies are P(pattern|n_e), no extra exposure or global scale.

---

## Step 5: Background (this config: bp_br_template)

- **background_source == "bp_br_template"**:  
  `B_pat[i] = run.background_Bp[i] + run.background_Br[i]` (no exposure or scaling in code).  
  So **B_pat** are fixed counts per pattern; they must already correspond to the **same exposure** as the one used for signal (config exposure). If they were derived with a different exposure, that would be an external mismatch, not a second scaling inside the app.

- Flat background (when used): `flat_rate_per_eV = flat_bkg_norm_per_kg_year / 1000` (config is per keV → divide by 1000 for per eV). Then `FoldToNe(*dRdE_flat, summary.exposure_kg_year, ...)` — exposure applied once. Dark current: from BackgroundBuilder with `norm_scale` (1.0) and N_active_pixels × N_exposures; no exposure in the sense of kg·year (it’s per-pixel, per-exposure).

---

## Step 6: Profile likelihood

- **Inputs**: `SetData(data)` = D_pat (counts per pattern); `SetBpBr(Bp, Br)`; for each grid point, **S_pat** = FoldNeToPatternRates(S_true, …) with **S_true = FoldToNe(dRdE_sig, summary.exposure_kg_year, …)**.
- **Model**: μ_i = S_pat[i] + θ·Br[i] (with Bp+Br template). S_pat is **counts**; θ is a scale on the reference component. No extra scaling of S_pat or exposure inside ProfileLikelihood.

---

## Step 7: Data (D_pat)

- **load_data(run.data_path, …)**: reads histogram **D_pat** from ROOT (or CSV). Values are **counts per pattern**, same order as pattern_roi. **No scaling** applied when loading.
- **exposure_kg_year** in the data file is only used for the **comparison warning** and for writing `exposure_kg_year_from_data_file` in the output; the scan **always** uses **config exposure** (summary.exposure_kg_year) for S_pat and B_pat.

---

## Summary table (scaling / normalization)

| Step | What | Scaling applied? |
|------|------|-------------------|
| 1 | exposure_kg_year | livetime_days × duty_cycle × mass_kg / 365.25 (once) |
| 2 | dR/dE from CSV | None; CSV in events/(kg·year·eV) |
| 3 | FoldToNe | counts = rate × **exposure_kg_year** × dE (once) |
| 4 | FoldNeToPatternRates | counts → counts with ε(pattern\|n_e); no extra scale |
| 5 | B_pat (bp_br_template) | None in code; user must supply for same exposure |
| 5b | Flat bkg | norm_per_kg_year_keV/1000 → per eV; then FoldToNe with exposure once |
| 6 | ProfileLikelihood | μ = S_pat + θ·Br; no extra scaling |
| 7 | Data D_pat | Loaded as-is (counts) |

**Conclusion:** There is **no additional or double scaling** of signal or exposure in this pipeline. The only exposure factor is **exposure_kg_year** in FoldToNe. If limits still disagree with a reference, check: (1) same exposure and mass, (2) Bp/Br and data at that exposure, (3) same rate units (events/(kg·year·eV)), (4) same CL/UL definition and test statistic.
