# Flow of `scan_dmelectron_ne_pydme_minuit.json` in the app

Walkthrough of how each part of the JSON is used when running  
`./build/ccdarksens_scan_dmelectron_pattern configs/scan_dmelectron_ne_pydme_minuit.json`.

---

## 1. Parse config

- **ConfigManager** reads the JSON and fills internal structs:
  - `run` → run options (label, outdir, data_path, use_profile_likelihood, background_source, theta_lo/hi, etc.)
  - `detector` → geometry (rows, cols, pixel_size_um, thickness_mm, mass_kg, target_element, density)
  - `experiment` → mode, livetime_days, duty_cycle, binning (ne_min, ne_max), roi_bins, pattern_roi, observable_bins
  - `backgrounds` → dark_current (lambda, norm_scale), timing, flat_background, background_efficiency_csv
  - `response` → pattern_mc (efficiency_csv, pattern_image, pattern_classifier, pcd)
  - `model` → type, material, mediator, rates_dir, filename_template, Emin/Emax, nbins, grid
- **jroot** keeps the raw JSON for a few keys that ConfigManager doesn’t expose (e.g. `response.pattern_mc.efficiency_csv_reference`).

---

## 2. Experiment setup

- **ExperimentSetup(setup)** uses `experiment` + `detector.mass_kg`:
  - **Exposure:** `exposure_kg_year = livetime_days * duty_cycle * mass_kg / 365.25`  
    (from `experiment.livetime_days`, `duty_cycle`, `detector.mass_kg`).
  - **summary.binning:** `ne_min` = 1, `ne_max` = 5 from `experiment.binning`.
  - **summary.roi_bins:** [1, 2, 3, 4, 5] from `experiment.roi_bins`.
  - **summary.observable_bins:** `"ne"` from `experiment.observable_bins`.
- So **use_pattern_bins = false** (observable is n_e, not pattern).

---

## 3. Detector response and efficiencies

- **ChargeIonization**, **ChargeTransport**, **PatternMC** are built from `detector` + `response.pattern_mc` (A_um2, b_umInv, alpha, sigma_readout_e, half_window_pix, row_length, etc.).
- **pattern_eff_map:** filled from `response.pattern_mc.efficiency_csv` (and optionally `efficiency_csv_reference`).  
  Format: (pattern_id, n_e) → efficiency. Used to build ε(n_e) for the pipeline.
- **PatternMC** gets accepted labels from `experiment.pattern_roi`; for n_e config `pattern_roi` is [], so a default single-pixel label is used. The important part for n_e is the **efficiency CSV**, which gives ε(ne) for ne = 1..5.
- **PrecomputeEpsilon(ne_min_bkg, ne_max, Ee_ref, pattern_eff_map)** builds a histogram **h_eps_ne** = ε(n_e) used by the pipeline and by **BackgroundBuilder** when `!use_pattern_bins` (so background is folded with the same ε(n_e)).

---

## 4. Background (DC + flat, no Bp/Br)

- **BackgroundBuilder** is created with `detector.rows/cols/active_fraction` and `ne_min_bkg`, `ne_max`.
- **Timing:** from `backgrounds.timing` (exposure_time_s, n_exposures) and `experiment.livetime_days`, `duty_cycle`.
- **Dark current:** from `backgrounds.dark_current` (lambda_e_per_pix_per_year, norm_scale).
- Because **use_pattern_bins == false**, **bld.SetPatternEfficiency(pe)** is called so DC (and flat) are folded with ε(n_e).
- **B_dc_ne** = `bld.BuildBkgAsimov()` (DC in n_e with ε applied).
- **Flat background:** `model` is used to build a flat dR/dE spectrum; `backgrounds.flat_background` gives norm (norm_per_kg_year_keV → flat_rate_per_eV). It is folded through the pipeline with ε(n_e): **B_flat_ne** = `pipe.Apply(dRdE_flat, exposure_kg_year, ne_min_bkg, ne_max, Ee_ref_eV)`.
- **B_tot** = B_dc_ne + B_flat_ne.  
  So background comes only from **backgrounds** (DC + flat) and **response** efficiencies; **run.background_source** is `"dc_flat"`, so Bp/Br are never used (only in pattern-bins path).

---

## 5. n_e observable and profile likelihood (your config)

- Because **use_pattern_bins** is false, the **else** branch runs (n_e observable).
- **B_ne:** for each `ne` in **summary.roi_bins** [1,2,3,4,5], take `B_tot->GetBinContent(ne)` → vector of 5 numbers.
- **run.use_profile_likelihood** is true → **profile_pl** is created.
  - **Data:** `run.data_path` is "" → **Asimov:** `profile_pl->SetData(B_ne)` (data = background).
  - **Background model:** B = θ·B_ne → `profile_pl->SetBpBr(zeros, B_ne)`.
  - **profile_S_null** = vector of 5 zeros (signal = 0).
  - Constraints from **run:** constrain_prior_strength, constrain_use_gamma_sign, constrain_use_tau_weighted, constrain_n_bins, theta_lo, theta_hi, pydme_style_ul.

---

## 6. Grid and scan loop

- **model.grid** is read from **jroot["model"]["grid"]**:  
  **mchi_list** and **sigma_list** from `grid.mchi_MeV` and `grid.sigma_e_cm2` (expand_axis for "values" or "logspace").
- For each **mchi** in mchi_list:
  - **n_e profile path:** nll_null = minimize NLL over θ at S = 0 (profile_S_null); reserve nll_values.
  - For each **sigma** in sigma_list:
    - **DMElectronModel** is configured from **model** (material, mediator, rates_dir, filename_template, Emin_eV, Emax_eV, nbins) + current mchi, sigma.
    - **dRdE_sig** = load rate from CSV (rates_dir + filename_template).
    - **S_true** = ion->FoldToNe(dRdE_sig, exposure_kg_year, ne_min, ne_max).
    - **S_obs** = pipe.Apply(dRdE_sig, exposure_kg_year, ne_min_bkg, ne_max, Ee_ref_eV) (signal in n_e with ε(n_e)).
    - **S_for_pl:** for each ne in roi_bins, S_for_pl[i] = S_obs->GetBinContent(ne).
    - **nll** = profile_pl->MinimizeOverScale(S_for_pl, theta_lo, theta_hi) (or Minuit); push into nll_values.
  - After the sigma loop: **nll_min** = min(nll_values); **q_μ** = 2(nll − nll_min); **UL** = smallest σ where q_μ ≥ target_q; fill **h_q**, **h_upper_limit**, **h_q0**.

---

## 7. Output

- **run.outdir** → create directory; output ROOT file: `outdir/scan_dmelectron_pattern.root`.
- Written: B_tot (if not bp_br_template), h_q, h_q0, h_upper_limit, exposure_kg_year.  
  No D_pat/D_ne when data_path is empty (Asimov).  
  **run.label** is only for your bookkeeping; the app doesn’t encode it in the file name.

---

## 8. Summary table (JSON → usage)

| JSON section        | Used for |
|---------------------|----------|
| **run**             | outdir, data_path (empty ⇒ Asimov), use_profile_likelihood, background_source (dc_flat ⇒ no Bp/Br), theta_lo/hi, constrain_*, profile_minimizer, pydme_style_ul, cl, verbosity |
| **detector**        | rows, cols, pixel_size_um, thickness_mm, active_fraction, mass_kg (exposure), target_element, density_g_cm3 |
| **experiment**      | mode (asimov), livetime_days, duty_cycle → exposure; binning (ne_min=1, ne_max=5); roi_bins [1..5]; observable_bins "ne"; pattern_roi [] |
| **response**        | pattern_mc (efficiency_csv → ε(n_e), pattern_image, classifier), pcd |
| **backgrounds**     | dark_current (lambda, norm_scale), timing, flat_background (norm, E range, nbins), background_efficiency_csv ("" ⇒ no migration CSV) |
| **model**           | material, mediator, rates_dir, filename_template, Emin/Emax, nbins; **model.grid** → mchi_MeV and sigma_e_cm2 lists for the scan |

For your Asimov n_e config: **data** = B_ne, **background** = DC + flat with same efficiencies as for data (from efficiency_csv), **observable** = 5 bins (1e–5e), **profile** over θ with B = θ·B_ne.
