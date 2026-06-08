# Data units and using real data in the framework

## Units (must match everywhere)

| Quantity | Units | Where |
|----------|--------|--------|
| **D_pat, S_pat, B_pat** | **Counts** (dimensionless) per pattern bin | Same order as `pattern_roi`. Expected counts for the **exposure** below. |
| **exposure_kg_year** | **kg·year** | Config: `livetime_days * duty_cycle * detector_mass_kg / 365`. Data ROOT: from `texp` (days) and `Nusedpix` × mass per pixel. |
| **texp** in CSV | **days** | Per-row exposure time (e.g. LBC/pydme style). |
| **Nusedpix** in CSV | **pixels** | Number of (unmasked) pixels for that row. |
| **mass_per_pixel_kg** (converter default) | **kg** | Si, 15 µm × 15 µm × 0.67 mm, ρ = 2.33 g/cm³ → 3.51e-10 kg. |

So **S_pat** and **B_pat** are expected **counts** for the exposure defined in the config. **Data** (D_pat) must be observed **counts** for the **same** exposure so the likelihood is consistent.

## Converter: CSV → ROOT

- **App:** `ccdarksens_csv_to_root_data`
- **Output:** ROOT file with `D_pat` (TH1D, counts per pattern), and when CSV has `texp` and `Nusedpix`, also `exposure_kg_year` and `mass_kg` (TParameter).
- **Exposure formula:** `exposure_kg_year = sum(texp * Nusedpix * mass_per_pixel_kg) / 365` (same 365 as `ExperimentSetup`).

## Using the data in the scan or example app

1. Convert the data CSV to ROOT (optional but recommended for exposure check):
   ```bash
   ./build/ccdarksens_csv_to_root_data data/Final_Combined_Image_Data.csv build/data_pattern.root configs/ccdarksens_scan_dmelectron_pattern.json
   ```
2. Set **run.data_path** in the config to either:
   - the **ROOT** file path (e.g. `build/data_pattern.root`), or  
   - a **CSV** path (one row of counts in `pattern_roi` order).
3. **Match exposure:** Set **experiment.livetime_days** and detector **mass** so that the config’s `exposure_kg_year` matches the data. If you use the ROOT file from the converter, it writes `exposure_kg_year`; the scan/example app will **warn** when it differs from the config by more than 1%.
4. Run with **use_profile_likelihood: true** and your **data_path**; the app loads D_pat (from ROOT) or counts (from CSV) and uses them as the observed counts in the likelihood.

## Summary

- **Counts:** D_pat, S_pat, B_pat are all in **counts** per pattern.
- **Exposure:** Same **kg·year** for data and for S_pat/B_pat; match config to data (or data to config) and use the exposure warning as a check.
- **No full scan** is run here; this only ensures units and data loading are ready when you enable **data_path** and run.
