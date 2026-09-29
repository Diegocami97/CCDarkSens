# QCDark2 n_e Exposure Projections Guide

A step-by-step walkthrough for collaboration members who want **QCDark2** DM-electron **sensitivity projections** in **n_e space** with **flat dark current + flat d.r.u.** backgrounds, and to compare limits at different exposures (0.5, 1, and 2 kg·year).

This complements the pattern-count study guide ([`QCDark2_Pattern_Counts_Guide.md`](QCDark2_Pattern_Counts_Guide.md)) and the QEDark LBC reproduction guide ([`LBC_QEDark_Reproduction_Guide.md`](LBC_QEDark_Reproduction_Guide.md)). Here the focus is **Si_comp dielectric rates**, **Asimov mode**, and **exposure scaling** via detector mass.

---

## 1. What this workflow does

| Item | Value |
|------|-------|
| Signal model | QCDark2 (Si **Si_comp** composite dielectric, scissor gap 1.2 eV) |
| Observable | n_e bins: 1, 2, 3, 4, 5 |
| Background | DC rate + flat d.r.u., folded to B(n_e); nuisance **θ** scales all of B |
| Data | **Asimov** — expected counts S = B at each (mχ, σₑ) point |
| Statistic | Profile likelihood, 90% CL upper limit on σₑ(mχ) |
| Mediator | **Heavy** (rates from `Si_comp_long_scan`) |

**Example configs** (`configs/examples/`):

| Exposure (kg·year) | Scan JSON |
|--------------------|-----------|
| 0.5 | `qcdark2_lbc_ne_flatbkg_proj_0p5kgy.json` |
| 1.0 | `qcdark2_lbc_ne_flatbkg_proj_1p0kgy.json` |
| 2.0 | `qcdark2_lbc_ne_flatbkg_proj_2p0kgy.json` |

**Rate generation** (one-time, if tables are missing):

| Step | Config |
|------|--------|
| Si_comp dense grid | `qcdark2_generate_si_comp_dense.json` |

**Dense grid:** 800 masses (0.2–1000 MeV) × 300 cross sections (10⁻⁴⁶–10⁻²⁶ cm²).

---

## 2. Prerequisites

### Software

- C++17, CMake ≥ 3.18, ROOT (with Minuit2), nlohmann/json
- Python 3.8+ with NumPy and the `ccdarkphys` package (for QCDark2 rate generation, if needed)

### Build

```bash
cmake -B build -S .
cmake --build build -j8 --target ccdarksens_scan_dmelectron_pattern ccdarksens_plot_dmelectron_limit
```

### Input files

| File | Purpose |
|------|---------|
| `data/qcdark2_rates/Si/heavy/Si_comp_long_scan/` | Pre-generated dR/dE CSVs (may already exist in the repo) |
| `data/Efficiencies_patterns_Nsims1000000_DCTrue_alpha1.csv` | P(pattern \| nₑ) — used to build ε(n_e) for folding |
| QCDark2 `Si_comp.h5` | Only needed if you regenerate rates, from the QCDark2 package's `dielectric_functions/composite/Si_comp.h5` |

---

## 3. Pipeline overview

```
Step 1   Ensure QCDark2 Si_comp rate CSVs exist   (Python, one-time)
           ↓
Step 2   Pick exposure (mass_kg / livetime)       (edit scan JSON)
           ↓
Step 3   Run n_e-space Asimov scan                (C++)
           ↓
Step 4   Plot and compare exposure curves         (C++)
```

---

## 4. Step 1 — QCDark2 Si_comp rate tables

Rates are **dR/dE** CSVs (events / kg / year / eV), one file per (mχ, σₑ) on the grid. The scan configs point to:

```
data/qcdark2_rates/Si/heavy/Si_comp_long_scan/
```

If that directory is already populated (800×300 grid), **skip this step**.

To generate or refresh tables:

1. Obtain the production composite dielectric `Si_comp.h5` from the QCDark2 package (`dielectric_functions/composite/Si_comp.h5`).
2. Edit `epsilon_h5` in the generate config if your copy lives elsewhere.
3. Run:

```bash
python3 utils/qcdark2_generate_grid.py configs/examples/qcdark2_generate_si_comp_dense.json
```

Output: `data/qcdark2_rates/Si/heavy/Si_comp_long_scan/`

Generation is parallelized and skips existing files (`skip_existing: true`). A full 800×300 grid is large — expect long runtimes on first generation.

---

## 5. Step 2 — Set exposure and backgrounds

### Ready-made exposures

Use one of the three example scan JSONs. Each sets **`detector.mass_kg`** with **`livetime_days = 365.25`** and **`duty_cycle = 1.0`**, so:

```
exposure_kg_year = mass_kg × livetime_days × duty_cycle / 365.25  →  mass_kg
```

| Config | `detector.mass_kg` | Exposure |
|--------|-------------------|----------|
| `qcdark2_lbc_ne_flatbkg_proj_0p5kgy.json` | 0.5 | 0.5 kg·year |
| `qcdark2_lbc_ne_flatbkg_proj_1p0kgy.json` | 1.0 | 1.0 kg·year |
| `qcdark2_lbc_ne_flatbkg_proj_2p0kgy.json` | 2.0 | 2.0 kg·year |

### Custom exposure

Copy one of the example JSONs and edit:

```json
"_comment_change_exposure": ">>> EDIT HERE: detector.mass_kg (primary), or livetime_days / duty_cycle",
"detector": {
  "mass_kg": 1.0
},
"experiment": {
  "livetime_days": 365.25,
  "duty_cycle": 1.0
},
"run": {
  "outdir": "outputs/my_custom_exposure",
  "label": "my_custom_exposure"
}
```

Also change **`run.outdir`** (and optionally **`run.label`**) so scans do not overwrite each other.

### n_e analysis space

These configs use:

```json
"experiment": {
  "observable_bins": "n_e",
  "roi_bins": [1, 2, 3, 4, 5]
}
```

Anything other than `"pattern"` selects **n_e-space** analysis in the scan app.

### Background model

```json
"run": {
  "background_source": "dc_flat_migration",
  "background_model": "scale"
}
```

- **`dc_flat_migration`** — builds B(n_e) from dark current + flat d.r.u. in `backgrounds`, folded through the response pipeline.
- **`scale`** — likelihood uses B(θ) = θ × B_template with a tau-weighted prior (strength 98, θ ∈ [0.5, 10]).

Default background parameters in the example configs:

```json
"backgrounds": {
  "dark_current": {
    "lambda_e_per_pix_per_year": 0.00041
  },
  "flat_background": {
    "norm_per_kg_year_keV": 1.0,
    "Emin_eV": 2,
    "Emax_eV": 20,
    "nbins": 200
  },
  "timing": {
    "exposure_time_s": 1800
  }
}
```

Signal and background both scale with exposure; the scan runs in **Asimov** mode (`experiment.mode = "asimov"`), so no `data_path` is required.

---

## 6. Step 3 — Run the scan

```bash
build/ccdarksens_scan_dmelectron_pattern configs/examples/qcdark2_lbc_ne_flatbkg_proj_0p5kgy.json
build/ccdarksens_scan_dmelectron_pattern configs/examples/qcdark2_lbc_ne_flatbkg_proj_1p0kgy.json
build/ccdarksens_scan_dmelectron_pattern configs/examples/qcdark2_lbc_ne_flatbkg_proj_2p0kgy.json
```

Outputs:

- `outputs/scan_pattern_proj_qcdark2_fullgrid_0p5kgy/scan_dmelectron_pattern.root`
- `outputs/scan_pattern_proj_qcdark2_fullgrid_1p0kgy/scan_dmelectron_pattern.root`
- `outputs/scan_pattern_proj_qcdark2_fullgrid_2p0kgy/scan_dmelectron_pattern.root`

At startup the scan logs exposure, n_e range, and background template when `verbosity ≥ 1`.

---

## 7. Step 4 — Plot limit curves

Single exposure:

```bash
build/ccdarksens_plot_dmelectron_limit --batch \
  --out-pdf outplots/qcdark2_ne_proj_1kgy.pdf \
  outputs/scan_pattern_proj_qcdark2_fullgrid_1p0kgy/scan_dmelectron_pattern.root \
  "QCDark2 n_e 1 kg·year" 1.642374415149816 heavy
```

Overlay three exposures:

```bash
build/ccdarksens_plot_dmelectron_limit --batch \
  --out-pdf outplots/qcdark2_ne_proj_exposure_compare.pdf \
  outputs/scan_pattern_proj_qcdark2_fullgrid_0p5kgy/scan_dmelectron_pattern.root "0.5 kg·year" \
  outputs/scan_pattern_proj_qcdark2_fullgrid_1p0kgy/scan_dmelectron_pattern.root "1.0 kg·year" \
  outputs/scan_pattern_proj_qcdark2_fullgrid_2p0kgy/scan_dmelectron_pattern.root "2.0 kg·year" \
  1.642374415149816 heavy
```

Use `1.642374415149816` as the 90% CL threshold for pydme-style PLR parity. The last argument (`heavy`) selects literature overlays if you add `--show-damic`.
