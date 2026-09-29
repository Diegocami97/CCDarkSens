# CCDarkSens — Beginner's Guide

A practical introduction for newcomers. After reading this you will understand what the code does, how to choose between the two analysis spaces, and how to run a sensitivity calculation from scratch.

---

## 1. What is CCDarkSens?

CCDarkSens is a C++17 framework for computing **dark matter sensitivity** with the DAMIC-M CCD detector. It takes a DM particle physics model as input and returns a 90% CL upper limit on the DM–electron cross section σ_e as a function of DM mass mχ.

The full chain is:

```
DM rate tables  dR/dE  [events / kg·year·eV]
        ↓  ChargeIonization: fold over P(n_e | E)
  S_true(n_e)  [counts per n_e bin]
        ↓  Detector response: diffusion + readout noise + EfficiencyMC
Observable space  ←─── YOUR FIRST CHOICE (see §2)
        ↓
  Background model  ←─── YOUR SECOND CHOICE (see §5)
        ↓
  Profile likelihood ratio → upper limit on σ_e(mχ)
```

---

## 2. The Two Observable Spaces — Choose One

This is the most important configuration decision. Everything else follows from it.

---

### Space A — n_e bins

The observable is **how many electrons were deposited**: bins at n_e = 1, 2, 3, 4, 5.

```
ε(n_e)   = P(event with n_e electrons passes the selection)   [scalar per bin]
S(n_e)   = dRdE-folded rate × ε(n_e) × exposure_kg_year      [counts]
B(n_e)   = dark current + flat d.r.u., both folded with ε(n_e)[counts]
```

**When to use it:** Simpler analysis, good for quick projections and cross-checks. Does not use the spatial shape of the cluster.

Set in config:
```json
"experiment": { "observable_bins": "n_e",      "roi_bins": [1,2,3,4,5] },
"response":   { "analysis_space": "ne" }
```

The scan app treats anything **other than** `"pattern"` as n_e mode. Some older configs use `"ne"` (no underscore); both work.

---

### Space B — Pattern bins

The observable is **which charge cluster shape was found**: the six patterns `{11}`, `{21}`, `{111}`, `{31}`, `{22}`, `{211}`.

A pattern is a sorted tuple of pixel charges (in electrons) forming an isolated cluster — `{2,1}` means a 2e pixel adjacent to a 1e pixel.

```
P(p | n_e)  = full probability matrix from EfficiencyMC          [matrix]
S_p         = Σ_{n_e} S_true(n_e) × P(p | n_e)                  [counts per pattern]
B_p         = B^rc_p (dark current coincidences) + θ × B^rad_p   [counts per pattern]
```

The pattern shape carries more information than n_e bins: a DM scatter deep in the bulk produces a diffuse cloud spread across pixels, while dark current produces compact single-pixel events. This is the observable used in the published SRDM result.

Set in config:
```json
"experiment": { "observable_bins": "pattern", "pattern_roi": [11,21,111,31,22,211] },
"response":   { "analysis_space": "pattern" }
```

---

### Comparison

| | n_e bins | Pattern bins |
|---|---|---|
| Observables | n_e = 1..5 | {11}, {21}, {111}, {31}, {22}, {211} |
| Efficiency | ε(n_e): one number per bin | P(p\|n_e): matrix from EfficiencyMC |
| Background | DC rate + flat d.r.u., folded | B^rc_p + θ·B^rad_p (see §5) |
| Reference | — | DAMIC-M SRDM paper Table I |

---

### Mediator and rate backend

| Backend | Physics | Typical `mediator` values | Rate directory example |
|---------|---------|---------------------------|------------------------|
| **QEDark** | Free-electron / band-gap model | `heavy`, `massless` (ultra-light) | `data/qedark_rates/Si/heavy/long_scan/` |
| **QCDark2** | Crystal dielectric ε(ω,q) | `heavy`, `light` | `data/qcdark2_rates/Si/heavy/Si_comp_long_scan/` |
| **SRDM CSV** | Pre-tabulated rates | (in CSV path) | project-specific |

Choose the backend in the **rate-generation** JSON (`utils/qedark_generate_grid.py` or `utils/qcdark2_generate_grid.py`) and point `model.rates_dir` in the scan config to the output folder. Filenames follow `dRdE_{material}_{mediator}_m{mchi}_s{sigma}.csv` (units: events/kg/year/eV).

---

## 3. Prerequisites, repo layout, and building

### Software

| Component | Purpose |
|-----------|---------|
| C++17 compiler, CMake ≥ 3.18 | Build CCDarkSens |
| [ROOT](https://root.cern/) (Core, Hist, RIO, Minuit2) | Histograms, I/O, minimisation |
| [nlohmann/json](https://github.com/nlohmann/json) | Config parsing |
| Python 3.8+ with NumPy | Rate generation (`ccdarkphys` package) |

### Repo map

```
apps/               → C++ source for executables (scan, plot, one-point, band, …)
configs/            → JSON scan and rate-generation configs
configs/examples/   → collaboration-ready worked examples + README
data/               → rates, efficiencies, ionization (p100K) tables, reference curves
outputs/            → scan ROOT files (created at run time)
utils/              → rate grid scripts, audits, plotting helpers
python/ccdarkphys/  → QEDark / QCDark2 rate physics
docs/               → deep dives, reproduction guides, audits
build/              → compiled binaries (after cmake)
```

### Build

```bash
cmake -B build -S .
cmake --build build -j8 --target ccdarksens_scan_dmelectron_pattern ccdarksens_plot_dmelectron_limit ccdarksens_example_one_point_pattern
```

All binaries land in `build/`. Add other targets as needed (`ccdarksens_band`, `ccdarksens_scan_srdm_pattern_csv`, …).

---

## 4. Detector Parameters

Every config starts with a `"detector"` block. These numbers fully define the CCD geometry and set the exposure.

```json
"detector": {
  "rows":             1300,       // number of pixel rows in the CCD
  "cols":             6300,       // number of pixel columns
  "pixel_size_um":    15.0,       // pixel pitch in micrometres
  "thickness_mm":     0.67,       // depletion depth (sets diffusion range)
  "active_fraction":  0.98,       // fraction of pixels not masked
  "target_element":   "Si",
  "density_g_cm3":    2.329,
  "mass_kg":          0.5         // fiducial detector mass
}
```

**How exposure is computed:**

```
exposure_kg_year = livetime_days × duty_cycle × mass_kg / 365.25
```

This single number converts every DM rate table (in events/kg/year/eV) into expected counts. It is computed in `ExperimentSetup` and passed to `ChargeIonization::FoldToNe`. It is applied **exactly once** per spectrum; nothing downstream multiplies by exposure again.

### Which mass / exposure for which analysis?

Different configs intentionally use different fiducial masses. **Always check the exposure printed at scan startup** (`[scan-pattern] Exposure = … kg·year`).

| Analysis | Typical `mass_kg` | `livetime_days` | `exposure_kg_year` | Notes |
|----------|-------------------|-----------------|---------------------|-------|
| **SRDM / Table I reference** | 0.5 | 85.6 | ≈ 0.117 | Published LBC dataset geometry |
| **QEDark LBC repro** (`configs/examples/lbc_qedark_*.json`) | 0.01523 | 85.36 | ≈ 0.0036 | Matches pydme export mass; `data_path` may carry its own exposure metadata |
| **n_e projections** (`qcdark2_lbc_ne_flatbkg_proj_*.json`) | 0.5 / 1 / 2 | 365.25 | = `mass_kg` | Asimov sensitivity studies at 1 year livetime |
| **Band-gap one-point dumps** | 0.5 | 365.25 | 0.5 | Quick spectrum comparisons |

If limits disagree with a reference, the first check is whether **exposure_kg_year** matches between your config, the data file, and the comparison curve.

The pixel mass used for background scaling is:
```
mass_per_pixel_kg = pixel_size_um² × thickness_mm × density_g_cm3 × unit_factors
                  ≈ 3.51 × 10⁻¹⁰ kg   (15 µm × 15 µm × 0.67 mm, Si)
```

---

## 5. Background Model

### 5.1 Asimov vs real data

By default the framework runs in **Asimov mode**: the observed counts are set equal to the expected background, `D_p = B_p`. This gives the **median expected sensitivity** — how well the experiment would do on average.

```json
"run": { "data_path": "" }   // empty = Asimov
```

To use real data, provide a path to a CSV or ROOT file containing the observed pattern counts in `pattern_roi` order:
```json
"run": { "data_path": "data/Final_Combined_Image_Data.csv" }
```

### 5.2 Background source — `background_source`

Controls how `B_p` (or `B(n_e)`) is formed:

| Value | Meaning | Typical use |
|---|---|---|
| `"dc_flat_migration"` | Compute B from dark current rate + flat d.r.u., then fold through signal efficiency. | n_e space or projections |
| `"bp_br_template"` | Use pre-computed Bp and Br vectors from config directly. | Pattern space, reproducing published results |

### 5.3 Flat dark current (DC) background

When `background_source = "dc_flat_migration"`, the dark current background in n_e is computed by `BackgroundBuilder`:

```json
"backgrounds": {
  "dark_current": {
    "lambda_e_per_pix_per_year": 0.0365,  // DC rate in e-/pixel/year
    "norm_scale": 1.0                      // optional rescaling
  },
  "timing": {
    "exposure_time_s": 1800,   // duration of one readout (sets n_e pileup)
    "n_exposures": null        // null = infer from livetime_days
  }
}
```

`BackgroundBuilder` converts the per-pixel rate into expected counts per n_e bin using:
```
λ_per_exp = lambda_e_per_pix_per_year × (exposure_time_s / 86400 / 365.25)
B_dc(n_e) = N_pixels_active × Poisson(n_e | λ_per_exp)
```
then optionally folds with `ε(n_e)` if in n_e-space mode.

### 5.4 Flat radiogenic background (d.r.u.)

A flat differential rate spectrum (events/kg/year/keV, i.e. d.r.u.) is treated identically to a signal: it is folded through `ChargeIonization` and then through the detector response pipeline to produce `B_flat(n_e)` or `B_flat_p`.

```json
"backgrounds": {
  "flat_background": {
    "norm_per_kg_year_keV": 1.0,  // flat rate in d.r.u. (events/kg/year/keV)
    "Emin_eV": 0,
    "Emax_eV": 20,
    "nbins":   200
  }
}
```

The total Asimov background is `B_tot = B_dc + B_flat`.

### 5.5 Pattern-space template (Bp/Br)

When `background_source = "bp_br_template"`, the framework skips the DC/flat construction and uses pre-computed vectors directly:

```json
"run": {
  "background_Bp": [141.4, 0.111, 0.042, 0.019, 2.5e-5, 5.8e-5],
  "background_Br": [0.039, 0.039, 0.016, 0.052, 0.011, 0.035]
}
```

These vectors are indexed in the same order as `pattern_roi = [11, 21, 111, 31, 22, 211]`.

- **Bp** = B^rc_p = expected random-coincidence background from dark current. Derived via:
  `B^rc_p = N_total × Σ_q Π Poisson(q_i|λ_i) × B[p|q]`
  where N_total ≈ 1.85 × 10⁹ pixels and B[p|q] is the background confusion matrix.

- **Br** = B^rad_p = expected radiogenic background, derived from external Geant4 simulations. Not recomputed inside CCDarkSens.

### 5.6 Background model — `background_model`

Controls the nuisance parameter structure in the likelihood:

| Value | Meaning |
|---|---|
| `"scale"` | `B(scale) = scale × B_template`. One nuisance scales the whole background. |
| `"Bp_theta_Br"` | `B(θ) = Bp + θ·Br`. Separates fixed random-coincidence from scaled radiogenic. Matches pydme. |

With `"Bp_theta_Br"` and **`constrain_use_tau_weighted: true`** (recommended for pydme parity), the prior is a tau-weighted Gamma form summed over pattern bins — see [`Framework_Architecture.md`](Framework_Architecture.md) §6.1 for the exact NLL construction. Set **`constrain_n_bins: 1`** when using a single pooled constraint. θ is bounded to [`theta_lo`, `theta_hi`] (typically 0.5–10) during minimisation.

With **`background_model: "scale"`** (common in n_e projections), one nuisance θ scales the entire B(n_e) template: B(θ) = θ × B_template.

---

## 6. Signal Efficiency

### From EfficiencyMC (default)

`EfficiencyMC` runs a pixel-level MC to compute `P(pattern p | n_e)`:

1. Sample absorption depth z uniformly in [0, thickness]
2. Compute diffusion spread `σ_xy(z) = √(−A·ln(1 − b·z)) · α` via `ChargeTransport`
3. Distribute n_e electrons as a 2D Gaussian cloud → bin into pixels
4. Add readout noise N(0, σ_ro) to every pixel
5. Run `PatternClassifier` on the 3-row neighbourhood → record which pattern is found
6. Repeat N_trials times → normalised counts = P(p | n_e)

Diffusion parameters come from `response.efficiency_mc` (legacy alias: `pattern_mc` — both are accepted):
```json
"efficiency_mc": {
  "A_um2":          803.25,   // diffusion amplitude (µm²)
  "b_umInv":        0.00065,  // depth factor (µm⁻¹)
  "alpha":          1.0,      // dimensionless scale
  "beta_per_keV":   0.0,      // energy dependence (set 0 for standard Si)
  "sigma_readout_e":0.16,     // readout noise (e-)
  "n_events_per_ne":200000    // MC trials per n_e value
}
```

### From a pre-computed CSV (faster)

If you already have the efficiency table from a previous run or from the collaboration reference, skip the MC entirely:

```json
"efficiency_mc": {
  "efficiency_csv": "data/Efficiencies_patterns_Nsims1000000_DCTrue_alpha1.csv"
}
```

The CSV format is `pattern_id, n_e, efficiency` (one row per (pattern, n_e) pair).

### Pattern classifier thresholds

```json
"pattern_classifier": {
  "Qmin_e":   0.60,   // minimum charge to call a pixel occupied (e-)
  "thr_M":    3.5,    // single-pixel pattern threshold (in units of σ_ro)
  "thr_MN":   4.0,    // two-pixel pattern threshold
  "thr_MNL":  5.5     // three-pixel pattern threshold
}
```

---

## 7. Key Apps

### 7.0 Which app should I use?

| Goal | App | Typical config |
|------|-----|----------------|
| Sanity-check one (mχ, σ) point | `ccdarksens_example_one_point_pattern` | any scan JSON |
| **QEDark / QCDark2 DM-e limits (pattern or n_e)** | **`ccdarksens_scan_dmelectron_pattern`** | `configs/examples/lbc_qedark_heavy_mediator.json` |
| SRDM limits from pre-made CSV rates | `ccdarksens_scan_srdm_pattern_csv` | SRDM-specific configs |
| Sensitivity bands (median ± quantiles) | `ccdarksens_band` | `configs/band_gap_one_point_*.json` |
| Plot limit curves | `ccdarksens_plot_dmelectron_limit` | scan ROOT output |

For current collaboration work (QEDark LBC repro, QCDark2 studies), **`ccdarksens_scan_dmelectron_pattern`** is the primary scan binary.

### 7.1 Single-point diagnostic — `ccdarksens_example_one_point_pattern`

Runs the full pipeline at one (mχ, σ) point. Use this first to verify your config before a full scan.

```bash
build/ccdarksens_example_one_point_pattern configs/band_gap_one_point_gap0p1_D-equal.json
# Or override the mass point:
build/ccdarksens_example_one_point_pattern configs/my_config.json 0.5 1e-35
```

What it prints:
```
S_pat[11] = ...      expected signal counts in pattern {11}
B_pat[11] = ...      expected background (Asimov: ≈ 141.4 for published config)
D_pat[11] = ...      observed (= B in Asimov mode)
CLs upper limit: σ_e < X.XXe-YY cm²
```

It also writes `efficiency_per_pattern.csv` (P(p|n_e) table) and a ROOT file to `outdir/`.

### 7.2 Full grid scan — `ccdarksens_scan_dmelectron_pattern` (primary)

Scans a 2D (mχ, σ) grid and finds the 90% CL upper limit on σₑ at each mass. Works in **pattern** or **n_e** space depending on `experiment.observable_bins`.

```bash
build/ccdarksens_scan_dmelectron_pattern configs/examples/lbc_qedark_heavy_mediator.json
build/ccdarksens_scan_dmelectron_pattern configs/examples/qcdark2_lbc_ne_flatbkg_proj_1p0kgy.json
```

The grid is defined in the JSON (dense-grid configs use `logspace`):
```json
"model": {
  "grid": {
    "mchi_MeV": {
      "logspace": { "start": 0.2, "stop": 1000.0, "num": 800, "endpoint": true }
    },
    "sigma_e_cm2": {
      "logspace": { "start_exp": -46, "stop_exp": -26, "num": 300, "endpoint": true }
    },
    "format": { "mchi": ".6f", "sigma": ".1e" }
  }
}
```

Smaller smoke grids can use explicit lists: `"values": [0.3, 0.5, 1.0]` or `"log10_range": [-42, -35], "nsteps": 50`.

**ROOT output** (in `run.outdir/scan_dmelectron_pattern.root`):

| Object | Role |
|--------|------|
| `q_mchi_sigma_pattern` | 2D histogram of PLR test statistic q(σ) at each grid point |
| `upper_limit_sigma_e_mchi_graph` | **Primary UL curve** — σₑ limit vs mχ from q-grid crossing |
| `upper_limit_sigma_e_mchi` | 1D histogram form of the same limit (legacy) |
| `upper_limit_sigma_e_mchi_pydme_bisection` | Optional pydme-style bisection diagnostic (when enabled) |
| `D_pat` / `D_ne` | Observed counts (Asimov: equals background) |

The plot app reads **`upper_limit_sigma_e_mchi_graph`** by default. Use `--from-qhist` to rebuild the limit from `q_mchi_sigma_pattern` instead.

### 7.2b Alternate — `ccdarksens_scan_srdm_pattern_csv`

Older path for SRDM-style scans that read pre-tabulated CSV rate files directly. Use **`ccdarksens_scan_dmelectron_pattern`** for QEDark/QCDark2 work unless you have a legacy SRDM config.

```bash
build/ccdarksens_scan_srdm_pattern_csv configs/your_srdm_config.json
```

### 7.3 Sensitivity band — `ccdarksens_band`

Wraps the scan to compute median sensitivity bands (±1σ, ±2σ) using toy MC or Asimov. Used for projection plots.

```bash
build/ccdarksens_band configs/band_gap_one_point_gap0p1_D-equal.json
```

See `apps/README_band.md` for options. Output: ROOT file with median + quantile curves.

### 7.4 Limit plot — `ccdarksens_plot_dmelectron_limit`

Takes the ROOT output from the scan or band and produces an exclusion plot.

```bash
build/ccdarksens_plot_dmelectron_limit --batch \
  --out-pdf outplots/my_limit.pdf \
  outputs/lbc_qedark_heavy/scan_dmelectron_pattern.root \
  "QEDark heavy" 1.642374415149816 heavy
```

Useful flags:

| Flag | Purpose |
|------|---------|
| `--batch` | Non-interactive PDF output |
| `--show-damic` | Overlay canonical pydme paper export (ScienceRun2024) |
| `--draw-both` | Draw stored UL and pydme bisection / q-map diagnostic |
| `--from-qhist` | Rebuild UL from `q_mchi_sigma_pattern` instead of stored graph |

The threshold **`1.642374415149816`** is the 90% CL q-value for pydme-style PLR parity. Last argument (`heavy` / `light`) selects literature curve family when using `--show-damic`.

Audit scripts for reproduction checks:

```bash
python3 utils/audit_qedark_ul_vs_pydme.py
python3 utils/plot_qedark_all_references.py
```

See [`Framework_Architecture.md`](Framework_Architecture.md) §6 for the UL-crossing methodology these scripts check against.

### 7.5 Rate table generation (step 0)

Before any scan, generate **dR/dE** CSVs for your (mχ, σ) grid.

**QEDark** (example: heavy mediator, dense grid):

```bash
python3 utils/qedark_generate_grid.py configs/examples/qedark_generate_heavy_dense.json
```

Output: `data/qedark_rates/Si/heavy/long_scan/dRdE_Si_heavy_m{mchi}_s{sigma}.csv`

**QCDark2** (example: Si_comp dielectric):

```bash
python3 utils/qcdark2_generate_grid.py configs/examples/qcdark2_generate_si_comp_dense.json
```

Output: `data/qcdark2_rates/Si/heavy/Si_comp_long_scan/`

Set `model.rates_dir` and `filename_template` in the scan JSON to match. Generation skips existing files when `skip_existing: true`. See worked guides in `configs/examples/README.md`.

---

## 8. The JSON Config — Full Reference

### `run` block

```json
"run": {
  "label":                    "my_run",
  "outdir":                   "outputs/my_run",
  "cl":                       0.9,
  "test_stat":                "PLR",
  "background_source":        "bp_br_template",
  "background_model":         "Bp_theta_Br",
  "background_Bp":            [141.4, ...],
  "background_Br":            [0.039, ...],
  "observed_counts":          [144, 0, 0, 1, 0, 0],
  "use_profile_likelihood":   true,
  "profile_minimizer":        "pydme",
  "constrain_prior_strength": 98,
  "constrain_use_tau_weighted": true,
  "constrain_use_gamma_sign": false,
  "constrain_n_bins":         1,
  "theta_lo":                 0.5,
  "theta_hi":                 10.0,
  "pydme_style_ul":           true,
  "smooth_ul_envelope":       false,
  "data_path":                ""
}
```

`observed_counts` overrides `data_path` when non-empty (one entry per pattern in `pattern_roi` order). For n_e Asimov projections, leave both empty.

### `detector` block

```json
"detector": {
  "rows":             1300,    // CCD rows
  "cols":             6300,    // CCD columns
  "pixel_size_um":    15.0,    // pixel pitch (µm)
  "thickness_mm":     0.67,    // depletion depth (mm) — sets diffusion range
  "active_fraction":  0.98,    // fraction of pixels not masked by quality cuts
  "target_element":   "Si",
  "density_g_cm3":    2.329,
  "mass_kg":          0.5      // fiducial mass — sets exposure with livetime_days
}
```

### `experiment` block

```json
"experiment": {
  "mode":            "asimov",   // always "asimov" for now
  "livetime_days":   365.25,     // exposure time → exposure_kg_year = livetime × duty × mass / 365.25
  "duty_cycle":      1.0,
  "observable_bins": "pattern",  // "pattern" | "n_e"  ← KEY CHOICE
  "pattern_roi":     [11, 21, 111, 31, 22, 211],  // patterns used in the fit
  "roi_bins":        [1, 2, 3, 4, 5],              // used in ne-space mode
  "binning": { "ne_min": 1, "ne_max": 5 }
}
```

### `response` block

```json
"response": {
  "mode":            "pattern",
  "analysis_space":  "pattern",   // "pattern" | "ne"  ← KEY CHOICE
  "efficiency_mc": {
    "efficiency_csv":    "",       // leave empty to use EfficiencyMC; or path to pre-computed CSV
    "n_events_per_ne":   200000,
    "sigma_readout_e":   0.16,
    "A_um2":             803.25,
    "b_umInv":           0.00065,
    "alpha":             1.0,
    "beta_per_keV":      0.0
  },
  "pattern_classifier": {
    "Qmin_e":  0.60,
    "thr_M":   3.5,
    "thr_MN":  4.0,
    "thr_MNL": 5.5
  },
  "charge_ionization": {
    "table_csv":   "data/p100K_gap0p1_eh0p1.csv",  // P(n_e | E) table
    "band_gap_eV": 0.1,
    "eh_pair_eV":  0.1
  }
}
```

### `response.cluster_fit_mc` block (WIMP-nucleon channel only — `analysis_space: "cluster_energy"`)

Not used by the DM-electron/dark-photon/Migdal channels above; this is the WIMP-nucleus SI channel's continuous-energy reconstruction (see `docs/ClusterFitMC_Design.md` for the full derivation and validation).

```json
"response": {
  "analysis_space":  "cluster_energy",
  "cluster_fit_mc": {
    "window_nx": 15, "window_ny": 15,           // pixel window size for the ΔLL fit
    "pixel_size_um": 15.0,
    "sigma_readout_e": 0.16,                     // switch to 1.8 for a literal PhysRevD.94.082006 reproduction
    "diffusion": {                               // same convention as efficiency_mc's diffusion block
      "A_um2": 803.25, "b_umInv": 0.00065, "alpha": 1.0, "beta_per_keV": 0.0, "thickness_um": 670.0
    },
    "fit_method": "nelder_mead",                 // "nelder_mead" | "minuit2"
    "sigma_xy_lo_px": 0.1, "sigma_xy_hi_px": 2.0, // search bounds on the fitted width
    "n_toys": 100000,                            // Phase 3: noise-tail calibration toy count
    "target_tail_prob": 1e-3,
    "ne_trials_per_point": 5000,                 // Phase 4: forward-sim trials per E_true grid point
    "sigma_xy_fid_min_px": 0.35, "sigma_xy_fid_max_px": 1.22,  // paper's fiducial (surface-rejection) cut
    "eh_pair_eV": 3.77,                          // PhysRevD.94.082006's stated value
    "fano_factor": 0.133,                        // the paper's own measured value; treat as a systematic to sweep
    "Etrue_min_eV": 10.0, "Etrue_max_eV": 300.0, "Etrue_npoints": 12,
    "Ereco_min_eV": 0.0, "Ereco_max_eV": 400.0, "Ereco_nbins": 40
  }
}
```

### `backgrounds` block (for `dc_flat_migration` source)

```json
"backgrounds": {
  "dark_current": {
    "lambda_e_per_pix_per_year": 0.0365,  // DC rate in e-/pixel/year
    "norm_scale": 1.0
  },
  "timing": {
    "exposure_time_s": 1800,   // single-readout duration (for pileup)
    "n_exposures":     null    // null = infer from livetime_days / exposure_time_s
  },
  "flat_background": {
    "norm_per_kg_year_keV": 1.0,  // flat radiogenic rate in d.r.u. (events/kg/year/keV)
    "Emin_eV": 0,
    "Emax_eV": 20,
    "nbins":   200
  }
}
```

### `model` block

```json
"model": {
  "type":              "dm_electron",
  "material":          "Si",
  "mediator":          "heavy",              // "heavy" | "light" | "ultralight"
  "rates_dir":         "data/qcdark2_rates/Si/heavy/...",
  "filename_template": "dRdE_{material}_{mediator}_m{mchi_MeV}_s{sigma_e_cm2}.csv",
  "Emin_eV": 0, "Emax_eV": 20, "nbins": 200,
  "grid": {
    "mchi_MeV": {
      "logspace": { "start": 0.2, "stop": 1000.0, "num": 800, "endpoint": true }
    },
    "sigma_e_cm2": {
      "logspace": { "start_exp": -46, "stop_exp": -26, "num": 300, "endpoint": true }
    },
    "format": { "mchi": ".6f", "sigma": ".1e" }
  }
}
```

Rate tables must be in **events/kg/year/eV**. The filename template uses `{key}` substitution; `mchi_MeV` is formatted according to `format.mchi`. Legacy configs may use `"values": [...]` or `"log10_range"` instead of `"logspace"`.

---

## 9. Reproducing published results

### 9.1 SRDM LBC (Table I, image CSV)

Use the SRDM config with the published Bp/Br values and the actual image data:

```json
"run": {
  "background_source": "bp_br_template",
  "background_model":  "Bp_theta_Br",
  "background_Bp": [141.4, 0.111, 0.042, 0.019, 2.5e-5, 5.8e-5],
  "background_Br": [0.039, 0.039, 0.016, 0.052, 0.011, 0.035],
  "data_path": "data/Final_Combined_Image_Data.csv",
  "use_profile_likelihood": true,
  "constrain_prior_strength": 98,
  "pydme_style_ul": true
},
"experiment": {
  "observable_bins": "pattern",
  "pattern_roi": [11, 21, 111, 31, 22, 211],
  "livetime_days": 85.57
}
```

Cross-check against Table I:

| Pattern | D_p (observed) | Bp (B^rc) | Br (B^rad) |
|---------|---------------|-----------|------------|
| {11}    | 144           | 141.4     | 0.039      |
| {21}    | 1             | 0.111     | 0.039      |
| {111}   | 0             | 0.042     | 0.016      |
| {31}    | 0             | 0.019     | 0.052      |
| {22}    | 0             | 2.5×10⁻⁵  | 0.011      |
| {211}   | 0             | 5.8×10⁻⁵  | 0.035      |

Data file: `data/Final_Combined_Image_Data.csv`. This path reproduces the **SRDM pattern analysis** from the collaboration paper tables.

### 9.2 QEDark LBC (2025 pydme reproduction)

Separate from §9.1 — uses **QEDark rate tables**, pydme-style likelihood, and the pydme pattern export:

| Item | Value |
|------|-------|
| Data | `build/data_pattern.root` or `run.observed_counts` |
| Rates | `data/qedark_rates/Si/heavy/long_scan/` |
| Configs | `configs/examples/lbc_qedark_heavy_mediator.json`, `lbc_qedark_light_mediator.json` |
| Paper reference | `collab_frameworks/pydme/figures/ScienceRun2024-figures/data/ScienceRun2024_results-Pattern/` |

Full walkthrough: [`LBC_QEDark_Reproduction_Guide.md`](LBC_QEDark_Reproduction_Guide.md).

---

## 10. Typical Workflow

> **Worked examples** (configs in `configs/examples/`):
> - LBC pydme + QEDark reproduction: [`LBC_QEDark_Reproduction_Guide.md`](LBC_QEDark_Reproduction_Guide.md)
> - QCDark2 pattern-count sensitivity: [`QCDark2_Pattern_Counts_Guide.md`](QCDark2_Pattern_Counts_Guide.md)
> - QCDark2 n_e exposure projections (flat DC + d.r.u.): [`QCDark2_NE_Exposure_Projections_Guide.md`](QCDark2_NE_Exposure_Projections_Guide.md)

1. **Generate rate tables** — `python3 utils/qedark_generate_grid.py` or `qcdark2_generate_grid.py` (see §7.5). Point `model.rates_dir` at the output folder.

2. **Choose observable space** — n_e bins for quick projections; pattern bins for full LBC-style analysis (§2).

3. **Copy and edit a config** — start from `configs/examples/` or the closest config in `configs/`. Update `rates_dir`, `grid`, `livetime_days`, `mass_kg`, and `background_source`.

4. **Sanity-check one point** — `ccdarksens_example_one_point_pattern your_config.json`. Verify S, B, D, and the preliminary limit.

5. **Scan the grid** — `ccdarksens_scan_dmelectron_pattern your_config.json`.

6. **Plot and validate** — `ccdarksens_plot_dmelectron_limit` with `--show-damic` for paper comparison; audit scripts in `utils/` for quantitative checks.

**QEDark heavy end-to-end** (copy-paste):

```bash
python3 utils/qedark_generate_grid.py configs/examples/qedark_generate_heavy_dense.json
build/ccdarksens_scan_dmelectron_pattern configs/examples/lbc_qedark_heavy_mediator.json
build/ccdarksens_plot_dmelectron_limit --batch --show-damic \
  --out-pdf outplots/lbc_qedark_heavy.pdf \
  outputs/lbc_qedark_heavy/scan_dmelectron_pattern.root \
  "QEDark heavy" 1.642374415149816 heavy
```

---

## 11. Troubleshooting checklist

If your limit looks wrong, check in this order:

1. **Exposure** — does the scan log `Exposure = X kg·year` match your data file and reference curve? (§4 table)
2. **Bin order** — `pattern_roi`, `background_Bp`, `background_Br`, and `observed_counts` must use the same order: `[11, 21, 111, 31, 22, 211]`.
3. **Rate units** — tables must be events/kg/year/eV; wrong units shift limits by orders of magnitude.
4. **Rate path** — `model.rates_dir` + `filename_template` must match generated CSV names on disk.
5. **Ionization table** — `response.charge_ionization.table_csv` and `band_gap_eV` / `eh_pair_eV` must match your pheno scenario.
6. **θ constraint** — for pydme parity use `constrain_use_tau_weighted: true`, `constrain_n_bins: 1`, `profile_minimizer: "pydme"`. θ pinned at `theta_lo` usually means the prior settings are wrong.
7. **One-point cross-check** — run `ccdarksens_example_one_point_pattern` at the same (mχ, σ) and compare S/B/D before debugging the full grid.
8. **Plot source** — default plot reads `upper_limit_sigma_e_mchi_graph`; use `--from-qhist` if you suspect the stored graph is stale.

---

## 12. Further reading

| Topic | Document |
|-------|----------|
| Full pipeline architecture (config → response → stats, module by module) | [`Framework_Architecture.md`](Framework_Architecture.md) |
| QEDark LBC reproduction (step-by-step) | [`LBC_QEDark_Reproduction_Guide.md`](LBC_QEDark_Reproduction_Guide.md) |
| QCDark2 pattern-count sensitivity | [`QCDark2_Pattern_Counts_Guide.md`](QCDark2_Pattern_Counts_Guide.md) |
| QCDark2 n_e exposure projections | [`QCDark2_NE_Exposure_Projections_Guide.md`](QCDark2_NE_Exposure_Projections_Guide.md) |
| UL methodology, likelihood/θ minimization | [`Framework_Architecture.md`](Framework_Architecture.md) §6 |
| WIMP-nucleon cluster-fit design | [`ClusterFitMC_Design.md`](ClusterFitMC_Design.md) |
| Example config index | [`configs/examples/README.md`](../configs/examples/README.md) |

---

## 13. Glossary

| Term | Meaning |
|------|---------|
| n_e | Number of electrons deposited by a DM scatter |
| pattern | Sorted pixel charge tuple identifying a cluster shape, e.g. {2,1} |
| ε(n_e) | Signal efficiency in n_e space: P(event passes selection \| n_e) |
| P(p\|n_e) | Signal efficiency in pattern space: probability of identifying pattern p |
| σ_ro | Readout noise (e-), typically 0.16 e- |
| λ_dc | Dark current rate in e-/pixel/year |
| d.r.u. | Differential rate unit: events/kg/year/keV (flat radiogenic background) |
| exposure_kg_year | livetime_days × duty_cycle × mass_kg / 365.25 |
| B^rc_p | Random-coincidence background in pattern p (from dark current) |
| B^rad_p | Radiogenic background in pattern p (from Geant4 simulations, hardcoded as Br) |
| Bp, Br | Total expected B^rc and B^rad over the full dataset |
| B[p\|q] | Confusion matrix: P(identified as p \| injected cluster q) |
| θ | Nuisance parameter scaling Br in the likelihood; profiled with prior strength 98 |
| Asimov | Pseudo-dataset where D_p = B_p; gives median expected sensitivity |
| PLR | Profile Likelihood Ratio test statistic q(σ) = 2(NLL(σ) − NLL_min); configs use `test_stat: "PLR"` |
| q threshold | 1.642… for 90% CL upper limit (pydme parity); passed to plot app |
| CLs | Modified frequentist CL: CL_s = CL_{s+b}/CL_b (toy MC path; PLR grid scan is default) |
