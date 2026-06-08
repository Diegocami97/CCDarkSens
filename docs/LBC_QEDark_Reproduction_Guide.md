# LBC QEDark Reproduction Guide

A step-by-step walkthrough for collaboration members who want to reproduce the **DAMIC-M 2025 LBC pattern-space DM-electron limits** using **QEDark** rate tables and the **pydme-style** profile likelihood in CCDarkSens.

This guide assumes you have read the overview in [`Beginners_Guide.md`](Beginners_Guide.md) and focuses on one concrete analysis path: **pattern bins + Bp/Br background + QEDark rates**.

---

## 1. What you will reproduce

| Item | Value |
|------|-------|
| Observable | Six pattern bins: `{11}`, `{21}`, `{111}`, `{31}`, `{22}`, `{211}` |
| Signal model | QEDark DM-electron rates (Si, band gap 1.2 eV, Eh 3.8 eV) |
| Background | `B = Bp + θ·Br` with tau-weighted Gamma prior (strength 98) |
| Data | Observed pattern counts from pydme export (`build/data_pattern.root`) |
| Statistic | Profile likelihood ratio (PLR), 90% CL upper limit on σₑ(mχ) |
| Mediators | **Heavy** and **ultra-light (massless)** — same grid, different rate tables |

**Example configs** (dense grid, ready to run):

| Mediator | Rate generation | Scan |
|----------|-----------------|------|
| Heavy | [`configs/examples/qedark_generate_heavy_dense.json`](../configs/examples/qedark_generate_heavy_dense.json) | [`configs/examples/lbc_qedark_heavy_mediator.json`](../configs/examples/lbc_qedark_heavy_mediator.json) |
| Ultra-light | [`configs/examples/qedark_generate_light_dense.json`](../configs/examples/qedark_generate_light_dense.json) | [`configs/examples/lbc_qedark_light_mediator.json`](../configs/examples/lbc_qedark_light_mediator.json) |

Both scans use the same **dense (mχ, σₑ) grid**:

- **800** DM masses, log-spaced from 0.2 to 1000 MeV
- **300** cross sections, log-spaced from 10⁻⁴⁶ to 10⁻²⁶ cm²

At each mass point the scan finds the 90% CL upper limit by profiling θ and crossing the q(σ) curve.

---

## 2. Prerequisites

### Software

- C++17 compiler, CMake ≥ 3.18
- [ROOT](https://root.cern/) (Core, Hist, RIO, Minuit2)
- [nlohmann/json](https://github.com/nlohmann/json)
- Python 3.8+ with NumPy (for rate generation)

### Build CCDarkSens

From the repository root:

```bash
cmake -B build -S .
cmake --build build -j8 --target ccdarksens_scan_dmelectron_pattern ccdarksens_plot_dmelectron_limit
```

Binaries are written to `build/`.

### Input data files

These must exist before running a scan:

| File | Purpose |
|------|---------|
| `build/data_pattern.root` | Observed pattern counts + exposure metadata (pydme LBC export) |
| `data/Efficiencies_patterns_Nsims1000000_DCTrue_alpha1.csv` | Precomputed P(pattern \| nₑ) efficiencies |

If `build/data_pattern.root` is missing, obtain it from the pydme LBC pattern analysis export or ask the analysis contact.

---

## 3. Pipeline overview

```
Step 1   Generate QEDark rate CSVs     (Python, one-time per halo/grid)
           ↓
Step 2   Run pattern-space scan       (C++, finds UL at each mχ)
           ↓
Step 3   Plot limit curve             (C++, overlay pydme reference)
```

---

## 4. Step 1 — Generate QEDark rate tables

Rate tables are differential spectra **dR/dE** in units of events / (kg·year·eV), one CSV per (mχ, σₑ) point on the grid.

### Heavy mediator

```bash
python3 utils/qedark_generate_grid.py configs/examples/qedark_generate_heavy_dense.json
```

Output directory: `data/qedark_rates/Si/heavy/long_scan/`

Filename pattern: `dRdE_Si_heavy_m{mchi}_s{sigma}.csv`

### Ultra-light (massless) mediator

```bash
python3 utils/qedark_generate_grid.py configs/examples/qedark_generate_light_dense.json
```

Output directory: `data/qedark_rates/Si/ultralight/long_scan/`

Filename pattern: `dRdE_Si_massless_m{mchi}_s{sigma}.csv`

### Notes

- **Grid size:** 800 × 300 = **240 000** CSV files per mediator. Generation is parallelized (`options.parallel` in the JSON). With `skip_existing: true`, re-running skips files already on disk.
- **Halo:** Default is Standard Halo Model with v_E = 263 km/s. Change `halo` in the rate-gen JSON if you need a different astrophysics setup.
- **Detector response in rates:** Band gap 1.2 eV, Eh pair 3.8 eV, 0.1 eV binning — matching the QEDark notebook defaults.
- **Verify:** After generation, check that files exist, e.g.:

  ```bash
  ls data/qedark_rates/Si/heavy/long_scan/dRdE_Si_heavy_m1.000000_s1.0e-30.csv
  ls data/qedark_rates/Si/ultralight/long_scan/dRdE_Si_massless_m1.000000_s1.0e-30.csv
  ```

The scan config `model.rates_dir` and `model.grid` must match the rate-generation JSON. The example configs are already aligned.

---

## 5. Step 2 — Run the scan and upper-limit extraction

The scan app folds each rate table through ionization and pattern efficiencies, builds the pydme likelihood, profiles the background scale θ at each σₑ trial, and stores q(mχ, σₑ) plus the 90% CL upper limit.

### Heavy mediator

```bash
build/ccdarksens_scan_dmelectron_pattern configs/examples/lbc_qedark_heavy_mediator.json
```

### Ultra-light mediator

```bash
build/ccdarksens_scan_dmelectron_pattern configs/examples/lbc_qedark_light_mediator.json
```

### Outputs

| Path | Content |
|------|---------|
| `outputs/lbc_qedark_heavy/scan_dmelectron_pattern.root` | Heavy-mediator results |
| `outputs/lbc_qedark_light/scan_dmelectron_pattern.root` | Ultra-light results |

Key objects inside each ROOT file:

| Object | Description |
|--------|-------------|
| `q_mchi_sigma_pattern` | TH2D of test statistic q vs (mχ, σₑ) |
| `upper_limit_sigma_e_mchi_graph` | **Primary** 90% CL UL curve (TGraph) |
| `upper_limit_sigma_e_mchi_pydme_bisection` | Optional pydme-style bisection diagnostic |
| `D_pat` | Observed pattern counts |
| `exposure_kg_year` | Exposure used in the likelihood |

### What the scan does at each mass point

1. Load dR/dE CSV for all σₑ on the grid (signal scales linearly with σₑ).
2. Fold through charge ionization → S_true(nₑ).
3. Apply pattern efficiencies ε(pattern \| nₑ) from the CSV (or EfficiencyMC if no CSV).
4. Build background `B = Bp + θ·Br` per pattern bin.
5. Load observed counts from `data_path` (or use Asimov: data = B).
6. For each σₑ: minimize θ, compute q(σ) = 2ΔlnL.
7. Find σ_UL where q crosses **1.642** (90% CL, 1 dof).

### Runtime expectation

A full **800-mass** scan is compute-intensive (profiling + 300 σₑ points per mass). Expect **many hours** on a single machine. Use `run.verbosity` and the log in `outputs/.../run.log` to monitor progress. For a quick sanity check, temporarily reduce `model.grid.mchi_MeV.num` to e.g. 5 in a copy of the config.

### Key config fields (pydme matching)

```json
"run": {
  "use_profile_likelihood": true,
  "data_path": "build/data_pattern.root",
  "background_source": "bp_br_template",
  "background_model": "Bp_theta_Br",
  "constrain_prior_strength": 98,
  "constrain_use_tau_weighted": true,
  "profile_minimizer": "pydme",
  "pydme_style_ul": true
}
```

---

## 6. Step 3 — Plot the limit curves

Use `ccdarksens_plot_dmelectron_limit` to draw the stored UL and optionally overlay the canonical pydme paper export.

### Heavy mediator

```bash
build/ccdarksens_plot_dmelectron_limit --batch --draw-both --show-damic \
  --out-pdf outplots/lbc_qedark_heavy.pdf \
  --out-root outplots/lbc_qedark_heavy.root \
  outputs/lbc_qedark_heavy/scan_dmelectron_pattern.root \
  "LBC QEDark heavy" 1.642374415149816 heavy
```

### Ultra-light mediator

```bash
build/ccdarksens_plot_dmelectron_limit --batch --draw-both --show-damic \
  --out-pdf outplots/lbc_qedark_light.pdf \
  --out-root outplots/lbc_qedark_light.root \
  outputs/lbc_qedark_light/scan_dmelectron_pattern.root \
  "LBC QEDark light" 1.642374415149816 light
```

### Plot flags explained

| Flag | Effect |
|------|--------|
| `--batch` | Save and exit (no interactive ROOT window) |
| `--draw-both` | Draw primary UL (`upper_limit_sigma_e_mchi_graph`) and pydme bisection diagnostic (if present) |
| `--show-damic` | Overlay canonical pydme paper-export reference curve |
| `--out-pdf` | Output plot path |
| `1.642374415149816` | 90% CL threshold (optional; default is 2.71 for a different convention — use this value for pydme parity) |
| `heavy` / `light` | Selects which pydme reference file to load |

### Compare both mediators on one plot

```bash
build/ccdarksens_plot_dmelectron_limit --batch --show-damic \
  --out-pdf outplots/lbc_qedark_both.pdf \
  outputs/lbc_qedark_heavy/scan_dmelectron_pattern.root "heavy" \
  outputs/lbc_qedark_light/scan_dmelectron_pattern.root "light" \
  1.642374415149816 heavy
```

(Use `heavy` or `light` as the last argument to pick the reference overlay; both scan curves are still drawn.)

### Export CSV

```bash
build/ccdarksens_plot_dmelectron_limit --batch \
  --out-csv outplots/lbc_qedark_heavy.csv \
  outputs/lbc_qedark_heavy/scan_dmelectron_pattern.root \
  "LBC QEDark heavy" 1.642374415149816 heavy
```

