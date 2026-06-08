# QCDark2 Pattern-Count Sensitivity Guide

A step-by-step walkthrough for collaboration members who want to run **QCDark2** DM-electron limits in **pattern space** and study how the result changes when you vary the **observed count in one or more pattern bins**.

This complements the QEDark LBC reproduction guide ([`LBC_QEDark_Reproduction_Guide.md`](LBC_QEDark_Reproduction_Guide.md)). Here the focus is **QCDark2 crystal dielectric rates** and **toy / what-if data vectors** set directly in the scan JSON.

---

## 1. What this workflow does

| Item | Value |
|------|-------|
| Signal model | QCDark2 (Si dielectric ε(ω,q), scissor gap 1.2 eV) |
| Observable | Six pattern bins: `{11}`, `{21}`, `{111}`, `{31}`, `{22}`, `{211}` |
| Background | `B = Bp + θ·Br` with tau-weighted constraint (strength 98) |
| Data | **`run.observed_counts`** in the JSON (overrides `data_path`) |
| Statistic | Profile likelihood, 90% CL upper limit on σₑ(mχ) |
| Mediators | **Heavy** and **light** — same dense grid, different rate tables |

**Example configs** (`configs/examples/`):

| Mediator | Rate generation | Scan |
|----------|-----------------|------|
| Heavy | `qcdark2_generate_heavy_dense.json` | `qcdark2_pattern_counts_heavy_mediator.json` |
| Light | `qcdark2_generate_light_dense.json` | `qcdark2_pattern_counts_light_mediator.json` |

**Dense grid (both mediators):** 800 masses (0.2–1000 MeV) × 300 cross sections (10⁻⁴⁶–10⁻²⁶ cm²).

---

## 2. Prerequisites

### Software

- C++17, CMake ≥ 3.18, ROOT (with Minuit2), nlohmann/json
- Python 3.8+ with NumPy and the `ccdarkphys` package (for QCDark2 rate generation)
- QCDark2 dielectric HDF5: `data/qcdark2_epsilon/Si/Si_fast_gap1p2.h5` (generate via [`qcdark2_dielectric_workflow.md`](qcdark2_dielectric_workflow.md) if missing)

### Build

```bash
cmake -B build -S .
cmake --build build -j8 --target ccdarksens_scan_dmelectron_pattern ccdarksens_plot_dmelectron_limit
```

### Other inputs

| File | Purpose |
|------|---------|
| `data/Efficiencies_patterns_Nsims1000000_DCTrue_alpha1.csv` | P(pattern \| nₑ) efficiencies |
| `data/qcdark2_epsilon/Si/Si_fast_gap1p2.h5` | Dielectric function for rate generation |

---

## 3. Pipeline overview

```
Step 1   Generate QCDark2 rate CSVs     (Python)
           ↓
Step 2   Edit observed_counts in JSON   (one value per pattern bin)
           ↓
Step 3   Run pattern-space scan         (C++)
           ↓
Step 4   Plot and compare curves        (C++)
```

---

## 4. Step 1 — Generate QCDark2 rate tables

Rates are **dR/dE** CSVs (events / kg / year / eV), one file per (mχ, σₑ) on the grid.

### Heavy mediator

```bash
python3 utils/qcdark2_generate_grid.py configs/examples/qcdark2_generate_heavy_dense.json
```

Output: `data/qcdark2_rates/Si/heavy/examples_dense_long_scan/`

### Light mediator

```bash
python3 utils/qcdark2_generate_grid.py configs/examples/qcdark2_generate_light_dense.json
```

Output: `data/qcdark2_rates/Si/light/examples_dense_long_scan/`

Generation is parallelized and skips existing files (`skip_existing: true`). A full 800×300 grid is large — expect long runtimes on first generation.

---

## 5. Step 2 — Set observed pattern counts

Open the scan JSON for your mediator. The fields you edit are at the top of the `run` block:

```json
"_comment_change_counts": ">>> EDIT HERE: run.observed_counts ...",
"run": {
  "data_path": "",
  "observed_counts": [144, 0, 0, 1, 0, 0],
  ...
}
```

### Pattern bin order

Counts must follow **`experiment.pattern_roi`** order:

| Index | Pattern | LBC nominal count |
|-------|---------|-------------------|
| 0 | `{11}` | 144 |
| 1 | `{21}` | 0 |
| 2 | `{111}` | 0 |
| 3 | `{31}` | 1 |
| 4 | `{22}` | 0 |
| 5 | `{211}` | 0 |

### Examples

- Add 10 events in `{11}`: `[154, 0, 0, 1, 0, 0]`
- Remove the single `{31}` event: `[144, 0, 0, 0, 0, 0]`
- Inject 3 events in `{211}`: `[144, 0, 0, 1, 0, 3]`

### Rules

- **`observed_counts` overrides `data_path`** when the array is non-empty.
- Leave **`data_path` empty** (`""`) when using in-config counts (recommended for this study).
- Array length must equal `pattern_roi` length (6).
- For multiple scenarios, **copy the JSON**, change `observed_counts` and `run.outdir`, then rerun.

---

## 6. Step 3 — Run the scan

```bash
build/ccdarksens_scan_dmelectron_pattern configs/examples/qcdark2_pattern_counts_heavy_mediator.json
build/ccdarksens_scan_dmelectron_pattern configs/examples/qcdark2_pattern_counts_light_mediator.json
```

Outputs:

- `outputs/qcdark2_pattern_counts_heavy/scan_dmelectron_pattern.root`
- `outputs/qcdark2_pattern_counts_light/scan_dmelectron_pattern.root`

The scan logs each pattern and its D value at startup when `verbosity ≥ 1`.

---

## 7. Step 4 — Plot limit curves

Compare baseline vs perturbed counts by plotting multiple ROOT files:

```bash
build/ccdarksens_plot_dmelectron_limit --batch \
  --out-pdf outplots/qcdark2_counts_heavy_baseline.pdf \
  outputs/qcdark2_pattern_counts_heavy/scan_dmelectron_pattern.root \
  "QCDark2 heavy baseline" 1.642374415149816 heavy
```

Overlay two scenarios (e.g. baseline vs perturbed `{11}`):

```bash
build/ccdarksens_plot_dmelectron_limit --batch \
  --out-pdf outplots/qcdark2_counts_compare.pdf \
  outputs/qcdark2_pattern_counts_heavy/scan_dmelectron_pattern.root "baseline" \
  outputs/qcdark2_pattern_counts_heavy_pat11_154/scan_dmelectron_pattern.root "pat11=154" \
  1.642374415149816 heavy
```

Use `1.642374415149816` as the 90% CL threshold for pydme-style PLR parity. The last argument (`heavy` or `light`) selects literature overlays if you add `--show-damic`.
