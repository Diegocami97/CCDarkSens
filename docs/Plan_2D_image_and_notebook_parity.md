# Plan: 2D binned image and notebook parity

Goal: add 2D binned image generation (like the notebook) and the remaining unimplemented features, building on current CCDarkSens (ChargeTransport, PixelSimulator, PatternClassifier, PatternMC).

---

## Notebook parity: what’s done vs remaining

**Already in C++ (matches notebook):** Diffusion (σ_xy from A, b, z, E), readout noise, same Pm/Pmn/Pmnl and thresholds (thr_M, thr_MN, thr_MNL), PatternClassifier `ScanRow` (1D), PatternMC `BuildPatternTable` → P(pattern|n_e), efficiency CSV write, `response.pattern_image` config (nrows_binned, ncols, row_binning, etc.).

**Not yet in C++ (to match notebook):**

| # | Notebook feature | C++ equivalent (plan) | Effort |
|---|-------------------|------------------------|--------|
| 1 | `generate_image_E(ne, Ee)` → 2D (3×50) binned image | PatternImageGenerator | Medium |
| 2 | `plot_image(image)` | TH2D + SaveAs PNG | Small |
| 3 | `scan_image(image)` (2D → pattern + charge) | Middle row → ScanRow (then optional full 2D + isolation) | Small |
| 4 | `Ploting_simulation_specific_E` (charge histos per n_e) | Charge histos per n_e → ROOT/PNG | Small |
| 5 | `sanity_checks` (q1 vs q2 etc.) | CSV/TTree “sanity” output | Small |
| 6 | `simulate_cluster(b,c,d)` | simulate_cluster(b,c,d) + optional classify | Small |
| 7 | `pattern_simulation_background` + Background_efficiencies.csv | New app or mode + CSV | Medium |
| 8 | `Generate_efficiencies` with E grid (ε vs E) | E-dependent ε CSV (total then per-pattern) | Small–Medium |
| 9 | Dark current in pattern simulation | **Done:** PatternMC `include_dc_pileup` + config | Small |

**Status:** All items above are now implemented (PatternImageGenerator, example app 2D/plot/scan/charge/sanity/simulate_cluster, ccdarksens_pattern_background_eff, epsilon_vs_E CSV, dark current). **Rough total (was):** **2 Medium + 7 Small** (or 2 Med + 6 Small + 1 Small–Medium). In person-weeks: about **1–2 weeks** for an experienced developer (assuming Phases 1–3 first, then 4–9 as needed).

---

## Current state (brief)

- **ChargeTransport**: `SampleDepthUm()`, `SigmaXYUm(z, E)`, `SampleCloudXY(x0, y0, sigma_xy_um, n_e, xs, ys)` — diffusion physics in μm.
- **PixelSimulator**: 2D grid `nx × ny`, `DepositElectron(x_um, y_um)` (maps to pixel with center at origin), `AddDarkCurrent()`, `AddReadoutNoise()`, `PixelCharges()` row-major. Modes: LocalPatch, RowSegment, FullCCD.
- **PatternClassifier**: `ScanRow(row)` on 1D row → list of PatternResult (label, charge, etc.); same Pm/Pmn/Pmnl and thresholds as notebook.
- **PatternMC**: 1D row (row_length), BuildPatternTable → P(pattern|n_e); optional PrintExampleRows.

Notebook: 2D image (3×50) from (300×50) with 100:1 row binning, diffusion from random z and σ_xy(z,E), readout (+ optional DC).

---

## 1. 2D binned image generator

**Objective:** Generate a single 2D binned image (e.g. 3×50) matching notebook’s `generate_image_E(ne, Ee)`.

**Option A – New class `PatternImageGenerator` (recommended)**  
Lives next to PatternMC, uses existing ChargeTransport and a dedicated 2D grid + binning.

- **Config from JSON:** Image size and binning are read from the same config file (e.g. `configs/ccdarksens_scan_dmelectron_pattern.json`) under a dedicated block so they are not confused with the 1D pattern MC settings.
  - **Current JSON:** `response.pattern_mc` has `row_length`, `rows_bin`, `cols_bin` — those are for the **1D row** (PatternMC), not the 2D notebook-style image. So today the 2D image size and binning are **not** in the config.
  - **Add** a block under `response`, e.g. `pattern_image`, with the 2D image parameters. Then the 2D generator and any app use this block.

- **Suggested keys** (under `response.pattern_image` in the JSON):
  - `nrows_binned` (e.g. 3) — number of rows after binning
  - `ncols` (e.g. 50) — number of columns
  - `row_binning` (e.g. 100) → internal grid has `(nrows_binned * row_binning) × ncols` pixels before binning
  - `pixel_size_um`, `sigma_readout_e`, `lambda_dc` (optional), `rng_seed`
  - Diffusion can be inherited from `pattern_mc` (A_um2, b_umInv, alpha, beta_per_keV) and detector thickness, or duplicated in `pattern_image` if desired.

- **Example** addition to `configs/ccdarksens_scan_dmelectron_pattern.json`:
```json
"response": {
  "pattern_mc": { ... },
  "pattern_image": {
    "nrows_binned": 3,
    "ncols": 50,
    "row_binning": 100,
    "pixel_size_um": 15.0,
    "sigma_readout_e": 0.21,
    "lambda_dc": 0.0,
    "rng_seed": 987654321
  }
}
```
  ConfigManager (or the app) reads `pattern_image` when building the 2D generator; if the block is absent, 2D image generation can be skipped or use defaults.
- **Logic:**
  1. Internal grid: `ny_raw = nrows_binned * row_binning`, `nx = ncols` (e.g. 300×50).
  2. Pick random z ∈ [0, thickness], compute σ_xy = SigmaXYUm(z, Ee). Pick center (x0, y0) in pixel coords (e.g. middle of grid with margin).
  3. Sample ne electrons: Gaussian in x,y with σ_xy (in μm), convert to pixel indices (using pixel_size_um and grid origin). Clip to [0, nx-1], [0, ny_raw-1], add 1.0 to those pixels.
  4. Row binning: sum every `row_binning` rows → shape (nrows_binned, ncols).
  5. Add readout noise (Gaussian) to each pixel; optionally add Poisson dark current.
- **Output:** 2D array (e.g. `std::vector<std::vector<double>>` or store in a small “image” struct). Optionally also expose as ROOT `TH2D` for plotting.

**Option B – New PixelSimulator mode**  
Add e.g. `PixelSimMode::Image2D` with `nx`, `ny`, `row_binning`. PixelSimulator would handle a large 2D grid and a “binning” step. More invasive and mixes two concerns (pixel response vs. notebook-style geometry).

**Recommendation:** Option A. Keep PixelSimulator as-is; implement 2D + binning in `PatternImageGenerator` and call ChargeTransport for (z, σ_xy, cloud).

**Files to add:**
- `include/ccdarksens/response/PatternImageGenerator.hh`
- `src/response/PatternImageGenerator.cc`
- Config: extend `ConfigManager` / JSON for `pattern_image` (or pass from app).

**Note:** ChargeTransport currently expects `SampleCloudXY(..., sigma_xy_um, ...)` in μm. For 2D image in pixel coords, either: (i) generate cloud in μm with a chosen center (e.g. center of (ncols/2, ny_raw/2) in μm) and pass to a PixelSimulator configured as a large 2D grid (nx=ncols, ny=ny_raw) with appropriate origin, or (ii) implement pixel-domain sampling in PatternImageGenerator using `SigmaXYUm(z,E)/pixel_size_um` as σ in pixels and add electrons to a 2D buffer, then bin. (i) reuses DepositElectron; (ii) is simpler for “exact” notebook layout. Prefer (ii) for minimal changes and clear notebook match.

---

## 2. Plot 2D image (plot_image)

**Objective:** Visualize the 2D binned image (e.g. save as PNG).

- In the **example app** (or a small tool): after generating one 2D image (from §1), fill a `TH2D` (bins = nrows_binned × ncols), draw with colz, save canvas to PNG (e.g. `example_2d_image.png`).
- Alternatively: write image to a simple CSV (row, col, value) and document “plot with Python/matplotlib” for flexibility.
- **Recommendation:** Use ROOT `TH2D` + `TCanvas::SaveAs` in the same app that generates the image so we have a direct C++ “plot_image” without Python.

---

## 3. Scan 2D image for patterns (scan_image)

**Objective:** Run pattern finding on the 2D binned image (not only on a 1D row).

- Notebook: scans 2D (3×50), uses middle row and neighbors for isolation (above/below empty).
- **Options:**
  - **A)** Add `PatternClassifier::ScanImage2D(ny, nx, image_row_major)` that loops over the middle row (or each row), extracts 5-pixel windows, checks above/below for isolation, calls existing 1D classification logic. Then we have full 2D scan in C++.
  - **B)** Only “middle row” path: extract middle row from 2D image, call existing `ScanRow(middle_row)`. No isolation check; matches notebook only for the middle-row-only case.
- **Recommendation:** Start with **B** (middle row → ScanRow) for minimal code and to get 2D image + one-row pattern ID working. Add **A** later if you need isolation and full 2D scan.

---

## 4. Charge distribution per n_e (Ploting_simulation_specific_E)

**Objective:** For each n_e, histogram total charge (or pattern charge) over Nsims events.

- **Where:** PatternMC or a small app. For each n_e, run Nsims trials; for each trial, get the chosen pattern’s total charge (or sum of pixel charges in the row). Fill a `TH1D` per n_e.
- **Output:** Write ROOT file with one TH1D per n_e (e.g. `charge_ne1`, `charge_ne2`, …), or one TH2D (n_e vs charge). Optionally also save PNG per n_e (charge distribution + qmin/qmax lines) via TCanvas.
- **Config:** Reuse same physics (ChargeTransport, PixelSimulator, row_length). Can be a flag in the example app (e.g. `verbosity` or `output_charge_histos`) or a dedicated small app.

---

## 5. Sanity checks (sanity_checks)

**Objective:** Scatter plots (e.g. q1 vs q2 for pattern 11) to verify pattern ID.

- Run N events per n_e; for each event store (pattern label, charge_1, charge_2, …). For pattern (1,1) we have two charges.
- **Output:** Write a TTree (pattern_id, q1, q2, …) or CSV, then either:
  - **C++:** Fill `TGraph` or `TH2D` for a few pattern types and save PNGs; or
  - **Python:** Document “load CSV/TTree and run the notebook’s sanity_checks”.
- **Recommendation:** Add optional mode in example app that writes a CSV or TTree “sanity” file (event, ne, pattern_id, q1, q2, q3). Use Python/notebook to plot initially; add ROOT scatter plots later if desired.

---

## 6. simulate_cluster (ideal cluster + noise)

**Objective:** Build an ideal 3×5 cluster with charges (b,c,d) in the middle row, add readout noise, return 2D image (and optionally classify).

- **Implementation:** New function or small class:
  - Allocate 3×5 (or 3×ncols) buffer; set [1][1]=b, [1][2]=c, [1][3]=d (middle row).
  - Add Gaussian(sigma_readout_e) to each pixel.
  - Return as 2D array (or TH2D). For classification: extract middle row (5 pixels), call `ScanRow(middle_row)`.
- **Place:** Could live in `PatternImageGenerator` (e.g. `SimulateCluster(b,c,d)`) or in a small utility used by the example app. Expose in example app when e.g. `debug.simulate_cluster: true` with (b,c,d) from config or fixed (1,1,0).

---

## 7. pattern_simulation_background + Background_efficiencies.csv

**Objective:** For each “ideal” pattern (1), (1,1), (2,1), …, simulate_cluster Nsims times, run classifier, count how often each pattern is identified, write Background_efficiencies.csv.

- **Implementation:**
  - Enumerate ideal patterns (e.g. (1), (1,1), (2,1), (1,2), (1,1,1), …) up to max electrons (e.g. 5).
  - For each ideal pattern (b,c,d), call simulate_cluster(b,c,d) Nsims times; each time extract middle row, ScanRow; count identified pattern IDs.
  - Efficiency = count / Nsims per (ideal pattern, identified pattern). Save to CSV (e.g. iden_pat, eff_(1), eff_(1,1), …) like the notebook.
- **Where:** New small app (e.g. `ccdarksens_pattern_background_eff`) or a mode in the example app. Depends on §6 (simulate_cluster).

---

## 8. E-dependent efficiency CSV (Generate_efficiencies with E grid)

**Objective:** Export ε(pattern, n_e, E) to CSV for use in plots/notebook.

- **Current:** `PatternMC::PrecomputeEpsilonVsEnergy(E_grid, ne_min, ne_max)` fills internal `epsilon_Ene_[iE][ne - ne_min]` (total ε per n_e per E, not per-pattern).
- **Gap:** No per-pattern table at each E; no CSV export.
- **Options:**
  - **A)** For each E in grid, call BuildPatternTable(ne_min, ne_max, E), then write rows (E, ne, pattern_id, P(pattern|n_e)) to CSV. Expensive but full parity.
  - **B)** Export current epsilon_Ene_ (total ε only) as CSV (E, ne, epsilon_total) for plotting.
- **Recommendation:** Add **B** first (small change in app or PatternMC). Then add **A** if you need per-pattern vs E (new method or loop in app).

---

## 9. Dark current in PatternMC

**Objective:** Optionally add dark current when building the pattern table.

- **Current:** In `PatternMC::BuildPatternTable`, `AddDarkCurrent()` is commented out.
- **Change:** Add config flag (e.g. in `pattern_mc`: `include_dark_current: true`) and set `PixelSimulator::SetDarkCurrent(lambda_dc)` when building the row. In BuildPatternTable, call `pixSim.AddDarkCurrent()` when the flag is set. Ensure `lambda_dc` is read from config (already in PixelSimulatorConfig).

---

## Implementation order (suggested)

| Phase | Item | Depends on | Effort (rough) |
|-------|------|------------|------------------|
| 1 | PatternImageGenerator (2D binned image) | — | Medium |
| 2 | Plot 2D image (TH2D + SaveAs PNG) | §1 | Small |
| 3 | Scan 2D image (middle row → ScanRow) | §1, §2 | Small |
| 4 | Charge distribution histos per n_e | — | Small |
| 5 | Sanity-check output (CSV/TTree) | — | Small |
| 6 | simulate_cluster(b,c,d) | — | Small |
| 7 | pattern_simulation_background + CSV | §6 | Medium |
| 8 | E-dependent ε CSV (total, then per-pattern if needed) | — | Small / Medium |
| 9 | Dark current in PatternMC | — | Small |

Do **1 → 2 → 3** first so you have “generate 2D image + plot + pattern from 2D” in C++. Then add **4, 5, 6, 7, 8, 9** in any order that fits your priorities.

---

## Config / apps

- **Config:** Add optional blocks `response.pattern_image` (for §1) and e.g. `response.pattern_mc.include_dark_current` (for §9). Use existing `run.verbosity` and `debug.*` where it helps (e.g. `debug.plot_2d_image`, `debug.simulate_cluster`).
- **Apps:** Implement §1–§3 and §2 in the **example one-point app** (when `observable_bins == "pattern"` and a flag): generate one 2D image, plot it, optionally run middle-row ScanRow and print result. Add **charge histos** and **sanity CSV** as optional outputs in the same app or a tiny second app. **Background efficiencies** (§7) are best in a dedicated small app that uses §6.

---

## Summary

- **2D binned image:** New `PatternImageGenerator` (configurable 2D grid + row binning, same diffusion/readout/DC as notebook).
- **Plot image:** Fill `TH2D` from that image, draw, `TCanvas::SaveAs` PNG.
- **Scan 2D:** Start with “middle row → ScanRow”; later optionally add full 2D scan with isolation.
- **Charge distros, sanity, simulate_cluster, background eff, E CSV, DC:** All achievable with the above and small additions to PatternMC and apps; order and app split as in the table.

This keeps the existing pipeline unchanged and adds the notebook-parity features in a modular way.
