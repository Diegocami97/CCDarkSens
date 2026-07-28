# CCDarkSens Refactoring Changelog

Summary of all structural changes made to the codebase during the refactoring sessions. Intended as a reference for updating Cursor context or onboarding collaborators.

---

## Custom Dielectric CSV Path (dark photon backend)

**Motivation:** Published optical conductivity / dielectric tensor data (e.g. from DFT codes like WIEN2k or digitized from papers) can now be fed directly into the dark photon absorption pipeline without adding a new hardcoded material. This enables rapid exploration of any narrow-gap material for which Im[ε_ii(ω)] is available.

**Files modified:**
- `python/ccdarkphys/darkphoton/entry.py` — added `_load_eps_csv`, `_loss_function`, `_pol_avg_loss`, `_write_custom_darkelf_files` helpers; extended `compute_dRdE` with `eps_csv` and `density_g_cm3` kwargs; updated `_header_lines` and `_cli()`.
- `utils/darkphoton_generate_grid.py` — threads `eps_csv` and `density_g_cm3` from JSON config through task dicts to `compute_dRdE`.

**Files added:**
- `configs/darkphoton_generate_custom_eps_template.json` — annotated template showing all required fields for the custom path.

**CSV format** (one `#`-prefixed header line, then comma-separated rows):
```
# omega_eV, re_eps_xx, im_eps_xx[, re_eps_yy, im_eps_yy, re_eps_zz, im_eps_zz]
0.01, 12.3, 0.0
...
```
yy and zz columns are optional; absent axes fall back to xx (isotropic). The polarization-averaged loss function W_avg = (1/3)(W_xx + W_yy + W_zz) is encoded into a temporary darkelf material directory and cleaned up after each call. `density_g_cm3` and `band_gap_eV` are required when `eps_csv` is set.

---

## Phase 1 — Dead Code Removal

### Files deleted
| File | Reason |
|------|--------|
| `include/ccdarksens/response/ClusterMC.hh` | ClusterMC backend abandoned in favour of PatternMC (now EfficiencyMC) |
| `src/response/ClusterMC.cc` | Same |
| `apps/ccdarksens_sensitivity_cluster_mc.cc` | ClusterMC-only app, never used |
| `apps/ccdarksens_check_dmelectron_chain_cluster_mc.cc` | Same |
| `apps/ccdarksens_check_dmelectron_chain_cluster_mc_old.cc` | Same |

### Files modified
- **`CMakeLists.txt`** — removed all three deleted app targets and `ClusterMC.cc` from `ccdarksens_core` source list.
- **`include/ccdarksens/io/ConfigManager.hh`** — removed `ClusterMCJSON` struct and its member `cmc` from `ResponseJSON`.

  > **Note:** `ClusterMCJSON` and the `cmc` member were subsequently re-introduced as a legacy-compat shim: `ConfigManager.cc` copies `cmc` fields into `emc` when a legacy `"cluster_mc"` JSON key is found, so old configs do not silently break. The struct exists only in the config layer; the corresponding C++ class is gone.

---

## Phase 2 — Extract Shared Diffusion Physics

**Problem:** The formula `σ_xy(z,E) = √(−A·ln(1 − b·z)) · (α + β·E_keV)` was duplicated in `Diffusion.cc` and `ChargeTransport.cc`.

### Files added
- **`include/ccdarksens/response/DiffusionPhysics.hh`** — single inline free function `ComputeSigmaXYUm(z, E, A, b, alpha, beta)`. Returns `0.0` (not NaN) when argument is out of range.

### Files modified
- **`src/response/Diffusion.cc`** — `sigma_xy_um_()` now delegates to `ComputeSigmaXYUm`.
- **`src/response/ChargeTransport.cc`** — `SigmaXYUm()` now delegates to `ComputeSigmaXYUm`.

---

## Phase 3 — Rename PatternMC → EfficiencyMC

**Rationale:** `PatternMC` handles both pattern-space and n_e-space acceptance; `EfficiencyMC` better describes its role.

### Files renamed
| Old | New |
|-----|-----|
| `include/ccdarksens/response/PatternMC.hh` | `EfficiencyMC.hh` |
| `src/response/PatternMC.cc` | `EfficiencyMC.cc` |

### Symbol renames (all files)
| Old | New |
|-----|-----|
| `PatternMC` (class) | `EfficiencyMC` |
| `PatternMCConfig` (struct) | `EfficiencyMCConfig` |
| `PatternMCJSON` (config struct) | `EfficiencyMCJSON` |
| `ResponseJSON::pmc` (member) | `ResponseJSON::emc` |
| `response_.pmc` (ConfigManager.cc) | `response_.emc` |
| `cfg.response().pmc` (all apps) | `cfg.response().emc` |
| local alias `pmj` (all apps) | `emj` |

### Files modified
- `include/ccdarksens/io/ConfigManager.hh` — struct renamed, member renamed, comments updated.
- `src/io/ConfigManager.cc` — all `response_.pmc` → `response_.emc`; comments updated.
- `src/response/DetectorResponsePipeline.cc` — `SetPatternMC()` → `SetEfficiencyMC()`, member `pmc_` → `emc_`.
- `include/ccdarksens/response/DetectorResponsePipeline.hh` — same.
- `CMakeLists.txt` — source filename updated.
- All apps that include PatternMC: `scan_dmelectron_pattern.cc`, `scan_dmelectron_grid_dualspace.cc`, `scan_dmelectron_pattern_Edep.cc`, `example_one_point_pattern.cc`, `xcheck_ne_pcd.cc`, `pattern_background_eff.cc`.

### JSON config key
The **JSON key** read from config files remains `"efficiency_mc"` (it was already updated in a prior step). No config file changes needed.

### Follow-up — complete the `pmc_` → `emc_` rename (2026-06-08)
The pipeline member rename was finished (the declaration in `DetectorResponsePipeline.hh` and the `ApplyEDependent` null-guard still used `pmc_`). For full consistency, the app-local variables were also renamed:

| Old | New |
|-----|-----|
| `EfficiencyMCConfig pmc` | `EfficiencyMCConfig emc_cfg` |
| `auto pmc_ptr` (shared_ptr) | `auto emc_ptr` |
| `EfficiencyMC patternMC` (instance) | `EfficiencyMC emc` |
| local alias `pmj` | `emj` |

Public API methods that legitimately refer to charge patterns (`SetPatternImageGenerator`, `GetPatternTable`, `PrecomputeEpsilonWithPatternEff`) were intentionally left unchanged.

---

## Phase 4 — Consolidate `safe_log`

**Problem:** Identical `safe_log` helper defined in anonymous namespaces in two stats source files.

### Files added
- **`include/ccdarksens/stats/StatsUtils.hh`** — `namespace ccdarksens::stats { inline double safe_log(double x) }`.

### Files modified
- `src/stats/PoissonAsimovPLR.cc` — removed local definition, added include.
- `src/stats/ProfileLikelihood.cc` — removed local definition, added include.

---

## Phase 5 — Extract Shared App Utilities

**Problem:** `MakeFlatEfficiency`, `SumROI`, and grid-axis expansion (`linspace`/`logspace`/`values`) copy-pasted across ≥4 apps each.

### Files added
- **`include/ccdarksens/utils/AppUtils.hh`** — three free functions:
  - `MakeFlatEfficiency(ne_min, ne_max)` → `std::unique_ptr<TH1D>`
  - `SumROI(h, roi_bins, ne_min)` → `double`
  - `ExpandAxis(spec, kind)` → `std::vector<double>`

### Files modified
All apps that previously had local copies now `#include "ccdarksens/utils/AppUtils.hh"` and call the shared versions:
- `ccdarksens_sensitivity.cc`
- `ccdarksens_sensitivity_fast.cc`
- `ccdarksens_scan_dmelectron_pattern.cc`
- `ccdarksens_scan_dmelectron_grid_dualspace.cc`
- `ccdarksens_scan_dmelectron_pattern_Edep.cc`

---

## Phase 6 — Deduplicate PatternImageJSON / PatternImageConfig

**Problem:** `PatternImageJSON` (in `ConfigManager.hh`) and `PatternImageConfig` (in `PatternImageGenerator.hh`) were nearly identical structs kept manually in sync.

### Change
- **`PatternImageJSON` removed** from `ConfigManager.hh`.
- `ConfigManager.cc` now directly populates `PatternImageConfig` during JSON parsing (the runtime type, no intermediate struct).
- `include_dark_current` field added to `PatternImageConfig` with default `false` to cover the one field that was previously only in `PatternImageJSON`.

---

## Phase 7 — Minor API Cleanup

### 7a. Stale "PatternMC" doc comments
All remaining `PatternMC` references in comments and string literals across:
- `include/ccdarksens/response/EfficiencyMC.hh` — file header, struct doc block, class doc block
- `src/response/EfficiencyMC.cc` — file header, `[PatternMC]` log prefixes, ROOT histogram titles
- `include/ccdarksens/response/PixelSimulator.hh` — one comment
- `include/ccdarksens/response/PatternClassifier.hh` — one comment
- `include/ccdarksens/io/ConfigManager.hh` — two inline comments
- `src/io/ConfigManager.cc` — two block comments
- `src/response/DetectorResponsePipeline.cc` — one error message string

### 7b. Expose PCD kernel parameters via config
`PCDJSON` struct in `ConfigManager.hh` previously had no `sigma_res_e`, `Dqmin`, `Dqmax` fields, so apps using `cfg.response().pcd.sigma_res_e` failed to compile.

- **`include/ccdarksens/io/ConfigManager.hh`** — added three fields to `PCDJSON`:
  ```cpp
  double sigma_res_e = 0.21;
  double Dqmin       = 0.5;
  double Dqmax       = 0.5;
  ```
- **`src/io/ConfigManager.cc`** — added parsing of the three new fields from the `"pcd"` JSON block.

---

## New Feature — B^rc_p Computation

### Files added
- **`apps/ccdarksens_compute_brc.cc`** — standalone app that computes the random-coincidence background:
  ```
  B^rc_p = N_total × Σ_q [ Π_i Poisson(q_i | λ_i) ] × B[p|q]
  ```
  Reads `Background_efficiencies.csv` and `Final_Combined_Image_Data.csv`. Accepts θ parameters on the command line. Prints the full result table plus a per-injected-pattern breakdown for the dominant pattern {11}.

  Usage: `build/ccdarksens_compute_brc [theta1] [theta2] [theta3]`
  Default θ = (3, 2, 2) → λ = (3×10⁻⁴, 2×10⁻⁴, 2×10⁻⁴) e⁻/pixel.

- **`CMakeLists.txt`** — registered `ccdarksens_compute_brc` target.

---

## New Documentation

### Files added
- **`docs/Beginners_Guide.md`** — end-to-end introduction covering: the two observable spaces (n_e bins vs pattern bins), detector parameters and exposure calculation, background model options (DC + flat d.r.u. vs Bp/Br template), Asimov mode, signal efficiency via EfficiencyMC or CSV, all key apps, full JSON config reference, LBC reproduction recipe, and a glossary.

---

## Summary Table

| What changed | Old name / location | New name / location |
|---|---|---|
| Efficiency MC class | `PatternMC` | `EfficiencyMC` |
| Efficiency MC config struct | `PatternMCConfig` | `EfficiencyMCConfig` |
| Efficiency MC JSON struct | `PatternMCJSON` | `EfficiencyMCJSON` |
| Config member accessor | `.pmc` | `.emc` |
| App local JSON alias | `pmj` | `emj` |
| Diffusion formula | duplicated in Diffusion.cc + ChargeTransport.cc | `DiffusionPhysics.hh::ComputeSigmaXYUm` |
| `safe_log` helper | duplicated in two stats .cc files | `StatsUtils.hh::ccdarksens::stats::safe_log` |
| App utilities (flat ε, SumROI, ExpandAxis) | duplicated in ≥4 apps | `AppUtils.hh::ccdarksens::utils::` |
| PatternImageJSON | separate struct in ConfigManager.hh | merged into PatternImageConfig |
| PCD kernel params | hardcoded in pipeline / missing from JSON | `PCDJSON::sigma_res_e`, `Dqmin`, `Dqmax` |
| ClusterMC backend | `ClusterMC.hh/.cc` + 3 apps | **deleted** |
| B^rc_p computation | not in C++ | `ccdarksens_compute_brc` app |
