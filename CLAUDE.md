# CCDarkSens — Claude Code Context

This file is read automatically at the start of every Claude Code session. It records architectural decisions and naming conventions so they are not accidentally undone.

---

## What this repo is

CCDarkSens is a C++17 dark matter sensitivity framework for DAMIC-M. It computes upper limits on the DM–electron cross section σ_e(mχ) using a profile likelihood ratio over either **n_e bins** or **pattern bins** as the observable. See `docs/Beginners_Guide.md` for a full introduction.

---

## Key naming conventions (do not revert)

| Concept | Correct name | Old name (deleted/renamed — do not use) |
|---|---|---|
| Efficiency MC class | `EfficiencyMC` | `PatternMC` — class deleted, file renamed |
| Efficiency MC config struct | `EfficiencyMCConfig` | `PatternMCConfig` |
| Efficiency MC JSON struct | `EfficiencyMCJSON` | `PatternMCJSON` |
| Config member accessor | `.emc` | `.pmc` |
| App local JSON alias | `emj` | `pmj` |
| Diffusion formula | `DiffusionPhysics.hh::ComputeSigmaXYUm` | duplicated in Diffusion.cc + ChargeTransport.cc |
| Shared stats helper | `StatsUtils.hh::ccdarksens::stats::safe_log` | duplicated local lambdas |
| Shared app utilities | `AppUtils.hh::ccdarksens::utils::` | copy-pasted per-app |

---

## Files that no longer exist (do not recreate)

- `include/ccdarksens/response/ClusterMC.hh` — deleted (ClusterMC backend abandoned)
- `src/response/ClusterMC.cc` — deleted
- `apps/ccdarksens_sensitivity_cluster_mc.cc` — deleted
- `apps/ccdarksens_check_dmelectron_chain_cluster_mc.cc` — deleted
- `apps/ccdarksens_check_dmelectron_chain_cluster_mc_old.cc` — deleted

`ClusterMCJSON` and the `cmc` member **still exist** in `ConfigManager` as a legacy read shim (copies into `emc` so old JSON configs don't silently break), but the corresponding C++ class is gone.

---

## Important shared headers (use these, don't duplicate)

- **`include/ccdarksens/response/DiffusionPhysics.hh`** — `ComputeSigmaXYUm(z, E, A, b, alpha, beta)`. Use this instead of reimplementing the formula.
- **`include/ccdarksens/stats/StatsUtils.hh`** — `ccdarksens::stats::safe_log(x)`. Use this in any stats code.
- **`include/ccdarksens/utils/AppUtils.hh`** — `MakeFlatEfficiency`, `SumROI`, `ExpandAxis`. Use these in apps instead of copy-pasting.

---

## Observable space — the primary config choice

Two observable spaces are supported. The choice is made in the JSON config:

```json
// n_e bins:
"experiment": { "observable_bins": "ne",      "roi_bins": [1,2,3,4,5] },
"response":   { "analysis_space":  "ne" }

// Pattern bins (SRDM published result):
"experiment": { "observable_bins": "pattern", "pattern_roi": [11,21,111,31,22,211] },
"response":   { "analysis_space":  "pattern" }
```

---

## Background model — the secondary config choice

```json
// Compute from dark current + flat d.r.u. (use for projections / n_e space):
"run": { "background_source": "dc_flat_migration", "background_model": "scale" }

// Use pre-computed Bp/Br vectors (use for pattern space / reproducing LBC result):
"run": {
  "background_source": "bp_br_template",
  "background_model":  "Bp_theta_Br",
  "background_Bp": [141.4, 0.111, 0.042, 0.019, 2.5e-5, 5.8e-5],
  "background_Br": [0.039, 0.039, 0.016, 0.052, 0.011, 0.035]
}
```

---

## Validated physics (do not change defaults without re-running validation)

The following values were cross-checked against the collaboration Python reference and pass at < 3σ pull:

| Quantity | Value | Validated by |
|---|---|---|
| `sigma_readout_e` | 0.16 e⁻ | `ccdarksens_validate_pattern_efficiency` |
| `thr_M` | 3.5 × σ_ro | same |
| `thr_MN` | 4.0 × σ_ro | same |
| `thr_MNL` | 5.5 × σ_ro | same |
| `Qmin_e` | 0.60 e⁻ | same |
| B[{11}\|{11}] | ≈ 0.908 | `ccdarksens_validate_background_efficiency` |
| B^rc_{11} diagonal term | 141.4 at θ₁=2.90 | `ccdarksens_compute_brc` |

---

## Key apps (production use)

| App | Purpose |
|-----|---------|
| `ccdarksens_example_one_point_pattern` | Single (mχ, σ) diagnostic run — run this first |
| `ccdarksens_scan_srdm_pattern_csv` | Full (mχ, σ) grid scan |
| `ccdarksens_band` | Sensitivity band (median ± 1σ, 2σ) |
| `ccdarksens_plot_dmelectron_limit` | Exclusion limit plot from ROOT output |
| `ccdarksens_compute_brc` | Compute B^rc_p from θ parameters and confusion matrix |
| `ccdarksens_validate_pattern_efficiency` | Cross-check EfficiencyMC vs reference CSV |
| `ccdarksens_validate_background_efficiency` | Cross-check B[p\|q] matrix vs reference CSV |

---

## Documentation index

| File | Contents |
|------|----------|
| `docs/Beginners_Guide.md` | Full introduction: concepts, all config fields, worked examples, glossary |
| `docs/Refactoring_Changelog.md` | Complete record of all structural changes made during refactoring |
| `docs/Exposure_usage_audit.md` | How `exposure_kg_year` flows through the pipeline (applied exactly once) |
| `docs/Pydme_CCDarkSens_crosscheck.md` | Crosscheck of likelihood and background against pydme |
| `docs/App_flow_walkthrough.md` | Step-by-step walkthrough of scan and example apps |
