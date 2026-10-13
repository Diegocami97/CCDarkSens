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
| Diffusion formula | `DiffusionPhysics.hh::ComputeSigmaXYUm` | previously duplicated in Diffusion.cc + ChargeTransport.cc — **resolved**, both now delegate to the shared header |
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

**Do not confuse this with `ClusterFitMC`** (`include/ccdarksens/response/ClusterFitMC.hh`) — a new, unrelated module for the WIMP-nucleon channel's per-event ΔLL reconstruction (see `docs/ClusterFitMC_Design.md`). Similar name, different purpose, not a revival of the deleted backend.

---

## Important shared headers (use these, don't duplicate)

- **`include/ccdarksens/response/DiffusionPhysics.hh`** — `ComputeSigmaXYUm(z, E, A, b, alpha, beta)`. Use this instead of reimplementing the formula.
- **`include/ccdarksens/stats/StatsUtils.hh`** — `ccdarksens::stats::safe_log(x)`. Use this in any stats code.
- **`include/ccdarksens/utils/AppUtils.hh`** — `MakeFlatEfficiency`, `SumROI`, `ExpandAxis`. Use these in apps instead of copy-pasting.

---

## Efficiency MC — where P(pattern | n_e) really comes from

Production scans do **not** run the MC: every pattern-/n_e-space config in `configs/examples/` sets `response.efficiency_mc.efficiency_csv`, which loads the collaboration table `data/Efficiencies_patterns_Nsims1000000_DCTrue_alpha1.csv`. That table was made with the reference 2D simulation (3×50 binned image, 100-row binning, Poisson dark-current pileup λ = 0.00041 per binned pixel). The MC code only runs when no CSV is given.

| Path | Config trigger | Status |
|---|---|---|
| **2D image** | `response.efficiency_mc.use_2d_image_efficiency: true` | Reproduces the reference table (`ccdarksens_validate_pattern_efficiency`, RMS pull ≈ 1σ at 1e6 trials). Set `response.pattern_image.lambda_dc` (0.00041 to match the table; the default 0 switches pileup off). `PatternImageConfig::randomize_center` defaults to `true` — the reference draws a random sub-pixel cloud centre, and a fixed centre sits on a pixel corner and gives a wrong table. |
| **1D row** | `use_2d_image_efficiency: false` and no CSV | Approximate single-row model, **not** validated against the reference (it does not reproduce the shipped table). Needs `pix_cfg.collapse_y = true` for the cloud to stay in the row. Diagnostic use only; do not enable `include_dc_pileup` here. |

**Two different signal tables exist — the published analysis uses the second one.** `data/efficiencies_paolo.csv` summed over the six search patterns gives (38.4, 65.0, 78.7, 86.3)% for a (2, 3, 4, 5) e⁻ deposit, exactly what the PRL (arXiv:2503.14617) quotes, (38, 65, 79, 86)%. The pydme reference table above gives (37.5, 63.9, 74.7, 24.7)%. The file is a merged table: its single-pixel rows follow the α = 1 notebook run, its six search-pattern rows come from the published selection with pydme's diffusion α = 0.83 and the header's 1.5e-4 pile-up.

The difference is the pattern selection, now selectable with `response.pattern_classifier.selection`:
- `"tree"` (default): the notebook decision tree of `efficiencies.py` (horizontal isolation, 311/221/32/41 classes). Reproduces the pydme tables (`ccdarksens_pattern_signal_eff` with `configs/examples/pattern_signal_eff_tree_pydme.json`: RMS pull 1.13 against the reference table).
- `"chain"`: the selection of the DAMIC-M analysis, written from the data-selection note and the PRL — a pair passing p₁₁<4 is promoted to 2+1 (p₂₁<4) and then 3+1 (p₃₁<4); 1+1 and 2+1 extend to 1+1+1 / 2+1+1 with a third 1e⁻ pixel on either side (p<5.5); the highest total charge wins; only the rows above and below must be clean. Options: `chain_max_pair` (3, or 4 to add 4+1), `chain_allow_22`, `chain_horizontal_isolation`, `chain_first_only`. With `configs/examples/pattern_signal_eff_chain.json` (α = 0.83, pile-up 1.5e-4) the six-pattern sums are (38.7, 65.5, 79.3, 86.3)% against Paolo's (38.4, 65.0, 78.7, 86.3)% and the rows at n_e = 2 and 3 agree within statistics (RMS pull 1.2, 1.1; 3.5 at n_e = 4). The n_e = 5 rows still differ (RMS pull 17: {111} 0.015 against 0.028, {31} 0.526 against 0.489), as does P({11} | 1 e⁻) (Paolo has exactly 0); energy-dependent diffusion, pile-up, the charge cut and horizontal isolation do not explain them, so the exact published rules are not known beyond the note. The effect on the limit is small: a 25-mass scan with `efficiencies_paolo.csv` against the chain table differs by 0.7% at the median and 1.2% at most (the pydme tree table differs by up to 7.6% at high mass).

Check any regenerated table with `ccdarksens_validate_pattern_efficiency --new` (or `ccdarksens_validate_background_efficiency --new`): pooled two-sample binomial pulls. A table written by `ccdarksens_pattern_signal_eff` has only the patterns the chosen selection can output.

The background identification matrix B[p|q] is a separate table, `data/Background_efficiencies.csv`, regenerated with `ccdarksens_pattern_background_eff` and checked with `ccdarksens_validate_background_efficiency`. Efficiency is a detector-geometry property; dark current otherwise enters only as a background rate.

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
| `ccdarksens_scan_dmelectron_pattern` | The DM-electron reference scan app (pattern- and n_e-space, every `run.*` statistical mode). **Do not modify** — kept as the trusted reference/consumer for `ccdarksens_scan_generic`'s parity gates; see `docs/GenericScanApp_Design.md`. |
| `ccdarksens_scan_generic` | **Preferred for new work.** Single config-driven scan app for any channel (dm_electron, dark_photon, migdal, wimp_nucleon) × any analysis space (pattern, n_e, cluster_energy), including the 2D `pydme` minimizer, `Bp_theta_Br` background, `single_bin_likelihood`, `smooth_ul_envelope`, and `run.mode="threshold_toys"` for `ccdarksens_band`'s toy-MC calibration. Verified bit-exact against `ccdarksens_scan_dmelectron_pattern` on real production configs — see `docs/GenericScanApp_Design.md` §6-§9. |
| `ccdarksens_scan_srdm_pattern_csv` | Secondary scan app reading pre-folded `S_pat` CSVs (pattern-space only); also the original home of the `threshold_toys` toy-MC algorithm, since ported into `ccdarksens_scan_generic`. |
| `ccdarksens_band` | Sensitivity band (median ± 1σ, 2σ); works unmodified against either scan app above via `band.scan_binary`. |
| `ccdarksens_plot_dmelectron_limit` | Exclusion limit plot from ROOT output (the reference plotter; do not modify) |
| `ccdarksens_plot_limit` | **Preferred plotter.** Limit curves from scan ROOT files for every channel; finds literature curves from the repo root (any working directory), `--out-csv` exports the curves, `--ratio-to <n\|label>` prints window medians of the ratio of each curve to a reference (every curve is drawn at one confidence level, `--cl`, default 0.90 = one-sided asymptotic q = 1.6424, rebuilt from each file's q-map and labelled on the plot; another threshold, e.g. the positional `2.71`, is accepted but flagged and labelled with its own C.L.; masses with no crossing are excluded from the ratio and counted); the legacy Migdal canvas is only written with `--legacy-migdal-canvas`. |
| `ccdarksens_validate_response_factory` | Diff `ResponseFactory`/`BackgroundFactory` output against reference-app ground truth (both `background_source` modes) |
| `ccdarksens_compute_brc` | Compute B^rc_p from θ parameters (positional) or from a rates config (`--config`, same file as `background_pattern_rates.json`) with optional per-segment rates (`--segments`), print a `background_Bp` line and write it with `--out`. |
| `utils/pattern_pvalues.py`, `utils/compare_tables.py`, `utils/background_pattern_counts.py`, `utils/compare_limit_curves.py` | **Independent Python cross-checks**, not part of the workflow: each re-implements one stage and must agree with the app (`ccdarksens_classify_pair`, `ccdarksens_validate_* --new`, `ccdarksens_compute_brc --config`, `ccdarksens_plot_limit --draw-both --ratio-to`). Usage and agreeing numbers: §14 of the walkthrough. |
| `ccdarksens_classify_pair` | Hand check of one pair/triple of pixels: the p_mn statistics of the note and the label the real `PatternClassifier` returns for `chain` and `tree` |
| `ccdarksens_pattern_signal_eff` | Generate P(pattern \| n_e) with the `tree` or `chain` selection, in the reference CSV layout (generation lives in `PatternTables.cc`, shared with the validators) |
| `ccdarksens_pattern_background_eff` | Generate the background matrix B[p\|q] (reference layout, rows = identified pattern, columns = `eff_<injected>`); honours `pattern_classifier.selection` |
| `ccdarksens_validate_pattern_efficiency` | Cross-check the signal table vs a reference CSV: built-in simulation, `--config cfg.json` (same generation code as `ccdarksens_pattern_signal_eff`) or `--new t.csv` (compare existing tables, pooled binomial pulls) |
| `ccdarksens_validate_background_efficiency` | Cross-check B[p\|q] vs a reference CSV: built-in simulation, `--config cfg.json` or `--new m.csv` |
| `ccdarksens_validate_cluster_fit_engine` | WIMP-nucleon channel: inject known signals, cross-check Nelder-Mead vs Minuit2 fit recovery |
| `ccdarksens_calibrate_noise_tail` | WIMP-nucleon channel: pure-noise ΔLL tail calibration, convergence/stability checks |
| `ccdarksens_build_cluster_fit_kernel` | WIMP-nucleon channel: build K[E_true,E_reco] kernel, efficiency-curve sanity checks |
| `ccdarksens_example_one_point_cluster` | WIMP-nucleon channel: full config-driven single-point run (calibrate → build kernel → fold → PLR) |
| `ccdarksens_validate_wimp_nucleon_paper_repro` | WIMP-nucleon channel: config-driven calibration+kernel reproduction check against PhysRevD.94.082006's own Figs. 6/9 (paper detector parameters, not the modern-projection defaults) |

---

## Documentation index

| File | Contents |
|------|----------|
| `docs/Framework_Architecture.md` | Standalone architecture reference: how a config file turns into an exclusion curve, layer by layer (ConfigManager → Detector/ExperimentSetup → ModelFactory → ResponseFactory/ResponseFold → BackgroundFactory → ProfileLikelihood/ScanUtils → `ccdarksens_scan_generic` orchestration). Grounded directly in current source; covers the `q_threshold` convention subtlety (scan's own asymptotic value vs. the plotting app's hardcoded default) in detail. |
| `docs/Beginners_Guide.md` | Full introduction: concepts, all config fields, worked examples, glossary |
| `docs/DM_Signal_Models_Physics_Reference.md` | Consolidated physics reference: DM-e (QEDark/QCDark2), dark photon absorption, Migdal effect, charge ionization (Klein's formula), and the shared PLR pipeline — Si and SrCd₂Sb₂ cases |
| `docs/ClusterFitMC_Design.md` | WIMP-nucleus SI channel Phase 3/4 (noise-tail ΔLL calibration, per-event cluster reconstruction): pixel-shape derivation, closed-form ΔLL, minimizer strategy, validation results — built up slice by slice |
| `docs/Pipeline_Raw_Data_to_Limit.md` | The full analysis from raw CCD images to the exclusion limit (calibration, six masks, pattern selection, efficiency tables, B^rc and B^rad, signal, likelihood), what the framework has at each stage, a reproduction ledger against the PRL/note numbers, the pitfalls found, and the recipe for adding a new pattern bin. Seven Mermaid diagrams. |
| `docs/Student_Examples_Background_Pattern_Counts.md` | Short pointer: Bp from per-pixel k-electron rates is `ccdarksens_compute_brc`; the full walkthrough is `docs/Student_Examples_Pattern_Pipeline_Walkthrough.md`. Keeps the account of why the published `Bp` cannot be reproduced exactly (paper Eq. 4 sums per image and 10-column group; those occupancies are not in the repo). |
| `docs/Student_Examples_Pattern_Pipeline_Walkthrough.md` | The pattern pipeline one command at a time, from the rate tables to the limit: generate the QEDark rates (the one Python tool, `utils/qedark_generate_grid.py`), classify a pair, generate the signal table and the background matrix, check them, compute Bp, wire a scan config, scan, plot; usage, expected numbers and failures of every app (`ccdarksens_classify_pair`, `_pattern_signal_eff`, `_pattern_background_eff`, `_validate_*`, `_compute_brc`, `_scan_generic`, `_plot_limit`). |
| `docs/Student_Examples_Row_Cut_Recovered_Events.md` | Worked example from one student finding: an event removed by the row-cut mask, how it classifies by hand (`ccdarksens_classify_pair`), how to attribute count changes, what a mask change does to N_pix / B^rc / exposure, and the effect on the heavy-mediator LBC limit of changed observed counts (configs `dm_electron_heavy_lbc_rowcut_{A,B,C,D,E}.json`; D and E also scale Bp and exposure by 1+f; compare with `ccdarksens_plot_limit --draw-both --ratio-to`). Includes the pitfall that a limit at the sigma-grid floor (1e-46) means q_0 is above threshold, not a limit. |
| `docs/GenericScanApp_Design.md` | `ccdarksens_scan_generic`: the `ResponseFold`/`ResponseFactory`/`BackgroundFactory` abstractions, and how every DM-electron reference-app statistical mode (2D `pydme` minimizer, `Bp_theta_Br`, `single_bin_likelihood`, `smooth_ul_envelope`, `threshold_toys`) was ported and verified bit-exact — built up slice by slice |
