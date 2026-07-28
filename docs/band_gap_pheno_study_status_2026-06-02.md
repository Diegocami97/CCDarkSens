# Band-Gap Pheno Study Status (2026-06-02)

## Scope and Goal

This study evaluates how Si DM-electron sensitivity changes when varying:

- **Scissor gap** in QCDark2 (bulk rates, `dR/dE`)
- **Ionization mapping** via scaled p100K tables (`P(n_e|E)`)

across six gap points:

- `0.1, 0.3, 0.5, 0.7, 0.9, 1.2` eV

and two ionization scenarios per gap:

- **D-equal**: `eh_pair_eV = gap`
- **B-thresh**: `eh_pair_eV = 3.8`

The current analysis focus is **electron-bin space** with ROI up to 5 electrons:

- `observable_bins = "ne"`
- `roi_bins = [1,2,3,4,5]`

---

## Method / Approach

## Physics layers (decoupled by construction)

- **Layer A (rates)**: QCDark2 dielectric + scissor -> `dR/dE`
- **Layer B (ionization)**: scaled p100K CSVs -> `P(n_e|E)` -> folded `S_true(n_e)`

QCDark2 mediator setting changes Layer A only; p100K tables are reused between heavy/light mediator studies.

## Execution structure

1. Generate/validate all required inputs (epsilon, rates, p100K)
2. One-point sanity scans (fixed `mchi`, `sigma`) and before/after spectra QA
3. Full Phase C scans (80x30 grid) for all cases
4. Overlay limit curves by scenario

---

## Inputs and Config Coverage

## Available epsilon files (QCDark2)

- `data/qcdark2_epsilon/Si/Si_fast_gap0p1.h5`
- `data/qcdark2_epsilon/Si/Si_fast_gap0p3.h5`
- `data/qcdark2_epsilon/Si/Si_fast_gap0p5.h5`
- `data/qcdark2_epsilon/Si/Si_fast_gap0p7.h5`
- `data/qcdark2_epsilon/Si/Si_fast_gap0p9.h5`
- `data/qcdark2_epsilon/Si/Si_fast_gap1p2.h5`

## Ionization tables

All D-equal/B-thresh tables for six gaps exist under:

- `data/p100K_gap*_eh*.csv`

## Full-scan grid

All Phase C runs use:

- `mchi`: logspace 0.2 -> 1000 MeV, 80 points
- `sigma_e`: logspace 1e-46 -> 1e-26 cm^2, 30 points

---

## Heavy Mediator Results (Completed)

Full heavy-mediator coupled scans are complete for all 12 cases.

## Limit overlays

- `outplots/band_gap_pheno/step5_limits/limit_sweep_B_thresh.pdf`
- `outplots/band_gap_pheno/step5_limits/limit_sweep_B_thresh.csv`
- `outplots/band_gap_pheno/step5_limits/limit_sweep_D_equal.pdf`
- `outplots/band_gap_pheno/step5_limits/limit_sweep_D_equal.csv`

## One-point and diagnostic outputs

Main one-point outputs and overlays are under:

- `outplots/band_gap_one_point_spectra/`

including:

- Per-case before/after: `before_after__*.pdf`
- Per-gap D-equal vs B-thresh: `compare_D-equal_vs_B-thresh__gap*.pdf`
- Scenario-only overlays: `compare_*_D-equal_all_gaps.pdf`, `compare_*_B-thresh_all_gaps.pdf`
- Folded low-`n_e` comparisons: `compare_ionization_*_ne1to5.pdf`
- Integral table: `integrals_summary.txt`

Raw p100K table comparisons are under:

- `outplots/band_gap_one_point_spectra/p100K_Pne/`

---

## Light Mediator Results (Completed)

The same full exercise was repeated with **QCDark2 light mediator** (`mediator = "light"`), same 80x30 grid and same reference point strategy.

## Rate grids generated

Light rates per scissor gap were generated under:

- `data/qcdark2_rates/Si/light/Si_fast_gap0p1/`
- `data/qcdark2_rates/Si/light/Si_fast_gap0p3/`
- `data/qcdark2_rates/Si/light/Si_fast_gap0p5/`
- `data/qcdark2_rates/Si/light/Si_fast_gap0p7/`
- `data/qcdark2_rates/Si/light/Si_fast_gap0p9/`
- `data/qcdark2_rates/Si/light/Si_fast_gap1p2/`

## Light Phase C scan outputs

- `outputs/scan_band_gap_light_{gap}_eh{eh}/scan_dmelectron_pattern.root` for all 12 cases

## Light limit overlays

- `outplots/band_gap_pheno/step5_limits/limit_sweep_light_B_thresh.pdf`
- `outplots/band_gap_pheno/step5_limits/limit_sweep_light_B_thresh.csv`
- `outplots/band_gap_pheno/step5_limits/limit_sweep_light_D_equal.pdf`
- `outplots/band_gap_pheno/step5_limits/limit_sweep_light_D_equal.csv`

---

## Key Observations So Far

## 1) D-equal vs B-thresh behavior is materially different in ROI `n_e <= 5`

- **B-thresh** tends to keep signal in low-electron bins.
- **D-equal** (especially small gap, e.g. 0.1) shifts significant signal to higher `n_e`, reducing ROI occupancy despite larger total rates.

## 2) Total rate and ROI signal are not equivalent

For D-equal, smaller gap can increase `dR/dE` integrals while still weakening low-`n_e` bins.

## 3) This effect is currently treated as model-intrinsic

Given current anchored p100K scaling (`E_gap` and `eh_pair_eV` choices), the observed redistribution is expected from the model assumptions.

---

## Approaches Implemented in Code

Automation and reproducibility scripts added/used:

- `utils/run_band_gap_one_point_spectra.py`
- `utils/plot_band_gap_one_point_spectra_compare.cc`
- `utils/plot_band_gap_one_point_p100K.py`
- `utils/run_band_gap_phase_c.py`
- `utils/run_band_gap_phase_c_batch.sh`
- `utils/gen_qcdark2_generate_Si_light_gap_configs.py`
- `utils/gen_band_gap_light_pheno_scan_configs.py`
- `utils/run_band_gap_light_phase_c_batch.sh`

---

## Current Concerns / Risks

## A) Interpretation risk: model-choice sensitivity

Differences between D-equal and B-thresh in low-`n_e` ROI can dominate conclusions. Limits should always be interpreted per ionization scenario, not gap-only.

## B) Numerical mapping artifacts at extreme D-equal (small gap)

At very small gap/`eh` (notably `0.1`), strong energy remapping can undersample or distort specific `P(n_e|E)` features on coarse grids. This is a known caveat in the current pheno mapping.

## C) Plot hygiene

Some legacy plots remain in output folders (older all-case overlays or smoke artifacts). The preferred references are the scenario-separated and per-gap overlays listed above.

---

## Recommended Next Decisions (Before Model Changes)

1. **Review heavy and light overlay limits side-by-side** for B-thresh and D-equal.
2. **Lock baseline claim** on the scenario that matches near-term physics priority for ROI `n_e <= 5` (likely B-thresh baseline + D-equal as bracket).
3. If needed, then open a dedicated follow-up to revise ionization pheno assumptions (do not mix with this completed baseline production set).

---

## Repro Commands (Core)

Heavy/light full scans were run via batch scripts. Useful entry points:

```bash
python3 utils/gen_qcdark2_generate_Si_light_gap_configs.py
python3 utils/gen_band_gap_light_pheno_scan_configs.py
bash utils/run_band_gap_light_phase_c_batch.sh
```

```bash
python3 utils/run_band_gap_phase_c.py scan --tier B-thresh
python3 utils/run_band_gap_phase_c.py scan --tier D-equal
python3 utils/run_band_gap_phase_c.py plot-limits --tier B-thresh
python3 utils/run_band_gap_phase_c.py plot-limits --tier D-equal
```

