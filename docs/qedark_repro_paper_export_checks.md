# QEdark heavy repro — paper export checks (2026-06-03)

Companion to [qedark_ul_audit_pydme_srdm.md](qedark_ul_audit_pydme_srdm.md).

## Canonical reference

Use the **paper figure export** (25 points):

`Software/pydme/figures/ScienceRun2024-figures/data/ScienceRun2024_results-1/DAMIC-M_2025_QEDark_DMe_heavymediator.txt`

Do **not** use `DAMIC-M_this_work_QEDark_hm.csv` for paper parity (different curve at low mass).

## Commands

```bash
# Comparison plot (stored UL + --from-qhist dashed)
python3 utils/plot_qedark_all_references.py

# Audit vs paper export (default)
python3 utils/audit_qedark_ul_vs_pydme.py

# Halo smoke test (rate scaling + optional few-mass scan)
PYTHONPATH=python python3 utils/halo_smoke_test_qedark.py --run-scan
```

## Results summary

### Stored UL vs paper export (full grid)

| Region | scan / paper |
|--------|----------------|
| m ≈ 2 MeV | 0.83 |
| 1.2–500 MeV median | 1.00 |
| 5–500 MeV median | 1.06 |

### `--from-qhist` vs paper export

| Region | qhist / paper |
|--------|----------------|
| m ≈ 2 MeV | 0.81 |
| 1.2–500 MeV median | 0.92 |

qhist makes the scan **more sensitive** (~14% on median). Keep stored UL as primary.

### Halo (v_E = 253.7 vs 263 km/s)

At 2 MeV, weakening rates by ~11% raises the UL by ~12% → scan/paper moves **0.83 → ~0.93**. Does not close the full 1–5 MeV offset.

## Paper-25 mass grid test

Script: [`utils/run_paper25_qedark_scan.py`](../utils/run_paper25_qedark_scan.py)  
Config: [`configs/scan_dmelectron_pattern_data_qedark_paper25.json`](../configs/scan_dmelectron_pattern_data_qedark_paper25.json)

Scans exactly the **25 paper-export masses** (rates symlinked from nearest `long_scan` point).

| Comparison | m ≤ 5 MeV median | 5–500 MeV median |
|------------|------------------|------------------|
| paper25 / paper | **0.83** | **1.07** |
| fullgrid / paper (interp.) | **0.75** | **1.04** |
| paper25 vs fullgrid @ same m | **≈ 0** (identical) | **≈ 0** |

**Conclusion:** Using the same 25 masses does **not** materially improve agreement. At 2 MeV, paper25/paper = fullgrid/paper = **0.83**. The low-m offset is not primarily from 800-vs-25 mass interpolation — it is in the UL/rate/statistics pipeline itself. Sub-MeV points (0.75–0.85 MeV) show a larger gap (~0.3–0.5).

Plot: `outplots/qedark_repro/paper25_vs_paper_export.pdf`

```bash
python3 utils/run_paper25_qedark_scan.py --force-rates --run-scan
python3 utils/run_paper25_qedark_scan.py --skip-rates --compare-only
```

## Baxter full grid (in progress / complete)

| Item | Path |
|------|------|
| Rates (v_E=253.67) | `data/qedark_rates/Si/heavy/long_scan_baxter_vE253p7/` |
| Mass ratios cache | `data/qedark_rates/Si/heavy/baxter_vE253p7_mass_ratios.json` |
| Scan config | `configs/scan_dmelectron_pattern_data_qedark_fullgrid_baxter.json` |
| Scan log | `outputs/scan_pattern_data_qedark_fullgrid_baxter/run.log` |
| Scan ROOT | `outputs/scan_pattern_data_qedark_fullgrid_baxter/scan_dmelectron_pattern.root` |

```bash
# Build rates (once; ~20 min)
PYTHONPATH=python python3 utils/build_baxter_halo_rates_fullgrid.py --workers 12

# Run scan (~800 masses; check log)
build/ccdarksens_scan_dmelectron_pattern configs/scan_dmelectron_pattern_data_qedark_fullgrid_baxter.json

# Compare when done
python3 utils/compare_baxter_fullgrid.py
python3 utils/plot_qedark_all_references.py --root outputs/scan_pattern_data_qedark_fullgrid_baxter/scan_dmelectron_pattern.root
```


- Paper export has **25** mass points; full grid has **800** — interpolation differences at m ≲ 5 MeV.
- Signal table provenance (pydme pre-folded pattern signals vs CCDarkSens pipeline).
- Confirm v_E used when the paper export was generated.
