# Band-gap Phase C — coupled pheno limit scans

## Quick start

```bash
cd /Users/diegovenegasvargas/Documents/CCDarkSens
source /path/to/root/bin/thisroot.sh

# 1) Regenerate all 12 scan JSONs + manifest
python3 utils/gen_band_gap_pheno_scan_configs.py
python3 utils/update_band_gap_pheno_manifest.py

# 2) Smoke test (5x3 grid, ~minutes)
python3 utils/run_band_gap_phase_c.py smoke

# 3) Production (12 x 80x30 grids — long; run in screen/tmux)
bash utils/run_band_gap_phase_c_batch.sh
# or one tier / one case:
python3 utils/run_band_gap_phase_c.py scan --tier B-thresh
python3 utils/run_band_gap_phase_c.py scan --case gap0p1_B-thresh

# 4) Limit overlays (after ROOT files exist)
python3 utils/run_band_gap_phase_c.py plot-limits --tier B-thresh
python3 utils/run_band_gap_phase_c.py plot-limits --tier D-equal
```

## Cases (12)

| Gap ε | D-equal (`eh = ε`) | B-thresh (`eh = 3.8`) |
|-------|--------------------|------------------------|
| 0.1–1.2 eV | `configs/scan_band_gap_pheno_*_eh{gap}.json` | `*_eh3p8.json` |

Outputs: `outputs/scan_band_gap_{gap}_eh{eh}/scan_dmelectron_pattern.root`

Limits: `outplots/band_gap_pheno/step5_limits/limit_sweep_B_thresh.pdf` (and D-equal)

## Settings

- `observable_bins: "ne"`, `roi_bins: [1,2,3,4,5]`
- Profile likelihood, Bp/Br background, 0.5 kg·yr
- Pheno model unchanged (anchored p100K); see one-point QA under `outputs/band_gap_one_point_spectra/`

## Logs

`outputs/phase_c_logs/scan_B-thresh.log`, `scan_D-equal.log`, etc.
