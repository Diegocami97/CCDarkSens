#!/usr/bin/env bash
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  run_band_gap_pheno_all_figures.sh -- I regenerate all band-gap
#  phenomenology validation figures under outplots/band_gap_pheno/ (dR/dE,
#  P(n_e|E), S(n_e), scan diagnostics, limits and the 2D surface).
# ============================================================================

# Regenerate all band-gap pheno validation figures under outplots/band_gap_pheno/
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

source "${ROOT_INSTALL:-}/bin/thisroot.sh" 2>/dev/null || true

mkdir -p outplots/band_gap_pheno/{step1_dRdE,step2_p100K,step3_Sne,step4_scan_diag,step5_limits,2d_surface}

echo "======== Build plotter ========"
cmake --build build -j --target ccdarksens_plot_dmelectron_limit

RATES_HEAVY="data/qcdark2_rates/Si/heavy"
CURVES=(
  "gap0p1=${RATES_HEAVY}/Si_fast_gap0p1"
  "gap0p3=${RATES_HEAVY}/Si_fast_gap0p3"
  "gap0p5=${RATES_HEAVY}/Si_fast_gap0p5"
  "gap0p7=${RATES_HEAVY}/Si_fast_gap0p7"
  "gap0p9=${RATES_HEAVY}/Si_fast_gap0p9"
  "gap1p2=${RATES_HEAVY}/Si_fast_gap1p2"
)

echo "======== Step 1: dR/dE overlays ========"
python3 utils/plot_band_gap_pheno_dRdE.py \
  --mchi-MeV 1.0 --sigma 1.0e-35 \
  --out outplots/band_gap_pheno/step1_dRdE/dRdE_m1p0_heavy.pdf \
  --curve "${CURVES[@]}" --Emin 0 --Emax 5

python3 utils/plot_band_gap_pheno_dRdE.py \
  --mchi-MeV 10.0 --sigma 1.0e-35 \
  --out outplots/band_gap_pheno/step1_dRdE/dRdE_m10_all_gaps.pdf \
  --curve "${CURVES[@]}" --Emin 0 --Emax 5

echo "======== Step 2: p100K scaling QA ========"
python3 utils/eval_p100K_scaling.py --outdir outplots/band_gap_pheno/step2_p100K
python3 utils/plot_p100K_scaling_compare.py --outdir outplots/band_gap_pheno/step2_p100K
# Pne_compare_scenarios.pdf written directly by plot_p100K_scaling_compare.py (Fig. 3)

echo "======== Step 5: limit curves (Phase C tiers) ========"
python3 utils/run_band_gap_phase_c.py plot-limits --all

echo "======== Step 5: fixed E_gap -> scan epsilon_h ========"
python3 utils/plot_band_gap_fixed_gap_eh_sweep.py --all-gaps --all-mediators

echo "======== Step 5: fixed epsilon_h=3.8 -> scan E_gap ========"
python3 utils/plot_band_gap_fixed_eh_gap_sweep.py --eh 3.8 --all-mediators

echo "======== 2D sensitivity heatmaps (m_chi = 1, 10, 100 MeV) ========"
bash utils/run_band_gap_2d_heatmaps.sh

echo "======== Done ========"
echo "Figures under: outplots/band_gap_pheno/"
find outplots/band_gap_pheno -name '*.pdf' | sort
