#!/usr/bin/env bash
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: run_band_gap_light_phase_c_batch.sh
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  run_band_gap_light_phase_c_batch.sh -- I run the full Phase C of the band-
#  gap phenomenology with the light mediator (QCDark2): generate the light
#  rate configs and the 12 light scan configs, then for every gap the dR/dE
#  grids and both scans (B-thresh and D-equal), and finally overlay the two
#  limit curves.
# ============================================================================

# Phase C full run for band-gap pheno with light mediator (QCDark2).
#
# This script:
#  1) generates light rate configs (qcdark2_generate_Si_light_gap*.json),
#  2) generates 12 light pheno scan configs (scan_band_gap_light_pheno_*.json),
#  3) for each gap: generates light dR/dE grids, then runs the two scans (B-thresh and D-equal),
#  4) finally overlays limit curves for B-thresh and D-equal.

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

source "${ROOT_INSTALL:-}/bin/thisroot.sh" 2>/dev/null || true

LOGDIR="$ROOT/outputs/phase_c_logs"

# Set CCDARK_QCDARK2_DIR to the local QCDark2 checkout.
QCDARK2_DIR="${CCDARK_QCDARK2_DIR:-}"
QCDARK2_VENV_PY="${QCDARK2_DIR:+$QCDARK2_DIR/.venv/bin/python3}"
QCDARK2_VENV_PY="${QCDARK2_VENV_PY:-python3}"

echo "[phase-c-light] generating configs..."
python3 utils/gen_qcdark2_generate_Si_light_gap_configs.py
python3 utils/gen_band_gap_light_pheno_scan_configs.py

BUILD_SCAN="$ROOT/build/ccdarksens_scan_dmelectron_pattern"
BUILD_PLOT="$ROOT/build/ccdarksens_plot_dmelectron_limit"

if [[ ! -f "$BUILD_SCAN" || ! -f "$BUILD_PLOT" ]]; then
  echo "[phase-c-light] building binaries..."
  cmake --build build -j --target ccdarksens_scan_dmelectron_pattern ccdarksens_plot_dmelectron_limit
fi

declare -a GAP_SHORTS=("0p1" "0p3" "0p5" "0p7" "0p9" "1p2")
declare -a GAP_VALS=("0.1" "0.3" "0.5" "0.7" "0.9" "1.2")

GAP_TO_TAG="gap" # tag = gap${gap_short}

# ----------------------------------------------------------------------------
# run_scan
#   Run the scan binary on a config and keep the log (arguments: config, tier label, log file).
# ----------------------------------------------------------------------------
run_scan() {
  local cfg="$1"
  local tierlabel="$2"
  local log="$3"
  echo "[phase-c-light] scan: $tierlabel cfg=$(basename "$cfg")"
  "$BUILD_SCAN" "$cfg" 2>&1 | tee "$log"
}

# ----------------------------------------------------------------------------
# run_rates
#   Generate a QCDark2 rate grid from a config with the QCDark2 Python environment and keep the log (arguments: config, log file).
# ----------------------------------------------------------------------------
run_rates() {
  local cfg="$1"
  local log="$2"
  echo "[phase-c-light] rates: $(basename "$cfg")"
  "$QCDARK2_VENV_PY" utils/qcdark2_generate_grid.py "$cfg" 2>&1 | tee "$log"
}

for idx in "${!GAP_SHORTS[@]}"; do
  gs="${GAP_SHORTS[$idx]}"
  gv="${GAP_VALS[$idx]}"
  tag="${GAP_TO_TAG}${gs}" # gap0p7

  rates_cfg="configs/qcdark2_generate_Si_light_${tag}.json"
  rates_log="$LOGDIR/rates_${tag}.log"

  # Generate dR/dE tables for this scissor gap & mediator.
  run_rates "$rates_cfg" "$rates_log"

  # B-thresh: eh = 3.8 -> 3p8
  cfg_b="configs/scan_band_gap_light_pheno_${gs}_eh3p8.json"
  run_scan "$cfg_b" "${gv} eV B-thresh" "$LOGDIR/scan_B_${tag}.log"

  # D-equal: eh = gap -> eh tag = gs
  cfg_d="configs/scan_band_gap_light_pheno_${gs}_eh${gs}.json"
  run_scan "$cfg_d" "${gv} eV D-equal" "$LOGDIR/scan_D_${tag}.log"
done

echo "[phase-c-light] plotting limit overlays (all mediators + tiers)..."
python3 utils/run_band_gap_phase_c.py plot-limits --all

echo "[phase-c-light] done."

