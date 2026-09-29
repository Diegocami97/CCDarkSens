#!/usr/bin/env bash
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: run_band_gap_phase_c_batch.sh
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  run_band_gap_phase_c_batch.sh -- I run the Phase C production scans one
#  after the other; the logs go to outputs/phase_c_logs/.
# ============================================================================

# Run Phase C production scans sequentially. Logs under outputs/phase_c_logs/
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
LOGDIR="$ROOT/outputs/phase_c_logs"
mkdir -p "$LOGDIR"

source "${ROOT_INSTALL:-}/bin/thisroot.sh" 2>/dev/null || true

# ----------------------------------------------------------------------------
# run_tier
#   Run the Phase C scans of one tier (B-thresh or D-equal) and keep the log.
# ----------------------------------------------------------------------------
run_tier() {
  local tier="$1"
  echo "======== Phase C tier: $tier ========"
  python3 utils/run_band_gap_phase_c.py scan --tier "$tier" 2>&1 | tee "$LOGDIR/scan_${tier}.log"
}

run_tier "B-thresh"
run_tier "D-equal"

echo "======== Phase C limit plots (heavy only) ========"
python3 utils/run_band_gap_phase_c.py plot-limits --tier B-thresh 2>&1 | tee "$LOGDIR/plot_B-thresh.log"
python3 utils/run_band_gap_phase_c.py plot-limits --tier D-equal 2>&1 | tee "$LOGDIR/plot_D-equal.log"

echo "Phase C batch complete."
