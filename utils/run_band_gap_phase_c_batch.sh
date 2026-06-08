#!/usr/bin/env bash
# Run Phase C production scans sequentially. Logs under outputs/phase_c_logs/
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
LOGDIR="$ROOT/outputs/phase_c_logs"
mkdir -p "$LOGDIR"

source "${ROOT_INSTALL:-}/bin/thisroot.sh" 2>/dev/null || true

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
