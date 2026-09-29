#!/usr/bin/env bash
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: run_band_gap_2d_batch.sh
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  run_band_gap_2d_batch.sh -- I run the 2D (E_gap, eh) band-gap scans for
#  the heavy and then the light mediator, skipping scans whose ROOT outputs
#  already exist.
# ============================================================================

# Run 2D (gap,eh) scans for heavy then light mediator.
# Skips scans whose ROOT outputs already exist.

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

source "${ROOT_INSTALL:-}/bin/thisroot.sh" 2>/dev/null || true

LOGDIR="$ROOT/outputs/phase_c_logs"
mkdir -p "$LOGDIR"
LOG="$LOGDIR/scan_2d_grid.log"
SUMMARY="$LOGDIR/scan_2d_grid_summary.tsv"

BUILD_SCAN="$ROOT/build/ccdarksens_scan_dmelectron_pattern"
if [[ ! -f "$BUILD_SCAN" ]]; then
  echo "[build] ccdarksens_scan_dmelectron_pattern"
  cmake --build build -j --target ccdarksens_scan_dmelectron_pattern
fi

echo -e "mediator\tgap\teh\tstatus\tconfig\toutroot" > "$SUMMARY"

# ----------------------------------------------------------------------------
# run_one
#   Run one 2D-grid scan (mediator, gap, eps_h): skip it if its ROOT output already exists or its config is missing, otherwise run the scan binary; every outcome is appended to the summary file.
# ----------------------------------------------------------------------------
run_one() {
  local mediator="$1"
  local gap="$2"
  local eh="$3"
  local gtag="${gap/./p}"   # 0.1 -> 0p1
  local etag="${eh/./p}"    # 1.0 -> 1p0 (passed already with one decimal)

  local cfg="configs/scan_band_gap_2d_${mediator}_${gtag}_eh${etag}.json"
  local out="outputs/scan_band_gap_2d_${mediator}_${gtag}_eh${etag}/scan_dmelectron_pattern.root"

  if [[ -f "$out" ]]; then
    echo "[skip] ${mediator} gap=${gap} eh=${eh} (exists)"
    echo -e "${mediator}\t${gap}\t${eh}\tskipped_exists\t${cfg}\t${out}" >> "$SUMMARY"
    return 0
  fi
  if [[ ! -f "$cfg" ]]; then
    echo "[missing-config] ${cfg}"
    echo -e "${mediator}\t${gap}\t${eh}\tmissing_config\t${cfg}\t${out}" >> "$SUMMARY"
    return 0
  fi

  echo "[run] ${mediator} gap=${gap} eh=${eh}"
  if "$BUILD_SCAN" "$cfg" 2>&1 | tee -a "$LOG"; then
    echo -e "${mediator}\t${gap}\t${eh}\tsuccess\t${cfg}\t${out}" >> "$SUMMARY"
  else
    echo -e "${mediator}\t${gap}\t${eh}\tfailed\t${cfg}\t${out}" >> "$SUMMARY"
  fi
}

# Fixed grids (eh as string preserves 1.0 formatting in tags).
GAPS=(0.1 0.3 0.5 0.7 0.9 1.2)
EHS=(0.5 1.0 1.5 2.0 2.5 3.8)

for mediator in heavy light; do
  echo "======== mediator=${mediator} ========" | tee -a "$LOG"
  for gap in "${GAPS[@]}"; do
    for eh in "${EHS[@]}"; do
      # valid only if eh >= gap
      awk "BEGIN{exit !($eh >= $gap)}" || continue
      # skip existing complete sets (B-thresh row + on-grid D-equal points)
      awk "BEGIN{exit !($eh == 3.8)}" && continue
      if [[ "$gap" == "$eh" ]]; then
        if [[ "$gap" == "0.5" || "$gap" == "0.7" || "$gap" == "0.9" || "$gap" == "1.2" ]]; then
          continue
        fi
      fi
      run_one "$mediator" "$gap" "$eh"
    done
  done
done

echo "[done] 2D batch scan complete. Summary: $SUMMARY"

