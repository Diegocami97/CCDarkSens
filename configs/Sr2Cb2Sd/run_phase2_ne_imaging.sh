#!/usr/bin/env bash
# Sr2Cb2Sd phase 2: n_e imaging one-point scans + signal spectrum figures.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"

echo "=== write ne_imaging JSON configs ==="
python3 configs/Sr2Cb2Sd/write_ne_imaging_configs.py

echo "=== one-point scans (10 per mediator batch; 20 total) ==="
for cfg in configs/Sr2Cb2Sd/ne_imaging_*.json; do
  echo "--- $cfg ---"
  build/ccdarksens_scan_dmelectron_pattern "$cfg"
done

echo "=== plot n_e signal spectra ==="
python3 utils/plot_Sr2Cb2Sd_ne_spectra.py

echo "Done. Figures under outplots/Sr2Cb2Sd/ne_signal_spectrum_{heavy,light}.pdf"
