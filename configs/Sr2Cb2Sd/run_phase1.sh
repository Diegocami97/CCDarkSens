#!/usr/bin/env bash
# Sr2Cb2Sd phase 1: Si_fast dielectric (equal footing), dense rates, limit scans.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"

# QCDark2 checkout (override if yours lives elsewhere)
export CCDARK_QCDARK2_DIR="${CCDARK_QCDARK2_DIR:-/Users/diegovenegasvargas/Documents/Software/QCDark2}"

echo "=== write JSON configs ==="
python3 configs/Sr2Cb2Sd/write_configs.py

echo "=== Step 1: Si_fast dielectrics (0.556, 0.603; 1.2 if missing) ==="
python3 utils/qcdark2_regenerate_epsilon.py \
  --template configs/qcdark2/Si_fast_scissor.in \
  --scissor 0.556 0.603 \
  2>&1 | tee data/qcdark2_epsilon/Si/scissor_sweep_Sr2Cb2Sd_0p556_0p603.log

if [[ ! -f data/qcdark2_epsilon/Si/Si_fast_gap1p2.h5 ]]; then
  python3 utils/qcdark2_regenerate_epsilon.py \
    --template configs/qcdark2/Si_fast_scissor.in \
    --scissor 1.2 \
    2>&1 | tee data/qcdark2_epsilon/Si/scissor_sweep_Sr2Cb2Sd_gap1p2.log
fi

echo "=== Step 2: p100K ionization tables ==="
python3 utils/build_p100K_scaled.py --band-gap-eV 0.556 --eh-pair-eV 2.06 \
  --scenario "Sr2Cb2Sd indirect Klein"
python3 utils/build_p100K_scaled.py --band-gap-eV 0.603 --eh-pair-eV 2.19 \
  --scenario "Sr2Cb2Sd direct Klein"
if [[ ! -f data/p100K_gap1p2_eh3p8.csv ]]; then
  python3 utils/build_p100K_scaled.py --band-gap-eV 1.2 --eh-pair-eV 3.8 \
    --scenario "Sr2Cb2Sd Si ref B-thresh"
fi

echo "=== Step 3: dense rate grids → data/Sr2Cb2Sd/rates/ ==="
for cfg in configs/Sr2Cb2Sd/qcdark2_rates_*.json; do
  echo "--- $cfg ---"
  python3 utils/qcdark2_generate_grid.py "$cfg"
done

echo "=== Step 4: limit scans → outputs/Sr2Cb2Sd/ ==="
for cfg in configs/Sr2Cb2Sd/scan_*_1kgy.json configs/Sr2Cb2Sd/scan_oscura_*.json; do
  echo "--- $cfg ---"
  build/ccdarksens_scan_dmelectron_pattern "$cfg"
done

echo "Done. Results under outputs/Sr2Cb2Sd/"
