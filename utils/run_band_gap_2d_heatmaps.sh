#!/usr/bin/env bash
# Extract sigma_UL on the 2D (E_gap, epsilon_h) grid and plot gain heatmaps.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

MCHI_LIST="${MCHI_LIST:-1 10 100}"

for m in $MCHI_LIST; do
  for med in heavy light; do
    echo "[2d] extract $med m_chi=${m} MeV"
    python3 utils/extract_band_gap_2d_sensitivity.py --mediator "$med" --mchi-MeV "$m"
  done
done

echo "[2d] plot heatmaps"
python3 utils/plot_band_gap_2d_heatmap.py --mchi-MeV $MCHI_LIST

echo "[2d] done -> outplots/band_gap_pheno/2d_surface/heatmap_*_mchi*MeV.pdf"
