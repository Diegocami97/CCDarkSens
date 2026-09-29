# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: band_gap_scan_paths.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  band_gap_scan_paths.py -- Resolve band-gap pheno / 2D-grid scan ROOT file
#  paths
# ============================================================================
"""Resolve band-gap pheno / 2D-grid scan ROOT paths."""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# Same grids as 2D heatmaps / Phase C pheno.
GAP_GRID = [0.1, 0.3, 0.5, 0.7, 0.9, 1.2]
EH_GRID = [0.5, 1.0, 1.5, 2.0, 2.5, 3.8]

# Si reference cell (heatmap star): used as overlay on limit sweeps.
SI_REF_GAP_EV = 1.2
SI_REF_EH_EV = 3.8


# ----------------------------------------------------------------------------
# si_reference_root
#   Path of the silicon-reference scan ROOT file for a mediator.
# ----------------------------------------------------------------------------
def si_reference_root(mediator: str) -> Path:
    return scan_root_path(mediator, SI_REF_GAP_EV, SI_REF_EH_EV)


# ----------------------------------------------------------------------------
# ev_tag
#   Energy formatted for file names: one decimal with '.' replaced by 'p' (1.2 -> "1p2").
# ----------------------------------------------------------------------------
def ev_tag(x: float) -> str:
    return f"{x:.1f}".replace(".", "p")


def scan_root_path(mediator: str, gap: float, eh: float) -> Path:
    """Path to scan_dmelectron_pattern.root for (mediator, E_gap, epsilon_h)."""
    gtag = ev_tag(gap)
    etag = ev_tag(eh)

    if abs(eh - 3.8) < 1e-12:
        if mediator == "heavy":
            return ROOT / f"outputs/scan_band_gap_{gtag}_eh3p8/scan_dmelectron_pattern.root"
        return ROOT / f"outputs/scan_band_gap_light_{gtag}_eh3p8/scan_dmelectron_pattern.root"

    if abs(eh - gap) < 1e-12 and gap in {0.5, 0.7, 0.9, 1.2}:
        if mediator == "heavy":
            return ROOT / f"outputs/scan_band_gap_{gtag}_eh{etag}/scan_dmelectron_pattern.root"
        return ROOT / f"outputs/scan_band_gap_light_{gtag}_eh{etag}/scan_dmelectron_pattern.root"

    return ROOT / f"outputs/scan_band_gap_2d_{mediator}_{gtag}_eh{etag}/scan_dmelectron_pattern.root"
