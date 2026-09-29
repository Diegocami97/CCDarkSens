# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: band_gap_limit_plot.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  band_gap_limit_plot.py -- Helper functions for band-gap limit overlay
#  figures
# ============================================================================
"""Helpers for band-gap limit overlay plots (ccdarksens_plot_dmelectron_limit)."""

from __future__ import annotations

from pathlib import Path

from band_gap_plot_labels import si_reference_label_root
from band_gap_scan_paths import SI_REF_GAP_EV, SI_REF_EH_EV, si_reference_root

ROOT = Path(__file__).resolve().parents[1]


# ----------------------------------------------------------------------------
# is_si_reference_cell
#   True if (gap, eh) is the silicon reference point.
# ----------------------------------------------------------------------------
def is_si_reference_cell(gap_eV: float, eh_eV: float) -> bool:
    return abs(gap_eV - SI_REF_GAP_EV) < 1e-9 and abs(eh_eV - SI_REF_EH_EV) < 1e-9


def reference_scan_cli_args(mediator: str, *, enabled: bool = True) -> list[str]:
    """
    CLI tokens for --reference-scan (solid Si reference curve on limit plots).

    Skipped if the reference ROOT file is missing.
    """
    if not enabled:
        return []
    ref = si_reference_root(mediator)
    if not ref.is_file():
        print(f"WARN: Si reference scan missing: {ref}")
        return []
    rel = ref.relative_to(ROOT).as_posix()
    return ["--reference-scan", rel, si_reference_label_root()]
