# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: band_gap_klein_plot.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  band_gap_klein_plot.py -- Shared Klein-tier helpers for band-gap
#  ne_imaging figure scripts.
# ============================================================================

"""Shared Klein-tier helpers for band-gap ne_imaging figure scripts."""
from __future__ import annotations

DEFAULT_KLEIN_GAPS = [0.1, 0.3, 0.5, 0.7, 0.9]
SI_EGAP_EV = 1.2
SI_EH_EV = 3.8
COLOR_SI_REF = "#D4A017"
KLEIN_GAP_PALETTE = [
    "#C0392B",
    "#D35400",
    "#7D3C98",
    "#2471A3",
    "#148F77",
    "#884EA0",
    "#1ABC9C",
]


# ----------------------------------------------------------------------------
# klein_eh
#   Klein-formula electron-hole pair energy for a band gap: eps_h = 2.8*E_gap + 0.5 eV, rounded to 2 decimals.
# ----------------------------------------------------------------------------
def klein_eh(gap_ev: float) -> float:
    return round(2.8 * gap_ev + 0.5, 2)


# ----------------------------------------------------------------------------
# parse_float_list
#   Parse a comma-separated string into floats (blank items skipped).
# ----------------------------------------------------------------------------
def parse_float_list(csv: str) -> list[float]:
    return [float(tok.strip()) for tok in csv.split(",") if tok.strip()]


def filter_klein_gaps(gaps: list[float]) -> list[float]:
    """Drop E_gap=1.2 — covered by the Si reference row."""
    out: list[float] = []
    for g in gaps:
        if abs(g - SI_EGAP_EV) < 1e-6:
            print(f"WARN: skipping E_gap={g:g} eV (same as Si ref; use Si row)")
            continue
        out.append(g)
    return out


# ----------------------------------------------------------------------------
# gap_color
#   Colour for the index-th band gap, cycling through the palette.
# ----------------------------------------------------------------------------
def gap_color(index: int) -> str:
    return KLEIN_GAP_PALETTE[index % len(KLEIN_GAP_PALETTE)]


# ----------------------------------------------------------------------------
# eh_label
#   epsilon_h formatted for a label: two decimals with trailing zeros and a trailing dot removed.
# ----------------------------------------------------------------------------
def eh_label(eh: float) -> str:
    return f"{eh:.2f}".rstrip("0").rstrip(".")


# ----------------------------------------------------------------------------
# klein_mathtext_label
#   Matplotlib math-text label "(E_gap, eps_h) = (gap, eh)"; eh defaults to the Klein value for that gap.
# ----------------------------------------------------------------------------
def klein_mathtext_label(gap_ev: float, eh: float | None = None) -> str:
    eh_v = klein_eh(gap_ev) if eh is None else eh
    return rf"$(E_{{\mathrm{{gap}}}},\ \varepsilon_h)=({gap_ev:g},\ {eh_label(eh_v)})$"


def gap_scan_id(gap_ev: float) -> str:
    """Subdir name under band_gap_ne_imaging_klein/, e.g. klein_gap0p1."""
    return "klein_" + f"gap{gap_ev:.1f}".replace(".", "p")


def er_file_tag(er_ev: float) -> str:
    """Filename token, e.g. Er4, Er2, Er3p5."""
    return "Er" + f"{er_ev:g}".replace(".", "p")
