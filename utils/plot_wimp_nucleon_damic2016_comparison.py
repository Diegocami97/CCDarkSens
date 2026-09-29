#!/usr/bin/env python3
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: plot_wimp_nucleon_damic2016_comparison.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: plot_wimp_nucleon_damic2016_comparison.py
#  Overlays this repo's WIMP-nucleon SI exclusion reproduction against the
#  real DAMIC 2016 (0.6 kg-day, PhysRevD.94.082006 / arXiv:1607.07410)
#  literature curves. Styled to match ccdarksens_plot_dmelectron_limit.cc's
#  actual gStyle/TCanvas settings (inward ticks all sides, thick lines,
#  borderless legend) -- that app has no flag to overlay an arbitrary
#  literature CSV, so this is a standalone matplotlib companion, not a
#  replacement.
#
#  Usage:
#    python3 utils/plot_wimp_nucleon_damic2016_comparison.py \
#      <ccdarksens_curve.csv> <out.png>
#
#  <ccdarksens_curve.csv> comes from:
#    ./build/ccdarksens_plot_dmelectron_limit <scan.root> "label" --migdal \
#      --batch --out-csv <ccdarksens_curve.csv> --out-pdf /tmp/_throwaway.pdf \
#      --out-root /tmp/_throwaway.root
#
#  See docs/ClusterFitMC_Design.md Sec. 6 for the full writeup and the
#  reasoning behind the literature-file provenance (data/previous_limits/WIMP/
#  DAMIC_2016_SNOLAB_{observed_radomir,expected_nominal,expected_lower,
#  expected_upper}.csv).
# ============================================================================

import csv
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

REPO_ROOT = Path(__file__).resolve().parents[1]
LITERATURE_DIR = REPO_ROOT / "data" / "previous_limits" / "WIMP"

# --- Mimic ccdarksens_plot_dmelectron_limit.cc's ROOT style ---
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Helvetica", "Arial", "DejaVu Sans"]
rcParams["axes.linewidth"] = 1.2
rcParams["xtick.direction"] = "in"
rcParams["ytick.direction"] = "in"
rcParams["xtick.top"] = True
rcParams["ytick.right"] = True
rcParams["xtick.minor.visible"] = True
rcParams["ytick.minor.visible"] = True
rcParams["xtick.major.size"] = 7
rcParams["ytick.major.size"] = 7
rcParams["xtick.minor.size"] = 4
rcParams["ytick.minor.size"] = 4
rcParams["axes.labelsize"] = 15
rcParams["xtick.labelsize"] = 12
rcParams["ytick.labelsize"] = 12
rcParams["legend.frameon"] = False
rcParams["legend.fontsize"] = 9.5


# ----------------------------------------------------------------------------
# load_two_col
#   Read a two-column comma-separated literature curve (comments skipped) as two arrays.
# ----------------------------------------------------------------------------
def load_two_col(path):
    m, s = [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            a, b = line.split(",")
            m.append(float(a))
            s.append(float(b))
    return np.array(m), np.array(s)


# ----------------------------------------------------------------------------
# load_ccdarksens_curve
#   Read the CCDarkSens limit CSV (mchi_MeV, sigma) and convert the mass from MeV to GeV.
# ----------------------------------------------------------------------------
def load_ccdarksens_curve(path):
    m, s = [], []
    with open(path) as f:
        r = csv.reader(f)
        for row in r:
            if not row or row[0].startswith("#") or row[0] == "mchi_MeV":
                continue
            m.append(float(row[0]) / 1000.0)  # MeV -> GeV
            s.append(float(row[1]))
    return np.array(m), np.array(s)


# ----------------------------------------------------------------------------
# main
#   Overlay the CCDarkSens WIMP-nucleon limit on the DAMIC 2016 observed and expected curves (points at the top edge of the scan grid are treated as not real limits) and save the figure. Usage: <program> <ccdarksens_curve.csv> <out.png>.
# ----------------------------------------------------------------------------
def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <ccdarksens_curve.csv> <out.png>", file=sys.stderr)
        return 1
    curve_path, out_path = sys.argv[1], sys.argv[2]

    ours_m, ours_s = load_ccdarksens_curve(curve_path)
    grid_edge = ours_s >= 0.9999 * ours_s.max()
    real = ~grid_edge

    m_rad, s_rad = load_two_col(LITERATURE_DIR / "DAMIC_2016_SNOLAB_observed_radomir.csv")
    m_nom, s_nom = load_two_col(LITERATURE_DIR / "DAMIC_2016_SNOLAB_expected_nominal.csv")
    m_lo, s_lo = load_two_col(LITERATURE_DIR / "DAMIC_2016_SNOLAB_expected_lower.csv")
    m_hi, s_hi = load_two_col(LITERATURE_DIR / "DAMIC_2016_SNOLAB_expected_upper.csv")

    fig, ax = plt.subplots(figsize=(9, 7))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    m_band = np.geomspace(max(m_lo.min(), m_hi.min()), min(m_lo.max(), m_hi.max()), 200)
    lo_i = np.exp(np.interp(np.log(m_band), np.log(m_lo), np.log(s_lo)))
    hi_i = np.exp(np.interp(np.log(m_band), np.log(m_hi), np.log(s_hi)))
    ax.fill_between(m_band, lo_i, hi_i, color="gray", alpha=0.30, lw=0,
                     label=r"DAMIC 2016 expected $\pm1\sigma$ band")
    ax.plot(m_nom, s_nom, color="dimgray", lw=1.6, ls=":",
            label="DAMIC 2016 expected median (nominal)")
    ax.plot(m_rad, s_rad, color="black", lw=3.0,
            label="DAMIC at SNOLAB 2016, observed 90% CL\n(radomir's curve)")
    ax.plot(ours_m[real], ours_s[real], color=(0.90, 0.10, 0.10), lw=3.0, ls="--",
            label="CCDarkSens repro (this work)")
    if grid_edge.any():
        ax.scatter(ours_m[grid_edge], ours_s[grid_edge], color=(0.90, 0.10, 0.10),
                   marker="x", s=45, lw=2.0, label="grid edge (no sensitivity found)")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(max(0.45, ours_m.min() * 0.9), ours_m.max() * 1.1)
    ax.set_ylim(1e-41, max(1e-33, ours_s[real].max() * 3 if real.any() else 1e-33))
    ax.set_xlabel(r"$m_{\chi}$  [GeV/$c^2$]")
    ax.set_ylabel(r"$\bar{\sigma}_n$  [cm$^2$]")
    ax.set_title("WIMP-nucleon SI exclusion — DAMIC 2016 (0.6 kg-day) reproduction", fontsize=13, pad=12)
    ax.legend(loc="upper right")
    fig.subplots_adjust(left=0.13, right=0.96, bottom=0.11, top=0.92)
    fig.savefig(out_path, dpi=170)
    print(f"saved {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
