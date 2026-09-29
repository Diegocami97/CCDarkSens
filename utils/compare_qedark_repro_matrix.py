#!/usr/bin/env python3
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: compare_qedark_repro_matrix.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  compare_qedark_repro_matrix.py -- Compare QEdark reproduction matrix few-
#  mass scans vs DAMIC-M reference
# ============================================================================
"""Compare qedark repro matrix few-mass scans vs DAMIC-M reference curves."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]

RUNS = [
    ("tau + config exp", "outputs/qedark_repro_matrix/fewmass_tau_config_exp/scan_dmelectron_pattern.root", "heavy"),
    ("pydme multibin + config exp", "outputs/qedark_repro_matrix/fewmass_pydme_multibin_config_exp/scan_dmelectron_pattern.root", "heavy"),
    ("tau + data exp (wrong B scale)", "outputs/qedark_repro_matrix/fewmass_tau_data_exp/scan_dmelectron_pattern.root", "heavy"),
    ("light tau + config exp", "outputs/qedark_repro_matrix/fewmass_light_tau_config_exp/scan_dmelectron_pattern.root", "light"),
]

REFS = {
    "heavy": ROOT / "data/previous_limits/heavy_mediator/DAMIC-M_this_work_QEDark_hm.csv",
    "light": ROOT / "data/previous_limits/light_mediator/DAMIC-M_this_work_QEDark_ulm_.csv",
}


# ----------------------------------------------------------------------------
# load_ref
#   Read a reference limit curve (two comma-separated columns) as an array.
# ----------------------------------------------------------------------------
def load_ref(path: Path):
    return np.genfromtxt(path, delimiter=",")


# ----------------------------------------------------------------------------
# load_ul
#   Upper-limit histogram of a scan ROOT file as (mass centres, limits), keeping the physical entries (positive and below 0.9e-26).
# ----------------------------------------------------------------------------
def load_ul(root_path: Path):
    import uproot

    f = uproot.open(root_path)
    ul = f["upper_limit_sigma_e_mchi"]
    xe = ul.axis().edges()
    xc = 0.5 * (xe[:-1] + xe[1:])
    v = ul.values()
    m = (v > 0) & (v < 0.9e-26)
    return xc[m], v[m]


# ----------------------------------------------------------------------------
# main
#   Plot the heavy- and light-mediator QEDark reproduction limits together with the reference curves, print the comparison and save the figure and table to --outdir.
# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=str(ROOT / "outplots/qedark_repro"))
    args = ap.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    rows = []
    fig, ax = plt.subplots(figsize=(9, 6))

    for med in ("heavy", "light"):
        ref = load_ref(REFS[med])
        rm, rs = ref[:, 0], ref[:, 1]
        ax.plot(rm, rs, "k-" if med == "heavy" else "k--", lw=2, label=f"ref ({med})")

    colors = ["C0", "C1", "C2", "C3"]
    for (label, rpath, med), col in zip(RUNS, colors):
        p = ROOT / rpath
        if not p.exists():
            print("skip missing", p)
            continue
        mc, sg = load_ul(p)
        ref = load_ref(REFS[med])
        rm, rs = ref[:, 0], ref[:, 1]
        ri = np.interp(mc, rm, rs)
        ratio = sg / ri
        med_ratio = np.median(ratio[(mc >= 5) & (mc <= 1000)])
        rows.append((label, med_ratio, np.median(ratio)))
        ax.plot(mc, sg, "-", color=col, lw=1.2, label=f"{label} (med={med_ratio:.3f})")
        print(f"{label:30s}  median(5-1000)={med_ratio:.4f}  median(all)={np.median(ratio):.4f}")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$m_\chi$ [MeV]")
    ax.set_ylabel(r"$\sigma_e$ [cm$^2$] (90% CL)")
    ax.set_title("QEdark repro matrix (few-mass) vs DAMIC-M 2025")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    pdf = outdir / "compare_repro_matrix_fewmass.pdf"
    fig.savefig(pdf, dpi=150)
    print("wrote", pdf)

    # CSV summary
    csv = outdir / "repro_matrix_summary.csv"
    with csv.open("w") as f:
        f.write("variant,median_ratio_5_1000_MeV,median_ratio_all_m\n")
        for label, mr, ma in rows:
            f.write(f"{label},{mr},{ma}\n")
    print("wrote", csv)


if __name__ == "__main__":
    main()
