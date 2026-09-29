#!/usr/bin/env python3
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: compare_baxter_fullgrid.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  compare_baxter_fullgrid.py -- Compare Baxter fullgrid scan vs paper export
#  and v_E=263 baseline
# ============================================================================
"""Compare Baxter fullgrid scan vs paper export and vE=263 baseline."""
from __future__ import annotations

import argparse
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
_pydme_ref_dir = os.environ.get("PYDME_REF_DIR", "")
PAPER_EXPORT = (
    Path(_pydme_ref_dir)
    / "ScienceRun2024-figures/data/ScienceRun2024_results-1/"
      "DAMIC-M_2025_QEDark_DMe_heavymediator.txt"
) if _pydme_ref_dir else None


# ----------------------------------------------------------------------------
# load_paper
#   Read the paper-export limit curve, converting masses above 1e4 from eV to MeV; returns the sorted masses and cross sections.
# ----------------------------------------------------------------------------
def load_paper(path: Path):
    rows = []
    for ln in path.read_text().splitlines():
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        parts = ln.replace("\t", ",").split(",")
        if parts[0].lower().startswith("mass"):
            continue
        m = float(parts[0])
        rows.append((m / 1e6 if m > 1e4 else m, float(parts[1])))
    arr = np.asarray(rows)
    order = np.argsort(arr[:, 0])
    return arr[order, 0], arr[order, 1]


# ----------------------------------------------------------------------------
# load_scan
#   Upper-limit curve of a scan ROOT file (graph if present, else the histogram), keeping physical points: positive, below 0.9e-26 and at masses of at least 0.5 MeV.
# ----------------------------------------------------------------------------
def load_scan(root: Path):
    import uproot

    f = uproot.open(root)
    if "upper_limit_sigma_e_mchi_graph" in f:
        xc, v = f["upper_limit_sigma_e_mchi_graph"].values()
        xc = np.asarray(xc, dtype=float)
        v = np.asarray(v, dtype=float)
    else:
        ul = f["upper_limit_sigma_e_mchi"]
        xe = ul.axis().edges()
        xc = 0.5 * (xe[:-1] + xe[1:])
        v = ul.values()
    ok = (v > 0) & (v < 0.9e-26) & (xc >= 0.5)
    return xc[ok], v[ok]


# ----------------------------------------------------------------------------
# main
#   Compare the full-grid scan with the Baxter halo, the baseline scan and the paper curve and write the comparison plots to --outdir.
# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--baxter-root",
        default=str(ROOT / "outputs/scan_pattern_data_qedark_fullgrid_baxter/scan_dmelectron_pattern.root"),
    )
    ap.add_argument(
        "--baseline-root",
        default=str(ROOT / "outputs/scan_pattern_data_qedark_fullgrid/scan_dmelectron_pattern.root"),
    )
    ap.add_argument("--outdir", default=str(ROOT / "outplots/qedark_repro"))
    args = ap.parse_args()

    pm, ps = load_paper(PAPER_EXPORT)
    mb, sb = load_scan(Path(args.baxter_root))
    m0, s0 = load_scan(Path(args.baseline_root))

    rp = np.interp(mb, pm, ps)
    r_bax = sb / rp
    r_base = np.interp(mb, m0, s0) / rp

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("=== Baxter fullgrid vs paper export ===")
    for label, m, r in [
        ("m≈1 MeV", 1.0, None),
        ("m≈2 MeV", 2.0, None),
        ("m≈5 MeV", 5.0, None),
    ]:
        i = int(np.argmin(np.abs(mb - m)))
        print(
            f"  {label}: baxter/pap={r_bax[i]:.3f}  baseline/pap={r_base[i]:.3f}  "
            f"baxter/baseline={sb[i]/np.interp(mb[i],m0,s0):.3f}"
        )
    m12 = (mb >= 1.2) & (mb <= 500)
    m5 = (mb >= 5) & (mb <= 500)
    print(f"  median baxter/paper (1.2-500): {np.median(r_bax[m12]):.3f}")
    print(f"  median baseline/paper (1.2-500): {np.median(r_base[m12]):.3f}")
    print(f"  median baxter/paper (5-500): {np.median(r_bax[m5]):.3f}")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True, gridspec_kw={"height_ratios": [2, 1]})
    ax1.loglog(pm, ps, "k-", lw=2, label="paper export")
    ax1.loglog(m0, s0, color="#1f77b4", lw=1.2, alpha=0.7, label="fullgrid vE=263")
    ax1.loglog(mb, sb, color="#d62728", lw=2, label="fullgrid Baxter vE=253.7")
    ax1.set_ylabel(r"$\bar{\sigma}_e$ UL [cm$^2$]")
    ax1.set_title("Heavy QEdark — Baxter 2021 halo full grid")
    ax1.grid(True, which="both", alpha=0.3)
    ax1.legend(fontsize=8)

    ax2.semilogx(mb, r_base, color="#1f77b4", lw=1.2, label="baseline / paper")
    ax2.semilogx(mb, r_bax, color="#d62728", lw=1.5, label="Baxter / paper")
    ax2.axhline(1.0, color="k", lw=0.8, alpha=0.5)
    ax2.set_xlabel(r"$m_\chi$ [MeV]")
    ax2.set_ylabel("scan / paper")
    ax2.set_ylim(0, 2.5)
    ax2.grid(True, which="both", alpha=0.3)
    ax2.legend(fontsize=8)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        p = outdir / f"compare_fullgrid_baxter_vs_paper.{ext}"
        fig.savefig(p, dpi=150)
        print("wrote", p)

    csv = outdir / "compare_fullgrid_baxter_summary.csv"
    with csv.open("w") as f:
        f.write("m_mev,baxter_ul,baseline_ul,paper_ul,baxter_over_paper,baseline_over_paper\n")
        for i in range(len(mb)):
            pap = rp[i]
            base = np.interp(mb[i], m0, s0)
            f.write(
                f"{mb[i]},{sb[i]},{base},{pap},{sb[i]/pap},{base/pap}\n"
            )
    print("wrote", csv)


if __name__ == "__main__":
    main()
