#!/usr/bin/env python3
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: audit_qedark_ul_vs_pydme.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  audit_qedark_ul_vs_pydme.py -- Audit CCDarkSens full-grid UL vs DAMIC-M
#  reference and pydme conventions.
# ============================================================================

"""Audit CCDarkSens full-grid UL vs DAMIC-M reference and pydme conventions."""
from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

# Reference file from the pydme package (optional — set PYDME_REF_DIR to enable).
# If unset, comparison panels are skipped gracefully.
_pydme_ref_dir = os.environ.get("PYDME_REF_DIR", "")
PAPER_EXPORT = (
    Path(_pydme_ref_dir)
    / "ScienceRun2024-figures/data/ScienceRun2024_results-1/"
      "DAMIC-M_2025_QEDark_DMe_heavymediator.txt"
) if _pydme_ref_dir else None


# ----------------------------------------------------------------------------
# load_paper_ref
#   Read the reference limit curve (mass, sigma) from a text/CSV file, converting masses above 1e4 from eV to MeV, and return it sorted by mass.
# ----------------------------------------------------------------------------
def load_paper_ref(path: Path):
    rows = []
    for ln in path.read_text().splitlines():
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        parts = ln.replace("\t", ",").split(",")
        if parts[0].lower().startswith("mass"):
            continue
        m_raw, s = float(parts[0]), float(parts[1])
        m_mev = m_raw / 1e6 if m_raw > 1e4 else m_raw
        rows.append((m_mev, s))
    arr = np.asarray(rows)
    return arr[np.argsort(arr[:, 0])]


# ----------------------------------------------------------------------------
# load_ul
#   Open a scan ROOT file with uproot and return the upper-limit curve (graph if present, else the histogram bin centres), q0, D_pat and the q histogram, plus the open file.
# ----------------------------------------------------------------------------
def load_ul(root_path: Path):
    import uproot

    f = uproot.open(root_path)
    if "upper_limit_sigma_e_mchi_graph" in f:
        xc, v = f["upper_limit_sigma_e_mchi_graph"].values()
        xc = np.asarray(xc, dtype=float)
        v = np.asarray(v, dtype=float)
    else:
        ul = f["upper_limit_sigma_e_mchi"]
        xe = ul.axis().edges()
        xc = 0.5 * (xe[:-1] + xe[1:])
        v = ul.values()
    q0 = f["q0_mchi"].values() if "q0_mchi" in f else None
    d_pat = f["D_pat"].values() if "D_pat" in f else None
    hq = f.get("q_mchi_sigma_pattern")
    return xc, v, q0, d_pat, hq, f


def ul_from_qhist(hq, ix: int, target_q: float) -> float:
    """Same crossing logic as ccdarksens_plot_dmelectron_limit.cc."""
    ye = hq.axis(1).edges()
    yc = 0.5 * (ye[:-1] + ye[1:])
    z = hq.values()[ix, :]
    for iy in range(len(yc) - 1):
        q1, q2 = z[iy], z[iy + 1]
        s1, s2 = yc[iy], yc[iy + 1]
        if q1 < target_q <= q2 and q2 > q1 and s1 > 0 and s2 > 0:
            log_s1, log_s2 = np.log10(s1), np.log10(s2)
            t = (target_q - q1) / (q2 - q1)
            return 10 ** (log_s1 + t * (log_s2 - log_s1))
    return np.nan


def constraint_pydme_sum(theta: float, br_row: float, n_gamma: int, strength: float = 98.0) -> float:
    """One pattern dataset: sum over gamma of -t*Br + strength*ln(t*Br), Br = br_row/len(gamma)."""
    br = br_row / n_gamma
    tbr = theta * br
    if tbr <= 0:
        return np.inf
    term = -tbr + strength * np.log(tbr)
    return n_gamma * term


def constraint_ccdarksens_tau(theta: float, bp, br, strength: float = 98.0) -> float:
    """CCDarkSens tau-weighted branch (n_bins=1), summed over pattern bins."""
    bp, br = np.asarray(bp), np.asarray(br)
    n = len(br)
    sbr = br.sum()
    if sbr <= 0:
        return np.inf
    tau = strength / sbr
    out = 0.0
    for i in range(n):
        x = theta * tau * br[i]
        if x <= 0:
            return np.inf
        nrc = strength * br[i] / sbr
        out += x - nrc * np.log(x)
    return out


# ----------------------------------------------------------------------------
# main
#   Compare the stored upper limit of a scan with the reference curve, print the median ratios over several mass ranges (and the same from the q histogram when available).
# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--root",
        default=str(ROOT / "outputs/scan_pattern_data_qedark_fullgrid/scan_dmelectron_pattern.root"),
    )
    ap.add_argument(
        "--ref",
        default=str(PAPER_EXPORT) if PAPER_EXPORT else "",
        help="Path to reference limit curve (DAMIC-M ScienceRun2024). "
             "Set PYDME_REF_DIR env var or pass explicitly. Skipped if empty.",
    )
    ap.add_argument(
        "--ref-this-work",
        action="store_true",
        help="Also compare against dense this_work CSV",
    )
    args = ap.parse_args()

    from scipy.stats import norm

    target_q = norm.ppf(0.9) ** 2

    xc, v, q0, d_pat, hq, _ = load_ul(Path(args.root))
    ref_path = Path(args.ref)
    if ref_path.suffix == ".csv" and "this_work" in ref_path.name:
        ref = np.loadtxt(ref_path, delimiter=",")
    else:
        ref = load_paper_ref(ref_path)
    rm, rs = ref[:, 0], ref[:, 1]

    m = (v > 0) & (v < 0.9e-26) & (xc >= 1)
    mc, sg = xc[m], v[m]
    ri = np.interp(mc, rm, rs)
    ratio = sg / ri

    print("=== Data ===")
    print("D_pat:", d_pat, " sum=", float(np.sum(d_pat)) if d_pat is not None else "n/a")
    print("target_q (90% CL):", target_q)

    print(f"\n=== Scan vs reference ({ref_path.name}) ===")
    print(f"median stored/paper (m>=1 MeV): {np.median(ratio):.4f}")
    print(f"median stored/paper (1.2-500 MeV): {np.median(ratio[(mc >= 1.2) & (mc <= 500)]):.4f}")
    print(f"median stored/paper (5-500 MeV): {np.median(ratio[(mc >= 5) & (mc <= 500)]):.4f}")

    if hq is not None:
        sg_q = np.array([ul_from_qhist(hq, int(np.argmin(np.abs(xc - m))), target_q) for m in mc])
        ratio_q = sg_q / ri
        print(f"median qhist/paper (1.2-500 MeV): {np.nanmedian(ratio_q[(mc >= 1.2) & (mc <= 500)]):.4f}")
        print(f"median qhist/stored: {np.nanmedian(sg_q / sg):.4f}")

        print("\n=== Stored UL vs q-map crossing (--from-qhist style) ===")
        for t in [1.0, 2.0, 10.0, 50.0, 100.0]:
            ix = int(np.argmin(np.abs(xc - t)))
            ul_st = v[ix] if v[ix] > 0 else np.nan
            ul_q = ul_from_qhist(hq, ix, target_q)
            ref_t = np.interp(xc[ix], rm, rs)
            if ul_q > 0:
                print(
                    f"  m~{xc[ix]:.1f} MeV  stored/pap={ul_st/ref_t:.3f}  "
                    f"qhist/pap={ul_q/ref_t:.3f}  qhist/stored={ul_st/ul_q:.3f}"
                )
            else:
                print(f"  m~{xc[ix]:.1f} MeV  no q crossing")

    if args.ref_this_work:
        tw = ROOT / "data/previous_limits/heavy_mediator/DAMIC-M_this_work_QEDark_hm.csv"
        if tw.exists():
            twref = np.loadtxt(tw, delimiter=",")
            rti = np.interp(mc, twref[:, 0], twref[:, 1])
            print(f"\n=== vs this_work CSV (not paper figure) ===")
            print(f"median stored/this_work (1.2-500 MeV): {np.median((sg / rti)[(mc >= 1.2) & (mc <= 500)]):.4f}")

    print("\n=== Constraint at theta=1 (Bp+Br template) ===")
    bp = [141.4, 0.111, 0.042, 0.019, 2.5e-5, 5.8e-5]
    br = [0.039, 0.039, 0.016, 0.052, 0.011, 0.035]
    for ng in [1, 100, 4450]:
        pydme_c = sum(constraint_pydme_sum(1.0, br[i], ng) for i in range(6))
        print(f"  pydme-style sum over 6 patterns, n_gamma={ng}: {pydme_c:.2f}")
    print(f"  CCDarkSens tau-weighted (n_bins=1): {constraint_ccdarksens_tau(1.0, bp, br):.2f}")

    print("\n=== Representative masses ===")
    print("mchi[MeV]  scan_UL   ref_UL   ratio   q0")
    for t in [2, 10, 50, 100, 200, 500, 1000]:
        ix = int(np.argmin(np.abs(mc - t)))
        qv = q0[m][ix] if q0 is not None else float("nan")
        print(f"{mc[ix]:8.1f}  {sg[ix]:.3e}  {ri[ix]:.3e}  {sg[ix]/ri[ix]:.3f}  {qv:.3f}")


if __name__ == "__main__":
    main()
