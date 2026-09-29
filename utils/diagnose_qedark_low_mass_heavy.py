#!/usr/bin/env python3
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: diagnose_qedark_low_mass_heavy.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  diagnose_qedark_low_mass_heavy.py -- Diagnose low-mass heavy-mediator
#  QEdark upper-limit deviations vs DAMIC-M
# ============================================================================
"""Low-mass heavy-mediator QEdark UL diagnostic vs DAMIC-M reference."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


# ----------------------------------------------------------------------------
# load_scan
#   Upper-limit histogram of a scan ROOT file as (mass centres, limits) plus the q0 histogram when it exists.
# ----------------------------------------------------------------------------
def load_scan(root_path: Path):
    import uproot

    f = uproot.open(root_path)
    h = f["upper_limit_sigma_e_mchi"]
    m = h.axis(0).centers()
    s = h.values()
    q0 = f["q0_mchi"].values() if "q0_mchi" in f else None
    return m, s, q0


# ----------------------------------------------------------------------------
# main
#   Diagnose the low-mass heavy-mediator QEDark limit: compare the scan result with the reference curve and write the diagnostic plots to --outdir.
# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--root",
        default=str(ROOT / "outputs/scan_pattern_data_qedark_fullgrid/scan_dmelectron_pattern.root"),
    )
    ap.add_argument(
        "--ref",
        default=str(ROOT / "data/previous_limits/heavy_mediator/DAMIC-M_this_work_QEDark_hm.csv"),
    )
    ap.add_argument(
        "--outdir",
        default=str(ROOT / "outplots/qedark_repro"),
    )
    args = ap.parse_args()

    m, s, q0 = load_scan(Path(args.root))
    ref = np.loadtxt(args.ref, delimiter=",")
    rm, rs = ref[:, 0], ref[:, 1]

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Valid comparison window (reference starts at ~1.03 MeV)
    valid = (m >= rm.min()) & (s > 0) & (s < 0.99e-26)
    mv, sv = m[valid], s[valid]
    rv = np.interp(mv, rm, rs)
    ratio = sv / rv

    pinned = (s >= 0.99e-26) & (m < rm.min())
    print("=== Low-mass heavy QEdark diagnostic ===")
    print(f"Reference mass range: {rm.min():.4f} – {rm.max():.1f} MeV")
    print(f"Grid-pinned (m < {rm.min():.2f} MeV): {pinned.sum()} masses (ignore in physics comparison)")
    for lo, hi in [(1, 5), (5, 10), (10, 20), (20, 50), (50, 500)]:
        mask = (mv >= lo) & (mv <= hi)
        if mask.any():
            print(f"  {lo:3d}–{hi:3d} MeV: median scan/ref = {np.median(ratio[mask]):.3f}  (n={mask.sum()})")

    print("\nRepresentative masses (interpolated ref):")
    for t in [1.03, 2, 3, 5, 7, 10, 14, 20, 50, 100, 200]:
        si = np.interp(t, m, s)
        ri = np.interp(t, rm, rs)
        print(f"  m={t:6.2f} MeV  scan={si:.3e}  ref={ri:.3e}  ratio={si/ri:.3f}")

    # Implied signal scale ~ ratio (UL linear in sigma at fixed q)
    print("\nImplied signal overshoot (1/ratio) if stats match:")
    for t in [1.03, 2, 5, 10]:
        r = np.interp(t, mv, ratio) if t >= mv.min() else np.nan
        if r > 0:
            print(f"  m={t:.2f} MeV: scan signal ~ {1/r:.2f}× reference")

    # Plot 1: UL curves
    fig, axes = plt.subplots(2, 1, figsize=(8, 8), sharex=True)

    ax = axes[0]
    ax.loglog(mv, sv, "b-", lw=1.2, label="CCDarkSens scan")
    ax.loglog(rm, rs, "k--", lw=1.5, label="Reference (QEDark hm)")
    ax.axvline(rm.min(), color="gray", ls=":", lw=1, label=f"ref start ({rm.min():.2f} MeV)")
    ax.set_ylabel(r"90% CL $\sigma_e$ [cm$^2$]")
    ax.set_title("Heavy mediator: low-mass UL comparison")
    ax.legend(loc="best")
    ax.grid(True, which="both", alpha=0.3)

    ax = axes[1]
    ax.semilogx(mv, ratio, "r-", lw=1.2)
    ax.axhline(1.0, color="k", ls="--", lw=1)
    ax.set_xlabel(r"$m_\chi$ [MeV]")
    ax.set_ylabel("scan / reference")
    ax.set_ylim(0, 1.5)
    ax.grid(True, which="both", alpha=0.3)

    fig.tight_layout()
    pdf = outdir / "low_mass_heavy_comparison.pdf"
    fig.savefig(pdf)
    plt.close(fig)

    # Plot 2: ratio at reference mass points only
    fig2, ax2 = plt.subplots(figsize=(8, 4))
    scan_at_ref = np.interp(rm, m, s)
    r_at_ref = scan_at_ref / rs
    ax2.semilogx(rm, r_at_ref, "o-", ms=3, lw=1)
    ax2.axhline(1.0, color="k", ls="--")
    ax2.set_xlabel(r"$m_\chi$ [MeV]")
    ax2.set_ylabel("scan / reference")
    ax2.set_title("Ratio at reference CSV mass points")
    ax2.set_ylim(0, 1.5)
    ax2.grid(True, which="both", alpha=0.3)
    fig2.tight_layout()
    pdf2 = outdir / "low_mass_heavy_ratio_at_ref_masses.pdf"
    fig2.savefig(pdf2)
    plt.close(fig2)

    csv_out = outdir / "low_mass_heavy_diagnostic.csv"
    np.savetxt(
        csv_out,
        np.column_stack([mv, sv, rv, ratio]),
        delimiter=",",
        header="mchi_MeV,scan_ul,ref_ul_interp,ratio",
        comments="",
    )
    print(f"\nWrote {pdf}, {pdf2}, {csv_out}")


if __name__ == "__main__":
    main()
