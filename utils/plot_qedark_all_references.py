#!/usr/bin/env python3
"""Overlay CCDarkSens QEdark heavy pattern scan vs reference curves (paper export primary)."""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[1]

_pydme_ref_dir = os.environ.get("PYDME_REF_DIR", "")
PAPER_EXPORT = (
    Path(_pydme_ref_dir)
    / "ScienceRun2024-figures/data/ScienceRun2024_results-1/"
      "DAMIC-M_2025_QEDark_DMe_heavymediator.txt"
) if _pydme_ref_dir else None

TARGET_Q = norm.ppf(0.9) ** 2


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


def _load_ul_from_root(f):
    """Prefer exact-mass TGraph; fall back to TH1D bin centers."""
    if "upper_limit_sigma_e_mchi_graph" in f:
        xc, v = f["upper_limit_sigma_e_mchi_graph"].values()
        return np.asarray(xc, dtype=float), np.asarray(v, dtype=float)
    ul = f["upper_limit_sigma_e_mchi"]
    xe = ul.axis().edges()
    xc = 0.5 * (xe[:-1] + xe[1:])
    v = ul.values()
    return xc, v


def load_scan_ul(root_path: Path):
    import uproot

    f = uproot.open(root_path)
    xc, v = _load_ul_from_root(f)
    ok = (v > 0) & (v < 0.9e-26) & (xc >= 0.5)
    return xc[ok], v[ok]


def load_scan_qhist(root_path: Path, mc: np.ndarray):
    import uproot

    f = uproot.open(root_path)
    hq = f.get("q_mchi_sigma_pattern")
    if hq is None:
        return None, None
    # q-map rows follow full scan mass order; mc may be a filtered subset — map by mass.
    xc_full, _ = _load_ul_from_root(f)
    sig_q = []
    for mchi in mc:
        ix = int(np.argmin(np.abs(xc_full - mchi)))
        sig_q.append(ul_from_qhist(hq, ix, TARGET_Q))
    return mc, np.asarray(sig_q)


def load_this_work(path: Path):
    arr = np.loadtxt(path, delimiter=",")
    return arr[:, 0], arr[:, 1]


def load_mass_ev_txt(path: Path):
    rows = []
    for ln in path.read_text().splitlines():
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        sep = "\t" if "\t" in ln else ","
        parts = ln.split(sep)
        if parts[0].lower().startswith("mass"):
            continue
        try:
            m_raw, s = float(parts[0]), float(parts[1])
        except ValueError:
            continue
        m_mev = m_raw / 1e6 if m_raw > 1e4 else m_raw
        rows.append((m_mev, s))
    arr = np.asarray(rows)
    order = np.argsort(arr[:, 0])
    return arr[order, 0], arr[order, 1]


def load_pydme_daily_mod(path: Path):
    import pandas as pd

    df = pd.read_csv(path)
    m = df["mass_MeV"].astype(float).values
    s = df["upper_limit"].astype(float).values
    return m, s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--root",
        default=str(ROOT / "outputs/scan_pattern_data_qedark_fullgrid/scan_dmelectron_pattern.root"),
    )
    ap.add_argument("--outdir", default=str(ROOT / "outplots/qedark_repro"))
    ap.add_argument(
        "--include-this-work",
        action="store_true",
        help="Include dense this_work CSV (not paper figure export)",
    )
    args = ap.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    refs = []

    def add(label, path, loader, style, **plot_kw):
        p = Path(path)
        if not p.exists():
            print(f"skip missing: {label} ({p})")
            return
        m, s = loader(p)
        refs.append({"label": label, "m": m, "s": s, "style": style, "plot_kw": plot_kw})
        print(f"loaded {label}: {len(m)} points, m=[{m.min():.3g}, {m.max():.3g}] MeV")

    add(
        "paper export (canonical)",
        PAPER_EXPORT,
        load_mass_ev_txt,
        "ref",
        color="black",
        ls="-",
        lw=2.2,
        zorder=5,
    )
    if args.include_this_work:
        add(
            "this_work (dense CSV, not paper fig)",
            ROOT / "data/previous_limits/heavy_mediator/DAMIC-M_this_work_QEDark_hm.csv",
            load_this_work,
            "ref",
            color="#888888",
            ls="--",
            lw=1.2,
            alpha=0.7,
        )
    add(
        "repo 2025 txt (hybrid low-m + paper)",
        ROOT / "data/previous_limits/heavy_mediator/DAMIC-M_2025_QEDark_DMe_heavymediator.txt",
        load_mass_ev_txt,
        "ref",
        color="#9467bd",
        ls="-.",
        lw=1.2,
        alpha=0.85,
    )
    add(
        "collab pydme Pattern export",
        ROOT
        / "collab_frameworks/pydme/analysis/DailyModulation/LBC-Sep2024/paper_figures/data/"
        "LBC_results/ScienceRun2024_results-Pattern/DAMIC-M_2025_QEDark_DMe_heavymediator.txt",
        load_mass_ev_txt,
        "ref",
        color="#2ca02c",
        ls=":",
        lw=1.4,
        alpha=0.85,
    )
    pydme_csv = (
        PAPER_EXPORT.parent.parent
        / "LBC_results/upperlimit_LBC_Sep2024_1e-rate_temp_corr_FDM_n0_qedark.csv"
    )
    add(
        "pydme LBC DailyModulation (FDM_n0 qedark)",
        pydme_csv,
        load_pydme_daily_mod,
        "ref",
        color="#ff7f0e",
        ls=(0, (3, 2)),
        lw=1.0,
        alpha=0.75,
    )

    root_path = Path(args.root)
    if not root_path.exists():
        raise SystemExit(f"scan ROOT missing: {root_path}")
    mc, sg_st = load_scan_ul(root_path)
    _, sg_qh = load_scan_qhist(root_path, mc)
    have_qhist = sg_qh is not None and np.any(np.isfinite(sg_qh))
    print(f"loaded CCDarkSens scan (stored UL): {len(mc)} points")
    if have_qhist:
        print(f"loaded q-map crossing UL (--from-qhist style): {np.sum(np.isfinite(sg_qh))} points")

    fig, (ax_ul, ax_ratio) = plt.subplots(
        2, 1, figsize=(9, 8), sharex=True, gridspec_kw={"height_ratios": [2.2, 1], "hspace": 0.08}
    )

    for r in refs:
        ax_ul.plot(r["m"], r["s"], label=r["label"], **r["plot_kw"])

    ax_ul.plot(
        mc,
        sg_st,
        color="#1f77b4",
        lw=2.2,
        label="CCDarkSens scan (stored UL)",
        zorder=10,
    )
    if have_qhist:
        ax_ul.plot(
            mc,
            sg_qh,
            color="#1f77b4",
            lw=1.4,
            ls="--",
            alpha=0.85,
            label="CCDarkSens scan (--from-qhist)",
            zorder=9,
        )

    ax_ul.set_yscale("log")
    ax_ul.set_ylabel(r"$\bar{\sigma}_e$ upper limit [cm$^2$] @ 90% CL")
    ax_ul.set_title("Heavy mediator QEdark — CCDarkSens vs references (paper export primary)")
    ax_ul.grid(True, which="both", alpha=0.25)
    ax_ul.legend(fontsize=7, loc="upper left")

    paper = next((r for r in refs if "canonical" in r["label"]), None)
    thisw = next((r for r in refs if "this_work" in r["label"]), None)

    if paper is not None:
        rp = np.interp(mc, paper["m"], paper["s"])
        ax_ratio.plot(mc, sg_st / rp, color="#1f77b4", lw=1.5, label="stored / paper export")
        if have_qhist:
            ax_ratio.plot(
                mc,
                sg_qh / rp,
                color="#1f77b4",
                lw=1.2,
                ls="--",
                alpha=0.85,
                label="qhist / paper export",
            )
    if thisw is not None:
        rt = np.interp(mc, thisw["m"], thisw["s"])
        ax_ratio.plot(mc, sg_st / rt, color="#888888", lw=1.0, ls=":", label="stored / this_work")

    ax_ratio.axhline(1.0, color="black", lw=0.8, alpha=0.5)
    ax_ratio.set_xscale("log")
    ax_ratio.set_xlabel(r"$m_\chi$ [MeV]")
    ax_ratio.set_ylabel("scan / reference")
    ax_ratio.set_ylim(0, 2.5)
    ax_ratio.grid(True, which="both", alpha=0.25)
    ax_ratio.legend(fontsize=7.5, loc="upper right")

    fig.tight_layout()
    pdf = outdir / "compare_scan_all_references_heavy_qedark.pdf"
    png = outdir / "compare_scan_all_references_heavy_qedark.png"
    fig.savefig(pdf, dpi=150)
    fig.savefig(png, dpi=150)
    print("wrote", pdf)
    print("wrote", png)

    csv_path = outdir / "compare_scan_all_references_summary.csv"
    with csv_path.open("w") as f:
        f.write(
            "reference,median_stored_ratio_1.2_500_MeV,median_qhist_ratio_1.2_500_MeV,"
            "median_stored_ratio_5_500_MeV,m2_MeV_stored_ratio,m2_MeV_qhist_ratio\n"
        )
        m12 = (mc >= 1.2) & (mc <= 500)
        m5 = (mc >= 5) & (mc <= 500)
        i2 = int(np.argmin(np.abs(mc - 2.0)))
        for r in refs:
            ri = np.interp(mc, r["m"], r["s"])
            rat_st = sg_st / ri
            rat_qh = sg_qh / ri if have_qhist else np.full_like(rat_st, np.nan)
            f.write(
                f"{r['label']},"
                f"{np.median(rat_st[m12]):.4f},"
                f"{np.nanmedian(rat_qh[m12]):.4f},"
                f"{np.median(rat_st[m5]):.4f},"
                f"{rat_st[i2]:.4f},"
                f"{rat_qh[i2]:.4f}\n"
            )
    print("wrote", csv_path)


if __name__ == "__main__":
    main()
