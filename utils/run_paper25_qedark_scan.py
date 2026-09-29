#!/usr/bin/env python3
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: run_paper25_qedark_scan.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  run_paper25_qedark_scan.py -- Prepare and run CCDarkSens QEdark scan on
#  the 25 paper-export mass points.
# ============================================================================

"""Prepare and run CCDarkSens QEdark scan on the 25 paper-export mass points."""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
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
SRC_RATES = ROOT / "data/qedark_rates/Si/heavy/long_scan"
PAPER_RATES = ROOT / "data/qedark_rates/Si/heavy/paper25"
CONFIG_PATH = ROOT / "configs/scan_dmelectron_pattern_data_qedark_paper25.json"
OUTDIR = ROOT / "outputs/scan_pattern_data_qedark_paper25"
ROOT_OUT = OUTDIR / "scan_dmelectron_pattern.root"
FULLGRID_ROOT = ROOT / "outputs/scan_pattern_data_qedark_fullgrid/scan_dmelectron_pattern.root"


# ----------------------------------------------------------------------------
# load_paper_masses
#   Masses (converted from eV to MeV) and limits of the paper-export curve.
# ----------------------------------------------------------------------------
def load_paper_masses(path: Path):
    masses, limits = [], []
    for ln in path.read_text().splitlines():
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        parts = ln.replace("\t", ",").split(",")
        if parts[0].lower().startswith("mass"):
            continue
        m_ev, s = float(parts[0]), float(parts[1])
        masses.append(m_ev / 1e6)
        limits.append(s)
    arr = np.asarray(list(zip(masses, limits)))
    order = np.argsort(arr[:, 0])
    return arr[order, 0], arr[order, 1]


# ----------------------------------------------------------------------------
# fmt_mchi
#   Mass formatted with 6 decimals, matching the rate-file names.
# ----------------------------------------------------------------------------
def fmt_mchi(m: float) -> str:
    return f"{m:.6f}"


# ----------------------------------------------------------------------------
# available_masses
#   Sorted masses found in the heavy-mediator silicon rate-file names of a directory.
# ----------------------------------------------------------------------------
def available_masses(rates_dir: Path):
    pat = re.compile(r"_m([0-9.]+)_s")
    masses = set()
    for p in rates_dir.glob("dRdE_Si_heavy_m*_s*.csv"):
        m = pat.search(p.name)
        if m:
            masses.add(float(m.group(1)))
    return sorted(masses)


# ----------------------------------------------------------------------------
# nearest_mass
#   The available mass closest to m.
# ----------------------------------------------------------------------------
def nearest_mass(m: float, avail: list[float]) -> float:
    return min(avail, key=lambda x: abs(x - m))


# ----------------------------------------------------------------------------
# prepare_rates
#   Build the rate directory for the paper masses: for each mass I copy the files of the nearest available mass, renamed to the requested mass. An existing directory is kept unless force is set.
# ----------------------------------------------------------------------------
def prepare_rates(masses: np.ndarray, src_dir: Path, dst_dir: Path, force: bool = False):
    avail = available_masses(src_dir)
    if dst_dir.exists() and not force:
        print(f"rates dir exists: {dst_dir} (use --force-rates to rebuild)")
        return
    if dst_dir.exists():
        shutil.rmtree(dst_dir)
    dst_dir.mkdir(parents=True)
    pat = re.compile(r"(_m)([0-9.]+)(_s)")

    mapping = []
    n_files = 0
    for m in masses:
        m_src = nearest_mass(float(m), avail)
        rel = abs(m_src - m) / m if m > 0 else 0.0
        mapping.append((float(m), m_src, rel))
        src_glob = list(src_dir.glob(f"dRdE_Si_heavy_m{fmt_mchi(m_src)}_s*.csv"))
        if not src_glob:
            raise SystemExit(f"no rate files for source mass {m_src} (paper {m})")
        for src in src_glob:
            dst_name = pat.sub(rf"\g<1>{fmt_mchi(m)}\g<3>", src.name)
            dst = dst_dir / dst_name
            try:
                os.link(src, dst)
            except OSError:
                shutil.copy2(src, dst)
            n_files += 1
    print(f"prepared {n_files} rate files for {len(masses)} masses → {dst_dir}")
    for m, m_src, rel in mapping:
        flag = "ok" if rel < 0.001 else "near"
        print(f"  paper {m:12.6f} MeV ← long_scan {m_src:.6f}  Δ={rel*100:.3f}%  [{flag}]")
    return mapping


# ----------------------------------------------------------------------------
# write_config
#   Write the scan config for the paper mass points from the full-grid config, with the label, output directory, rates directory and mass list replaced.
# ----------------------------------------------------------------------------
def write_config(masses: list[float], path: Path):
    base = json.loads(
        (ROOT / "configs/scan_dmelectron_pattern_data_qedark_fullgrid.json").read_text()
    )
    base["_comment"] = "LBC pattern analysis on the 25 paper-export mass points (heavy QEdark)."
    base["run"]["label"] = "scan_pattern_data_qedark_paper25"
    base["run"]["outdir"] = "outputs/scan_pattern_data_qedark_paper25"
    base["model"]["rates_dir"] = "data/qedark_rates/Si/heavy/paper25"
    base["model"]["grid"]["mchi_MeV"] = {"values": masses}
    path.write_text(json.dumps(base, indent=2) + "\n")
    print(f"wrote {path}")


# ----------------------------------------------------------------------------
# run_scan
#   Run the scan binary from the build directory on a config; returns 1 if the binary does not exist.
# ----------------------------------------------------------------------------
def run_scan(build_dir: Path, config: Path) -> int:
    exe = build_dir / "ccdarksens_scan_dmelectron_pattern"
    if not exe.exists():
        print(f"ERROR: {exe} not found")
        return 1
    cmd = [str(exe), str(config)]
    print("run:", " ".join(cmd))
    return subprocess.call(cmd, cwd=ROOT)


def load_scan_by_index(root_path: Path):
    """Load UL values in grid order (bin index = mass index for explicit values grid)."""
    import uproot

    f = uproot.open(root_path)
    ul = f["upper_limit_sigma_e_mchi"]
    v = ul.values()
    ok = (v > 0) & (v < 0.9e-26)
    return v[ok]


# ----------------------------------------------------------------------------
# compare_and_plot
#   Compare the scan limits at the paper masses (and the full grid, if given) with the paper curve, print the ratios and save the plots to outdir.
# ----------------------------------------------------------------------------
def compare_and_plot(
    paper_m: np.ndarray,
    paper_s: np.ndarray,
    scan_s: np.ndarray,
    full_m: np.ndarray | None,
    full_s: np.ndarray | None,
    outdir: Path,
):
    outdir.mkdir(parents=True, exist_ok=True)
    n = min(len(paper_m), len(scan_s))
    paper_m = paper_m[:n]
    paper_s = paper_s[:n]
    scan_at_paper = scan_s[:n]

    full_at_paper = None
    if full_m is not None and full_s is not None:
        full_at_paper = np.interp(paper_m, full_m, full_s)

    ratio_p25 = scan_at_paper / paper_s
    ratio_fg = full_at_paper / paper_s if full_at_paper is not None else None

    print("\n=== Point-by-point vs paper export ===")
    print(f"{'m [MeV]':>12} {'paper':>12} {'paper25':>12} {'fullgrid':>12} {'p25/pap':>10} {'fg/pap':>10}")
    for i, m in enumerate(paper_m):
        fg = full_at_paper[i] if full_at_paper is not None else float("nan")
        print(
            f"{m:12.4f} {paper_s[i]:12.3e} {scan_at_paper[i]:12.3e} {fg:12.3e} "
            f"{ratio_p25[i]:10.3f} {fg/paper_s[i] if full_at_paper is not None else float('nan'):10.3f}"
        )

    m_low = paper_m <= 5.0
    m_mid = (paper_m >= 5.0) & (paper_m <= 500.0)
    print(f"\nMedian p25/paper (m<=5 MeV):   {np.nanmedian(ratio_p25[m_low]):.3f}")
    print(f"Median p25/paper (5-500 MeV):  {np.nanmedian(ratio_p25[m_mid]):.3f}")
    if ratio_fg is not None:
        print(f"Median fullgrid/paper (m<=5): {np.nanmedian(ratio_fg[m_low]):.3f}")
        print(f"Median fullgrid/paper (5-500): {np.nanmedian(ratio_fg[m_mid]):.3f}")
        print(f"Median |p25-fg|/fg at paper m: {np.nanmedian(np.abs(ratio_p25 - ratio_fg)):.4f}")

    csv = outdir / "paper25_vs_paper_export.csv"
    with csv.open("w") as f:
        f.write("m_mev,paper_ul,paper25_ul,fullgrid_ul,paper25_over_paper,fullgrid_over_paper\n")
        for i, m in enumerate(paper_m):
            fg = full_at_paper[i] if full_at_paper is not None else ""
            fg_r = ratio_fg[i] if ratio_fg is not None else ""
            f.write(
                f"{m},{paper_s[i]},{scan_at_paper[i]},{fg},{ratio_p25[i]},{fg_r}\n"
            )
    print(f"wrote {csv}")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True, gridspec_kw={"height_ratios": [2, 1]})
    ax1.loglog(paper_m, paper_s, "k-o", ms=4, lw=1.5, label="paper export")
    ax1.loglog(paper_m, scan_at_paper, "b-s", ms=4, lw=1.5, label="CCDarkSens paper25 scan")
    if full_at_paper is not None:
        ax1.loglog(paper_m, full_at_paper, "c^", ms=3, lw=1.0, alpha=0.8, label="fullgrid @ paper m")
    ax1.set_ylabel(r"$\bar{\sigma}_e$ UL [cm$^2$]")
    ax1.set_title("Heavy QEdark — 25 paper mass points")
    ax1.grid(True, which="both", alpha=0.3)
    ax1.legend(fontsize=8)

    ax2.semilogx(paper_m, ratio_p25, "b-o", ms=4, label="paper25 / paper")
    if ratio_fg is not None:
        ax2.semilogx(paper_m, ratio_fg, "c^", ms=3, alpha=0.8, label="fullgrid / paper")
    ax2.axhline(1.0, color="k", lw=0.8, alpha=0.5)
    ax2.set_xlabel(r"$m_\chi$ [MeV]")
    ax2.set_ylabel("scan / paper")
    ax2.set_ylim(0, 2.5)
    ax2.grid(True, which="both", alpha=0.3)
    ax2.legend(fontsize=8)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        p = outdir / f"paper25_vs_paper_export.{ext}"
        fig.savefig(p, dpi=150)
        print(f"wrote {p}")


# ----------------------------------------------------------------------------
# main
#   Command line: prepare the rate files and the config for the 25 paper mass points, optionally run the scan (--run-scan), and compare the result with the paper (--compare-only skips the rest).
# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build-dir", default=str(ROOT / "build"))
    ap.add_argument("--run-scan", action="store_true")
    ap.add_argument("--force-rates", action="store_true")
    ap.add_argument("--skip-rates", action="store_true")
    ap.add_argument("--compare-only", action="store_true")
    args = ap.parse_args()

    paper_m, paper_s = load_paper_masses(PAPER_EXPORT)
    masses = [float(m) for m in paper_m]
    print(f"paper export: {len(masses)} masses, m=[{masses[0]:.4g}, {masses[-1]:.4g}] MeV")

    if not args.skip_rates and not args.compare_only:
        prepare_rates(paper_m, SRC_RATES, PAPER_RATES, force=args.force_rates)
        write_config(masses, CONFIG_PATH)

    if args.run_scan and not args.compare_only:
        rc = run_scan(Path(args.build_dir), CONFIG_PATH)
        if rc != 0:
            sys.exit(rc)

    if not ROOT_OUT.exists():
        print(f"scan output missing: {ROOT_OUT}")
        if not args.compare_only:
            print("Re-run with --run-scan")
        sys.exit(1)

    sm = load_scan_by_index(ROOT_OUT)
    fm = fs = None
    if FULLGRID_ROOT.exists():
        import uproot

        fg = uproot.open(FULLGRID_ROOT)
        ful = fg["upper_limit_sigma_e_mchi"]
        fxe = ful.axis().edges()
        fxc = 0.5 * (fxe[:-1] + fxe[1:])
        fv = ful.values()
        ok = (fv > 0) & (fv < 0.9e-26)
        fm, fs = fxc[ok], fv[ok]
    compare_and_plot(paper_m, paper_s, sm, fm, fs, ROOT / "outplots/qedark_repro")


if __name__ == "__main__":
    main()
