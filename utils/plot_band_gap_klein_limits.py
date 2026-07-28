#!/usr/bin/env python3
# ============================================================================
#  CCDarkSens — plot_band_gap_klein_limits
#  Overlay DM-e limit curves for the Klein-tier band-gap pheno scan.
#  One curve per E_gap (Klein), plus Si reference (1.2, 3.8) eV.
#
#  Author: Diego Venegas-Vargas
# ============================================================================
"""
Generate limit-overlay PDFs (heavy + light) for the Klein ladder + Si ref:

  python3 utils/plot_band_gap_klein_limits.py
  python3 utils/plot_band_gap_klein_limits.py --dc-tier dc10x
  python3 utils/plot_band_gap_klein_limits.py --dc-tier dc100x --mediator light

Outputs (outplots/band_gap_pheno/step5_limits_refix/ or step6_dc_sensitivity/):
  limit_sweep_klein.pdf              baseline heavy
  limit_sweep_light_klein.pdf        baseline light
  limit_sweep_klein_dc10x.pdf        10x DC
  limit_sweep_light_klein_dc100x.pdf 100x DC
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "utils"))
from band_gap_plot_labels import si_reference_label_root  # noqa: E402

BUILD_PLOT = ROOT / "build" / "ccdarksens_plot_dmelectron_limit"
OUT_STEP5 = ROOT / "outplots" / "band_gap_pheno" / "step5_limits_refix"
OUT_STEP6 = ROOT / "outplots" / "band_gap_pheno" / "step6_dc_sensitivity"

GAPS = [0.1, 0.3, 0.5, 0.7, 0.9, 1.2]
KLEIN_EH = {g: round(2.8 * g + 0.5, 2) for g in GAPS}

DC_LAMBDA = {
    "baseline": 0.0365,
    "dc10x": 0.365,
    "dc100x": 3.65,
}
EXPOSURE_KG_YR = 0.5
DAYS_PER_YEAR = 365.25


def ev_tag(x: float) -> str:
    return ("%.2f" % x).replace(".", "p")


def gap_tag(g: float) -> str:
    return ("%.1f" % g).replace(".", "p")


def _dc_suffix(dc_tier: str) -> str:
    return "" if dc_tier == "baseline" else f"_{dc_tier}"


def klein_scan_root(mediator: str, g: float, dc_tier: str) -> Path:
    eh = KLEIN_EH[g]
    gtag = gap_tag(g)
    etag = ev_tag(eh)
    prefix = "scan_band_gap_light_klein" if mediator == "light" else "scan_band_gap_klein"
    suffix = _dc_suffix(dc_tier)
    return (
        ROOT
        / "outputs"
        / f"{prefix}_{gtag}_eh{etag}{suffix}"
        / "scan_dmelectron_pattern.root"
    )


def si_ref_scan_root(mediator: str, dc_tier: str) -> Path:
    suffix = _dc_suffix(dc_tier)
    if mediator == "light":
        sub = f"scan_band_gap_light_1p2_eh3p8{suffix}"
    else:
        sub = f"scan_band_gap_1p2_eh3p8{suffix}"
    return ROOT / "outputs" / "refix_eff_roi" / sub / "scan_dmelectron_pattern.root"


def klein_curve_label(g: float) -> str:
    """ROOT TLatex legend text — no $ (ROOT ExpandFileName treats $ as env vars)."""
    eh = KLEIN_EH[g]
    eh_s = f"{eh:.2f}".rstrip("0").rstrip(".")
    return f"E_{{gap}} = {g:g} eV, #varepsilon_{{h}} = {eh_s} eV"


def format_dc_per_pix_per_day(lam_e_per_pix_per_year: float) -> str:
    lam_day = lam_e_per_pix_per_year / DAYS_PER_YEAR
    if lam_day < 0.01:
        return f"{lam_day:.2e}"
    return f"{lam_day:g}"


def plot_title(dc_tier: str) -> str:
    lam_day = format_dc_per_pix_per_day(DC_LAMBDA[dc_tier])
    return (
        f"Low band gap study "
        f"({EXPOSURE_KG_YR:g} kg#cdotyr, DC = {lam_day} e^{{-}}/pix/day)"
    )


def plot_one(mediator: str, dc_tier: str, dry_run: bool = False) -> int:
    paths: list[str] = []
    labels: list[str] = []

    si_root = si_ref_scan_root(mediator, dc_tier)
    if si_root.is_file():
        paths.append(si_root.relative_to(ROOT).as_posix())
        labels.append(si_reference_label_root())
    else:
        print(f"WARN: missing Si ref {si_root}")

    for g in GAPS:
        root = klein_scan_root(mediator, g, dc_tier)
        if not root.is_file():
            print(f"WARN: missing {root}")
            continue
        paths.append(root.relative_to(ROOT).as_posix())
        labels.append(klein_curve_label(g))

    if len(paths) < 2:
        print(
            f"ERROR: need >= 2 ROOT files for {mediator} Klein ({dc_tier})",
            file=sys.stderr,
        )
        return 1

    out_dir = OUT_STEP5 if dc_tier == "baseline" else OUT_STEP6
    out_dir.mkdir(parents=True, exist_ok=True)
    med_slug = "" if mediator == "heavy" else "light_"
    tier_slug = "" if dc_tier == "baseline" else f"_{dc_tier}"
    out_pdf = out_dir / f"limit_sweep_{med_slug}klein{tier_slug}.pdf"
    out_csv = out_dir / f"limit_sweep_{med_slug}klein{tier_slug}.csv"
    title = plot_title(dc_tier)

    cmd = [str(BUILD_PLOT)]
    for p, lab in zip(paths, labels):
        cmd.extend([p, lab])
    cmd.extend(
        [
            mediator,
            "--from-qhist",
            "--plain-legend",
            "--batch",
            "--title",
            title,
            "--out-pdf",
            str(out_pdf),
            "--out-csv",
            str(out_csv),
        ]
    )

    print(f"[klein-limits] {mediator} {dc_tier} -> {out_pdf.name}  ({len(paths)} curves)")
    if dry_run:
        print("[run]", " ".join(cmd))
        return 0
    return subprocess.call(cmd, cwd=str(ROOT))


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--mediator", choices=["heavy", "light", "both"], default="both")
    ap.add_argument(
        "--dc-tier",
        choices=["baseline", "dc10x", "dc100x", "all"],
        default="baseline",
        help="Dark-current tier (default: baseline). Use 'all' for every tier.",
    )
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    if not BUILD_PLOT.is_file() and not args.dry_run:
        rc = subprocess.call(
            ["cmake", "--build", "build", "-j", "--target", "ccdarksens_plot_dmelectron_limit"],
            cwd=str(ROOT),
        )
        if rc != 0:
            return rc

    tiers = (
        ["baseline", "dc10x", "dc100x"]
        if args.dc_tier == "all"
        else [args.dc_tier]
    )
    meds = ["heavy", "light"] if args.mediator == "both" else [args.mediator]
    rc = 0
    for tier in tiers:
        for med in meds:
            if plot_one(med, tier, dry_run=args.dry_run) != 0:
                rc = 1
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
