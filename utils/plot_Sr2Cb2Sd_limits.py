#!/usr/bin/env python3
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: plot_Sr2Cb2Sd_limits.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  plot_Sr2Cb2Sd_limits.py -- Overlay DM-e limit curves for the Sr2Cb2Sd
#  phase-1 study (Si_fast epsilon).
# ============================================================================
"""
Six Sr2Cb2Sd limit figures (heavy + light):

  Baseline (2):
    limit_baseline_heavy.pdf   — DAMIC-M + SrCd indirect/direct + OSCURA
    limit_baseline_light.pdf

  DC study, indirect gap only (2):
    limit_dc_indirect_heavy.pdf  — DAMIC-M + 1x/100x/1000x DC + OSCURA
    limit_dc_indirect_light.pdf

  DC study, direct gap only (2):
    limit_dc_direct_heavy.pdf
    limit_dc_direct_light.pdf

  python3 utils/plot_Sr2Cb2Sd_limits.py
  python3 utils/plot_Sr2Cb2Sd_limits.py --which baseline
  python3 utils/plot_Sr2Cb2Sd_limits.py --which dc
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "utils"))
from band_gap_plot_labels import pheno_param_label_root  # noqa: E402

BUILD_PLOT = ROOT / "build" / "ccdarksens_plot_dmelectron_limit"
OUT_DIR = ROOT / "outplots" / "Sr2Cb2Sd"

SI_GAP_EV = 1.2
SI_EH_EV = 3.8
SRCD_INDIRECT_GAP_EV = 0.556
SRCD_INDIRECT_EH_EV = 2.06
SRCD_DIRECT_GAP_EV = 0.603
SRCD_DIRECT_EH_EV = 2.19


# ----------------------------------------------------------------------------
# _phys
#   ROOT label of the physical parameters (gap, eps_h) of a curve.
# ----------------------------------------------------------------------------
def _phys(gap_eV: float, eh_eV: float) -> str:
    return pheno_param_label_root(gap_eV, eh_eV)


# ----------------------------------------------------------------------------
# _label
#   Two-line ROOT legend entry: experiment name and exposure on the first line, the (gap, eps_h) label on the second.
# ----------------------------------------------------------------------------
def _label(name: str, exposure: str, gap_eV: float, eh_eV: float) -> str:
    # Two-line ROOT legend: avoids wide boxes and nested-brace TLatex glitches.
    return f"#splitline{{{name} ({exposure})}}{{{_phys(gap_eV, eh_eV)}}}"


# ----------------------------------------------------------------------------
# scan_root
#   Path of the scan ROOT file of a run tag under outputs/Sr2Cb2Sd/.
# ----------------------------------------------------------------------------
def scan_root(run_tag: str) -> Path:
    return ROOT / "outputs" / "Sr2Cb2Sd" / run_tag / "scan_dmelectron_pattern.root"


# ----------------------------------------------------------------------------
# run_plot
#   Plot several limit curves of one mediator with ccdarksens_plot_dmelectron_limit (missing scans are skipped; with dry_run only the command is printed). Returns the plotter's exit code.
# ----------------------------------------------------------------------------
def run_plot(
    *,
    mediator: str,
    curves: list[tuple[str, str]],
    out_stem: str,
    title: str,
    dry_run: bool,
) -> int:
    paths: list[str] = []
    labels: list[str] = []
    for tag, label in curves:
        p = scan_root(tag)
        if not p.is_file():
            print(f"WARN: missing {p}")
            continue
        paths.append(p.relative_to(ROOT).as_posix())
        labels.append(label)

    if len(paths) < 2:
        print(f"ERROR: need >= 2 curves for {out_stem} ({mediator})", file=sys.stderr)
        return 1

    out_pdf = OUT_DIR / f"{out_stem}.pdf"
    out_csv = OUT_DIR / f"{out_stem}.csv"
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    cmd = [
        str(BUILD_PLOT),
        "--batch",
        "--from-qhist",
        "--plain-legend",
        "--legend-right",
        "--title",
        title,
        "--out-pdf",
        str(out_pdf),
        "--out-csv",
        str(out_csv),
    ]
    for p, lab in zip(paths, labels):
        cmd.extend([p, lab])
    cmd.append(mediator)

    print(f"[Sr2Cb2Sd] {out_stem} ({mediator}) -> {out_pdf.name}  ({len(paths)} curves)")
    if dry_run:
        print("[run]", " ".join(cmd))
        return 0
    return subprocess.call(cmd, cwd=str(ROOT))


# ----------------------------------------------------------------------------
# damic_curve
#   Run tag and legend label of the DAMIC-M silicon reference curve (1 kg-yr).
# ----------------------------------------------------------------------------
def damic_curve(med: str) -> tuple[str, str]:
    return (
        f"Sr2Cb2Sd_si_ref_gap1p2_{med}_1kgy",
        _label("DAMIC-M", "1 kg-yr", SI_GAP_EV, SI_EH_EV),
    )


# ----------------------------------------------------------------------------
# oscura_curve
#   Run tag and legend label of the OSCURA silicon curve (30 kg-yr).
# ----------------------------------------------------------------------------
def oscura_curve(med: str) -> tuple[str, str]:
    return (
        f"Sr2Cb2Sd_oscura_gap1p2_{med}_30kgy",
        _label("OSCURA", "30 kg-yr", SI_GAP_EV, SI_EH_EV),
    )


# ----------------------------------------------------------------------------
# baseline_curves
#   Curves of the baseline figure: DAMIC-M silicon, the SrCd2Sb2 indirect-gap and direct-gap cases, and OSCURA.
# ----------------------------------------------------------------------------
def baseline_curves(med: str) -> list[tuple[str, str]]:
    return [
        damic_curve(med),
        (
            f"Sr2Cb2Sd_srcd_indirect_gap0p556_{med}_1kgy",
            _label(
                "SrCd_{2}Sb_{2} indirect",
                "1 kg-yr",
                SRCD_INDIRECT_GAP_EV,
                SRCD_INDIRECT_EH_EV,
            ),
        ),
        (
            f"Sr2Cb2Sd_srcd_direct_gap0p603_{med}_1kgy",
            _label(
                "SrCd_{2}Sb_{2} direct",
                "1 kg-yr",
                SRCD_DIRECT_GAP_EV,
                SRCD_DIRECT_EH_EV,
            ),
        ),
        oscura_curve(med),
    ]


# ----------------------------------------------------------------------------
# _dc_tier_label
#   Legend entry of a dark-current tier (dc times the baseline) with the (gap, eps_h) label.
# ----------------------------------------------------------------------------
def _dc_tier_label(dc: str, gap_eV: float, eh_eV: float) -> str:
    return f"#splitline{{DC {dc}#times}}{{{_phys(gap_eV, eh_eV)}}}"


# ----------------------------------------------------------------------------
# dc_indirect_curves
#   Curves of the dark-current study for the indirect gap: DAMIC-M reference plus the 1x and 100x dark-current tiers.
# ----------------------------------------------------------------------------
def dc_indirect_curves(med: str) -> list[tuple[str, str]]:
    g, eh = SRCD_INDIRECT_GAP_EV, SRCD_INDIRECT_EH_EV
    return [
        damic_curve(med),
        (
            f"Sr2Cb2Sd_srcd_indirect_gap0p556_{med}_1kgy",
            _dc_tier_label("1", g, eh),
        ),
        (
            f"Sr2Cb2Sd_srcd_indirect_gap0p556_{med}_1kgy_dc100x",
            _dc_tier_label("100", g, eh),
        ),
        (
            f"Sr2Cb2Sd_srcd_indirect_gap0p556_{med}_1kgy_dc1000x",
            _dc_tier_label("1000", g, eh),
        ),
        oscura_curve(med),
    ]


# ----------------------------------------------------------------------------
# dc_direct_curves
#   Curves of the dark-current study for the direct gap: DAMIC-M reference plus the 1x and 100x dark-current tiers.
# ----------------------------------------------------------------------------
def dc_direct_curves(med: str) -> list[tuple[str, str]]:
    g, eh = SRCD_DIRECT_GAP_EV, SRCD_DIRECT_EH_EV
    return [
        damic_curve(med),
        (
            f"Sr2Cb2Sd_srcd_direct_gap0p603_{med}_1kgy",
            _dc_tier_label("1", g, eh),
        ),
        (
            f"Sr2Cb2Sd_srcd_direct_gap0p603_{med}_1kgy_dc100x",
            _dc_tier_label("100", g, eh),
        ),
        (
            f"Sr2Cb2Sd_srcd_direct_gap0p603_{med}_1kgy_dc1000x",
            _dc_tier_label("1000", g, eh),
        ),
        oscura_curve(med),
    ]


# ----------------------------------------------------------------------------
# plot_baseline
#   Baseline figure for one mediator (1 kg-yr, DC = 1e-5 e-/pix/day).
# ----------------------------------------------------------------------------
def plot_baseline(mediator: str, dry_run: bool) -> int:
    title = (
        "Sr2Cb2Sd projection (1 kg-yr, DC = 10^{-5} e^{-}/pix/day)"
        if mediator == "heavy"
        else "Sr2Cb2Sd projection, light mediator (1 kg-yr, DC = 10^{-5} e^{-}/pix/day)"
    )
    return run_plot(
        mediator=mediator,
        curves=baseline_curves(mediator),
        out_stem=f"limit_baseline_{mediator}",
        title=title,
        dry_run=dry_run,
    )


# ----------------------------------------------------------------------------
# _dc_study_title
#   Title of a dark-current study figure.
# ----------------------------------------------------------------------------
def _dc_study_title(gap_branch: str, mediator: str) -> str:
    med_label = "heavy" if mediator == "heavy" else "light"
    return f"Sr2Cb2Sd DC study, {gap_branch} gap, {med_label} mediator (1 kg-yr)"


# ----------------------------------------------------------------------------
# plot_dc_indirect
#   Dark-current study figure for the indirect gap.
# ----------------------------------------------------------------------------
def plot_dc_indirect(mediator: str, dry_run: bool) -> int:
    return run_plot(
        mediator=mediator,
        curves=dc_indirect_curves(mediator),
        out_stem=f"limit_dc_indirect_{mediator}",
        title=_dc_study_title("indirect", mediator),
        dry_run=dry_run,
    )


# ----------------------------------------------------------------------------
# plot_dc_direct
#   Dark-current study figure for the direct gap.
# ----------------------------------------------------------------------------
def plot_dc_direct(mediator: str, dry_run: bool) -> int:
    return run_plot(
        mediator=mediator,
        curves=dc_direct_curves(mediator),
        out_stem=f"limit_dc_direct_{mediator}",
        title=_dc_study_title("direct", mediator),
        dry_run=dry_run,
    )


# ----------------------------------------------------------------------------
# main
#   Command line: choose the mediator (heavy, light or both) and which figures to make (all six, baseline, or the DC studies); --dry-run only prints the plot commands.
# ----------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--mediator", choices=["heavy", "light", "both"], default="both")
    ap.add_argument(
        "--which",
        choices=["all", "baseline", "dc", "dc-indirect", "dc-direct"],
        default="all",
        help="all = six figures; baseline = two; dc = four DC figures",
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

    meds = ["heavy", "light"] if args.mediator == "both" else [args.mediator]
    rc = 0
    for med in meds:
        if args.which in ("all", "baseline"):
            if plot_baseline(med, args.dry_run) != 0:
                rc = 1
        if args.which in ("all", "dc", "dc-indirect"):
            if plot_dc_indirect(med, args.dry_run) != 0:
                rc = 1
        if args.which in ("all", "dc", "dc-direct"):
            if plot_dc_direct(med, args.dry_run) != 0:
                rc = 1
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
