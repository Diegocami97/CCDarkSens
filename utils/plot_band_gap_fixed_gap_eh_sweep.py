#!/usr/bin/env python3
# ============================================================================
#  CCDarkSens — plot_band_gap_fixed_gap_eh_sweep
#  Plot σ_UL vs ε_h at fixed band gap for heavy and light mediators
#
#  Author: Diego Venegas-Vargas
# ============================================================================
"""
Overlay DM-e limits at fixed E_gap, scanning epsilon_h (2D-grid row).

Complement to plot_band_gap_fixed_eh_gap_sweep.py (fixed epsilon_h, scan E_gap).

  python3 utils/plot_band_gap_fixed_gap_eh_sweep.py --gap 0.1 --mediator light
  python3 utils/plot_band_gap_fixed_gap_eh_sweep.py --gap 0.3 0.5 0.7 0.9 --all-mediators
  python3 utils/plot_band_gap_fixed_gap_eh_sweep.py --all-gaps --all-mediators
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "utils"))

from band_gap_limit_plot import is_si_reference_cell, reference_scan_cli_args  # noqa: E402
from band_gap_plot_labels import pheno_param_label_root  # noqa: E402
from band_gap_scan_paths import (  # noqa: E402
    EH_GRID,
    GAP_GRID,
    SI_REF_EH_EV,
    SI_REF_GAP_EV,
    ev_tag,
    scan_root_path,
)

BUILD_PLOT = ROOT / "build" / "ccdarksens_plot_dmelectron_limit"
OUT_PLOTS = ROOT / "outplots" / "band_gap_pheno" / "step5_limits"
MEDIATORS = ("heavy", "light")


def eh_values(gap: float, include_dequal: bool) -> list[float]:
    ehs = list(EH_GRID)
    if include_dequal and gap not in ehs and gap >= 0:
        ehs.insert(0, gap)
    return sorted(ehs)


def sweep_title(mediator: str, gap: float) -> str:
    med = "light" if mediator == "light" else "heavy"
    return f"Band-gap pheno: E_{{gap}} = {gap:g} eV ({med} mediator)"


def plot_one(
    mediator: str, gap: float, include_dequal: bool, si_reference: bool, dry_run: bool
) -> int:
    paths: list[str] = []
    labels: list[str] = []
    for eh in eh_values(gap, include_dequal):
        if eh < gap - 1e-9:
            continue
        root = scan_root_path(mediator, gap, eh)
        if not root.is_file():
            print(f"WARN: skip missing {root}")
            continue
        if si_reference and is_si_reference_cell(gap, eh):
            continue  # drawn via --reference-scan (solid overlay)
        paths.append(root.relative_to(ROOT).as_posix())
        labels.append(pheno_param_label_root(gap, eh))

    if len(paths) < 2:
        print(f"ERROR: {mediator} gap={gap:g}: need >= 2 ROOT files", file=sys.stderr)
        return 1

    med_slug = "" if mediator == "heavy" else f"{mediator}_"
    out_stem = f"limit_sweep_{med_slug}gap{ev_tag(gap)}_eh_scan"
    out_pdf = OUT_PLOTS / f"{out_stem}.pdf"
    out_csv = OUT_PLOTS / f"{out_stem}.csv"
    OUT_PLOTS.mkdir(parents=True, exist_ok=True)

    print(f"[eh-sweep] {mediator} E_gap={gap:g} eV -> {out_pdf.name} ({len(paths)} curves)")
    cmd = [str(BUILD_PLOT)]
    for p, lab in zip(paths, labels):
        cmd.extend([p, lab])
    cmd.extend(reference_scan_cli_args(mediator, enabled=si_reference))
    cmd.extend(
        [
            mediator,
            "--from-qhist",
            "--batch",
            "--title",
            sweep_title(mediator, gap),
            "--out-pdf",
            str(out_pdf),
            "--out-csv",
            str(out_csv),
        ]
    )
    print("[run]", " ".join(cmd))
    if dry_run:
        return 0
    return subprocess.call(cmd, cwd=ROOT)


def gaps_for_args(args: argparse.Namespace) -> list[float]:
    if args.all_gaps:
        return list(GAP_GRID)
    return args.gap if args.gap else [0.1]


def main() -> int:
    gap_help = (
        f"Fixed E_gap [eV]; study grid: {', '.join(f'{g:g}' for g in GAP_GRID)}. "
        "Per gap, epsilon_h runs over {0.5,1,1.5,2,2.5,3.8} with eh>=E_gap "
        "(5 curves at 0.7/0.9, 4 at 1.2; use --include-dequal for eh=E_gap)."
    )
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--gap",
        type=float,
        nargs="*",
        metavar="eV",
        help=gap_help,
    )
    ap.add_argument(
        "--all-gaps",
        action="store_true",
        help=f"Same as --gap {' '.join(str(g) for g in GAP_GRID)}",
    )
    ap.add_argument("--mediator", choices=MEDIATORS, default="light")
    ap.add_argument("--all-mediators", action="store_true")
    ap.add_argument(
        "--include-dequal",
        action="store_true",
        help=f"Also include epsilon_h = E_gap (e.g. {0.1} eV at gap 0.1)",
    )
    ap.add_argument(
        "--no-si-reference",
        action="store_true",
        help=f"Do not overlay Si reference (E_gap={SI_REF_GAP_EV}, epsilon_h={SI_REF_EH_EV})",
    )
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    if not BUILD_PLOT.is_file() and not args.dry_run:
        rc = subprocess.call(
            ["cmake", "--build", "build", "-j", "--target", "ccdarksens_plot_dmelectron_limit"],
            cwd=ROOT,
        )
        if rc != 0:
            return rc

    mediators = list(MEDIATORS) if args.all_mediators else [args.mediator]
    gaps = gaps_for_args(args)
    rc = 0
    for gap in gaps:
        for med in mediators:
            if plot_one(
                med, gap, args.include_dequal, not args.no_si_reference, args.dry_run
            ) != 0:
                rc = 1
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
