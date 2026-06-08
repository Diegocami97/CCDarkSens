#!/usr/bin/env python3
# ============================================================================
#  CCDarkSens — plot_band_gap_fixed_eh_gap_sweep
#  Plot σ_UL vs band gap at fixed ε_h for heavy and light mediators
#
#  Author: Diego Venegas-Vargas
# ============================================================================
"""
Overlay DM-e limits at fixed epsilon_h, scanning E_gap (2D-grid column).

Complement to plot_band_gap_fixed_gap_eh_sweep.py (fixed E_gap, scan epsilon_h).
Same as run_band_gap_phase_c.py plot-limits --tier B-thresh when --eh 3.8.

  python3 utils/plot_band_gap_fixed_eh_gap_sweep.py --eh 3.8 --mediator light
  python3 utils/plot_band_gap_fixed_eh_gap_sweep.py --eh 3.8 --all-mediators
  python3 utils/plot_band_gap_fixed_eh_gap_sweep.py --eh 1.0 --gap 0.3 0.5 0.7
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
DEFAULT_EH = 3.8


def gaps_for_eh(eh: float, gap_list: list[float] | None) -> list[float]:
    """Gaps on the 2D grid with epsilon_h >= E_gap."""
    base = gap_list if gap_list is not None else list(GAP_GRID)
    return [g for g in base if eh >= g - 1e-9]


def sweep_title(mediator: str, eh: float) -> str:
    med = "light" if mediator == "light" else "heavy"
    return f"Band-gap pheno: #varepsilon_{{h}} = {eh:g} eV ({med} mediator)"


def plot_one(
    mediator: str,
    eh: float,
    gaps: list[float],
    si_reference: bool,
    dry_run: bool,
) -> int:
    paths: list[str] = []
    labels: list[str] = []
    for gap in gaps:
        root = scan_root_path(mediator, gap, eh)
        if not root.is_file():
            print(f"WARN: skip missing {root}")
            continue
        if si_reference and is_si_reference_cell(gap, eh):
            continue
        paths.append(root.relative_to(ROOT).as_posix())
        labels.append(pheno_param_label_root(gap, eh))

    if len(paths) < 2:
        print(f"ERROR: {mediator} eh={eh:g}: need >= 2 ROOT files", file=sys.stderr)
        return 1

    med_slug = "" if mediator == "heavy" else f"{mediator}_"
    out_stem = f"limit_sweep_{med_slug}eh{ev_tag(eh)}_gap_scan"
    out_pdf = OUT_PLOTS / f"{out_stem}.pdf"
    out_csv = OUT_PLOTS / f"{out_stem}.csv"
    OUT_PLOTS.mkdir(parents=True, exist_ok=True)

    print(f"[gap-sweep] {mediator} epsilon_h={eh:g} eV -> {out_pdf.name} ({len(paths)} curves)")
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
            sweep_title(mediator, eh),
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


def main() -> int:
    eh_choices = ", ".join(f"{e:g}" for e in EH_GRID)
    gap_choices = ", ".join(f"{g:g}" for g in GAP_GRID)
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--eh",
        type=float,
        default=DEFAULT_EH,
        help=f"Fixed epsilon_h [eV] (default {DEFAULT_EH}; grid: {eh_choices})",
    )
    ap.add_argument(
        "--gap",
        type=float,
        nargs="*",
        metavar="eV",
        help=f"E_gap values to include (default: all grid gaps with eh>=E_gap: {gap_choices})",
    )
    ap.add_argument(
        "--all-gaps",
        action="store_true",
        help=f"Use full gap grid: {gap_choices}",
    )
    ap.add_argument("--mediator", choices=MEDIATORS, default="light")
    ap.add_argument("--all-mediators", action="store_true")
    ap.add_argument(
        "--no-si-reference",
        action="store_true",
        help=f"Do not overlay Si reference (E_gap={SI_REF_GAP_EV}, epsilon_h={SI_REF_EH_EV}, gold)",
    )
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    gap_list = list(GAP_GRID) if args.all_gaps or not args.gap else list(args.gap)
    gaps = gaps_for_eh(args.eh, gap_list)
    if not gaps:
        print(f"ERROR: no valid gaps for epsilon_h={args.eh:g} eV", file=sys.stderr)
        return 1

    if not BUILD_PLOT.is_file() and not args.dry_run:
        rc = subprocess.call(
            ["cmake", "--build", "build", "-j", "--target", "ccdarksens_plot_dmelectron_limit"],
            cwd=ROOT,
        )
        if rc != 0:
            return rc

    mediators = list(MEDIATORS) if args.all_mediators else [args.mediator]
    si_ref = not args.no_si_reference
    rc = 0
    for med in mediators:
        if plot_one(med, args.eh, gaps, si_ref, args.dry_run) != 0:
            rc = 1
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
