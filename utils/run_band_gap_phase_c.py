#!/usr/bin/env python3
# ============================================================================
#  CCDarkSens — run_band_gap_phase_c
#  Orchestrate Phase C band-gap pheno limit scans (6 gaps × D-equal + B-thresh)
#
#  Author: Diego Venegas-Vargas
# ============================================================================
"""
Phase C: band-gap pheno limit scans (6 gaps x D-equal + B-thresh).

  python3 utils/run_band_gap_phase_c.py gen-configs
  python3 utils/run_band_gap_phase_c.py smoke
  python3 utils/run_band_gap_phase_c.py scan --tier B-thresh
  python3 utils/run_band_gap_phase_c.py scan --tier D-equal
  python3 utils/run_band_gap_phase_c.py scan --case gap0p1_B-thresh
  python3 utils/run_band_gap_phase_c.py plot-limits --tier B-thresh
  python3 utils/run_band_gap_phase_c.py plot-limits --all   # 4 PDFs: heavy/light x B/D

Fixed epsilon_h, scan E_gap (same as plot-limits --tier B-thresh for eh=3.8):
  python3 utils/plot_band_gap_fixed_eh_gap_sweep.py --eh 3.8 --all-mediators
"""

from __future__ import annotations

import argparse
import copy
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "utils"))
from band_gap_limit_plot import is_si_reference_cell, reference_scan_cli_args  # noqa: E402
from band_gap_plot_labels import limit_curve_label, limit_sweep_title  # noqa: E402

BUILD_SCAN = ROOT / "build" / "ccdarksens_scan_dmelectron_pattern"
BUILD_PLOT = ROOT / "build" / "ccdarksens_plot_dmelectron_limit"
GAPS = [0.1, 0.3, 0.5, 0.7, 0.9, 1.2]
EH_B = 3.8
OUT_PLOTS = ROOT / "outplots" / "band_gap_pheno" / "step5_limits"
MEDIATORS = ("heavy", "light")
PLOT_TIERS = ("B-thresh", "D-equal")


def gap_tag(g: float) -> str:
    return "gap1p2" if abs(g - 1.2) < 1e-9 else f"gap{g:.1f}".replace(".", "p")


def eh_tag(eh: float) -> str:
    return f"{eh:g}".replace(".", "p")


def case_id(gap: float, scenario: str) -> str:
    eh = gap if scenario == "D-equal" else EH_B
    gt = gap_tag(gap)
    gap_short = gt[3:] if gt.startswith("gap") else gt
    return f"{gap_short}_eh{eh_tag(eh)}"


def config_path(gap: float, scenario: str) -> Path:
    eh = gap if scenario == "D-equal" else EH_B
    gs = gap_tag(gap)[3:]
    return ROOT / "configs" / f"scan_band_gap_pheno_{gs}_eh{eh_tag(eh)}.json"


def list_cases(tier: str | None = None) -> list[tuple[float, str]]:
    out: list[tuple[float, str]] = []
    scenarios = ["B-thresh", "D-equal"] if tier is None else [tier.replace("_", "-")]
    for g in GAPS:
        for s in scenarios:
            if s in ("B-thresh", "D-equal"):
                out.append((g, s))
    return out


def run_cmd(cmd: list[str], dry_run: bool = False) -> int:
    print("[run]", " ".join(cmd))
    if dry_run:
        return 0
    return subprocess.call(cmd, cwd=ROOT)


def cmd_gen_configs(_: argparse.Namespace) -> int:
    return subprocess.call([sys.executable, "utils/gen_band_gap_pheno_scan_configs.py"], cwd=ROOT)


def cmd_smoke(args: argparse.Namespace) -> int:
    tpl = json.loads(config_path(0.1, "B-thresh").read_text(encoding="utf-8"))
    cfg = copy.deepcopy(tpl)
    cfg["_comment"] = "Phase C smoke: 3x3 grid, B-thresh gap0p1"
    cfg["run"]["label"] = "scan_band_gap_pheno_smoke"
    cfg["run"]["outdir"] = "outputs/scan_band_gap_smoke"
    cfg["model"]["grid"]["mchi_MeV"]["logspace"] = {
        "start": 5.0,
        "stop": 50.0,
        "num": 3,
        "endpoint": True,
    }
    cfg["model"]["grid"]["sigma_e_cm2"]["logspace"] = {
        "start_exp": -44,
        "stop_exp": -34,
        "num": 5,
        "endpoint": True,
    }
    cfg["model"]["grid"]["mchi_MeV"]["logspace"]["start"] = 8.0
    cfg["model"]["grid"]["mchi_MeV"]["logspace"]["stop"] = 20.0
    smoke_path = ROOT / "configs" / "scan_band_gap_pheno_smoke.json"
    smoke_path.write_text(json.dumps(cfg, indent=2) + "\n", encoding="utf-8")
    print(f"[ok] wrote {smoke_path}")

    if not BUILD_SCAN.is_file() and not args.dry_run:
        print("Building ccdarksens_scan_dmelectron_pattern...")
        rc = run_cmd(["cmake", "--build", "build", "-j", "--target", "ccdarksens_scan_dmelectron_pattern"])
        if rc != 0:
            return rc

    rc = run_cmd([str(BUILD_SCAN), "configs/scan_band_gap_pheno_smoke.json"], args.dry_run)
    if rc != 0:
        return rc

    root_file = ROOT / "outputs/scan_band_gap_smoke/scan_dmelectron_pattern.root"
    if not args.dry_run and not root_file.is_file():
        print(f"ERROR: missing {root_file}", file=sys.stderr)
        return 1

    if BUILD_PLOT.is_file() or args.dry_run:
        OUT_PLOTS.mkdir(parents=True, exist_ok=True)
        rc = run_cmd(
            [
                str(BUILD_PLOT),
                "outputs/scan_band_gap_smoke/scan_dmelectron_pattern.root",
                "smoke B-thresh 0.1",
                "heavy",
                "--from-qhist",
                "--batch",
                "--title",
                "Phase C smoke (3x3)",
                "--out-pdf",
                str(OUT_PLOTS / "limit_smoke.pdf"),
                "--out-csv",
                str(OUT_PLOTS / "limit_smoke.csv"),
            ],
            args.dry_run,
        )
        if rc != 0:
            print("WARN: plot step failed (build ccdarksens_plot_dmelectron_limit if needed)", file=sys.stderr)
    else:
        print("WARN: skip plot — build/ccdarksens_plot_dmelectron_limit not found")

    scan_ok = root_file.is_file() if not args.dry_run else True
    if not scan_ok:
        return 1
    print("\nSmoke scan OK (ROOT written).")
    if rc != 0:
        print("Smoke limit plot skipped or failed — OK if sigma grid too narrow; full grid uses -46:-26.")
    return 0


def cmd_scan(args: argparse.Namespace) -> int:
    if args.case:
        # gap0p1_B-thresh
        parts = args.case.split("_", 1)
        if len(parts) != 2:
            print("ERROR: --case format gap0p1_B-thresh", file=sys.stderr)
            return 1
        gt, scen = parts[0], parts[1].replace("_", "-")
        gap = next(g for g in GAPS if gap_tag(g) == gt)
        cases = [(gap, scen)]
    else:
        cases = list_cases(args.tier)

    if not BUILD_SCAN.is_file() and not args.dry_run:
        rc = run_cmd(["cmake", "--build", "build", "-j", "--target", "ccdarksens_scan_dmelectron_pattern"])
        if rc != 0:
            return rc

    failed = []
    for gap, scen in cases:
        cp = config_path(gap, scen)
        if not cp.is_file():
            print(f"ERROR: missing {cp}", file=sys.stderr)
            failed.append(str(cp))
            continue
        rc = run_cmd([str(BUILD_SCAN), cp.relative_to(ROOT).as_posix()], args.dry_run)
        if rc != 0:
            failed.append(case_id(gap, scen))

    if failed:
        print(f"Failed: {failed}", file=sys.stderr)
        return 1
    print(f"Completed {len(cases)} scan(s).")
    return 0


def scan_root_path(gap: float, scenario: str, mediator: str) -> Path:
    cid = case_id(gap, scenario)
    prefix = "scan_band_gap_light" if mediator == "light" else "scan_band_gap"
    return ROOT / "outputs" / f"{prefix}_{cid}" / "scan_dmelectron_pattern.root"


def plot_limits_combo(
    mediator: str, tier: str, *, si_reference: bool = True, dry_run: bool = False
) -> int:
    """One overlay PDF: six gaps for the given mediator and tier (B-thresh or D-equal)."""
    cases = list_cases(tier)
    paths: list[str] = []
    labels: list[str] = []
    for gap, scen in cases:
        root = scan_root_path(gap, scen, mediator)
        if not root.is_file():
            print(f"WARN: skip missing {root}")
            continue
        eh = gap if scen == "D-equal" else EH_B
        if si_reference and is_si_reference_cell(gap, eh):
            continue
        paths.append(root.relative_to(ROOT).as_posix())
        labels.append(limit_curve_label(gap, scen))

    if len(paths) < 2:
        print(f"ERROR: {mediator} / {tier}: need at least 2 completed ROOT files", file=sys.stderr)
        return 1

    tier_slug = tier.replace("-", "_")
    med_slug = "" if mediator == "heavy" else f"{mediator}_"
    OUT_PLOTS.mkdir(parents=True, exist_ok=True)
    out_pdf = OUT_PLOTS / f"limit_sweep_{med_slug}{tier_slug}.pdf"
    print(f"[plot-limits] {mediator} {tier} -> {out_pdf.name}")
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
            limit_sweep_title(mediator, tier=tier, root=True),
            "--out-pdf",
            str(out_pdf),
            "--out-csv",
            str(OUT_PLOTS / f"limit_sweep_{med_slug}{tier_slug}.csv"),
        ]
    )
    return run_cmd(cmd, dry_run)


def cmd_plot_limits(args: argparse.Namespace) -> int:
    if not BUILD_PLOT.is_file() and not args.dry_run:
        rc = run_cmd(["cmake", "--build", "build", "-j", "--target", "ccdarksens_plot_dmelectron_limit"])
        if rc != 0:
            return rc

    si_reference = not getattr(args, "no_si_reference", False)

    if args.all:
        combos = [(m, t) for m in MEDIATORS for t in PLOT_TIERS]
    else:
        if args.tier is None:
            print("ERROR: pass --tier B-thresh|D-equal or use --all", file=sys.stderr)
            return 1
        combos = [(args.mediator, args.tier)]

    rc = 0
    for mediator, tier in combos:
        if plot_limits_combo(mediator, tier, si_reference=si_reference, dry_run=args.dry_run) != 0:
            rc = 1
    return rc


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true")
    sub = ap.add_subparsers(dest="command", required=True)

    sub.add_parser("gen-configs", help="Write all 12 pheno scan JSONs")

    p_smoke = sub.add_parser("smoke", help="3x3 grid + optional limit plot")
    p_smoke.set_defaults(func=cmd_smoke)

    p_scan = sub.add_parser("scan", help="Run full grid scan(s)")
    p_scan.add_argument("--tier", choices=["B-thresh", "D-equal"], help="Run all cases in tier")
    p_scan.add_argument("--case", help="Single case e.g. gap0p1_B-thresh")
    p_scan.set_defaults(func=cmd_scan)

    p_plot = sub.add_parser("plot-limits", help="Overlay limit curves")
    p_plot.add_argument("--tier", choices=list(PLOT_TIERS), help="B-thresh or D-equal row")
    p_plot.add_argument(
        "--mediator",
        choices=list(MEDIATORS),
        default="heavy",
        help="heavy: outputs/scan_band_gap_*; light: outputs/scan_band_gap_light_*",
    )
    p_plot.add_argument(
        "--all",
        action="store_true",
        help="Plot all four overlays: heavy/light x B-thresh/D-equal",
    )
    p_plot.add_argument(
        "--no-si-reference",
        action="store_true",
        help="Do not overlay Si reference limit (E_gap=1.2 eV, epsilon_h=3.8 eV)",
    )
    p_plot.set_defaults(func=cmd_plot_limits)

    args = ap.parse_args()
    if args.command == "gen-configs":
        return cmd_gen_configs(args)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
