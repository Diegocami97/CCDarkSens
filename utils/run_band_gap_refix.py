#!/usr/bin/env python3
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: run_band_gap_refix.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  run_band_gap_refix.py -- Non-destructive regeneration of band-gap pheno
#  limit scans + B/D overlays after the efficiency/ROI fix (single-pixel
#  epsilon in n_e space).
# ============================================================================
"""Re-run the band-gap pheno scans with the corrected n_e-space efficiency and
regenerate the heavy/light x B-thresh/D-equal limit overlays, WITHOUT touching
the pre-fix outputs.

  Scans  -> outputs/refix_eff_roi/<original_dir_name>/scan_dmelectron_pattern.root
  Plots  -> outplots/band_gap_pheno/step5_limits_refix/limit_sweep_*.pdf

Stages:
  python3 utils/run_band_gap_refix.py scan   # run all 25 redirected scans
  python3 utils/run_band_gap_refix.py plot   # 4 B/D overlays from refix scans
  python3 utils/run_band_gap_refix.py all    # scan then plot

NOTE: the *_eh_scan / *_gap_scan sweep plots in step5_limits also depend on the
scan_band_gap_2d_* grid (also n_e-space, also affected by the fix); those are NOT
regenerated here — handle separately if needed.
"""
from __future__ import annotations

import argparse
import glob
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REFIX_BASE = "outputs/refix_eff_roi"  # relative to ROOT (kept untracked-friendly)
# Redirected configs live in configs/ (same dir as originals) because the scan
# binary resolves data/rate paths relative to the config's directory; a suffix
# keeps them distinct and the source glob excludes them.
REFIX_SUFFIX = "__refix"
REFIX_PLOTS = ROOT / "outplots" / "band_gap_pheno" / "step5_limits_refix"
BUILD_SCAN = ROOT / "build" / "ccdarksens_scan_dmelectron_pattern"

CONFIG_GLOBS = [
    "configs/scan_band_gap*pheno*eh*.json",  # 25 pheno/light (B/D overlays)
    "configs/scan_band_gap_2d_*.json",       # 50 2D grid cells (eh/gap sweeps)
]


def redirected_configs() -> list[Path]:
    """Write copies of every scan config with run.outdir redirected to refix."""
    out: list[Path] = []
    seen: set[str] = set()
    for pattern in CONFIG_GLOBS:
        for f in sorted(glob.glob(str(ROOT / pattern))):
            stem = Path(f).stem
            if REFIX_SUFFIX in stem or stem in seen:
                continue  # skip our own copies / duplicates across globs
            seen.add(stem)
            d = json.loads(Path(f).read_text(encoding="utf-8"))
            name = Path(d["run"]["outdir"]).name  # e.g. scan_band_gap_0p5_eh0p5
            d["run"]["outdir"] = f"{REFIX_BASE}/{name}"
            dst = ROOT / "configs" / f"{stem}{REFIX_SUFFIX}.json"
            dst.write_text(json.dumps(d, indent=2) + "\n", encoding="utf-8")
            out.append(dst)
    return out


# ----------------------------------------------------------------------------
# cmd_scan
#   Sub-command scan: re-run every redirected config into the refix output tree with the scan binary, skipping outputs that exist unless --force; returns 1 if the binary is missing or any scan fails.
# ----------------------------------------------------------------------------
def cmd_scan(args: argparse.Namespace) -> int:
    if not BUILD_SCAN.is_file():
        print(f"ERROR: missing {BUILD_SCAN}; build it first.", file=sys.stderr)
        return 1
    cfgs = redirected_configs()
    force = getattr(args, "force", False)
    print(f"[refix] {len(cfgs)} redirected configs -> {REFIX_BASE}/ (force={force})", flush=True)
    failed: list[str] = []
    skipped = 0
    for i, cp in enumerate(cfgs, 1):
        d = json.loads(cp.read_text(encoding="utf-8"))
        out_root = ROOT / d["run"]["outdir"] / "scan_dmelectron_pattern.root"
        if out_root.is_file() and not force:
            skipped += 1
            continue
        print(f"\n===== [{i}/{len(cfgs)}] {cp.name} =====", flush=True)
        rc = subprocess.call([str(BUILD_SCAN), cp.relative_to(ROOT).as_posix()], cwd=ROOT)
        if rc != 0:
            failed.append(cp.name)
    print(f"\n[refix] scans complete. ran={len(cfgs)-skipped-len(failed)} "
          f"skipped(existing)={skipped} failed={failed}", flush=True)
    return 1 if failed else 0


def _redirect_path(p: Path) -> Path:
    """Map a default outputs/<name>/file path to the refix base."""
    return ROOT / REFIX_BASE / p.parent.name / p.name


# ----------------------------------------------------------------------------
# cmd_plot
#   Sub-command plot: redirect the scan-path resolvers of the Phase C plotting scripts to the refix outputs and make the gap and eps_h sweep plots.
# ----------------------------------------------------------------------------
def cmd_plot(_args: argparse.Namespace) -> int:
    sys.path.insert(0, str(ROOT / "utils"))
    import band_gap_scan_paths as bsp  # noqa: E402
    import run_band_gap_phase_c as h  # noqa: E402
    import plot_band_gap_fixed_gap_eh_sweep as eh_sweep  # noqa: E402
    import plot_band_gap_fixed_eh_gap_sweep as gap_sweep  # noqa: E402

    # Redirect the reference-overlay resolver (used by si_reference_root).
    _orig_bsp = bsp.scan_root_path

    def _bsp_path(mediator, gap, eh):
        return _redirect_path(_orig_bsp(mediator, gap, eh))

    bsp.scan_root_path = _bsp_path

    # B/D overlays (run_band_gap_phase_c has its own scan_root_path + OUT_PLOTS).
    h.OUT_PLOTS = REFIX_PLOTS
    _orig_h = h.scan_root_path
    h.scan_root_path = lambda gap, scen, med: _redirect_path(_orig_h(gap, scen, med))

    # Sweep plotters import scan_root_path into their own namespace.
    for mod in (eh_sweep, gap_sweep):
        mod.OUT_PLOTS = REFIX_PLOTS
        mod.scan_root_path = _bsp_path

    REFIX_PLOTS.mkdir(parents=True, exist_ok=True)
    rc = 0

    # 1) heavy/light x B-thresh/D-equal overlays
    for med in ("heavy", "light"):
        for tier in ("B-thresh", "D-equal"):
            if h.plot_limits_combo(med, tier, si_reference=True, dry_run=False) != 0:
                rc = 1

    # 2) fixed-eh gap scan at eh=3.8 (eh3p8_gap_scan), heavy + light
    gaps_38 = gap_sweep.gaps_for_eh(3.8, list(gap_sweep.GAP_GRID))
    for med in ("heavy", "light"):
        if gap_sweep.plot_one(med, 3.8, gaps_38, True, False) != 0:
            rc = 1

    # 3) fixed-gap eh scan for every gap (gap{X}_eh_scan), heavy + light
    for med in ("heavy", "light"):
        for gap in eh_sweep.GAP_GRID:
            if eh_sweep.plot_one(med, gap, False, True, False) != 0:
                rc = 1

    print(f"\n[refix] all overlays written to {REFIX_PLOTS}", flush=True)
    return rc


# ----------------------------------------------------------------------------
# main
#   Command line with the sub-commands scan, plot and all.
# ----------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = ap.add_subparsers(dest="command", required=True)
    p_scan = sub.add_parser("scan")
    p_scan.add_argument("--force", action="store_true", help="re-run even if output exists")
    p_scan.set_defaults(func=cmd_scan)
    sub.add_parser("plot").set_defaults(func=cmd_plot)
    p_all = sub.add_parser("all")
    p_all.add_argument("--force", action="store_true")
    p_all.set_defaults(func=lambda a: cmd_scan(a) or cmd_plot(a))
    args = ap.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
