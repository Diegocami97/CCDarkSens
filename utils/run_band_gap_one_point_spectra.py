#!/usr/bin/env python3
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: run_band_gap_one_point_spectra.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  run_band_gap_one_point_spectra.py -- Generate and compare signal spectra
#  at one (E_gap, ε_h) diagnostic point
# ============================================================================
"""
Run one (m_chi, sigma_e) point per band-gap / ionization case and dump spectra.

Uses ccdarksens_scan_dmelectron_pattern with dump_point_spectra_root=true
and experiment.observable_bins = "ne" (electron / n_e space, not pattern).

ROOT histograms per case:
  dRdE__*        — before ionization (QCDark2 rate table)
  S_true_ne__*   — after ChargeIonization::FoldToNe
  S_obs_ne__*    — after pattern efficiency on n_e (still n_e bins)

Author: Diego Venegas-Vargas
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
from band_gap_plot_labels import pheno_param_label_mpl  # noqa: E402

SCAN_BIN = ROOT / "build" / "ccdarksens_scan_dmelectron_pattern"
MANIFEST = ROOT / "configs" / "band_gap_one_point_spectra_manifest.json"
TEMPLATE = ROOT / "configs" / "scan_dmelectron_band_gap_study_0p3.json"
BUILD_P100K = ROOT / "utils" / "build_p100K_scaled.py"


# ----------------------------------------------------------------------------
# gap_tag
#   File-name tag of a band gap, e.g. "gap0p7" (1.2 eV gives "gap1p2").
# ----------------------------------------------------------------------------
def gap_tag(g: float) -> str:
    return "gap1p2" if abs(g - 1.2) < 1e-9 else f"gap{g:.1f}".replace(".", "p")


# ----------------------------------------------------------------------------
# eh_tag
#   Electron-hole pair energy formatted for file names with '.' replaced by 'p'.
# ----------------------------------------------------------------------------
def eh_tag(eh: float) -> str:
    return f"{eh:g}".replace(".", "p")


# ----------------------------------------------------------------------------
# ionization_csv_path
#   Path of the scaled p100K ionization table of a (gap, eps_h) pair.
# ----------------------------------------------------------------------------
def ionization_csv_path(gap_eV: float, eh_eV: float) -> Path:
    return ROOT / "data" / f"p100K_{gap_tag(gap_eV)}_eh{eh_tag(eh_eV)}.csv"


# ----------------------------------------------------------------------------
# build_cases
#   Cases to run from the manifest: every gap in the D-equal (eps_h = gap) and B-thresh (fixed eps_h) scenarios, with their rate directories and ionization tables.
# ----------------------------------------------------------------------------
def build_cases(manifest: dict) -> list[dict]:
    gaps = [float(g) for g in manifest.get("gaps_eV", [])]
    scenarios = manifest.get("scenarios", ["D-equal", "B-thresh"])
    eh_thresh = float(manifest.get("eh_pair_B_thresh_eV", 3.8))
    cases: list[dict] = []

    for g in gaps:
        tag = gap_tag(g)
        rates_dir = f"data/qcdark2_rates/Si/heavy/Si_fast_{tag}"
        for scen in scenarios:
            if scen == "D-equal":
                eh = g
            elif scen == "B-thresh":
                eh = eh_thresh
            else:
                continue
            label = pheno_param_label_mpl(g, eh)
            ion = ionization_csv_path(g, eh)
            cases.append(
                {
                    "id": f"{tag}_{scen}",
                    "label": label,
                    "band_gap_eV": g,
                    "eh_pair_eV": eh,
                    "scenario": scen,
                    "rates_dir": rates_dir,
                    "ionization_csv": str(ion.relative_to(ROOT)),
                }
            )
    return cases


# ----------------------------------------------------------------------------
# ensure_p100k_tables
#   Build the p100K tables that are missing for the cases (only printed with dry_run); returns 0 on success.
# ----------------------------------------------------------------------------
def ensure_p100k_tables(cases: list[dict], dry_run: bool) -> int:
    missing = []
    for c in cases:
        p = ROOT / c["ionization_csv"]
        if not p.is_file():
            missing.append(c)
    if not missing:
        return 0
    print(f"Building {len(missing)} missing p100K table(s)...")
    for c in missing:
        cmd = [
            sys.executable,
            str(BUILD_P100K),
            "--band-gap-eV",
            str(c["band_gap_eV"]),
            "--eh-pair-eV",
            str(c["eh_pair_eV"]),
            "--scenario",
            c["scenario"],
        ]
        print(f"  {' '.join(cmd)}")
        if dry_run:
            continue
        rc = subprocess.call(cmd, cwd=ROOT)
        if rc != 0:
            return rc
        if not ionization_csv_path(c["band_gap_eV"], c["eh_pair_eV"]).is_file():
            print(f"ERROR: expected {c['ionization_csv']} after build", file=sys.stderr)
            return 1
    return 0


# ----------------------------------------------------------------------------
# make_config
#   One-point scan config of a case: n_e-space dump of dR/dE and S_true(n_e) at a single (mass, cross-section) point, from the template.
# ----------------------------------------------------------------------------
def make_config(base: dict, case: dict, point: dict, manifest: dict, outdir: Path) -> dict:
    cfg = copy.deepcopy(base)
    cfg["_comment"] = (
        f"One-point spectra dump (n_e space): {case['label']}. "
        "dRdE + S_true(n_e); no pattern-space observable."
    )
    obs = manifest.get("run_defaults", {}).get("observable_bins", "ne")
    cfg["run"] = {
        "label": case["id"],
        "outdir": str(outdir.relative_to(ROOT)),
        "n_toys": 0,
        "rng_seed": 12345,
        "verbosity": 1,
        "cl": 0.9,
        "test_stat": "PLR",
        "dump_point_spectra_root": True,
        "use_profile_likelihood": False,
        "data_path": "",
        "background_source": "dc_flat",
        "background_model": "scale",
        "profile_minimizer": "brent",
        "pydme_style_ul": False,
    }
    mchi = point["mchi_MeV"]
    sigma = point["sigma_e_cm2"]
    fmt = point.get("format", {"mchi": ".6f", "sigma": ".1e"})
    cfg["model"]["rates_dir"] = case["rates_dir"]
    cfg["model"]["grid"] = {
        "mchi_MeV": {"values": [mchi]},
        "sigma_e_cm2": {"values": [sigma]},
        "format": fmt,
    }
    cfg.setdefault("experiment", {})["observable_bins"] = obs
    resp = cfg.setdefault("response", {})
    ci = resp.setdefault("charge_ionization", {})
    ci["table_csv"] = case["ionization_csv"]
    ci["band_gap_eV"] = case["band_gap_eV"]
    ci["eh_pair_eV"] = case["eh_pair_eV"]
    ci["scenario"] = case["scenario"]
    return cfg


# ----------------------------------------------------------------------------
# run_case
#   Run the scan binary on one config (only printed with dry_run) and return its exit code.
# ----------------------------------------------------------------------------
def run_case(cfg_path: Path, dry_run: bool) -> int:
    cmd = [str(SCAN_BIN), str(cfg_path.relative_to(ROOT))]
    print(f"[run] {' '.join(cmd)}")
    if dry_run:
        return 0
    return subprocess.call(cmd, cwd=ROOT)


# ----------------------------------------------------------------------------
# plot_results
#   Run the ROOT comparison macro and then the p100K plotting script on the results directory.
# ----------------------------------------------------------------------------
def plot_results(results_dir: Path, dry_run: bool) -> int:
    macro = "utils/plot_band_gap_one_point_spectra_compare.cc"
    rel = results_dir.relative_to(ROOT).as_posix()
    cmd = ["root", "-l", "-b", "-q", f'{macro}("{rel}")']
    print(f"[plot] {' '.join(cmd)}")
    if dry_run:
        return 0
    rc = subprocess.call(cmd, cwd=ROOT)
    if rc != 0:
        return rc
    p100k = ["python3", "utils/plot_band_gap_one_point_p100K.py", "--manifest", str(MANIFEST)]
    print(f"[plot] {' '.join(p100k)}")
    return subprocess.call(p100k, cwd=ROOT)


# ----------------------------------------------------------------------------
# main
#   Command line: prepare the p100K tables and the one-point configs, run the scans and plot the results; options select cases, skip steps (--plot-only, --no-plot, --skip-p100k-build) or only print (--dry-run).
# ----------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", type=Path, default=MANIFEST)
    ap.add_argument("--template", type=Path, default=TEMPLATE)
    ap.add_argument("--cases", nargs="*", help="Subset of case ids (default: all)")
    ap.add_argument("--plot-only", action="store_true", help="Skip scans, only plot")
    ap.add_argument("--no-plot", action="store_true", help="Skip plot step")
    ap.add_argument("--skip-p100k-build", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    base = json.loads(args.template.read_text(encoding="utf-8"))
    point = manifest["reference_point"]
    out_base = ROOT / manifest["run_defaults"]["outdir_base"]

    all_cases = build_cases(manifest)
    if args.cases:
        ids = set(args.cases)
        all_cases = [c for c in all_cases if c["id"] in ids]
        missing = ids - {c["id"] for c in all_cases}
        if missing:
            print(f"WARNING: unknown case ids: {sorted(missing)}", file=sys.stderr)

    obs = manifest.get("run_defaults", {}).get("observable_bins", "ne")
    print(f"Reference point: m_chi={point['mchi_MeV']} MeV, sigma_e={point['sigma_e_cm2']} cm^2")
    print(f"Observable: {obs} (not pattern)")
    print(f"Cases ({len(all_cases)}):")
    for c in all_cases:
        print(f"  {c['id']:22s}  {c['label']}")

    if not args.plot_only:
        if not SCAN_BIN.is_file():
            print(f"ERROR: build scan binary first: {SCAN_BIN}", file=sys.stderr)
            return 1
        if not args.skip_p100k_build:
            rc = ensure_p100k_tables(all_cases, args.dry_run)
            if rc != 0:
                return rc

        for case in all_cases:
            case_out = out_base / case["id"]
            case_out.mkdir(parents=True, exist_ok=True)
            cfg_path = ROOT / "configs" / f"band_gap_one_point_{case['id']}.json"
            cfg = make_config(base, case, point, manifest, case_out)
            cfg["run"]["outdir"] = str(case_out.relative_to(ROOT))
            cfg_path.write_text(json.dumps(cfg, indent=2) + "\n", encoding="utf-8")
            rc = run_case(cfg_path, args.dry_run)
            if rc != 0:
                print(f"ERROR: scan failed for {case['id']} (exit {rc})", file=sys.stderr)
                return rc

    if not args.no_plot and not args.dry_run:
        rc = plot_results(out_base, False)
        if rc != 0:
            return rc

    print(f"\nDone. Results under {out_base}/")
    print("  Per case: <case_id>/scan_dmelectron_pattern.root")
    print("  Plots:    outplots/band_gap_one_point_spectra/")
    print("    before_after__<case_id>.pdf              — single scenario")
    print("    compare_D-equal_vs_B-thresh__gap*.pdf    — per gap, two curves")
    print("    compare_*_D-equal_all_gaps.pdf           — all gaps, D-equal only")
    print("    compare_*_B-thresh_all_gaps.pdf          — all gaps, B-thresh only")
    print("    compare_ionization_*_ne1to5.pdf          — folded S_true(n_e) from scans (not P tables)")
    print("    p100K_Pne/compare_p100K_*.pdf            — P(n_e|E) from data/p100K_gap*_eh*.csv")
    print("\nRebuild scan app after ConfigManager changes:")
    print("  cmake --build build -j --target ccdarksens_scan_dmelectron_pattern")
    return 0


if __name__ == "__main__":
    sys.exit(main())
