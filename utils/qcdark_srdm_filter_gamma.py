#!/usr/bin/env python3
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: qcdark_srdm_filter_gamma.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  qcdark_srdm_filter_gamma.py -- Filters QCDark SRDM pattern-signal CSV rows
#  to a target γ value and can emit matching sigma-grid JSON snippets.
# ============================================================================

"""
Filter QCDark SRDM pattern-space signal CSVs by a fixed gamma value.

Input format (example):
  gamma,xsec,S1,S2,...,S311

This script keeps only rows where:
  abs(gamma - gamma_target) <= gamma_tol

It writes a new CSV (optionally stripping the gamma column) and can also
emit a JSON snippet with the unique sorted xsec values so you can reuse
them as a sigma grid.
"""

import argparse
import csv
import json
import math
import os
import re
import glob
from typing import List, Optional


# ----------------------------------------------------------------------------
# _parse_mX_from_filename
#   DM mass parsed from a file name containing "mX<value>", or None.
# ----------------------------------------------------------------------------
def _parse_mX_from_filename(path: str) -> Optional[float]:
    # Matches e.g. .../pattern_signal_summed_mX0.010000_full_QCD.csv
    m = re.search(r"mX([0-9]+(?:\.[0-9]+)?)", os.path.basename(path))
    if not m:
        return None
    try:
        return float(m.group(1))
    except ValueError:
        return None


# ----------------------------------------------------------------------------
# main
#   Filter the SRDM pattern-signal CSVs to a target gamma value (within --gamma-tol), for a single file (--in-csv/--out-csv) or a batch (--in-glob/--out-dir).
# ----------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-csv", default="", help="Input CSV path (single-file mode)")
    ap.add_argument(
        "--in-glob",
        default="",
        help="Input glob for batch mode, e.g. 'data/qcdark_srdm/pattern_signal_summed_mX*_full_QCD.csv'.",
    )
    ap.add_argument("--gamma", required=True, type=float, help="Target gamma value")
    ap.add_argument("--gamma-tol", default=1e-12, type=float, help="Absolute tolerance on gamma match")
    ap.add_argument("--out-csv", default="", help="Output CSV path (single-file mode)")
    ap.add_argument("--out-dir", default="outputs/qcdark_srdm_filtered_gamma", help="Output dir for batch mode")
    ap.add_argument("--strip-gamma", action="store_true", help="Remove gamma column in output")
    ap.add_argument(
        "--out-sigma-json",
        default="",
        help="Optional JSON path to write a {\"sigma_e_cm2\": {\"values\": [...]}} snippet from unique xsec.",
    )
    ap.add_argument(
        "--write-sigma-json-per-file",
        action="store_true",
        help="When in batch mode, write one sigma grid json per input file into out-dir.",
    )
    args = ap.parse_args()

    gamma = float(args.gamma)
    gamma_dec = f"{gamma:.10f}".rstrip("0").rstrip(".")
    gamma_tag = gamma_dec.replace(".", "p").replace("-", "m")

    in_files: List[str]
    single_mode = bool(args.in_csv) or bool(args.out_csv)
    if args.in_glob:
        in_files = sorted(glob.glob(args.in_glob))
        if not in_files:
            raise RuntimeError(f"No files matched --in-glob={args.in_glob}")
        single_mode = False
    else:
        if not args.in_csv:
            raise RuntimeError("Provide either --in-csv (single mode) or --in-glob (batch mode).")
        in_files = [args.in_csv]
        single_mode = True

    if single_mode:
        in_csv = args.in_csv
        out_csv = args.out_csv
        if not out_csv:
            raise RuntimeError("In single-file mode you must provide --out-csv.")
        if not os.path.exists(in_csv):
            raise FileNotFoundError(in_csv)

    else:
        os.makedirs(args.out_dir, exist_ok=True)

    def process_one(in_csv: str, out_csv: str, sigma_json_path: str = "") -> None:
        if not os.path.exists(in_csv):
            raise FileNotFoundError(in_csv)

        with open(in_csv, "r", newline="") as f:
            reader = csv.DictReader(f)
            if reader.fieldnames is None:
                raise RuntimeError("Missing CSV header")

            fieldnames = list(reader.fieldnames)
            if "gamma" not in fieldnames or "xsec" not in fieldnames:
                raise RuntimeError(f"Expected columns 'gamma' and 'xsec'. Found: {fieldnames}")

            gamma_col = "gamma"
            xsec_col = "xsec"
            pattern_cols = [c for c in fieldnames if c not in (gamma_col, xsec_col)]

            kept_rows: List[dict] = []
            xsec_values = set()

            for row in reader:
                try:
                    g = float(row[gamma_col])
                    s = float(row[xsec_col])
                except Exception:
                    continue

                if abs(g - gamma) <= args.gamma_tol:
                    kept_rows.append(row)
                    xsec_values.add(s)

        if not kept_rows:
            raise RuntimeError(
                f"No rows matched gamma={gamma} within tol={args.gamma_tol}. "
                f"Check CSV uses the same gamma convention/precision."
            )

        kept_rows.sort(key=lambda r: float(r[xsec_col]))
        unique_xsec_sorted = sorted(xsec_values)

        if args.strip_gamma:
            out_fieldnames = [xsec_col] + pattern_cols
        else:
            out_fieldnames = [gamma_col, xsec_col] + pattern_cols

        os.makedirs(os.path.dirname(out_csv) or ".", exist_ok=True)
        with open(out_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=out_fieldnames)
            writer.writeheader()
            for row in kept_rows:
                if args.strip_gamma:
                    out_row = {xsec_col: row[xsec_col]}
                    for c in pattern_cols:
                        out_row[c] = row[c]
                    writer.writerow(out_row)
                else:
                    writer.writerow(row)

        mX = _parse_mX_from_filename(in_csv)
        print(f"[filter-gamma] input={in_csv}")
        print(f"[filter-gamma] parsed mX={mX if mX is not None else 'n/a'} (from filename)")
        print(f"[filter-gamma] kept_rows={len(kept_rows)}")
        print(f"[filter-gamma] unique xsec values={len(unique_xsec_sorted)}")
        print(f"[filter-gamma] wrote={out_csv}")

        if sigma_json_path:
            snippet = {"sigma_e_cm2": {"values": unique_xsec_sorted}}
            os.makedirs(os.path.dirname(sigma_json_path) or ".", exist_ok=True)
            with open(sigma_json_path, "w") as f:
                json.dump(snippet, f, indent=2)
            print(f"[filter-gamma] wrote sigma grid snippet={sigma_json_path}")

    if single_mode:
        # Optional sigma json path (single-mode uses user-provided --out-sigma-json)
        process_one(in_files[0], out_csv, sigma_json_path=args.out_sigma_json if args.out_sigma_json else "")
    else:
        print(f"[filter-gamma] batch mode: matched {len(in_files)} input files")
        for in_csv in in_files:
            base = os.path.basename(in_csv)
            if base.lower().endswith(".csv"):
                base_noext = base[:-4]
            else:
                base_noext = base

            out_name = base_noext + f"_gamma{gamma_tag}.csv"
            out_path = os.path.join(args.out_dir, out_name)

            sigma_json_path = ""
            if args.write_sigma_json_per_file:
                sigma_json_path = os.path.join(
                    args.out_dir,
                    base_noext + f"_sigma_grid_gamma{gamma_tag}.json",
                )

            process_one(in_csv, out_path, sigma_json_path=sigma_json_path)


if __name__ == "__main__":
    main()

