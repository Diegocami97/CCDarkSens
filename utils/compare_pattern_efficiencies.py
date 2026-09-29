#!/usr/bin/env python3
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: compare_pattern_efficiencies.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  compare_pattern_efficiencies.py -- Compare CCDarkSens pattern-efficiency
#  tables against LBC reference data
# ============================================================================
"""Compare pattern-efficiency tables and LBC data counts."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ROI = [11, 21, 111, 31, 22, 211]


# ----------------------------------------------------------------------------
# load_eff
#   Read a (pattern, ne, efficiency) CSV into a DataFrame with integer pattern and ne columns.
# ----------------------------------------------------------------------------
def load_eff(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, comment="#")
    df["pattern"] = df["pattern"].astype(int)
    df["ne"] = df["ne"].astype(int)
    return df


# ----------------------------------------------------------------------------
# main
#   Compare the reference pattern-efficiency tables (Paolo's and the 1M-simulation DCTrue one) and print their metadata and differences.
# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--paolo", default=str(ROOT / "data/efficiencies_paolo.csv"))
    ap.add_argument(
        "--dctrue",
        default=str(ROOT / "data/Efficiencies_patterns_Nsims1000000_DCTrue_alpha1.csv"),
    )
    ap.add_argument("--data", default=str(ROOT / "data/Final_Combined_Image_Data.csv"))
    args = ap.parse_args()

    pa = load_eff(Path(args.paolo))
    dc = load_eff(Path(args.dctrue))

    print("=== Simulation metadata (comment line) ===")
    for label, p in [("paolo", args.paolo), ("DCTrue 1M", args.dctrue)]:
        with open(p) as f:
            print(f"  {label}: {f.readline().strip()}")

    print("\n=== ROI: P(pattern | n_e), ratio paolo / DCTrue ===")
    print(f"{'pat':>4} {'ne':>2}  {'paolo':>10} {'DCTrue':>10} {'ratio':>8}")
    for p in ROI:
        for ne in range(1, 6):
            vp = pa.query("pattern == @p and ne == @ne")["Efficiency"].values
            vd = dc.query("pattern == @p and ne == @ne")["Efficiency"].values
            vp = float(vp[0]) if len(vp) else 0.0
            vd = float(vd[0]) if len(vd) else 0.0
            r = vp / vd if vd > 0 else (np.inf if vp > 0 else np.nan)
            print(f"{p:4d} {ne:2d}  {vp:10.6g} {vd:10.6g} {r:8.3f}")

    print("\n=== Σ_{ne=2..5} efficiency per ROI pattern ===")
    for name, df in [("paolo", pa), ("DCTrue", dc)]:
        print(name)
        for p in ROI:
            s = df[(df["pattern"] == p) & (df["ne"] >= 2) & (df["ne"] <= 5)]["Efficiency"].sum()
            print(f"  pat {p:3d}: {s:.4f}")

    print("\n=== Pattern 11 (dominant signal channel) ===")
    for name, df in [("paolo", pa), ("DCTrue", dc)]:
        s = df[(df["pattern"] == 11) & (df["ne"] >= 2) & (df["ne"] <= 5)]["Efficiency"].sum()
        print(f"  {name} sum ne=2..5: {s:.4f}")

    data = pd.read_csv(args.data)
    pat_cols = sorted(c for c in data.columns if "Count_Candidate" in c)
    print(f"\n=== {Path(args.data).name} ===")
    print(f"  images: {len(data)}")
    for c in pat_cols:
        pid = c.replace("Count_Candidate_", "")
        print(f"  pat {pid}: {int(data[c].sum())}")
    print(f"  total: {int(data[pat_cols].sum().sum())}")

    Vpix = 0.0015 * 0.0015 * 0.0669
    rho = 2.33
    Mpix_100 = Vpix * rho * 100
    exp_pydme_gday = (data["Nusedpix"] * data["texp"] * Mpix_100).sum()
    print(f"  pydme exposure (binning=100): {exp_pydme_gday:.1f} g·day")
    print("  CCDarkSens scan uses config exposure ~3.56e-3 kg·year, not this value.")


if __name__ == "__main__":
    main()
