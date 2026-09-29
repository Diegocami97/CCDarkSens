#!/usr/bin/env python3
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: compare_pydme_pattern_signal.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  compare_pydme_pattern_signal.py -- Compare CCDarkSens-style pattern signal
#  to pydme Verne/QEDark rate tables (m=2 MeV).
# ============================================================================

"""Compare CCDarkSens-style pattern signal to pydme Verne/QEDark rate tables (m=2 MeV)."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

ROOT = Path(__file__).resolve().parents[1]
ROI = [11, 21, 111, 31, 22, 211]
NE_FOLD = [1, 2, 3, 4, 5]


# ----------------------------------------------------------------------------
# load_eff
#   Read a (pattern, ne, efficiency) CSV into a DataFrame with integer pattern and ne columns.
# ----------------------------------------------------------------------------
def load_eff(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, comment="#")
    df["pattern"] = df["pattern"].astype(int)
    df["ne"] = df["ne"].astype(int)
    return df


def pydme_fold_pattern_rates(
    signal_path: Path, eff_path: Path, xsec: float, ref_xsec: float = 1e-30
) -> dict[int, np.ndarray]:
    """Replicate pydme qedark4dm.get_pattern_rates summed columns (events/g/day)."""
    eff = load_eff(eff_path)
    sig = pd.read_csv(signal_path)
    gamma = sig["gamma"].values
    xsec_col = sig["xsec"].values
    scale = xsec / ref_xsec if np.allclose(xsec_col, ref_xsec) else xsec / xsec_col

    out: dict[int, np.ndarray] = {}
    for p in sorted(eff["pattern"].unique()):
        sp = np.zeros(len(gamma))
        for ne in NE_FOLD:
            col = f"S{ne}"
            if col not in sig.columns:
                continue
            row = eff[(eff["pattern"] == p) & (eff["ne"] == ne)]
            if row.empty:
                continue
            e = float(row["Efficiency"].values[0])
            sp += sig[col].values * e
        out[int(p)] = sp * scale
    return out


def rate_to_counts_gday(rate_gday: np.ndarray, exposure_gday: float) -> float:
    """Mean rate (events/g/day) × exposure (g·day) → expected counts."""
    return float(np.mean(rate_gday)) * exposure_gday


# ----------------------------------------------------------------------------
# main
#   Compare the CCDarkSens-style pattern signal with the pydme (Verne/QEDark) signal and post-diffusion tables for one mass (default 2 MeV) and cross section, using the two reference efficiency tables.
# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sigma", type=float, default=1.0800523745162496e-36)
    ap.add_argument("--mass", type=float, default=2.0)
    args = ap.parse_args()

    mass = args.mass
    sigma = args.sigma
    pydme_root = ROOT / "collab_frameworks/pydme/pydme/data/rates_from_verne/FDM_n2/qedark"
    signal_2 = pydme_root / f"signal_mX{mass:.6f}_full_QED.csv"
    post_diff_2 = pydme_root / f"post_diffusion_signal_summed_mX{mass:.6f}_full_QED.csv"

    eff_paolo = ROOT / "data/efficiencies_paolo.csv"
    eff_dc = ROOT / "data/Efficiencies_patterns_Nsims1000000_DCTrue_alpha1.csv"

    # CCDarkSens config exposure (kg·year)
    exp_kg_y = 0.01523 * 85.356 / 365.0
    exp_gday = exp_kg_y * 1000.0 * 365.25

    ref = np.loadtxt(
        ROOT / "data/previous_limits/heavy_mediator/DAMIC-M_this_work_QEDark_hm.csv",
        delimiter=",",
    )

    print("=== pydme pattern signal comparison ===")
    print(f"m = {mass} MeV, sigma_e = {sigma:.3e} cm^2")
    print(f"CCDarkSens exposure: {exp_kg_y:.4e} kg·year = {exp_gday:.1f} g·day\n")

    if not signal_2.is_file():
        print(f"Missing {signal_2}")
        return

    for eff_name, eff_p in [("paolo eff", eff_paolo), ("DCTrue 1M eff", eff_dc)]:
        folded = pydme_fold_pattern_rates(signal_2, eff_p, sigma)
        print(f"--- pydme fold: {signal_2.name} + {eff_name} ---")
        print(f"  (rates in events/g/day, scaled from file xsec=1e-30 to {sigma:.3e})")
        for p in ROI:
            if p not in folded:
                continue
            r = folded[p]
            cnt = rate_to_counts_gday(r, exp_gday)
            print(f"  pat {p:3d}: mean rate={np.mean(r):.4e} ev/g/day  →  counts≈{cnt:.3f}")

    # Pre-folded pydme tables in repo
    for label, path in [
        ("pattern_signal_summed m=1 QCD (SRDM screening)", ROOT / "collab_frameworks/pydme/pydme/data/SRDM_rates/FDM_n2/qcdark/screening/pattern_signal_summed_mX1.000000_full_QCD.csv"),
        ("post_diffusion_summed m=2 QED (S1..S5 only)", post_diff_2),
    ]:
        if not path.is_file():
            print(f"\nSkip missing: {path}")
            continue
        df = pd.read_csv(path)
        print(f"\n--- {label} ---")
        print(f"  file: {path.name}")
        if "S11" in df.columns:
            s11 = df["S11"].values * (sigma / 1e-30)
            print(f"  S11: mean={np.mean(s11):.4e} ev/g/day  counts≈{rate_to_counts_gday(s11, exp_gday):.3f}")
            for p in ROI:
                col = f"S{p}"
                if col in df.columns:
                    r = df[col].values * (sigma / 1e-30)
                    print(f"  {col}: mean={np.mean(r):.4e}  counts≈{rate_to_counts_gday(r, exp_gday):.3f}")
        else:
            print("  columns:", [c for c in df.columns if c.startswith("S")][:8], "...")
            for ne in range(1, 6):
                col = f"S{ne}"
                if col in df.columns:
                    r = df[col].values * (sigma / 1e-30)
                    print(f"  {col}: mean={np.mean(r):.4e} ev/g/day")

    # Reference UL at mass
    rm, rs = ref[:, 0], ref[:, 1]
    print(f"\n=== Reference UL at m≈{mass} MeV ===")
    print(f"  sigma_e(ref) = {np.interp(mass, rm, rs):.3e}")

    print("\n=== CCDarkSens one-point (DCTrue eff, from prior run) ===")
    print("  pat 11 S_pat ≈ 13.8 counts at sigma ~ 1.08e-36 (fullgrid pipeline)")
    print("  Implied ref-UL scale (linear in sigma): S_pat × (ref_sigma/sigma)")


if __name__ == "__main__":
    main()
