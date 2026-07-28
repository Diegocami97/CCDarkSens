#!/usr/bin/env python3
"""Build stellar-only exclusion band CSV for dark-photon limit plots.

Combines Sun + red-giant stellar limits and removes parameter-space already
covered by the direct-detection envelope (DAMIC-M, SuperCDMS, XENONnT lo).

Output columns:
  mass_ev   lower edge of stellar-only fill (combined stellar envelope)
  eps_hi    upper edge; 0 = fill to top of plot (no DD cap at this mass)

Rows are omitted where the stellar envelope is already inside DD coverage
(y_stellar >= y_dd).
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

KEV_TO_EV = 1000.0
REPO = Path(__file__).resolve().parents[1]
DP = REPO / "data" / "previous_limits" / "dark_photon"


def load_txt(path: Path) -> list[tuple[float, float]]:
    pts: list[tuple[float, float]] = []
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.lower().startswith("mass"):
                continue
            parts = line.replace(",", " ").split()
            if len(parts) < 2:
                continue
            pts.append((float(parts[0]), float(parts[1])))
    return pts


def load_csv_cols(path: Path, x_col: int, y_col: int, x_scale: float = 1.0) -> list[tuple[float, float]]:
    pts: list[tuple[float, float]] = []
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = [p for p in line.split(",") if p]
            try:
                nums = [float(p) for p in parts]
            except ValueError:
                continue
            if len(nums) <= max(x_col, y_col):
                continue
            pts.append((nums[x_col] * x_scale, nums[y_col]))
    return pts


def load_xenon_lo(path: Path) -> list[tuple[float, float]]:
    pts: list[tuple[float, float]] = []
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "mass" in line.lower():
                continue
            parts = [p for p in line.split(",") if p]
            if len(parts) < 2:
                continue
            try:
                m = float(parts[0]) * KEV_TO_EV
                lo = float(parts[1])
            except ValueError:
                continue
            if m > 0 and lo > 0:
                pts.append((m, lo))
    return pts


def log_interp(pts: list[tuple[float, float]], m: float) -> float | None:
    if not pts:
        return None
    xs = np.array([p[0] for p in pts], dtype=float)
    ys = np.array([p[1] for p in pts], dtype=float)
    if m < xs.min() or m > xs.max():
        return None
    return float(10 ** np.interp(np.log10(m), np.log10(xs), np.log10(ys)))


def lower_envelope(curves: list[list[tuple[float, float]]], m: float) -> float | None:
    vals = [log_interp(c, m) for c in curves if c]
    vals = [v for v in vals if v is not None and v > 0]
    return min(vals) if vals else None


def build_masked_rows(m_min: float, m_max: float, n: int) -> list[tuple[float, float, float]]:
    sun = load_csv_cols(DP / "dark_photon_stellar_limits_SUN.csv", 1, 2, KEV_TO_EV)
    rg = load_csv_cols(DP / "dark_photon_stellar_limits_Red_Giant.csv", 1, 2, KEV_TO_EV)
    damic = load_txt(DP / "DAMIC-M_2025_DAMICmodel_HP.txt")
    scdms = load_txt(DP / "SCDMS_HP2024.txt")
    xenon = load_xenon_lo(DP / "dark_photon_xenonnt_hp_bracket.csv")

    dd_curves = [damic, scdms, xenon]
    rows: list[tuple[float, float, float]] = []

    for m in np.logspace(np.log10(m_min), np.log10(m_max), n):
        y_st = lower_envelope([sun, rg], m)
        if y_st is None:
            continue
        y_dd = lower_envelope(dd_curves, m)
        if y_dd is not None and y_st < y_dd:
            rows.append((float(m), float(y_st), float(y_dd)))
        else:
            # Full stellar fill to plot top; DD grey is drawn on top where applicable.
            rows.append((float(m), float(y_st), 0.0))

    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out",
        type=Path,
        default=DP / "dark_photon_stellar_limits_combined_masked.csv",
    )
    parser.add_argument("--m-min", type=float, default=0.01, help="eV")
    parser.add_argument("--m-max", type=float, default=30.0, help="eV")
    parser.add_argument("--n", type=int, default=2000)
    args = parser.parse_args()

    rows = build_masked_rows(args.m_min, args.m_max, args.n)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["mass_ev", "eps_lo", "eps_hi"])
        w.writerow(["# stellar-only band; eps_hi=0 => fill to plot top"])
        for m, lo, hi in rows:
            w.writerow([f"{m:.8g}", f"{lo:.8g}", f"{0.0 if hi == 0.0 else hi:.8g}"])

    print(f"Wrote {len(rows)} row(s) to {args.out}")


if __name__ == "__main__":
    main()
