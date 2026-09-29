#!/usr/bin/env python3
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: extract_band_gap_2d_sensitivity.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  extract_band_gap_2d_sensitivity.py -- Extract σ_UL on the 2D (E_gap, ε_h)
#  grid for heatmap generation
# ============================================================================
"""
Extract sigma_UL at a chosen m_chi for the 2D (gap, eh) grid from scan ROOT outputs.

Output (per mediator and m_chi):
  outplots/band_gap_pheno/2d_surface/sensitivity_2d_{mediator}_mchi{mchi_tag}.npy
  outplots/band_gap_pheno/2d_surface/sensitivity_2d_{mediator}_mchi{mchi_tag}.csv
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
import uproot

import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "utils"))
from band_gap_scan_paths import EH_GRID as EHS  # noqa: E402
from band_gap_scan_paths import ev_tag, scan_root_path as output_root_path  # noqa: E402

OUTDIR = ROOT / "outplots" / "band_gap_pheno" / "2d_surface"

GAPS = [0.1, 0.3, 0.5, 0.7, 0.9, 1.2]
Q_THRESHOLD = 2.71


# ----------------------------------------------------------------------------
# mchi_tag
#   Mass formatted for file names: '<int>MeV' for integers, otherwise with '.' replaced by 'p'.
# ----------------------------------------------------------------------------
def mchi_tag(mev: float) -> str:
    if abs(mev - round(mev)) < 1e-9:
        return f"{int(round(mev))}MeV"
    return f"{mev:g}MeV".replace(".", "p")


# ----------------------------------------------------------------------------
# sensitivity_paths
#   The .npy and .csv output paths of the 2D sensitivity map for a mediator and mass.
# ----------------------------------------------------------------------------
def sensitivity_paths(mediator: str, mchi_mev: float) -> tuple[Path, Path]:
    tag = mchi_tag(mchi_mev)
    base = OUTDIR / f"sensitivity_2d_{mediator}_mchi{tag}"
    return base.with_suffix(".npy"), base.with_suffix(".csv")


# ----------------------------------------------------------------------------
# is_valid
#   A (gap, eh) cell is physical only if eps_h >= E_gap.
# ----------------------------------------------------------------------------
def is_valid(gap: float, eh: float) -> bool:
    return eh >= gap


# ----------------------------------------------------------------------------
# ul_from_qhist
#   Upper limit at the requested mass read from the q(m_chi, sigma) histogram of a scan file; NaN if the histogram is missing.
# ----------------------------------------------------------------------------
def ul_from_qhist(root_path: Path, target_mchi_mev: float) -> float:
    with uproot.open(root_path) as f:
        h = None
        for name in ("q_mchi_sigma_pattern", "q_mchi_sigma"):
            if name in f:
                h = f[name]
                break
        if h is None:
            return np.nan

        vals, x_edges, y_edges = h.to_numpy()
        x_cent = 0.5 * (x_edges[:-1] + x_edges[1:])
        y_cent = 0.5 * (y_edges[:-1] + y_edges[1:])

        ix = int(np.argmin(np.abs(x_cent - target_mchi_mev)))
        q_line = vals[ix, :]  # q(sigma) for selected mass
        if not np.any(np.isfinite(q_line)):
            return np.nan

        # Find first threshold crossing in ascending sigma bins.
        valid = np.isfinite(q_line)
        y = y_cent[valid]
        q = q_line[valid]
        for i in range(len(q) - 1):
            q0, q1 = q[i], q[i + 1]
            if q0 < Q_THRESHOLD <= q1:
                # Interpolate linearly in log10(sigma).
                ls0 = np.log10(y[i])
                ls1 = np.log10(y[i + 1])
                if q1 == q0:
                    return float(y[i])
                t = (Q_THRESHOLD - q0) / (q1 - q0)
                return float(10 ** (ls0 + t * (ls1 - ls0)))

        # If always above threshold, use lowest sigma bin.
        if np.nanmin(q) >= Q_THRESHOLD:
            return float(np.nanmin(y))
        # If never reaches threshold, no UL in sampled range.
        return np.nan


# ----------------------------------------------------------------------------
# main
#   Fill the (eps_h, E_gap) grid with sigma_UL for one mediator and mass (NaN for unphysical cells) and save it as .npy and .csv.
# ----------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mediator", choices=["heavy", "light"], required=True)
    ap.add_argument("--mchi-MeV", type=float, default=1.0)
    args = ap.parse_args()

    OUTDIR.mkdir(parents=True, exist_ok=True)

    arr = np.full((len(EHS), len(GAPS)), np.nan, dtype=float)
    rows: list[tuple[float, float, str, float]] = []

    for iy, eh in enumerate(EHS):
        for ix, gap in enumerate(GAPS):
            if not is_valid(gap, eh):
                rows.append((gap, eh, "unphysical", np.nan))
                continue
            rp = output_root_path(args.mediator, gap, eh)
            if not rp.is_file():
                rows.append((gap, eh, "missing_root", np.nan))
                continue
            ul = ul_from_qhist(rp, args.mchi_MeV)
            arr[iy, ix] = ul
            rows.append((gap, eh, "ok" if np.isfinite(ul) else "no_crossing", float(ul)))

    npy_path, csv_path = sensitivity_paths(args.mediator, args.mchi_MeV)
    np.save(npy_path, arr)
    print(f"[ok] wrote {npy_path} (m_chi = {args.mchi_MeV:g} MeV)")

    col = f"sigma_ul_cm2_at_mchi_{mchi_tag(args.mchi_MeV)}"
    with open(csv_path, "w", newline="", encoding="utf-8") as fp:
        w = csv.writer(fp)
        w.writerow(["band_gap_eV", "eh_pair_eV", "status", col])
        for r in rows:
            w.writerow(r)
    print(f"[ok] wrote {csv_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

