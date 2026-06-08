#!/usr/bin/env python3
# ============================================================================
#  CCDarkSens — plot_band_gap_pheno_ratio
#  Plot σ_UL ratio vs reference from limit CSV (band-gap pheno Step 5)
#
#  Author: Diego Venegas-Vargas
# ============================================================================
"""Plot sigma_UL ratio vs reference from limit CSV (Step 5)."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def load_limit_csv(path: Path) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Parse ccdarksens_plot_dmelectron_limit --out-csv format."""
    curves: dict[str, list[tuple[float, float]]] = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.lower().startswith("label,"):
                continue
            parts = [p.strip() for p in line.split(",")]
            if len(parts) < 3:
                continue
            label, mchi, sigma = parts[0], float(parts[1]), float(parts[2])
            curves.setdefault(label, []).append((mchi, sigma))
    out: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for lab, pts in curves.items():
        arr = np.asarray(pts)
        order = np.argsort(arr[:, 0])
        out[lab] = (arr[order, 0], arr[order, 1])
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", type=Path, required=True)
    ap.add_argument("--reference", required=True, help="Label of reference curve in CSV")
    ap.add_argument("--compare", nargs="*", help="Labels to plot (default: all except reference)")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    curves = load_limit_csv(args.csv)
    if args.reference not in curves:
        raise SystemExit(f"reference {args.reference!r} not in CSV; have {list(curves)}")
    m_ref, s_ref = curves[args.reference]
    labels = args.compare or [k for k in curves if k != args.reference]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 5))
    for lab in labels:
        if lab not in curves:
            print(f"[warn] skip unknown label {lab!r}")
            continue
        m, s = curves[lab]
        # interpolate ref onto m grid of this curve
        s_ref_i = np.interp(m, m_ref, s_ref, left=np.nan, right=np.nan)
        ratio = s / s_ref_i
        ax.plot(m, ratio, lw=1.5, label=lab)
    ax.axhline(1.0, color="k", ls="--", lw=0.8)
    ax.set_xscale("log")
    ax.set_xlabel(r"$m_\chi$ [MeV]")
    ax.set_ylabel(r"$\sigma_{\mathrm{UL}} / \sigma_{\mathrm{UL}}^{\mathrm{ref}}$")
    ax.legend(loc="best", fontsize=9)
    ax.grid(True, alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(args.out)
    plt.close(fig)
    print(f"[ok] wrote {args.out}")


if __name__ == "__main__":
    main()
