#!/usr/bin/env python3
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: plot_band_gap_pheno_p100K.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  plot_band_gap_pheno_p100K.py -- Compare P(n_e|E) for reference vs scaled
#  p100K tables (band-gap pheno Step 2)
# ============================================================================
"""Compare P(n_e|E) for reference vs scaled p100K tables (Step 2)."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Reuse loader/plot from builder
sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_p100K_scaled import load_p100k_csv, plot_p100k_compare  # noqa: E402


# ----------------------------------------------------------------------------
# main
#   Plot P(n_e = --ne | E) of the reference table (--ref) and of the tables given with --compare LABEL=CSV, and save the figure to --out.
# ----------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", default="data/p100K_table.csv")
    ap.add_argument("--compare", nargs="+", required=True, metavar="LABEL=CSV")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--ne", type=int, default=1)
    args = ap.parse_args()

    E_ref, P_ref = load_p100k_csv(Path(args.ref))
    tables = []
    for spec in args.compare:
        if "=" not in spec:
            ap.error(f"expected LABEL=path, got {spec!r}")
        label, path = spec.split("=", 1)
        E, P = load_p100k_csv(Path(path))
        if E.shape != E_ref.shape or not (E == E_ref).all():
            print(f"[warn] {path} energy grid differs from ref; plotting on its own E axis")
        tables.append((label, E, P))
    plot_p100k_compare(E_ref, P_ref, tables, args.out, ne_plot=args.ne)


if __name__ == "__main__":
    main()
