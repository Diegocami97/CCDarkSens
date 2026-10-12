#!/usr/bin/env python3
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  compare_tables.py -- Entry-by-entry comparison of a regenerated signal table
#  P(pattern | n_e) or background matrix B[p|q] with a reference table, using a
#  pooled two-sample binomial test. An independent cross-check of the --new
#  mode of ccdarksens_validate_pattern_efficiency / _background_efficiency.
# ============================================================================

"""
Compare a regenerated efficiency table with a reference table.

    python3 utils/compare_tables.py signal --ref data/efficiencies_paolo.csv \\
        --new outputs/pattern_signal_eff_chain/Efficiencies_patterns.csv \\
        --n-ref 100000 --n-new 1000000 --patterns 1 2 3 4 5 11 21 22 31 111 211

    python3 utils/compare_tables.py background --ref data/Background_efficiencies.csv \\
        --new outputs/pattern_background_eff/Background_efficiencies.csv \\
        --n-ref 100000 --n-new 1000000

For every entry the pull is (new - ref) / sigma with sigma from a pooled two-sample
binomial test, so entries that are zero or have one count in the reference need no
floor. A faithful reproduction has an RMS pull close to 1 over the non-zero entries and a
count of entries beyond 3 sigma close to 0.27% of them. Several --new files (independent
runs) are averaged. The per-entry table is written to --out (CSV).
"""

import argparse
import math
import sys

import pandas as pd


def pooled_pull(p_new, n_new, p_ref, n_ref):
    k = p_new * n_new + p_ref * n_ref
    ph = k / (n_new + n_ref)
    var = ph * (1.0 - ph) * (1.0 / n_new + 1.0 / n_ref)
    return 0.0 if var <= 0.0 else (p_new - p_ref) / math.sqrt(var)


def load_signal(path):
    d = pd.read_csv(path, comment="#")
    return d.groupby(["pattern", "ne"]).Efficiency.sum()


def load_background(path):
    d = pd.read_csv(path).set_index("iden_pat")
    s = d.stack()
    s.index = [(int(i), int(c.split("_")[-1])) for i, c in s.index]
    return s


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("kind", choices=["signal", "background"])
    ap.add_argument("--ref", required=True)
    ap.add_argument("--new", required=True, nargs="+")
    ap.add_argument("--n-ref", type=float, required=True, help="trials behind each reference entry")
    ap.add_argument("--n-new", type=float, required=True, help="trials per --new file")
    ap.add_argument("--patterns", type=int, nargs="*", help="signal only: restrict to these pattern codes")
    ap.add_argument("--out", default=None, help="write the per-entry comparison here (CSV)")
    args = ap.parse_args()

    load = load_signal if args.kind == "signal" else load_background
    ref = load(args.ref)
    frames = [load(f) for f in args.new]
    new = pd.concat(frames, axis=1).fillna(0.0).mean(axis=1)
    n_new = args.n_new * len(args.new)

    m = pd.concat([ref, new], axis=1, keys=["ref", "new"]).fillna(0.0)
    if args.kind == "signal" and args.patterns:
        m = m[[idx[0] in set(args.patterns) for idx in m.index]]
    m["diff"] = m.new - m.ref
    m["pull"] = [pooled_pull(a, n_new, b, args.n_ref) for a, b in zip(m.new, m.ref)]
    nz = m[(m.ref > 0) | (m.new > 0)]
    rms = math.sqrt((nz.pull ** 2).mean())
    n3 = int((nz.pull.abs() > 3).sum())
    print(f"{args.kind}: {len(m)} entries, {len(nz)} non-zero in either table")
    print(f"  RMS pull (non-zero entries): {rms:.2f}   chi2/n: {rms ** 2:.2f}")
    print(f"  entries beyond 3 sigma: {n3} (expected {0.0027 * len(nz):.1f})   beyond 2 sigma: "
          f"{int((nz.pull.abs() > 2).sum())} (expected {0.0455 * len(nz):.1f})   max |pull|: {nz.pull.abs().max():.2f}")
    if args.kind == "signal":
        by_ne = nz.groupby(level=1)
        print("  per n_e:  sum(ref)  sum(new)  RMS pull")
        for ne, g in by_ne:
            print(f"    {int(ne):2d}   {g.ref.sum():8.4f} {g.new.sum():8.4f}  {math.sqrt((g.pull ** 2).mean()):8.2f}")
    worst = nz.reindex(nz["pull"].abs().sort_values(ascending=False).index).head(8)
    print("  largest pulls:")
    for idx, r in worst.iterrows():
        print(f"    {idx}: ref {r.ref:.5g}  new {r.new:.5g}  pull {r.pull:+.2f}")
    if args.out:
        m.sort_index().to_csv(args.out, float_format="%.6g")
        print(f"  wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
