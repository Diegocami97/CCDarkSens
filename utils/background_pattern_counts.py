#!/usr/bin/env python3
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  background_pattern_counts.py -- From per-pixel k-electron rates to the
#  random-coincidence background counts B^rc_p in the six pattern bins, using
#  the pattern-migration matrix B[p|q] from data/Background_efficiencies.csv.
#  An independent cross-check of ccdarksens_compute_brc --config.
# ============================================================================

"""
Compute the random-coincidence pattern background

    B^rc_p = sum_s N_s * sum_q [ prod_i Poisson(q_i | lambda_{k=q_i, s}) ] * B[p|q]

where q runs over the 26 injected ideal patterns in the columns of the
migration matrix (the empty pattern 0 has weight 1), q_i is the electron
count of the i-th pixel of q, and lambda_k is the mean number of k-electron
events per pixel readout (k = 1, 2, and 3 or more). The outer sum runs over
"segments" s (an image, a CCD, a group of columns) of N_s pixel readouts that
each have their own rates, as in Eq. (4) of the DAMIC-M SRDM paper. With one
segment this is the uniform-rate formula used by pydme's Background_pattern_.

Usage:
    python3 utils/background_pattern_counts.py configs/examples/background_pattern_rates.json
    python3 utils/background_pattern_counts.py CONFIG --lambda1 5.8e-4 --lambda2 4e-4

The command-line overrides replace the corresponding "rates" entries of the
JSON. If "pixels.segments_csv" is set, the rates come from that table instead
(columns n_pixels, lambda_1 and optionally lambda_2, lambda_3plus; missing
columns fall back to the "rates" block). The result is printed and written as
a ready-to-paste "background_Bp" list for a scan config.
"""

import argparse
import csv
import json
import math
import os
import sys


def poisson_pmf(k, lam):
    """P(k | lam). A zero mean puts all weight on k = 0."""
    if lam <= 0.0:
        return 1.0 if k == 0 else 0.0
    return math.exp(-lam + k * math.log(lam) - math.lgamma(k + 1))


def load_migration_matrix(path):
    """Return {identified_code: {injected_code: B[p|q]}} from the reference-layout CSV."""
    with open(path, newline="") as f:
        rows = list(csv.reader(f))
    injected = [int(h.split("_")[-1]) for h in rows[0][1:]]
    return {int(r[0]): {q: float(v) for q, v in zip(injected, r[1:])} for r in rows[1:]}


def load_n_pixels(path):
    """Sum of the Nusedpix column of the LBC image table."""
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        return sum(int(float(r["Nusedpix"])) for r in reader)


def load_segments(path, default):
    """Read per-segment pixel counts and rates; absent rate columns use `default`."""
    segs = []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            lam = dict(default)
            for k, col in ((1, "lambda_1"), (2, "lambda_2"), (3, "lambda_3plus")):
                if col in r and r[col] != "":
                    lam[k] = float(r[col])
            segs.append((float(r["n_pixels"]), lam))
    return segs


def pattern_weight(q_code, lam):
    """Poisson probability of the ideal pattern q; lam maps k = 1, 2, 3 to a mean."""
    if q_code == 0:
        return 1.0
    weight = 1.0
    for ch in str(q_code):
        k = int(ch)
        weight *= poisson_pmf(k, lam[min(k, 3)])
    return weight


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("config")
    ap.add_argument("--lambda1", type=float, help="override rates.lambda_1")
    ap.add_argument("--lambda2", type=float, help="override rates.lambda_2")
    ap.add_argument("--lambda3", type=float, help="override rates.lambda_3plus")
    ap.add_argument("--n-pixels", type=float, help="override pixels.n_pixel_readouts")
    args = ap.parse_args()

    with open(args.config) as f:
        cfg = json.load(f)

    rates = cfg["rates"]
    lam = {
        1: args.lambda1 if args.lambda1 is not None else rates["lambda_1"],
        2: args.lambda2 if args.lambda2 is not None else rates["lambda_2"],
        3: args.lambda3 if args.lambda3 is not None else rates["lambda_3plus"],
    }
    seg_path = cfg["pixels"].get("segments_csv")
    if seg_path:
        segments = load_segments(seg_path, lam)
        n_pix = sum(n for n, _ in segments)
    else:
        n_pix = args.n_pixels if args.n_pixels is not None else cfg["pixels"].get("n_pixel_readouts")
        if n_pix is None:
            n_pix = load_n_pixels(cfg["pixels"]["image_table_csv"])
        n_pix = float(n_pix)
        segments = [(n_pix, lam)]

    B = load_migration_matrix(cfg["migration_matrix_csv"])
    patterns = cfg["patterns"]
    breakdown = cfg.get("breakdown_pattern", patterns[0])

    print("=" * 66)
    print("  Random-coincidence pattern background  B^rc_p")
    print("=" * 66)
    if seg_path:
        print(f"  rates per segment from {seg_path} ({len(segments)} segments)")
        for k, name in ((1, "lambda_1"), (2, "lambda_2"), (3, "lambda_3+")):
            mean = sum(n * l[k] for n, l in segments) / n_pix
            print(f"  pixel-weighted mean {name:<9} = {mean:.4e} e-/pixel")
    else:
        print(f"  lambda_1  = {lam[1]:.4e} e-/pixel   (theta_1 = {lam[1] * 1e4:.4g})")
        print(f"  lambda_2  = {lam[2]:.4e} e-/pixel   (theta_2 = {lam[2] * 1e4:.4g})")
        print(f"  lambda_3+ = {lam[3]:.4e} e-/pixel   (theta_3 = {lam[3] * 1e4:.4g})")
    print(f"  N_pixel readouts = {n_pix:,.0f}")
    print(f"  migration matrix = {cfg['migration_matrix_csv']}")
    print("-" * 66)
    print(f"  {'pattern':<9}{'B^rc_p':>13}{'diagonal (q=p)':>17}{'misidentified':>17}")

    results = []
    for p in patterns:
        total = diag = 0.0
        for q, bpq in B[p].items():
            term = sum(n * pattern_weight(q, l) for n, l in segments) * bpq
            total += term
            if q == p:
                diag = term
        results.append(total)
        print(f"  {{{p}}}".ljust(11) + f"{total:>13.4e}{diag:>17.4e}{total - diag:>17.4e}")

    print("-" * 66)
    print(f"  Which injected patterns feed {{{breakdown}}} (contribution > 1e-3 counts):")
    print(f"  {'q':<8}{'P(q)':>13}{'B[p|q]':>13}{'counts':>13}")
    terms = []
    for q, bpq in B[breakdown].items():
        w = sum(n * pattern_weight(q, l) for n, l in segments) / n_pix  # pixel-averaged P(q)
        terms.append((n_pix * w * bpq, q, w, bpq))
    for counts, q, w, bpq in sorted(terms, reverse=True):
        if counts > 1e-3:
            print(f"  {q:<8}{w:>13.3e}{bpq:>13.3e}{counts:>13.4f}")

    snippet = {"background_Bp": [float(f"{r:.6g}") for r in results]}
    print("-" * 66)
    print("  Paste into a scan config under \"run\" (pattern order " + str(patterns) + "):")
    print("  " + json.dumps(snippet))
    print("=" * 66)

    out = cfg.get("output_json")
    if out:
        os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
        with open(out, "w") as f:
            json.dump({"patterns": patterns,
                       "rates_source": seg_path or {"1": lam[1], "2": lam[2], "3plus": lam[3]},
                       "n_pixel_readouts": n_pix, **snippet}, f, indent=2)
        print(f"  Wrote {out}")


if __name__ == "__main__":
    sys.exit(main())
