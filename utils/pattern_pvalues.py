#!/usr/bin/env python3
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  pattern_pvalues.py -- Hand check of the pattern-identification variables
#  p_mn, p_111 and p_211 of the DAMIC-M data-selection note (Eqs. 9-13) for one
#  candidate, printed with the promotion chain step by step. An independent
#  cross-check of ccdarksens_classify_pair, which runs the C++ classifier.
# ============================================================================

"""
For two consecutive pixels (q_j, q_j+1) of one CCD, the note defines

    p_mn = -log[ cdf(q_j ; m) * cdf(q_j+1 ; n) ]        cdf(q ; m) = Phi((q - m) / sigma)

where Phi is the standard normal CDF and sigma is the readout noise of that
CCD. The two pixels are a candidate for 1+1 if p_11 < 4, for 2+1 if
min(p_21, p_12) < 4, and for 3+1 if min(p_31, p_13) < 4. A third pixel next
to a 1+1 (or 2+1) candidate turns it into 1+1+1 (or 2+1+1) when
p_111 = p_11 - max(log cdf(q_left ; 1), log cdf(q_right ; 1)) is below 5.5
(p_211 likewise starting from p_21).

The noise can be given directly (--sigma) or through the one-electron bound
that the LBC plots print: bound = 3.75 sigma, so sigma = bound / 3.75.

Usage:
    python3 utils/pattern_pvalues.py 3.17 1.85 --bound 0.64
    python3 utils/pattern_pvalues.py 1.05 0.80 --sigma 0.16 --left 0.7
"""

import argparse
import math


def cdf(q, m, sigma):
    """Lower-tail Gaussian CDF of a pixel read as q when m electrons are present."""
    return 0.5 * (1.0 + math.erf((q - m) / (sigma * math.sqrt(2.0))))


def p_mn(q1, q2, m, n, sigma):
    """-log of the product of the two lower-tail CDFs (Eq. 9), capped for printing."""
    prob = cdf(q1, m, sigma) * cdf(q2, n, sigma)
    return -math.log(max(prob, 1e-300))


def p_pair(q1, q2, m, n, sigma):
    """min(p_mn, p_nm): the pair may be ordered either way along the row (Eqs. 10, 11)."""
    return min(p_mn(q1, q2, m, n, sigma), p_mn(q1, q2, n, m, sigma))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("q1", type=float, help="charge of the first pixel (e-)")
    ap.add_argument("q2", type=float, help="charge of the second pixel (e-)")
    grp = ap.add_mutually_exclusive_group(required=True)
    grp.add_argument("--sigma", type=float, help="readout noise of the CCD (e-)")
    grp.add_argument("--bound", type=float, help="one-electron bound of the CCD (= 3.75 sigma)")
    ap.add_argument("--left", type=float, default=None, help="charge of the pixel to the left (e-)")
    ap.add_argument("--right", type=float, default=None, help="charge of the pixel to the right (e-)")
    ap.add_argument("--thr-pair", type=float, default=4.0, help="threshold for pairs (default 4)")
    ap.add_argument("--thr-triple", type=float, default=5.5, help="threshold for triples (default 5.5)")
    ap.add_argument("--qmax", type=float, default=5.6, help="maximum total charge (default 5.6 e-)")
    a = ap.parse_args()

    sigma = a.sigma if a.sigma is not None else a.bound / 3.75
    print(f"sigma = {sigma:.4f} e-   pair = ({a.q1}, {a.q2})   total charge = {a.q1 + a.q2:.2f} e-")
    print()
    print("  hypothesis (m,n)   cdf(q1;m)   cdf(q2;n)     p_mn")
    for m, n in [(1, 1), (2, 1), (1, 2), (3, 1), (1, 3), (2, 2), (3, 2), (2, 3)]:
        c1, c2 = cdf(a.q1, m, sigma), cdf(a.q2, n, sigma)
        print(f"        ({m},{n})        {c1:9.4f}   {c2:9.4f}   {p_mn(a.q1, a.q2, m, n, sigma):8.2f}")

    p11 = p_mn(a.q1, a.q2, 1, 1, sigma)
    p21 = p_pair(a.q1, a.q2, 2, 1, sigma)
    p31 = p_pair(a.q1, a.q2, 3, 1, sigma)
    p22 = p_mn(a.q1, a.q2, 2, 2, sigma)
    print()
    print(f"  p_11 = {p11:.2f}   p_21 = min(p_21,p_12) = {p21:.2f}   p_31 = min(p_31,p_13) = {p31:.2f}   p_22 = {p22:.2f}")

    label = None
    if p11 < a.thr_pair:
        label = "11"
        print(f"  p_11 < {a.thr_pair}: candidate 1+1")
        if p21 < a.thr_pair:
            label = "21"
            print(f"  p_21 < {a.thr_pair}: promoted to 2+1")
            if p31 < a.thr_pair:
                label = "31"
                print(f"  p_31 < {a.thr_pair}: promoted to 3+1")
            elif p22 < a.thr_pair:
                label = "22"
                print(f"  p_31 fails, p_22 < {a.thr_pair}: labelled 2+2")
    else:
        print(f"  p_11 >= {a.thr_pair}: not a candidate")

    # third-pixel extension of 1+1 and 2+1 (Eqs. 12, 13)
    if label in ("11", "21") and (a.left is not None or a.right is not None):
        logs = [math.log(max(cdf(q, 1, sigma), 1e-300)) for q in (a.left, a.right) if q is not None]
        base = p11 if label == "11" else p21
        ptri = base - max(logs)
        name = "p_111" if label == "11" else "p_211"
        print(f"  {name} = {base:.2f} - max(log cdf) = {ptri:.2f}")
        if ptri < a.thr_triple:
            label = "111" if label == "11" else "211"
            print(f"  {name} < {a.thr_triple}: extended to {label}")

    print()
    if label is None:
        print("  result: no pattern")
    elif a.q1 + a.q2 >= a.qmax:
        print(f"  result: {{{label}}} but total charge {a.q1 + a.q2:.2f} >= {a.qmax}: excluded")
    else:
        print(f"  result: {{{label}}}   (total charge {a.q1 + a.q2:.2f} < {a.qmax})")


if __name__ == "__main__":
    main()
