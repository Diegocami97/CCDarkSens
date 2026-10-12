#!/usr/bin/env python3
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  compare_limit_curves.py -- Ratio of the upper-limit curves of several
#  scan_generic.root files to the first one, per mass window, using the
#  continuous bisection limit, with a count of the grid-crossing points that sit
#  on the floor of the sigma grid. An independent cross-check of
#  ccdarksens_plot_limit --draw-both --ratio-to.
# ============================================================================

"""
Why the bisection graph and not the stored grid limit: when the observed counts
are above the expectation, the profile-likelihood statistic at sigma -> 0 can
already exceed the threshold (q_0 > q_threshold). The grid-crossing limit then
has no crossing and is reported at the lowest sigma of the grid (1e-46 cm^2 in
the shipped configs): that value means "the data prefer a signal at this mass",
it is not an upper limit. The bisection limit starts at the best-fit sigma and
steps upward, so it is a valid one-sided limit in that case too.

Usage:
    python3 utils/compare_limit_curves.py BASE.root A.root [B.root ...] \
        [--labels base A B] [--plot out.png]

All files must come from the same mass grid. Needs uproot and numpy (and
matplotlib for --plot).
"""

import argparse

import numpy as np
import uproot

BISECTION = "upper_limit_sigma_e_mchi_pydme_bisection"
GRID = "upper_limit_sigma_e_mchi_graph"
WINDOWS = [(0.7, 3.0), (3.0, 10.0), (10.0, 30.0), (30.0, 100.0), (100.0, 1000.0)]


def read_graph(path, name):
    """Return (mass, sigma) arrays of a TGraph stored in the scan output."""
    g = uproot.open(path)[name]
    return np.array(g.member("fX")), np.array(g.member("fY"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+", help="scan_generic.root files; the first one is the reference")
    ap.add_argument("--labels", nargs="*", default=None, help="one label per file")
    ap.add_argument("--sigma-floor", type=float, default=1e-46,
                    help="lowest sigma of the scan grid (cm^2), the value reported when there is no crossing")
    ap.add_argument("--plot", default=None, help="write a PNG of the ratio curves")
    a = ap.parse_args()
    labels = a.labels if a.labels else [f"run{i}" for i in range(len(a.files))]
    if len(labels) != len(a.files):
        ap.error("--labels must have one entry per file")

    bis = [read_graph(f, BISECTION) for f in a.files]
    grid = [read_graph(f, GRID) for f in a.files]
    x0, y0 = bis[0]
    floor = a.sigma_floor * (1.0 + 1e-6)

    print(f"reference: {labels[0]}  ({len(x0)} bisection points, "
          f"{x0.min():.2f}-{x0.max():.0f} MeV)")
    print()
    print("grid-crossing points sitting on the sigma floor (q_0 above threshold, not a limit):")
    for lab, (gx, gy) in zip(labels, grid):
        print(f"  {lab:>12}: {int((gy <= floor).sum())} of {len(gy)}")
    print()
    head = f"{'window (MeV)':>16} {'n':>5}  " + "  ".join(f"{l:>20}" for l in labels[1:])
    print("ratio of the bisection limit to the reference, median (16%-84%) over the window:")
    print(head)
    for lo, hi in WINDOWS:
        s = (x0 >= lo) & (x0 < hi)
        cells = []
        for (x, y) in bis[1:]:
            if len(x) != len(x0) or not np.allclose(x, x0):
                raise SystemExit("mass grids differ between files")
            r = y[s] / y0[s]
            cells.append(f"{np.median(r):6.2f} ({np.percentile(r, 16):4.2f}-{np.percentile(r, 84):4.2f})")
        print(f"{lo:7.1f}-{hi:7.1f} {int(s.sum()):5d}  " + "  ".join(f"{c:>20}" for c in cells))

    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 7), sharex=True,
                                       gridspec_kw={"height_ratios": [2, 1]})
        for i, (lab, (x, y)) in enumerate(zip(labels, bis)):
            ax1.loglog(x, y, color=f"C{i}", label=lab)
        ax1.set_ylabel(r"$\sigma_e$ upper limit (cm$^2$)")
        ax1.legend()
        for i, (lab, (x, y)) in enumerate(zip(labels, bis)):
            if i > 0:
                ax2.semilogx(x, y / y0, color=f"C{i}", label=lab)
        ax2.axhline(1.0, color="k", lw=0.8)
        ax2.set_xlabel(r"$m_\chi$ (MeV)")
        ax2.set_ylabel("ratio to reference")
        fig.tight_layout()
        fig.savefig(a.plot, dpi=150)
        print(f"\nwrote {a.plot}")


if __name__ == "__main__":
    main()
