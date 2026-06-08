#!/usr/bin/env python3
# ============================================================================
#  CCDarkSens — plot_band_gap_2d_heatmap
#  Plot σ_UL gain heatmap on the 2D (E_gap, ε_h) grid
#
#  Author: Diego Venegas-Vargas
# ============================================================================
"""
Plot 2D sensitivity-gain heatmaps for heavy and light mediators at chosen m_chi.

Gain definition:
  gain = sigma_UL(Si_ref) / sigma_UL(cell)
where Si_ref is (gap,eh)=(1.2,3.8).

  python3 utils/plot_band_gap_2d_heatmap.py
  python3 utils/plot_band_gap_2d_heatmap.py --mchi-MeV 10 100
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "outplots" / "band_gap_pheno" / "2d_surface"

GAPS = np.array([0.1, 0.3, 0.5, 0.7, 0.9, 1.2], dtype=float)
EHS = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.8], dtype=float)


def is_valid(gap: float, eh: float) -> bool:
    return eh >= gap


def mchi_tag(mev: float) -> str:
    if abs(mev - round(mev)) < 1e-9:
        return f"{int(round(mev))}MeV"
    return f"{mev:g}MeV".replace(".", "p")


def sensitivity_npy_path(mediator: str, mchi_mev: float) -> Path:
    tagged = OUTDIR / f"sensitivity_2d_{mediator}_mchi{mchi_tag(mchi_mev)}.npy"
    if tagged.is_file():
        return tagged
    # Legacy 1 MeV naming (pre m_chi-tagged outputs).
    if abs(mchi_mev - 1.0) < 1e-9:
        legacy = OUTDIR / f"sensitivity_2d_{mediator}.npy"
        if legacy.is_file():
            return legacy
    return tagged


def edges(vals: np.ndarray) -> np.ndarray:
    e = np.zeros(vals.size + 1, dtype=float)
    e[1:-1] = 0.5 * (vals[:-1] + vals[1:])
    e[0] = vals[0] - 0.5 * (vals[1] - vals[0])
    e[-1] = vals[-1] + 0.5 * (vals[-1] - vals[-2])
    return e


def draw_heatmap(mediator: str, mchi_mev: float) -> None:
    npy_path = sensitivity_npy_path(mediator, mchi_mev)
    if not npy_path.is_file():
        raise FileNotFoundError(
            f"missing {npy_path}; run extract_band_gap_2d_sensitivity.py "
            f"--mediator {mediator} --mchi-MeV {mchi_mev:g}"
        )
    arr = np.load(npy_path)
    # arr indexing: [eh_idx, gap_idx]
    ref_ix = int(np.where(np.isclose(GAPS, 1.2))[0][0])
    ref_iy = int(np.where(np.isclose(EHS, 3.8))[0][0])
    sigma_ref = arr[ref_iy, ref_ix]
    gain = sigma_ref / arr

    # Mask invalid cells and non-finite values.
    valid_mask = np.array([[is_valid(g, e) for g in GAPS] for e in EHS], dtype=bool)
    gain_plot = np.where(valid_mask, gain, np.nan)
    gain_plot = np.where(np.isfinite(gain_plot) & (gain_plot > 0), gain_plot, np.nan)
    log_gain = np.log10(gain_plot)

    # Symmetric color window around gain=1 (log10=0).
    vmax = np.nanmax(np.abs(log_gain)) if np.any(np.isfinite(log_gain)) else 1.0
    vmax = max(vmax, 0.5)
    vmin = -vmax

    x_edges = edges(GAPS)
    y_edges = edges(EHS)

    fig, ax = plt.subplots(figsize=(8.0, 6.2))
    im = ax.pcolormesh(
        x_edges,
        y_edges,
        log_gain,
        cmap="coolwarm",
        shading="flat",
        vmin=vmin,
        vmax=vmax,
    )

    # Unphysical cells: gray + hatch.
    for iy, eh in enumerate(EHS):
        for ix, gap in enumerate(GAPS):
            if is_valid(gap, eh):
                continue
            rect = mpatches.Rectangle(
                (x_edges[ix], y_edges[iy]),
                x_edges[ix + 1] - x_edges[ix],
                y_edges[iy + 1] - y_edges[iy],
                facecolor="lightgray",
                edgecolor="gray",
                hatch="//",
                linewidth=0.6,
                alpha=0.9,
            )
            ax.add_patch(rect)

    # D-equal line (eh = gap).
    xx = np.linspace(GAPS.min(), GAPS.max(), 200)
    yy = xx
    ax.plot(xx, yy, "k--", lw=1.4, label=r"D-equal: $\varepsilon_h = E_{\mathrm{gap}}$")

    # Klein line: eh = 2.8*gap + 0.5
    yy_k = 2.8 * xx + 0.5
    in_view = (yy_k >= EHS.min()) & (yy_k <= EHS.max())
    ax.plot(xx[in_view], yy_k[in_view], "k-", lw=1.8, label=r"Klein: $\varepsilon_h = 2.8E_{\mathrm{gap}} + 0.5$")

    # Si reference marker.
    ax.plot([1.2], [3.8], marker="*", ms=14, color="gold", markeredgecolor="black", label="Si reference")

    ax.set_xlabel(r"$E_{\mathrm{gap}}$ [eV]")
    ax.set_ylabel(r"$\varepsilon_h$ [eV]")
    mchi_title = (
        str(int(round(mchi_mev))) if abs(mchi_mev - round(mchi_mev)) < 1e-9 else f"{mchi_mev:g}"
    )
    ax.set_title(
        f"{mediator.capitalize()} mediator: sensitivity gain at $m_\\chi = {mchi_title}$ MeV"
    )
    ax.set_xlim(x_edges[0], x_edges[-1])
    ax.set_ylim(y_edges[0], y_edges[-1])
    ax.grid(alpha=0.2)

    cbar = fig.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label(r"$\log_{10}\left(\sigma_{\mathrm{UL}}^{\mathrm{Si}} / \sigma_{\mathrm{UL}}\right)$")

    ax.legend(loc="upper left", fontsize=9, framealpha=0.95)

    fig.tight_layout()
    out = OUTDIR / f"heatmap_{mediator}_mchi{mchi_tag(mchi_mev)}.pdf"
    fig.savefig(out)
    plt.close(fig)
    print(f"[ok] wrote {out}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--mchi-MeV",
        type=float,
        nargs="+",
        default=[1.0, 10.0, 100.0],
        help="Reference m_chi values (MeV); default: 1 10 100",
    )
    ap.add_argument(
        "--mediator",
        choices=["heavy", "light", "both"],
        default="both",
    )
    args = ap.parse_args()

    OUTDIR.mkdir(parents=True, exist_ok=True)
    mediators = ["heavy", "light"] if args.mediator == "both" else [args.mediator]
    rc = 0
    for mchi in args.mchi_MeV:
        for med in mediators:
            try:
                draw_heatmap(med, mchi)
            except FileNotFoundError as exc:
                print(f"WARN: {exc}")
                rc = 1
    return rc


if __name__ == "__main__":
    raise SystemExit(main())

