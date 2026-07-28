#!/usr/bin/env python3
# ============================================================================
#  CCDarkSens — plot_cluster_visualization
#  Figure 3: 2D pixel cluster images — Si ref vs Klein ladder, vs depth
#
#  Author: Diego Venegas-Vargas
# ============================================================================
"""
Visualize pixel clusters using the CCDarkSens diffusion model (analytic).

  python3 utils/plot_cluster_visualization.py
  python3 utils/plot_cluster_visualization.py --Er 2,4,6
  python3 utils/plot_cluster_visualization.py --klein-gaps 0.1,0.5 --Er 4

Writes one PDF/PNG per E_r: cluster_visualization_Er{tag}.pdf
"""
from __future__ import annotations

import argparse
import sys
from math import erf, sqrt
from pathlib import Path

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import BoundaryNorm  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "utils"))
from band_gap_klein_plot import (  # noqa: E402
    COLOR_SI_REF,
    DEFAULT_KLEIN_GAPS,
    SI_EGAP_EV,
    SI_EH_EV,
    er_file_tag,
    filter_klein_gaps,
    gap_color,
    klein_eh,
    klein_mathtext_label,
    parse_float_list,
)

OUTDIR = ROOT / "outplots" / "band_gap_pheno" / "ne_imaging"
DEFAULT_ER_VALUES = [4.0, 2.0]

A_UM2 = 803.25
B_UMINV = 0.00065
ALPHA = 1.0
BETA_PER_KEV = 0.0
SIGMA_READOUT_E = 0.16
PIXEL_SIZE_UM = 15.0

GRID_PIX = 9
DEPTHS_UM = [50.0, 337.0, 620.0]
DEPTH_LABELS = ["shallow", "mid", "deep"]
RNG_SEED = 20260609


def build_materials(klein_gaps: list[float]) -> list[dict]:
    si = dict(
        Egap=SI_EGAP_EV,
        eh=SI_EH_EV,
        color=COLOR_SI_REF,
        is_si=True,
        mat=r"Si ref $(1.2,\ 3.8)$",
    )
    klein = []
    for i, g in enumerate(klein_gaps):
        eh = klein_eh(g)
        klein.append(
            dict(
                Egap=g,
                eh=eh,
                color=gap_color(i),
                is_si=False,
                mat=klein_mathtext_label(g, eh),
            )
        )
    return [si] + klein


def _scenario_ne(Er: float, Egap: float, eh: float) -> float:
    return max(Er - Egap, 0.0) / eh


def compute_sigma_xy_um(z_um: float, E_eV: float) -> float:
    E_keV = max(0.0, E_eV) * 1e-3
    inside = 1.0 - B_UMINV * z_um
    if inside <= 0.0:
        return 0.0
    return sqrt(-A_UM2 * np.log(inside)) * (ALPHA + BETA_PER_KEV * E_keV)


def _gauss_pixel_fractions(sigma_um: float) -> np.ndarray:
    half = GRID_PIX // 2
    edges_um = (np.arange(GRID_PIX + 1) - (half + 0.5)) * PIXEL_SIZE_UM
    if sigma_um <= 0.0:
        frac_1d = np.zeros(GRID_PIX)
        frac_1d[half] = 1.0
    else:
        cdf = 0.5 * (1.0 + np.array([erf(e / (sqrt(2.0) * sigma_um)) for e in edges_um]))
        frac_1d = np.diff(cdf)
    return np.outer(frac_1d, frac_1d)


def make_cluster_image(n_e: float, sigma_um: float, rng: np.random.Generator) -> np.ndarray:
    frac = _gauss_pixel_fractions(sigma_um)
    expected = n_e * frac
    noisy = expected + rng.normal(0.0, SIGMA_READOUT_E, size=expected.shape)
    return np.rint(noisy).astype(int)


def _render_figure(
    Er: float, out_path: Path, rng: np.random.Generator, materials: list[dict]
) -> int:
    nrows, ncols = len(materials), len(DEPTHS_UM)

    images: list[list[np.ndarray]] = []
    vmax = 1
    for mat in materials:
        n_e = _scenario_ne(Er, mat["Egap"], mat["eh"])
        row_imgs = []
        for z in DEPTHS_UM:
            sigma = compute_sigma_xy_um(z, Er)
            img = make_cluster_image(n_e, sigma, rng)
            vmax = max(vmax, int(img.max()))
            row_imgs.append(img)
        images.append(row_imgs)

    cmap = plt.get_cmap("Blues")
    bounds = np.arange(-0.5, vmax + 1.5, 1.0)
    norm = BoundaryNorm(bounds, cmap.N)

    fig, axes = plt.subplots(
        nrows, ncols, figsize=(4.0 * ncols, 3.4 * nrows), squeeze=False
    )
    left, right, top, bottom = 0.05, 0.86, 0.93, 0.03
    fig.subplots_adjust(left=left, right=right, top=top, bottom=bottom, hspace=0.55, wspace=0.08)

    im = None
    for r, mat in enumerate(materials):
        n_e = _scenario_ne(Er, mat["Egap"], mat["eh"])
        title_color = mat["color"]
        for c, z in enumerate(DEPTHS_UM):
            ax = axes[r][c]
            sigma = compute_sigma_xy_um(z, Er)
            img = images[r][c]
            im = ax.imshow(img, cmap=cmap, norm=norm, origin="lower")

            for spine in ax.spines.values():
                spine.set_edgecolor(title_color)
                spine.set_linewidth(2.0 if mat.get("is_si") else 1.2)

            ax.set_xticks(np.arange(-0.5, GRID_PIX, 1), minor=True)
            ax.set_yticks(np.arange(-0.5, GRID_PIX, 1), minor=True)
            ax.grid(which="minor", color="0.6", linewidth=0.5)
            ax.tick_params(which="both", length=0)
            ax.set_xticks([])
            ax.set_yticks([])

            if c == 0:
                ax.set_ylabel(mat["mat"], fontsize=8.5, color=title_color, labelpad=8)

            ax.set_title(
                rf"{DEPTH_LABELS[c]} — $z={z:.0f}\ \mu$m" + "\n"
                rf"$\langle n_e\rangle={n_e:.2f}$,  $\sigma_{{xy}}={sigma:.1f}\ \mu$m",
                fontsize=8,
                color=title_color,
            )

        row_top = axes[r][0].get_position().y1
        row_bot = axes[r][0].get_position().y0
        cax = fig.add_axes([right + 0.02, row_bot, 0.012, row_top - row_bot])
        cbar = fig.colorbar(im, cax=cax, ticks=np.arange(0, min(vmax, 12) + 1))
        cbar.set_label(r"charge [$e^-$]", fontsize=7)

    fig.suptitle(
        rf"Pixel clusters at $E_r={Er:g}$ eV — Si ref vs Klein "
        rf"($\varepsilon_h=2.8\,E_{{\mathrm{{gap}}}}+0.5$), "
        rf"pitch {PIXEL_SIZE_UM:.0f} $\mu$m, $\sigma_{{ro}}={SIGMA_READOUT_E}\,e^-$",
        fontsize=11,
        y=0.995,
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight", dpi=150)
    fig.savefig(out_path.with_suffix(".png"), bbox_inches="tight", dpi=150)
    plt.close(fig)

    print(f"\n=== E_r = {Er} eV ===")
    for mat in materials:
        n_e = _scenario_ne(Er, mat["Egap"], mat["eh"])
        print(f"  {mat['mat']:40s}  <n_e>={n_e:.2f}")
    print(f"Wrote {out_path}")
    print(f"Wrote {out_path.with_suffix('.png')}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--out-dir", type=Path, default=OUTDIR)
    ap.add_argument(
        "--klein-gaps",
        type=str,
        default=",".join(str(g) for g in DEFAULT_KLEIN_GAPS),
        metavar="CSV",
        help="Comma-separated E_gap [eV] for Klein rows (default: 0.1,0.3,0.5,0.7,0.9)",
    )
    ap.add_argument(
        "--Er",
        type=str,
        default=",".join(str(e) for e in DEFAULT_ER_VALUES),
        metavar="CSV",
        help="Comma-separated recoil energies [eV]; one figure per value (default: 4,2)",
    )
    args = ap.parse_args()

    klein_gaps = filter_klein_gaps(parse_float_list(args.klein_gaps))
    er_values = parse_float_list(args.Er)
    if not klein_gaps:
        print("ERROR: no Klein gaps after filtering", file=sys.stderr)
        return 1
    if not er_values:
        print("ERROR: --Er must list at least one energy", file=sys.stderr)
        return 1

    materials = build_materials(klein_gaps)
    rng = np.random.default_rng(RNG_SEED)
    rc = 0
    for Er in er_values:
        out = args.out_dir / f"cluster_visualization_{er_file_tag(Er)}.pdf"
        rc |= _render_figure(Er, out, rng, materials)
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
