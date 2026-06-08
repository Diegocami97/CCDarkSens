#!/usr/bin/env python3
# ============================================================================
#  CCDarkSens — plot_band_gap_one_point_p100K
#  Diagnostic spectra and limit plot for a single (E_gap, ε_h) point
#
#  Author: Diego Venegas-Vargas
# ============================================================================
"""
Compare P(n_e | E) ionization tables (scaled p100K CSVs) for band-gap one-point cases.

Reads data/p100K_gap*_eh*.csv directly — not folded S_true(n_e) from scans.

Author: Diego Venegas-Vargas
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import List, Tuple

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from band_gap_plot_labels import pheno_param_label_mpl  # noqa: E402
from build_p100K_scaled import load_p100k_csv  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "configs" / "band_gap_one_point_spectra_manifest.json"
DEFAULT_OUT = ROOT / "outplots/band_gap_one_point_spectra/p100K_Pne"


def gap_tag(g: float) -> str:
    return "gap1p2" if abs(g - 1.2) < 1e-9 else f"gap{g:.1f}".replace(".", "p")


def eh_tag(eh: float) -> str:
    return f"{eh:g}".replace(".", "p")


def ionization_csv_path(gap_eV: float, eh_eV: float) -> Path:
    return ROOT / "data" / f"p100K_{gap_tag(gap_eV)}_eh{eh_tag(eh_eV)}.csv"


def build_cases(manifest: dict) -> list[dict]:
    gaps = [float(g) for g in manifest.get("gaps_eV", [])]
    scenarios = manifest.get("scenarios", ["D-equal", "B-thresh"])
    eh_thresh = float(manifest.get("eh_pair_B_thresh_eV", 3.8))
    cases: list[dict] = []
    for g in gaps:
        tag = gap_tag(g)
        for scen in scenarios:
            if scen == "D-equal":
                eh = g
            elif scen == "B-thresh":
                eh = eh_thresh
            else:
                continue
            cases.append(
                {
                    "id": f"{tag}_{scen}",
                    "band_gap_eV": g,
                    "eh_pair_eV": eh,
                    "scenario": scen,
                    "ionization_csv": str(ionization_csv_path(g, eh).relative_to(ROOT)),
                }
            )
    return cases


def ion_legend(c: dict) -> str:
    return pheno_param_label_mpl(c["band_gap_eV"], c["eh_pair_eV"])


def short_legend(c: dict) -> str:
    return pheno_param_label_mpl(c["band_gap_eV"], c["eh_pair_eV"])


def load_table(c: dict) -> Tuple[np.ndarray, np.ndarray]:
    path = ROOT / c["ionization_csv"]
    if not path.is_file():
        raise FileNotFoundError(path)
    return load_p100k_csv(path)


def tab20(n: int) -> List[tuple]:
    cmap = plt.get_cmap("tab20")
    return [cmap(i % 20) for i in range(n)]


def plot_pne_vs_E_subplots(
    series: List[Tuple[str, np.ndarray, np.ndarray, dict]],
    out: Path,
    suptitle: str,
    Emax: float,
    ne_max: int = 5,
    ylog: bool = False,
) -> None:
    """One subplot per n_e: P(n_e | E) vs E for each table in series."""
    out.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(ne_max, 1, figsize=(10, 2.2 * ne_max), sharex=True, squeeze=False)

    for j in range(ne_max):
        ax = axes[j, 0]
        for i, (label, E, P, _meta) in enumerate(series):
            m = (E >= 0) & (E <= Emax)
            y = np.clip(P[j, m], 1e-12 if ylog else 0.0, 1.0)
            ax.plot(E[m], y, lw=1.8, label=label)
        ax.set_ylabel(rf"$P(n_e={j + 1} \mid E)$")
        if ylog:
            ax.set_yscale("log")
            ax.set_ylim(1e-6, 1.05)
        else:
            ax.set_ylim(0, 1.05)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, loc="upper right")

    axes[-1, 0].set_xlabel(r"Recoil energy $E$ [eV]")
    axes[0, 0].set_xlim(0, Emax)
    fig.suptitle(suptitle, fontsize=11)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"[ok] {out}")


def plot_pne_vs_ne_at_E(
    series: List[Tuple[str, np.ndarray, np.ndarray]],
    energies: List[float],
    out: Path,
    suptitle: str,
    ne_max: int = 5,
) -> None:
    """Bar/line: P(n_e | E_r) vs n_e at fixed recoil energies (table-native view)."""
    out.parent.mkdir(parents=True, exist_ok=True)
    nE = len(energies)
    fig, axes = plt.subplots(1, nE, figsize=(3.4 * nE, 4.2), squeeze=False)
    x = np.arange(1, ne_max + 1)

    for ax, e0 in zip(axes[0], energies):
        for label, E, P in series:
            row = np.array([np.interp(e0, E, P[k], left=0.0, right=0.0) for k in range(P.shape[0])])
            ax.plot(x, row[:ne_max], "o-", lw=1.6, ms=5, label=label)
        ax.set_xticks(x)
        ax.set_xlabel(r"$n_e$")
        ax.set_ylabel(r"$P(n_e \mid E)$")
        ax.set_title(rf"$E={e0:g}$ eV")
        ax.set_ylim(0, 1.05)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7)

    fig.suptitle(suptitle, fontsize=11)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"[ok] {out}")


def plot_heatmap(E: np.ndarray, P: np.ndarray, out: Path, title: str, Emax: float, ne_max: int = 5) -> None:
    """P(n_e, E) for n_e = 1..ne_max (direct table view)."""
    out.parent.mkdir(parents=True, exist_ok=True)
    m = (E >= 0) & (E <= Emax)
    Ee = E[m]
    Z = P[:ne_max, m]

    fig, ax = plt.subplots(figsize=(9, 4))
    im = ax.imshow(
        Z,
        aspect="auto",
        origin="lower",
        extent=[Ee[0], Ee[-1], 0.5, ne_max + 0.5],
        vmin=0,
        vmax=1,
        cmap="viridis",
    )
    ax.set_xlabel(r"Recoil energy $E$ [eV]")
    ax.set_ylabel(r"$n_e$")
    ax.set_yticks(range(1, ne_max + 1))
    ax.set_title(title)
    fig.colorbar(im, ax=ax, label=r"$P(n_e \mid E)$")
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"[ok] {out}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    ap.add_argument(
        "--Emax",
        type=float,
        default=50.0,
        help="Energy axis max [eV] (Si reference tables extend to 50 eV)",
    )
    ap.add_argument("--ne-max", type=int, default=5, help="Show n_e = 1 .. ne-max")
    ap.add_argument(
        "--E-slices",
        type=float,
        nargs="+",
        default=[0.15, 0.3, 0.5, 1.0, 2.0],
        help="Energies for P(n_e|E) vs n_e plots",
    )
    ap.add_argument("--skip-heatmap", action="store_true")
    args = ap.parse_args()

    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    cases = build_cases(manifest)
    outdir = args.outdir
    ne_max = args.ne_max

    by_gap: dict[float, list[dict]] = {}
    by_scen: dict[str, list[dict]] = {"D-equal": [], "B-thresh": []}
    for c in cases:
        by_gap.setdefault(c["band_gap_eV"], []).append(c)
        by_scen.setdefault(c["scenario"], []).append(c)
    for g in by_gap:
        by_gap[g].sort(key=lambda x: x["scenario"])
    for sc in by_scen:
        by_scen[sc].sort(key=lambda x: x["band_gap_eV"])

    # Per scissor gap: D-equal vs B-thresh P(n_e|E) tables
    for gap, gap_cases in sorted(by_gap.items()):
        if len(gap_cases) < 2:
            continue
        gtag = gap_tag(gap)
        series = []
        bar_series = []
        for c in gap_cases:
            E, P = load_table(c)
            series.append((short_legend(c), E, P, c))
            bar_series.append((short_legend(c), E, P))
        plot_pne_vs_E_subplots(
            series,
            outdir / f"compare_p100K_D-equal_vs_B-thresh__{gtag}.pdf",
            rf"$P(n_e \mid E)$ tables — $E_{{\mathrm{{gap}}}} = {gap:g}$ eV",
            args.Emax,
            ne_max=ne_max,
        )
        plot_pne_vs_ne_at_E(
            bar_series,
            args.E_slices,
            outdir / f"compare_p100K_P_vs_ne_at_E__{gtag}.pdf",
            rf"$P(n_e \mid E_r)$ at fixed $E_r$ — $E_{{\mathrm{{gap}}}} = {gap:g}$ eV",
            ne_max=ne_max,
        )

    # All D-equal / all B-thresh ionization tables (6 curves per n_e panel)
    for scen, scen_cases in by_scen.items():
        if not scen_cases:
            continue
        series = []
        for c in scen_cases:
            E, P = load_table(c)
            g = c["band_gap_eV"]
            series.append((pheno_param_label_mpl(g, c["eh_pair_eV"]), E, P, c))
        tag = scen.replace("-", "_")
        plot_pne_vs_E_subplots(
            series,
            outdir / f"compare_p100K_{tag}_all_gaps.pdf",
            rf"$P(n_e \mid E)$ — all {scen} tables ($n_e \leq {ne_max}$)",
            args.Emax,
            ne_max=ne_max,
        )

    # Every case: heatmap of table (optional QA)
    if not args.skip_heatmap:
        for c in cases:
            E, P = load_table(c)
            plot_heatmap(
                E,
                P,
                outdir / "heatmaps" / f"p100K_{c['id']}.pdf",
                ion_legend(c),
                args.Emax,
                ne_max=ne_max,
            )

    print(f"\nP(n_e|E) table plots under {outdir}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
