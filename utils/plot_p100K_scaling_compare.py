#!/usr/bin/env python3
# ============================================================================
#  CCDarkSens — plot_p100K_scaling_compare
#  Comparison plots for scaled p100K tables vs Si reference across band-gap scenarios
#
#  Author: Diego Venegas-Vargas
# ============================================================================
"""Comparison plots for scaled p100K tables vs Si reference (extended grid)."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import List, Tuple

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from band_gap_plot_labels import eh_for_scenario, pheno_param_label_mpl  # noqa: E402
from build_p100K_scaled import (  # noqa: E402
    DEFAULT_E_MIN_EV,
    REF_EH_EV,
    REF_GAP_EV,
    build_scaled_table,
    extend_energy_grid,
    load_p100k_csv,
)

REF_CSV = Path("data/p100K_table.csv")


def mean_ne(P: np.ndarray) -> np.ndarray:
    ne = np.arange(1, P.shape[0] + 1, dtype=float)
    return (P * ne[:, None]).sum(axis=0)


def load_manifest(path: Path) -> List[dict]:
    with open(path, encoding="utf-8") as f:
        man = json.load(f)
    out = [man["reference"]]
    out.extend(man.get("scenarios", []))
    return out


def extended_ref() -> Tuple[np.ndarray, np.ndarray, str]:
    """Reference on same extended E grid as scaled tables (identity map)."""
    E_raw, P_raw = load_p100k_csv(REF_CSV)
    E = extend_energy_grid(E_raw, E_min=DEFAULT_E_MIN_EV)
    P = build_scaled_table(E, P_raw, E_raw, REF_GAP_EV, REF_EH_EV)
    return E, P, "Si ref (1.2 / 3.8 eV)"


def style_axes(ax, Emax: float, ylab: str) -> None:
    ax.set_xlabel(r"Recoil energy $E$ [eV]")
    ax.set_ylabel(ylab)
    ax.set_xlim(0, Emax)
    ax.grid(True, alpha=0.3)


def plot_overlay(
    E_ref: np.ndarray,
    P_ref: np.ndarray,
    ref_label: str,
    scaled: List[Tuple[str, str, np.ndarray, np.ndarray]],
    out: Path,
    Emax: float,
    quantity: str,
) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(9, 5.5))
    m = E_ref <= Emax

    if quantity == "P1":
        ax.plot(E_ref[m], P_ref[0, m], "k-", lw=2.5, label=ref_label)
        ylab = r"$P(n_e=1 \mid E)$"
        for label, sc, E, P in scaled:
            mm = E <= Emax
            ax.plot(E[mm], P[0, mm], lw=1.4, label=f"{sc}")
    else:
        ax.plot(E_ref[m], mean_ne(P_ref)[m], "k-", lw=2.5, label=ref_label)
        ylab = r"$\langle n_e \rangle$"
        for label, sc, E, P in scaled:
            mm = E <= Emax
            ax.plot(E[mm], mean_ne(P)[mm], lw=1.4, label=f"{sc}")

    style_axes(ax, Emax, ylab)
    ax.legend(loc="best", fontsize=9)
    fig.suptitle(f"p100K comparison — {ylab}  ($E_\\mathrm{{max}}$={Emax} eV)", fontsize=11)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"[ok] {out}")


def plot_pne_compare_scenarios(
    E_ref: np.ndarray,
    P_ref: np.ndarray,
    ref_label: str,
    scaled: List[Tuple[str, str, np.ndarray, np.ndarray]],
    out: Path,
    Emax: float,
    scenarios: tuple[str, ...] = ("B-thresh", "D-equal"),
) -> None:
    """Figure 3 style: P(n_e=1|E_r) and P(n_e=2|E_r) vs Si ref for pheno tiers."""
    out.parent.mkdir(parents=True, exist_ok=True)
    picks = [(lbl, sc, E, P) for lbl, sc, E, P in scaled if sc in scenarios]
    if not picks:
        print(f"[warn] no curves for scenarios {scenarios}; skip {out}")
        return

    m_ref = E_ref <= Emax
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharex=True)

    for ax, ne_idx, ylab in zip(
        axes,
        (0, 1),
        (r"$P(n_e = 1 \mid E_r)$", r"$P(n_e = 2 \mid E_r)$"),
    ):
        ax.plot(E_ref[m_ref], P_ref[ne_idx, m_ref], "k-", lw=2.5, label=ref_label)
        for label, sc, E, P in picks:
            mm = E <= Emax
            if P.shape[0] <= ne_idx:
                continue
            # label like "0.1 eV gap, 0.1 eV eh (equal-scale)" -> parse gap from text
            gap_ev = float(label.split()[0])
            eh_ev = eh_for_scenario(gap_ev, sc)
            leg = pheno_param_label_mpl(gap_ev, eh_ev)
            ax.plot(E[mm], P[ne_idx, mm], lw=1.6, label=leg)
        ax.set_ylabel(ylab)
        ax.set_xlim(0, Emax)
        ax.set_ylim(0, 1.02)
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best", fontsize=8)

    axes[1].set_xlabel(r"Recoil energy $E_r$ [eV]")
    fig.suptitle(
        r"Scaled p100K: threshold shift ($B$-thresh) vs coupled scales ($D$-equal)",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"[ok] {out}")


def plot_grid(
    E_ref: np.ndarray,
    P_ref: np.ndarray,
    ref_label: str,
    scaled: List[Tuple[str, str, np.ndarray, np.ndarray]],
    out: Path,
    Emax: float,
) -> None:
    n = len(scaled)
    ncols = 2
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(11, 3.2 * nrows), squeeze=False)
    m_ref = E_ref <= Emax

    for idx, (label, sc, E, P) in enumerate(scaled):
        r, c = divmod(idx, ncols)
        ax = axes[r, c]
        mm = E <= Emax
        ax.plot(E_ref[m_ref], P_ref[0, m_ref], "k--", lw=1.2, alpha=0.7, label="Si ref")
        ax.plot(E_ref[m_ref], mean_ne(P_ref)[m_ref], "k:", lw=1.2, alpha=0.7)
        ax.plot(E[mm], P[0, mm], lw=1.8, color="C0", label=r"$P(n_e{=}1|E)$")
        ax.plot(E[mm], mean_ne(P)[mm], lw=1.8, color="C1", label=r"$\langle n_e\rangle$")
        ax.set_title(f"{sc}\n{label}", fontsize=9)
        style_axes(ax, Emax, "probability / mean")
        ax.legend(fontsize=7, loc="upper left")

    for idx in range(n, nrows * ncols):
        r, c = divmod(idx, ncols)
        axes[r, c].set_visible(False)

    fig.suptitle("Scaled p100K vs Si reference (dashed/black)", fontsize=12)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"[ok] {out}")


def plot_multiplicity_bars(
    scaled: List[Tuple[str, str, np.ndarray, np.ndarray]],
    out: Path,
    energies: List[float],
    ne_max: int = 6,
) -> None:
    """Bar charts of P(n_e|E) at fixed E for each scenario."""
    n = len(scaled)
    fig, axes = plt.subplots(1, n, figsize=(3.2 * n, 4), squeeze=False)
    x = np.arange(1, ne_max + 1)

    for ax, (label, sc, E, P) in zip(axes[0], scaled):
        for e0 in energies:
            row = np.zeros(P.shape[0])
            for k in range(P.shape[0]):
                row[k] = np.interp(e0, E, P[k], left=0.0, right=0.0)
            ax.bar(x + 0.12 * energies.index(e0), row[:ne_max], width=0.1, alpha=0.75, label=f"{e0} eV")
        ax.set_xticks(x)
        ax.set_xlabel(r"$n_e$")
        ax.set_ylabel(r"$P(n_e \mid E)$")
        ax.set_title(sc, fontsize=9)
        ax.set_ylim(0, 1.05)
        ax.legend(fontsize=7)
        ax.grid(True, axis="y", alpha=0.3)

    fig.suptitle("Multiplicity distributions at selected energies")
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"[ok] {out}")


def tab20_colors(n: int) -> List[tuple]:
    """Matplotlib tab20 RGB (matches utils/plot_p100K_table.cc palette)."""
    rgb = [
        (31, 119, 180), (255, 127, 14), (44, 160, 44), (214, 39, 40), (148, 103, 189),
        (140, 86, 75), (227, 119, 194), (127, 127, 127), (188, 189, 34), (23, 190, 207),
        (174, 199, 232), (255, 187, 120), (152, 223, 138), (255, 152, 150), (197, 176, 213),
        (196, 156, 148), (247, 182, 210), (199, 199, 199), (219, 219, 141), (158, 218, 229),
    ]
    out = []
    for i in range(n):
        r, g, b = rgb[i % len(rgb)]
        out.append((r / 255.0, g / 255.0, b / 255.0))
    return out


def plot_p100K_Pne_fan(
    E: np.ndarray,
    P: np.ndarray,
    title: str,
    out: Path,
    Emax: float = 50.0,
    ne_max: int = 20,
) -> None:
    """All P(n_e|E_r) curves on one axes (like plot_p100K_table.cc / PRD-style figure)."""
    out.parent.mkdir(parents=True, exist_ok=True)
    k = min(ne_max, P.shape[0])
    m = (E >= 0) & (E <= Emax)
    E_plot = E[m]

    fig, ax = plt.subplots(figsize=(11, 7.5))
    colors = tab20_colors(k)

    for j in range(k):
        ax.plot(
            E_plot,
            P[j, m],
            lw=1.8,
            color=colors[j],
            label=rf"$n_e = {j + 1}$",
        )

    ax.set_xlim(0, Emax)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel(r"Recoil energy $E_r$ [eV]")
    ax.set_ylabel(r"$P(n_e \mid E_r)$")
    ax.set_title(title, fontsize=11)
    ax.grid(True, alpha=0.25)
    ax.legend(
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        fontsize=8,
        frameon=False,
        ncol=1,
    )
    fig.tight_layout()
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"[ok] {out}")


def scenario_fan_basename(ent: dict) -> str:
    sc = ent.get("scenario", "case")
    gap = ent["band_gap_eV"]
    eh = ent["eh_pair_eV"]
    g = str(gap).replace(".", "p")
    h = str(eh).replace(".", "p")
    return f"Pne_fan_{sc}_gap{g}_eh{h}"


def plot_ratio_mean_ne(
    E_ref: np.ndarray,
    P_ref: np.ndarray,
    scaled: List[Tuple[str, str, np.ndarray, np.ndarray]],
    out: Path,
    Emax: float,
) -> None:
    fig, ax = plt.subplots(figsize=(9, 5))
    mn_ref = mean_ne(P_ref)
    m = E_ref <= Emax
    ax.axhline(1.0, color="k", ls="--", lw=0.8)

    for label, sc, E, P in scaled:
        mm = E <= Emax
        ratio = np.divide(
            mean_ne(P)[mm],
            np.interp(E[mm], E_ref, mn_ref, left=np.nan, right=np.nan),
            where=np.interp(E[mm], E_ref, mn_ref) > 0,
        )
        ax.plot(E[mm], ratio, lw=1.4, label=sc)

    style_axes(ax, Emax, r"$\langle n_e \rangle_\mathrm{pheno} / \langle n_e \rangle_\mathrm{ref}$")
    ax.legend(fontsize=9)
    fig.suptitle("Mean multiplicity ratio to Si reference (extended grid)")
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"[ok] {out}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", default="configs/band_gap_pheno_scenarios.json")
    ap.add_argument("--outdir", default="outplots/band_gap_pheno/step2_p100K")
    ap.add_argument(
        "--Emax-pheno",
        type=float,
        default=5.0,
        help="Zoom for pheno / turn-on region [eV]",
    )
    ap.add_argument(
        "--Emax-low",
        type=float,
        default=2.0,
        help="Tight zoom on turn-on [eV]",
    )
    ap.add_argument(
        "--Emax-high",
        type=float,
        default=50.0,
        help="Full table range (match original p100K up to 50 eV)",
    )
    ap.add_argument(
        "--ne-max",
        type=int,
        default=20,
        help="Max n_e curves in P(n_e|E_r) fan plots",
    )
    ap.add_argument(
        "--skip-fan",
        action="store_true",
        help="Skip per-case P(n_e|E_r) fan plots",
    )
    args = ap.parse_args()

    outdir = Path(args.outdir)
    E_ref, P_ref, ref_label = extended_ref()

    scaled: List[Tuple[str, str, np.ndarray, np.ndarray]] = []
    for ent in load_manifest(Path(args.manifest)):
        if ent.get("scenario") == "ref":
            continue
        csv = Path(ent["ionization_csv"])
        E, P = load_p100k_csv(csv)
        scaled.append((ent.get("label", csv.name), ent.get("scenario", "?"), E, P))

    fan_dir = outdir / "Pne_fan"
    if not args.skip_fan:
        for ent in load_manifest(Path(args.manifest)):
            csv = Path(ent["ionization_csv"])
            if not csv.exists():
                print(f"[warn] skip fan plot, missing {csv}")
                continue
            E, P = load_p100k_csv(csv)
            base = scenario_fan_basename(ent)
            label = ent.get("label", csv.name)
            sc = ent.get("scenario", "?")
            plot_p100K_Pne_fan(
                E,
                P,
                f"{label}\n({sc}, gap={ent['band_gap_eV']} eV, eh={ent['eh_pair_eV']} eV)",
                fan_dir / f"{base}.pdf",
                Emax=args.Emax_high,
                ne_max=args.ne_max,
            )
        # Extended Si ref on same 0.05–50 eV grid as scaled tables
        plot_p100K_Pne_fan(
            E_ref,
            P_ref,
            f"{ref_label}\n(extended grid, identity map)",
            fan_dir / "Pne_fan_ref_extended.pdf",
            Emax=args.Emax_high,
            ne_max=args.ne_max,
        )

    # Pheno zoom (0–5 eV) and turn-on (0–2 eV)
    plot_overlay(
        E_ref, P_ref, ref_label, scaled,
        outdir / "compare_Pne1_pheno_0to5eV.pdf", args.Emax_pheno, "P1",
    )
    plot_overlay(
        E_ref, P_ref, ref_label, scaled,
        outdir / "compare_Pne1_lowE_0to2eV.pdf", args.Emax_low, "P1",
    )
    plot_overlay(
        E_ref, P_ref, ref_label, scaled,
        outdir / "compare_mean_ne_pheno_0to5eV.pdf", args.Emax_pheno, "mean",
    )
    plot_overlay(
        E_ref, P_ref, ref_label, scaled,
        outdir / "compare_mean_ne_lowE_0to2eV.pdf", args.Emax_low, "mean",
    )
    # Full original table range (1.1–50 eV on ref; scaled from 0.05 eV)
    plot_overlay(
        E_ref, P_ref, ref_label, scaled,
        outdir / "compare_Pne1_0to50eV.pdf", args.Emax_high, "P1",
    )
    plot_overlay(
        E_ref, P_ref, ref_label, scaled,
        outdir / "compare_mean_ne_0to50eV.pdf", args.Emax_high, "mean",
    )
    plot_grid(E_ref, P_ref, ref_label, scaled, outdir / "compare_grid_scenarios.pdf", args.Emax_low)
    plot_pne_compare_scenarios(
        E_ref,
        P_ref,
        ref_label,
        scaled,
        outdir / "Pne_compare_scenarios.pdf",
        args.Emax_pheno,
    )
    plot_ratio_mean_ne(E_ref, P_ref, scaled, outdir / "compare_ratio_mean_ne_0to50eV.pdf", args.Emax_high)
    plot_multiplicity_bars(
        scaled, outdir / "compare_multiplicity_bars.pdf",
        energies=[0.1, 0.2, 0.5, 1.0, 2.0],
    )

    # Per-gap overlays (same scissor, different ionization)
    by_gap: dict = {}
    for ent in load_manifest(Path(args.manifest)):
        if ent.get("scenario") == "ref":
            continue
        g = float(ent["band_gap_eV"])
        by_gap.setdefault(g, []).append(ent)

    for gap, ents in sorted(by_gap.items()):
        gtag = str(gap).replace(".", "p")
        sub = []
        for ent in ents:
            csv = Path(ent["ionization_csv"])
            E, P = load_p100k_csv(csv)
            sub.append((ent.get("label", ""), ent.get("scenario", "?"), E, P))
        plot_overlay(
            E_ref, P_ref, ref_label, sub,
            outdir / f"compare_gap{gtag}_Pne1_lowE.pdf", args.Emax_low, "P1",
        )
        plot_overlay(
            E_ref, P_ref, ref_label, sub,
            outdir / f"compare_gap{gtag}_mean_ne_lowE.pdf", args.Emax_low, "mean",
        )


if __name__ == "__main__":
    main()
