#!/usr/bin/env python3
# ============================================================================
#  CCDarkSens — eval_p100K_scaling
#  Evaluate and QA p100K table scaling across band-gap scenarios (Step 2)
#
#  Author: Diego Venegas-Vargas
# ============================================================================
"""Evaluate and compare p100K scaling scenarios (Step 2 QA).

Loads reference + scaled tables from band_gap_pheno_scenarios.json (or CLI),
checks anchor points, prints summary tables, and writes comparison figures.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_p100K_scaled import (  # noqa: E402
    REF_EH_EV,
    REF_GAP_EV,
    load_p100k_csv,
    map_E_prime,
)

REF_CSV = Path("data/p100K_table.csv")


def interp_at(E: np.ndarray, P: np.ndarray, E0: float) -> np.ndarray:
    """Linear interp each n_e row at E0."""
    out = np.zeros(P.shape[0])
    for n in range(P.shape[0]):
        out[n] = np.interp(E0, E, P[n], left=0.0, right=0.0)
    return out


def mean_ne(E: np.ndarray, P: np.ndarray) -> np.ndarray:
    ne = np.arange(1, P.shape[0] + 1, dtype=float)
    return (P * ne[:, None]).sum(axis=0)


def prob_ge1(P: np.ndarray) -> np.ndarray:
    return P.sum(axis=0)


def load_manifest(path: Path) -> List[dict]:
    with open(path, encoding="utf-8") as f:
        man = json.load(f)
    entries = [man["reference"]]
    entries.extend(man.get("scenarios", []))
    return entries


def anchor_check(
    E_ref: np.ndarray,
    P_ref: np.ndarray,
    gap_new: float,
    eh_new: float,
    E: np.ndarray,
    P: np.ndarray,
) -> Tuple[float, float]:
    """At E=gap_new, compare P(n=1) to ref at E=1.2 eV."""
    pref = interp_at(E_ref, P_ref, REF_GAP_EV)
    pnew = interp_at(E, P, gap_new)
    return float(pref[0]), float(pnew[0])


def print_summary_table(
    label: str,
    scenario: str,
    gap: float,
    eh: float,
    E: np.ndarray,
    P: np.ndarray,
    energies: List[float],
) -> None:
    print(f"\n--- {label} [{scenario}]  gap={gap} eV  eh={eh} eV ---")
    hdr_Ep = "E'(E)"
    print(f"{'E [eV]':>8}  {hdr_Ep:>8}  {'P(n=1)':>8}  {'<n_e>':>8}  {'sum P':>8}")
    for E0 in energies:
        if E0 < gap - 1e-9:
            Ep = float("nan")
            p1, mn, sp = 0.0, 0.0, 0.0
        else:
            Ep = map_E_prime(E0, gap, REF_GAP_EV, eh, REF_EH_EV)
            row = interp_at(E, P, E0)
            p1 = row[0]
            mn = float((row * np.arange(1, len(row) + 1)).sum())
            sp = float(row.sum())
        ep_str = f"{Ep:8.2f}" if Ep == Ep else f"{'—':>8}"
        print(f"{E0:8.2f}  {ep_str}  {p1:8.4f}  {mn:8.3f}  {sp:8.4f}")


def plot_group(
    E_ref: np.ndarray,
    P_ref: np.ndarray,
    curves: List[Tuple[str, np.ndarray, np.ndarray, str]],
    out: Path,
    title: str,
    Emax: float,
) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    mask_ref = E_ref <= Emax

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    ax_p1, ax_mean, ax_sum, ax_multi = axes.ravel()

    ax_p1.plot(E_ref[mask_ref], P_ref[0, mask_ref], "k-", lw=2, label="Si ref")
    ax_mean.plot(E_ref[mask_ref], mean_ne(E_ref, P_ref)[mask_ref], "k-", lw=2, label="Si ref")
    ax_sum.plot(E_ref[mask_ref], prob_ge1(P_ref)[mask_ref], "k-", lw=2, label="Si ref")

    colors = plt.cm.tab10(np.linspace(0, 0.9, max(len(curves), 1)))
    for i, (lab, E, P, sc) in enumerate(curves):
        m = E <= Emax
        c = colors[i]
        ax_p1.plot(E[m], P[0, m], lw=1.5, color=c, label=f"{sc}: {lab}")
        ax_mean.plot(E[m], mean_ne(E, P)[m], lw=1.5, color=c, label=f"{sc}")
        ax_sum.plot(E[m], prob_ge1(P)[m], lw=1.5, color=c, label=f"{sc}")

    ne_show = min(5, P_ref.shape[0])
    x = np.arange(ne_show) + 1
    if curves:
        _, E_p, P_p, sc_p = curves[0]
        e_pts = [0.1, 0.2, 0.5, 1.0, 1.2, 2.0]
        for j, e0 in enumerate(e_pts):
            if e0 > Emax:
                continue
            row = interp_at(E_p, P_p, e0)
            ax_multi.bar(x + 0.15 * j, row[:ne_show], width=0.12, alpha=0.7, label=f"E={e0} eV")
        ax_multi.set_title(f"P(n_e|E) bars — {sc_p} @ selected E")
    ax_multi.set_xlabel(r"$n_e$")
    ax_multi.set_ylabel(r"$P(n_e|E)$")
    ax_multi.legend(fontsize=7, ncol=2)

    for ax, ylab in [
        (ax_p1, r"$P(n_e=1\mid E)$"),
        (ax_mean, r"$\langle n_e \rangle$"),
        (ax_sum, r"$\sum_n P(n_e\mid E)$"),
    ]:
        ax.set_xlabel(r"Recoil energy $E$ [eV]")
        ax.set_ylabel(ylab)
        ax.set_xlim(0, Emax)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7, loc="best")

    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"[ok] {out}")


def plot_all_overlay(
    E_ref: np.ndarray,
    P_ref: np.ndarray,
    all_curves: List[Tuple[str, np.ndarray, np.ndarray, str, float, float]],
    out: Path,
    Emax: float,
) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))
    m = E_ref <= Emax
    ax1.plot(E_ref[m], P_ref[0, m], "k-", lw=2.5, label="Si ref (1.2/3.8)")
    ax2.plot(E_ref[m], mean_ne(E_ref, P_ref)[m], "k-", lw=2.5, label="Si ref")

    for lab, E, P, sc, gap, eh in all_curves:
        if sc == "ref":
            continue
        mm = E <= Emax
        short = f"{sc} gap={gap:g} eh={eh:g}"
        ax1.plot(E[mm], P[0, mm], lw=1.3, label=short)
        ax2.plot(E[mm], mean_ne(E, P)[mm], lw=1.3, label=short)

    for ax, ylab in [(ax1, r"$P(n_e=1\mid E)$"), (ax2, r"$\langle n_e \rangle$")]:
        ax.set_xlabel(r"$E$ [eV]")
        ax.set_ylabel(ylab)
        ax.set_xlim(0, Emax)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7, loc="best")
    fig.suptitle("All p100K scaling scenarios vs Si reference")
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"[ok] {out}")


def main() -> None:
    ap = argparse.ArgumentParser(description="Evaluate p100K scaling cases.")
    ap.add_argument("--manifest", default="configs/band_gap_pheno_scenarios.json")
    ap.add_argument("--outdir", default="outplots/band_gap_pheno/step2_p100K")
    ap.add_argument("--Emax", type=float, default=5.0, help="Plot x-axis max [eV]")
    ap.add_argument(
        "--E-check",
        nargs="+",
        type=float,
        default=[0.05, 0.1, 0.15, 0.2, 0.5, 1.0, 1.2, 2.0, 3.0, 5.0],
    )
    args = ap.parse_args()

    E_ref, P_ref = load_p100k_csv(REF_CSV)
    entries = load_manifest(Path(args.manifest))
    outdir = Path(args.outdir)

    print("=" * 72)
    print("p100K scaling evaluation (anchored map)")
    print(f"  ref: E_gap={REF_GAP_EV} eV, eh={REF_EH_EV} eV  ->  {REF_CSV}")
    print("=" * 72)

    all_curves: List[Tuple[str, np.ndarray, np.ndarray, str, float, float]] = []
    by_gap: Dict[float, List[Tuple[str, np.ndarray, np.ndarray, str]]] = {}

    for ent in entries:
        csv = Path(ent["ionization_csv"])
        gap = float(ent["band_gap_eV"])
        eh = float(ent["eh_pair_eV"])
        sc = ent.get("scenario", "?")
        lab = ent.get("label", csv.name)

        E, P = load_p100k_csv(csv)
        all_curves.append((lab, E, P, sc, gap, eh))
        if sc != "ref":
            by_gap.setdefault(gap, []).append((lab, E, P, sc))

        if sc == "ref":
            print(f"\n[ref] {lab}")
            continue

        p_ref1, p_new1 = anchor_check(E_ref, P_ref, gap, eh, E, P)
        # Prefer exact grid row at gap if present
        at_gap = np.where(np.isclose(E, gap))[0]
        if len(at_gap):
            p_new1 = float(P[0, at_gap[0]])
        print(f"\n[{sc}] {lab}")
        print(f"  file: {csv}")
        print(f"  scale factor eh_ref/eh_new = {REF_EH_EV / eh:.4g}")
        ok = abs(p_ref1 - p_new1) < 0.02
        print(
            f"  anchor @ E={gap} eV: P(n=1|ref@1.2)={p_ref1:.4f}  "
            f"P(n=1|new)={p_new1:.4f}  {'OK' if ok else 'CHECK'}"
        )

        idx = np.where(P[0] > 1e-6)[0]
        turn_on = float(E[idx[0]]) if len(idx) else float("nan")
        print(f"  first P(n=1)>0 at E={turn_on:.3g} eV (target gap={gap} eV)")

        print_summary_table(lab, sc, gap, eh, E, P, args.E_check)

    for gap, curves in sorted(by_gap.items()):
        gtag = str(gap).replace(".", "p")
        plot_group(
            E_ref,
            P_ref,
            curves,
            outdir / f"eval_gap{gtag}_panel.pdf",
            title=f"p100K scaling @ scissor gap {gap} eV",
            Emax=args.Emax,
        )

    plot_all_overlay(E_ref, P_ref, all_curves, outdir / "eval_all_scenarios.pdf", args.Emax)
    plot_all_overlay(E_ref, P_ref, all_curves, outdir / "eval_all_scenarios_lowE.pdf", Emax=2.0)

    summary_path = outdir / "eval_summary.txt"
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("p100K scaling evaluation summary\n")
        f.write(f"ref: gap={REF_GAP_EV} eh={REF_EH_EV}\n\n")
        for ent in entries:
            if ent.get("scenario") == "ref":
                continue
            csv = Path(ent["ionization_csv"])
            gap = float(ent["band_gap_eV"])
            eh = float(ent["eh_pair_eV"])
            E, P = load_p100k_csv(csv)
            f.write(f"{ent.get('scenario')} | gap={gap} eh={eh} | {csv.name}\n")
            for e0 in args.E_check:
                if e0 < gap:
                    f.write(f"  E={e0:.2f}: below gap -> P=0\n")
                else:
                    row = interp_at(E, P, e0)
                    idx = min(int(np.searchsorted(E, e0)), len(E) - 1)
                    f.write(
                        f"  E={e0:.2f}: P1={row[0]:.4f} "
                        f"<ne>={mean_ne(E, P)[idx]:.3f}\n"
                    )
            f.write("\n")
    print(f"[ok] {summary_path}")


if __name__ == "__main__":
    main()
