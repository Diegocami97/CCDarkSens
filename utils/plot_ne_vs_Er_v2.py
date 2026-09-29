#!/usr/bin/env python3
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: plot_ne_vs_Er_v2.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  plot_ne_vs_Er_v2.py -- Figure 1: mean charge yield <n_e>(E_r) under two
#  gap/eh conventions
# ============================================================================
"""
Mean number of electron-hole pairs <n_e> vs recoil energy E_r for the band-gap
pheno study, shown side by side under two charge-yield conventions.

  python3 utils/plot_ne_vs_Er_v2.py
  python3 utils/plot_ne_vs_Er_v2.py --klein-gaps 0.1,0.3,0.5
  python3 utils/plot_ne_vs_Er_v2.py --y-max 15 --x-max 10

Output: outplots/band_gap_pheno/ne_imaging/ne_vs_Er_conventions.pdf
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "utils"))
from band_gap_klein_plot import (  # noqa: E402
    COLOR_SI_REF,
    DEFAULT_KLEIN_GAPS,
    SI_EGAP_EV,
    SI_EH_EV,
    filter_klein_gaps,
    gap_color,
    klein_eh,
    klein_mathtext_label,
    parse_float_list,
)

DEFAULT_OUT = ROOT / "outplots" / "band_gap_pheno" / "ne_imaging" / "ne_vs_Er_conventions.pdf"

X_MIN = 0.0
N_E_THRESHOLD = 1.0
ER_DAMIC_THRESHOLD_EV = 4.0


# ----------------------------------------------------------------------------
# build_scenarios
#   Curves to draw: the silicon reference (solid) and one Klein-formula case per gap (dashed), as (label, E_gap, eps_h, colour, line style).
# ----------------------------------------------------------------------------
def build_scenarios(klein_gaps: list[float]) -> list[tuple[str, float, float, str, str]]:
    scenarios: list[tuple[str, float, float, str, str]] = [
        (r"Si ref  $(1.2,\ 3.8)$", SI_EGAP_EV, SI_EH_EV, COLOR_SI_REF, "-"),
    ]
    for i, g in enumerate(klein_gaps):
        eh = klein_eh(g)
        scenarios.append((klein_mathtext_label(g, eh), g, eh, gap_color(i), "--"))
    return scenarios


# ----------------------------------------------------------------------------
# y_max_for_scenarios
#   Largest n_e reached at x_max by any scenario with either convention, used to set the y-axis range.
# ----------------------------------------------------------------------------
def y_max_for_scenarios(
    scenarios: list[tuple[str, float, float, str, str]], x_max: float
) -> float:
    e_max = np.array([x_max])
    peak = 0.0
    for _, eg, eh, _, _ in scenarios:
        for conv in (1, 2):
            if conv == 1:
                ne = (e_max - eg) / eh
            else:
                ne = e_max / eh
            ne = np.where(e_max >= eg, ne, 0.0)
            peak = max(peak, float(ne.max()))
    return float(np.ceil(peak * 1.15 + 0.5))


# ----------------------------------------------------------------------------
# mean_ne
#   Mean n_e versus recoil energy: convention 1 is (E - E_gap)/eps_h, convention 2 is E/eps_h; zero below the gap. Raises ValueError for any other convention.
# ----------------------------------------------------------------------------
def mean_ne(E_r: np.ndarray, E_gap: float, eps_h: float, convention: int) -> np.ndarray:
    if convention == 1:
        ne = (E_r - E_gap) / eps_h
    elif convention == 2:
        ne = E_r / eps_h
    else:
        raise ValueError(f"convention must be 1 or 2, got {convention}")
    ne = np.where(E_r >= E_gap, ne, 0.0)
    return np.clip(ne, 0.0, None)


# ----------------------------------------------------------------------------
# draw_panel
#   Draw one panel (one convention): the n_e(E) curves of all scenarios over shaded single-carrier and few-carrier bands.
# ----------------------------------------------------------------------------
def draw_panel(
    ax,
    convention: int,
    subtitle: str,
    scenarios: list[tuple[str, float, float, str, str]],
    x_max: float,
    y_max: float,
) -> None:
    E_r = np.linspace(X_MIN, x_max, 1600)
    bands = [
        (0.0, 1.5, "#7FB685", r"single-carrier ($n_e\sim1$)"),
        (1.5, 5.0, "#E8A87C", r"few-carrier ($n_e\sim2$–5)"),
        (5.0, y_max, "#9B8BB4", r"multi-carrier ($n_e\gtrsim5$)"),
    ]

    for y_lo, y_hi, color, _ in bands:
        ax.axhspan(y_lo, y_hi, color=color, alpha=0.10, zorder=0)

    for y_lo, y_hi, color, label in bands:
        y_mid = 0.5 * (max(y_lo, 0.0) + min(y_hi, y_max))
        ax.text(
            x_max * 0.995,
            y_mid,
            label,
            color=color,
            fontsize=8,
            ha="right",
            va="center",
            alpha=0.9,
            zorder=1,
        )

    for label, E_gap, eps_h, color, ls in scenarios:
        ne = mean_ne(E_r, E_gap, eps_h, convention)
        lw = 2.4 if label.startswith("Si ref") else 2.0
        ax.plot(E_r, ne, color=color, linestyle=ls, lw=lw, label=label, zorder=4)

        ne_ref = float(
            mean_ne(np.array([ER_DAMIC_THRESHOLD_EV]), E_gap, eps_h, convention)[0]
        )
        if ne_ref > 0.0 and ne_ref <= y_max:
            ax.plot(ER_DAMIC_THRESHOLD_EV, ne_ref, "o", color=color, ms=5, zorder=5)
            ax.annotate(
                rf"$\langle n_e\rangle={ne_ref:.1f}$",
                xy=(ER_DAMIC_THRESHOLD_EV, ne_ref),
                xytext=(-3, 4),
                textcoords="offset points",
                ha="right",
                va="bottom",
                color=color,
                fontsize=7,
                zorder=6,
            )

    ax.axvline(
        ER_DAMIC_THRESHOLD_EV,
        color="0.15",
        linestyle="-.",
        lw=1.4,
        label=rf"DAMIC-M $1e^-$ threshold ($\approx{ER_DAMIC_THRESHOLD_EV:.0f}$ eV)",
        zorder=7,
    )

    ax.set_xlim(X_MIN, x_max)
    ax.set_ylim(0.0, y_max)
    ax.set_xlabel(r"Recoil energy $E_r$ [eV]")
    ax.set_ylabel(r"Mean charge yield $\langle n_e\rangle$")
    ax.set_title(subtitle, fontsize=9, pad=8)
    ax.legend(loc="upper left", fontsize=8, framealpha=0.9)
    ax.grid(True, which="both", ls=":", alpha=0.3)


# ----------------------------------------------------------------------------
# main
#   Command line: make the two-convention n_e versus E_r figure for the chosen Klein gaps and save it to --out.
# ----------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument(
        "--klein-gaps",
        type=str,
        default=",".join(str(g) for g in DEFAULT_KLEIN_GAPS),
        metavar="CSV",
        help="Comma-separated E_gap values [eV] for Klein curves "
        f"(default: {','.join(str(g) for g in DEFAULT_KLEIN_GAPS)}). "
        "E_gap=1.2 is skipped (use Si ref).",
    )
    ap.add_argument("--x-max", type=float, default=8.0, help="E_r axis maximum [eV]")
    ap.add_argument(
        "--y-max",
        type=float,
        default=None,
        help="⟨n_e⟩ axis maximum (default: auto from curves at E_r=x-max)",
    )
    args = ap.parse_args()

    klein_gaps = filter_klein_gaps(parse_float_list(args.klein_gaps))
    if not klein_gaps:
        print("ERROR: no Klein gaps after filtering", file=sys.stderr)
        return 1

    scenarios = build_scenarios(klein_gaps)
    y_max = args.y_max if args.y_max is not None else y_max_for_scenarios(scenarios, args.x_max)

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12.0, 5.2), sharey=True)
    draw_panel(
        axA,
        convention=1,
        subtitle=r"Convention 1: $(E_r - E_{\mathrm{gap}})/\varepsilon_h$ "
        r"— gap energy consumed by first pair",
        scenarios=scenarios,
        x_max=args.x_max,
        y_max=y_max,
    )
    draw_panel(
        axB,
        convention=2,
        subtitle=r"Convention 2: $E_r/\varepsilon_h$ "
        r"— pair scale absorbs gap cost (consistent with PRD 102, 063026)",
        scenarios=scenarios,
        x_max=args.x_max,
        y_max=y_max,
    )

    fig.suptitle(
        r"Mean charge yield $\langle n_e\rangle$ vs recoil energy — band-gap pheno scenarios",
        fontsize=12,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out)
    fig.savefig(args.out.with_suffix(".png"), dpi=150)
    plt.close(fig)
    print(f"Klein gaps: {klein_gaps}")
    print(f"Wrote {args.out}")
    print(f"Wrote {args.out.with_suffix('.png')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
