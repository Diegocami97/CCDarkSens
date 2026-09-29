#!/usr/bin/env python3
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: plot_Sr2Cb2Sd_ne_spectra.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  plot_Sr2Cb2Sd_ne_spectra.py -- Phase 2: n_e-space signal spectra for
#  Sr2Cb2Sd (Option C backgrounds).
# ============================================================================
"""
Sr2Cb2Sd n_e imaging figures (2 PDFs: heavy + light mediator).

Each figure: S_true | S_obs panels with
  - 3 signal bars: DAMIC-M, SrCd indirect, SrCd direct (1 kg·yr, 1× DC)
  - 3 B_tot lines: DAMIC-M reference @ 1× / 100× / 1000× DC (Si p100K flat fold)

  python3 utils/plot_Sr2Cb2Sd_ne_spectra.py
  python3 configs/Sr2Cb2Sd/write_ne_imaging_configs.py
  configs/Sr2Cb2Sd/run_phase2_ne_imaging.sh
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import uproot  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUTBASE = ROOT / "outputs" / "Sr2Cb2Sd" / "ne_imaging"
OUTPLOTS = ROOT / "outplots" / "Sr2Cb2Sd"

NE_BINS = [1, 2, 3, 4, 5]
MCHI_MEV = 1.000194
SIGMA_E_CM2 = 1.1e-35
LOG_FLOOR = 1e-3

COLOR_DAMIC = "#D4A017"
COLOR_INDIRECT = "#E41A1C"
COLOR_DIRECT = "#4DAF4A"

SIGNAL_SCENARIOS = [
    ("damic", r"DAMIC-M $(1.2,\ 3.8)$", COLOR_DAMIC),
    ("srcd_indirect", r"$(E_{\mathrm{gap}},\ \varepsilon_h)=(0.556,\ 2.06)$", COLOR_INDIRECT),
    ("srcd_direct", r"$(E_{\mathrm{gap}},\ \varepsilon_h)=(0.603,\ 2.19)$", COLOR_DIRECT),
]

BKG_TIERS = [
    ("1x", 0.00365, "0.25", (4, 2)),
    ("dc100x", 0.365, "0.45", (2, 2)),
    ("dc1000x", 3.65, "0.65", (1, 1)),
]

MEDIATORS = [
    ("heavy", r"Heavy mediator ($F_{\mathrm{DM}}=1$)"),
    ("light", r"Light mediator ($F_{\mathrm{DM}}\propto 1/q^2$)"),
]


# ----------------------------------------------------------------------------
# scan_path
#   Path of the scan ROOT file of a scenario, mediator and dark-current tier.
# ----------------------------------------------------------------------------
def scan_path(scen_key: str, med: str, dc_suffix: str) -> Path:
    return OUTBASE / f"{scen_key}_{med}_{dc_suffix}" / "scan_dmelectron_pattern.root"


# ----------------------------------------------------------------------------
# _first_key
#   Name (without the cycle number) of the first object in a ROOT file whose name starts with prefix, or None.
# ----------------------------------------------------------------------------
def _first_key(f, prefix: str) -> str | None:
    for k in f.keys():
        kk = k.split(";")[0]
        if kk.startswith(prefix):
            return kk
    return None


# ----------------------------------------------------------------------------
# _values_by_ne
#   Histogram content at the bin nearest to each n_e of the list.
# ----------------------------------------------------------------------------
def _values_by_ne(hist, ne_list: list[int] = NE_BINS) -> np.ndarray:
    centers = hist.axis().centers()
    vals = hist.values()
    out = []
    for ne in ne_list:
        i = int(np.argmin(np.abs(centers - ne)))
        out.append(float(vals[i]))
    return np.array(out)


# ----------------------------------------------------------------------------
# load_scan
#   Read S_true(n_e), S_obs(n_e) and the total background B_tot(n_e) from a scan file (None entries and a warning if the file is missing).
# ----------------------------------------------------------------------------
def load_scan(scen_key: str, med: str, dc_suffix: str) -> dict:
    path = scan_path(scen_key, med, dc_suffix)
    out = {"S_true": None, "S_obs": None, "B_tot": None}
    if not path.is_file():
        print(f"  WARN: missing {path}")
        return out
    with uproot.open(path) as f:
        st_key = _first_key(f, "S_true_ne__")
        so_key = _first_key(f, "S_obs_ne__")
        if st_key:
            out["S_true"] = _values_by_ne(f[st_key])
        if so_key:
            out["S_obs"] = _values_by_ne(f[so_key])
        if "B_tot_ne" in [k.split(";")[0] for k in f.keys()]:
            out["B_tot"] = _values_by_ne(f["B_tot_ne"])
    return out


# ----------------------------------------------------------------------------
# dc_rate_text
#   Legend text of a dark-current rate in e-/pix/yr.
# ----------------------------------------------------------------------------
def dc_rate_text(lam_e_per_pix_per_year: float) -> str:
    return rf"DC $= {lam_e_per_pix_per_year:g}\ \mathrm{{e}}^-\!/\mathrm{{pix}}/\mathrm{{yr}}$"


# ----------------------------------------------------------------------------
# _draw_grouped_bars
#   Draw one group of bars per n_e, one bar per scenario, dropping values below the log floor.
# ----------------------------------------------------------------------------
def _draw_grouped_bars(ax, data_by_scen: list, key: str, scenarios) -> None:
    x = np.arange(len(NE_BINS))
    n = len(scenarios)
    width = 0.8 / n
    for j, (_, label, color) in enumerate(scenarios):
        vals = data_by_scen[j][key]
        if vals is None:
            continue
        plotted = np.where(vals > LOG_FLOOR, vals, np.nan)
        offs = x + (j - (n - 1) / 2.0) * width
        ax.bar(
            offs,
            plotted,
            width=width,
            color=color,
            edgecolor="black",
            linewidth=0.4,
            label=label,
            zorder=3,
        )


# ----------------------------------------------------------------------------
# _draw_backgrounds
#   Draw the background of each dark-current tier as short horizontal lines at every n_e bin.
# ----------------------------------------------------------------------------
def _draw_backgrounds(ax, bkg_data: list[tuple[np.ndarray | None, str, str, tuple]]) -> None:
    x = np.arange(len(NE_BINS))
    for b_tot, dc_label, color, _dashes in bkg_data:
        if b_tot is None:
            continue
        for i in range(len(NE_BINS)):
            b = b_tot[i]
            if b <= 0:
                continue
            ax.hlines(
                b,
                x[i] - 0.45,
                x[i] + 0.45,
                colors=color,
                linestyles="--",
                lw=1.6,
                zorder=4,
                label=rf"$B_{{\mathrm{{tot}}}}(n_e)$ ({dc_label} + flat, Si ion.)"
                if i == 0
                else None,
            )


# ----------------------------------------------------------------------------
# plot_mediator
#   Signal n_e spectra of the scenarios with the background tiers for one mediator, saved as PDF; returns 1 if there are no signal scans.
# ----------------------------------------------------------------------------
def plot_mediator(med_key: str, med_label: str, out_path: Path) -> int:
    signal_data = [
        load_scan(scen_key, med_key, "1x") for scen_key, _, _ in SIGNAL_SCENARIOS
    ]
    bkg_data = []
    for dc_suffix, lam, color, dashes in BKG_TIERS:
        row = load_scan("damic", med_key, dc_suffix)
        bkg_data.append((row["B_tot"], dc_rate_text(lam), color, dashes))

    if not any(d["S_true"] is not None for d in signal_data):
        print(f"ERROR: no signal scans for {med_key}", file=sys.stderr)
        return 1
    if not any(b[0] is not None for b in bkg_data):
        print(f"ERROR: no DAMIC B_tot scans for {med_key}", file=sys.stderr)
        return 1

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.2))

    for ax, key, subtitle in (
        (axes[0], "S_true", r"$S_{\mathrm{true}}(n_e)$ — ionization (FoldToNe)"),
        (axes[1], "S_obs", r"$S_{\mathrm{obs}}(n_e)=S_{\mathrm{true}}\cdot\varepsilon(n_e)$"),
    ):
        _draw_grouped_bars(ax, signal_data, key, SIGNAL_SCENARIOS)
        _draw_backgrounds(ax, bkg_data)
        ax.set_yscale("log")
        ax.set_xticks(np.arange(len(NE_BINS)))
        ax.set_xticklabels([str(ne) for ne in NE_BINS])
        ax.set_xlabel(r"$n_e$ [electrons]")
        ax.set_title(subtitle, fontsize=9)
        ax.grid(True, which="both", axis="y", ls=":", alpha=0.3)
        ax.legend(loc="upper right", fontsize=7.5, framealpha=0.92)

    axes[0].set_ylabel(med_label + "\nExpected counts")

    fig.suptitle(
        r"$n_e$-space signal spectrum — Sr$_2$Cd$_2$Sb$_2$ — "
        rf"$m_\chi={MCHI_MEV:.2f}$ MeV, $\bar\sigma_e={SIGMA_E_CM2:.1e}$ cm$^2$, 1 kg·yr",
        fontsize=12,
    )
    fig.text(
        0.5,
        0.005,
        r"Signals: material-specific ionization; backgrounds: DAMIC-M reference "
        r"($B_{\mathrm{tot}}$ at 1$\times$/100$\times$/1000$\times$ DC, flat via Si 1.2/3.8 eV). "
        r"Lower-gap material populates higher $n_e$ $\rightarrow$ discrimination against $n_e{=}1$ DC.",
        ha="center",
        va="bottom",
        fontsize=9,
        style="italic",
    )
    fig.tight_layout(rect=(0, 0.04, 1, 0.94))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    fig.savefig(out_path.with_suffix(".png"), dpi=150)
    plt.close(fig)
    print(f"Wrote {out_path}")
    return 0


# ----------------------------------------------------------------------------
# main
#   Command line: make the n_e spectrum figure for the heavy mediator, the light one, or both.
# ----------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mediator", choices=["heavy", "light", "both"], default="both")
    args = ap.parse_args()

    meds = MEDIATORS if args.mediator == "both" else [m for m in MEDIATORS if m[0] == args.mediator]
    rc = 0
    for med_key, med_label in meds:
        out = OUTPLOTS / f"ne_signal_spectrum_{med_key}.pdf"
        if plot_mediator(med_key, med_label, out) != 0:
            rc = 1
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
