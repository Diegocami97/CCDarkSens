#!/usr/bin/env python3
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: plot_band_gap_pheno_dRdE.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  plot_band_gap_pheno_dRdE.py -- Overlay dR/dE spectra from QCDark2 rate
#  CSVs for band-gap pheno Step 1
# ============================================================================
"""Overlay dR/dE from QCDark2 rate CSVs for band-gap pheno Step 1."""

from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


# ----------------------------------------------------------------------------
# RateCsvMeta
#   Metadata parsed from the comment header of a rate CSV: material, mediator, mass, cross section and the ionization table it was made with.
# ----------------------------------------------------------------------------
@dataclass(frozen=True)
class RateCsvMeta:
    path: Path
    material: str
    mediator: str
    mchi_MeV: float
    sigma_e_cm2: float
    table_path: str

    # ----------------------------------------------------------------------------
    # RateCsvMeta.from_csv
    #   Parse the '#' header lines of a rate CSV into a RateCsvMeta (missing entries stay empty or NaN).
    # ----------------------------------------------------------------------------
    @classmethod
    def from_csv(cls, path: Path) -> RateCsvMeta:
        material = mediator = table_path = ""
        mX_eV = sigma_e = float("nan")
        with open(path, encoding="utf-8-sig") as f:
            for line in f:
                if not line.startswith("#"):
                    break
                raw = line[1:].strip()
                for chunk in raw.split(","):
                    chunk = chunk.strip()
                    if chunk.startswith("material ="):
                        material = chunk.split("=", 1)[1].strip()
                    elif chunk.startswith("mediator ="):
                        mediator = chunk.split("=", 1)[1].strip()
                    elif chunk.startswith("table ="):
                        table_path = chunk.split("=", 1)[1].strip()
                if raw.startswith("mX (eV)"):
                    mX_eV = float(raw.split("=", 1)[1].strip())
                elif raw.startswith("sigma_e (cm^2)"):
                    sigma_e = float(raw.split("=", 1)[1].strip())
        mchi_MeV = mX_eV / 1e6 if math.isfinite(mX_eV) else float("nan")
        if not math.isfinite(mchi_MeV) or not math.isfinite(sigma_e):
            parsed = _parse_rate_path(path)
            if parsed is None:
                raise ValueError(f"could not parse metadata from {path}")
            mchi_MeV, sig_s = parsed
            sigma_e = float(sig_s.replace("p", "."))
        return cls(
            path=path,
            material=material or "Si",
            mediator=mediator or "heavy",
            mchi_MeV=mchi_MeV,
            sigma_e_cm2=sigma_e,
            table_path=table_path,
        )


# ----------------------------------------------------------------------------
# read_dRdE_csv
#   Read the (E, dR/dE) columns of a rate CSV, skipping comments and unparsable rows.
# ----------------------------------------------------------------------------
def read_dRdE_csv(path: Path) -> tuple[np.ndarray, np.ndarray]:
    E, R = [], []
    with open(path, encoding="utf-8-sig") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = [p.strip() for p in line.split(",")]
            if len(parts) < 2:
                continue
            try:
                E.append(float(parts[0]))
                R.append(float(parts[1]))
            except ValueError:
                continue
    return np.asarray(E), np.asarray(R)


# ----------------------------------------------------------------------------
# RateFileMatch
#   A rate file found for the requested point: its path, mass and cross-section string.
# ----------------------------------------------------------------------------
@dataclass(frozen=True)
class RateFileMatch:
    path: Path
    mchi_MeV: float
    sigma_str: str


_RATE_RE = re.compile(r"_m([0-9.]+)_s(.+)\.csv$", re.IGNORECASE)


# ----------------------------------------------------------------------------
# _parse_rate_path
#   Mass and cross-section string parsed from a rate file name, or None if the name does not match.
# ----------------------------------------------------------------------------
def _parse_rate_path(path: Path) -> tuple[float, str] | None:
    m = _RATE_RE.search(path.name)
    if not m:
        return None
    return float(m.group(1)), m.group(2)


def find_rate_file(rates_dir: Path, mchi_MeV: float, sigma: str) -> RateFileMatch:
    """Pick closest (m_chi, sigma_e) on the generated rate grid."""
    sigma_target = float(sigma)
    log_sig_t = math.log10(sigma_target)
    candidates = list(rates_dir.glob("dRdE_*.csv"))
    if not candidates:
        raise FileNotFoundError(f"no rate files in {rates_dir}")

    best_path: Path | None = None
    best_score = float("inf")
    best_mchi = mchi_MeV
    best_sig = sigma

    for p in candidates:
        parsed = _parse_rate_path(p)
        if parsed is None:
            continue
        mchi_f, sig_s = parsed
        dm_rel = abs(mchi_f - mchi_MeV) / max(mchi_MeV, 1e-6)
        try:
            dlog = abs(math.log10(float(sig_s.replace("p", "."))) - log_sig_t)
        except ValueError:
            dlog = 1e6
        score = dm_rel**2 + dlog**2
        if score < best_score:
            best_score = score
            best_path = p
            best_mchi = mchi_f
            best_sig = sig_s

    if best_path is None:
        raise FileNotFoundError(f"no parseable rate file in {rates_dir}")

    if best_score > 1e-6:
        print(
            f"[warn] {rates_dir.name}: requested mchi={mchi_MeV} MeV, sigma={sigma}; "
            f"using {best_path.name}"
        )
    return RateFileMatch(best_path, best_mchi, best_sig)


def gap_ev_from_tag(tag: str) -> float:
    """gap0p1 -> 0.1 eV (matches QCDark2 table tag)."""
    return float(tag.replace("gap", "").replace("p", "."))


# ----------------------------------------------------------------------------
# gap_legend_label
#   Legend label with the band gap of a curve and the name of the table it used.
# ----------------------------------------------------------------------------
def gap_legend_label(tag: str, meta: RateCsvMeta) -> str:
    gap_ev = gap_ev_from_tag(tag)
    table = Path(meta.table_path).name if meta.table_path else tag
    return rf"$E_{{\mathrm{{gap}}}} = {gap_ev:g}$ eV; {table}"


# ----------------------------------------------------------------------------
# format_mchi_latex
#   Mass as a matplotlib math-text string.
# ----------------------------------------------------------------------------
def format_mchi_latex(mchi_MeV: float) -> str:
    return rf"$m_\chi = {mchi_MeV:.9f}$ MeV"


# ----------------------------------------------------------------------------
# format_sigma_latex
#   Cross section as a matplotlib math-text string.
# ----------------------------------------------------------------------------
def format_sigma_latex(sigma_cm2: float) -> str:
    return rf"$\sigma_e = {sigma_cm2:.16e}$ cm$^2$"


# ----------------------------------------------------------------------------
# format_run_title
#   Plot title with the mass, cross section, material and mediator of a run.
# ----------------------------------------------------------------------------
def format_run_title(meta: RateCsvMeta) -> str:
    return (
        f"{format_mchi_latex(meta.mchi_MeV)}, "
        f"{format_sigma_latex(meta.sigma_e_cm2)} "
        f"({meta.material}, {meta.mediator} mediator; QCDark2)"
    )


# ----------------------------------------------------------------------------
# main
#   Overlay dR/dE of several band-gap curves (--curve LABEL=RATES_DIR) at one mass and cross section and save the figure to --out.
# ----------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mchi-MeV", type=float, default=10.0)
    ap.add_argument("--sigma", default="1.0e-40")
    ap.add_argument(
        "--curve",
        nargs="+",
        required=True,
        metavar="LABEL=RATES_DIR",
        help='e.g. "gap0p1=data/qcdark2_rates/.../Si_fast_gap0p1"',
    )
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--Emin", type=float, default=0.0)
    ap.add_argument("--Emax", type=float, default=10.0)
    args = ap.parse_args()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(9, 5))

    title_meta: RateCsvMeta | None = None

    for spec in args.curve:
        if "=" not in spec:
            ap.error(f"expected LABEL=dir, got {spec!r}")
        tag, rdir = spec.split("=", 1)
        match = find_rate_file(Path(rdir), args.mchi_MeV, args.sigma)
        meta = RateCsvMeta.from_csv(match.path)
        if title_meta is None:
            title_meta = meta
        elif (
            abs(meta.mchi_MeV - title_meta.mchi_MeV) > 1e-9
            or abs(meta.sigma_e_cm2 - title_meta.sigma_e_cm2) > 0
        ):
            print(
                f"[warn] {tag} uses different (mchi, sigma) than first curve: "
                f"{meta.path.name}"
            )
        E, R = read_dRdE_csv(match.path)
        mask = (E >= args.Emin) & (E <= args.Emax)
        pos = R[mask] > 0
        ax.plot(E[mask][pos], R[mask][pos], lw=1.5, label=gap_legend_label(tag, meta))
        print(
            f"[ok] {tag}: {match.path.name}  "
            f"mchi={meta.mchi_MeV:.9f} MeV  sigma_e={meta.sigma_e_cm2:.16e} cm^2  "
            f"table={Path(meta.table_path).name}"
        )

    ax.set_xlabel(r"Recoil energy $E$ [eV]")
    ax.set_ylabel(r"$dR/dE$ [events/(kg·year·eV)]")
    ax.set_xlim(args.Emin, args.Emax)
    ax.set_yscale("log")
    ax.legend(loc="best", fontsize=7)
    ax.grid(True, alpha=0.3, which="both")
    if title_meta is not None:
        ax.set_title(format_run_title(title_meta), fontsize=10)
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    plt.close(fig)
    print(f"[ok] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
