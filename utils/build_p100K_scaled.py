#!/usr/bin/env python3
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: build_p100K_scaled.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  build_p100K_scaled.py -- Build pheno-scaled P(n_e|E) tables from the Si
#  reference p100K_table.csv using the anchored energy map (see
#  docs/band_gap_pheno_p100K_scaling_explained.md).
# ============================================================================

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import List, Tuple

import numpy as np

REF_GAP_EV = 1.2
REF_EH_EV = 3.8
DEFAULT_REF = "data/p100K_table.csv"
DEFAULT_E_MIN_EV = 0.05


# ----------------------------------------------------------------------------
# infer_energy_step
#   Energy step of a grid: E[1] - E[0] rounded to 6 decimals (0.05 eV for a grid with fewer than 2 points).
# ----------------------------------------------------------------------------
def infer_energy_step(E: np.ndarray) -> float:
    if E.size >= 2:
        return float(np.round(E[1] - E[0], 6))
    return 0.05


def extend_energy_grid(
    E_ref: np.ndarray,
    E_min: float = DEFAULT_E_MIN_EV,
    step: float | None = None,
) -> np.ndarray:
    """Prepend a uniform low-E grid [E_min, …) matching ref step, then merge with E_ref."""
    if step is None:
        step = infer_energy_step(E_ref)
    if E_min >= E_ref[0]:
        return E_ref.copy()
    # np.arange upper bound exclusive → stop just below first ref point
    E_low = np.arange(E_min, E_ref[0] - 0.25 * step, step)
    return np.unique(np.concatenate([E_low, E_ref]))


def fmt_ev_tag(x: float) -> str:
    """0.1 -> 0p1, 3.8 -> 3p8 for filenames."""
    s = f"{x:g}"
    return s.replace(".", "p")


# ----------------------------------------------------------------------------
# default_out_path
#   Default output file data/p100K_gap<gap>_eh<eh>.csv for a (band gap, e-h pair energy) pair.
# ----------------------------------------------------------------------------
def default_out_path(gap_eV: float, eh_eV: float) -> Path:
    return Path(f"data/p100K_gap{fmt_ev_tag(gap_eV)}_eh{fmt_ev_tag(eh_eV)}.csv")


def load_p100k_csv(path: Path) -> Tuple[np.ndarray, np.ndarray]:
    """Return E [eV], P [n_ne, n_E] with columns P(n_e=1..N|E)."""
    rows: List[List[float]] = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = [p.strip() for p in line.split(",")]
            vals = [float(p) for p in parts if p]
            if len(vals) >= 2:
                rows.append(vals)
    if len(rows) < 2:
        raise ValueError(f"malformed p100K table: {path}")
    arr = np.asarray(rows, dtype=float)
    E = arr[:, 0]
    P = arr[:, 1:].T
    return E, P


# ----------------------------------------------------------------------------
# interp_linear_clamped
#   Linear interpolation of y(x) at xq, with the end values held outside the range and the result clipped to [0,1] because it is a probability.
# ----------------------------------------------------------------------------
def interp_linear_clamped(x: np.ndarray, y: np.ndarray, xq: float) -> float:
    if xq <= x[0]:
        return float(np.clip(y[0], 0.0, 1.0))
    if xq >= x[-1]:
        return float(np.clip(y[-1], 0.0, 1.0))
    i1 = int(np.searchsorted(x, xq, side="right"))
    i0 = i1 - 1
    t = (xq - x[i0]) / (x[i1] - x[i0])
    v = (1.0 - t) * y[i0] + t * y[i1]
    return float(np.clip(v, 0.0, 1.0))


# ----------------------------------------------------------------------------
# map_E_prime
#   Anchored energy map: the reference-table energy that corresponds to E for a new band gap and pair energy, E' = gap_ref + (E - gap_new)*(eh_ref/eh_new).
# ----------------------------------------------------------------------------
def map_E_prime(
    E: float,
    gap_new: float,
    gap_ref: float,
    eh_new: float,
    eh_ref: float,
) -> float:
    if eh_new <= 0:
        raise ValueError("eh_pair_eV must be > 0")
    return gap_ref + (E - gap_new) * (eh_ref / eh_new)


# ----------------------------------------------------------------------------
# build_scaled_table
#   P(n_e | E) of the new (gap, eh) pair: for each energy above the new gap I map E to E' and interpolate the reference probabilities there; below the gap the probabilities are zero.
# ----------------------------------------------------------------------------
def build_scaled_table(
    E_grid: np.ndarray,
    P_ref: np.ndarray,
    E_ref: np.ndarray,
    gap_new: float,
    eh_new: float,
    gap_ref: float = REF_GAP_EV,
    eh_ref: float = REF_EH_EV,
) -> np.ndarray:
    n_ne = P_ref.shape[0]
    P_new = np.zeros((n_ne, E_grid.size), dtype=float)
    for i, E in enumerate(E_grid):
        if E < gap_new:
            continue
        Ep = map_E_prime(E, gap_new, gap_ref, eh_new, eh_ref)
        for n in range(n_ne):
            P_new[n, i] = interp_linear_clamped(E_ref, P_ref[n], Ep)
    return P_new


# ----------------------------------------------------------------------------
# write_p100k_csv
#   Write a scaled ionization table as CSV with a header that records the reference file, gaps and pair energies, and the scenario.
# ----------------------------------------------------------------------------
def write_p100k_csv(
    out: Path,
    E: np.ndarray,
    P: np.ndarray,
    *,
    gap_new: float,
    eh_new: float,
    gap_ref: float,
    eh_ref: float,
    scenario: str = "",
    ref_path: str = DEFAULT_REF,
) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8") as f:
        f.write(
            f"# Scaled p100K (anchored pheno map)\n"
            f"# ref: {ref_path}  E_gap_ref={gap_ref} eV  eh_ref={eh_ref} eV\n"
            f"# new: E_gap_new={gap_new} eV  eh_new={eh_new} eV"
        )
        if scenario:
            f.write(f"  scenario={scenario}")
        f.write(
            f"\n# E grid extended to E_min={E[0]:g} eV (0.05 eV steps; ref lookup starts at 1.1 eV)\n"
            f"# E'(E) = {gap_ref} + (E - {gap_new}) * ({eh_ref}/{eh_new})\n"
        )
        for i in range(E.size):
            row = [f"{E[i]:g}"] + [f"{P[n, i]:g}" for n in range(P.shape[0])]
            f.write(",".join(row) + "\n")
    print(f"[ok] wrote {out}  ({E.size} rows, {P.shape[0]} n_e columns)")


# ----------------------------------------------------------------------------
# plot_p100k_compare
#   Plot P(n_e = ne_plot | E) of the reference table and of the scaled tables in one PDF.
# ----------------------------------------------------------------------------
def plot_p100k_compare(
    E_ref: np.ndarray,
    P_ref: np.ndarray,
    tables: List[Tuple[str, np.ndarray, np.ndarray]],
    out_pdf: Path,
    ne_plot: int = 1,
) -> None:
    import matplotlib.pyplot as plt

    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 5))
    idx = ne_plot - 1
    ax.plot(E_ref, P_ref[idx], "k-", lw=2, label="reference (Si p100K)")
    for label, E, P in tables:
        ax.plot(E, P[idx], lw=1.5, label=label)
    ax.set_xlabel(r"Recoil energy $E$ [eV]")
    ax.set_ylabel(rf"$P(n_e={ne_plot}\mid E)$")
    ax.set_xlim(left=0)
    ax.set_ylim(-0.02, 1.05)
    ax.legend(loc="best", fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)
    print(f"[ok] plot {out_pdf}")


# ----------------------------------------------------------------------------
# build_from_manifest
#   Build every scaled table listed in the scenario manifest JSON (reference and scenarios), extending the energy grid down to E_min, skipping existing files unless force is set, with an optional comparison plot.
# ----------------------------------------------------------------------------
def build_from_manifest(
    manifest_path: Path,
    plot: bool,
    E_min: float,
    force: bool,
) -> None:
    with open(manifest_path, encoding="utf-8") as f:
        man = json.load(f)
    entries = []
    if "reference" in man:
        entries.append(man["reference"])
    entries.extend(man.get("scenarios", []))
    ref_csv = Path(man.get("reference", {}).get("ionization_csv", DEFAULT_REF))
    E_ref_raw, P_ref = load_p100k_csv(ref_csv if ref_csv.exists() else Path(DEFAULT_REF))
    E_grid = extend_energy_grid(E_ref_raw, E_min=E_min)
    print(f"[grid] E_min={E_min} eV  step={infer_energy_step(E_ref_raw)} eV  "
          f"rows={E_grid.size} (was {E_ref_raw.size})")
    built: List[Tuple[str, np.ndarray, np.ndarray]] = []
    for ent in entries:
        csv = ent.get("ionization_csv")
        if not csv or csv == str(ref_csv) and ent.get("scenario") == "ref":
            continue
        gap = float(ent["band_gap_eV"])
        eh = float(ent["eh_pair_eV"])
        out = Path(csv)
        if out.exists() and not force:
            print(f"[skip] exists {out} (use --force to overwrite)")
            continue
        P_new = build_scaled_table(E_grid, P_ref, E_ref_raw, gap, eh)
        write_p100k_csv(
            out,
            E_grid,
            P_new,
            gap_new=gap,
            eh_new=eh,
            gap_ref=REF_GAP_EV,
            eh_ref=REF_EH_EV,
            scenario=ent.get("scenario", ""),
        )
        built.append((ent.get("label", out.name), E_grid, P_new))
    if plot and built:
        plot_p100k_compare(
            E_ref_raw,
            P_ref,
            built,
            Path("outplots/band_gap_pheno/step2_p100K/manifest_all_Pne1.pdf"),
        )


# ----------------------------------------------------------------------------
# main
#   Command line: build one scaled table from --band-gap-eV and --eh-pair-eV, or all tables of a manifest (--from-manifest); optional comparison plot (--plot).
# ----------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser(description="Build scaled p100K ionization CSV tables.")
    ap.add_argument("--band-gap-eV", type=float, help="New pheno band gap [eV]")
    ap.add_argument("--eh-pair-eV", type=float, help="New pheno e-h pair scale [eV]")
    ap.add_argument("--ref-csv", default=DEFAULT_REF, help="Reference p100K CSV")
    ap.add_argument("--ref-gap-eV", type=float, default=REF_GAP_EV)
    ap.add_argument("--ref-eh-eV", type=float, default=REF_EH_EV)
    ap.add_argument("--out", type=Path, help="Output CSV (default: data/p100K_gap*_eh*.csv)")
    ap.add_argument("--scenario", default="", help="Scenario tag for header comment")
    ap.add_argument("--plot", action="store_true", help="Write comparison PDF")
    ap.add_argument(
        "--from-manifest",
        type=Path,
        metavar="JSON",
        help="Build all ionization_csv paths in band_gap_pheno_scenarios.json",
    )
    ap.add_argument(
        "--E-min",
        type=float,
        default=DEFAULT_E_MIN_EV,
        help="Extend output energy grid down to this value [eV] (default: 0.05)",
    )
    ap.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing output CSVs",
    )
    args = ap.parse_args()

    if args.from_manifest:
        build_from_manifest(args.from_manifest, plot=args.plot, E_min=args.E_min, force=args.force)
        return

    if args.band_gap_eV is None or args.eh_pair_eV is None:
        ap.error("provide --band-gap-eV and --eh-pair-eV, or use --from-manifest")

    ref_path = Path(args.ref_csv)
    E_ref_raw, P_ref = load_p100k_csv(ref_path)
    E_grid = extend_energy_grid(E_ref_raw, E_min=args.E_min)
    gap_new = args.band_gap_eV
    eh_new = args.eh_pair_eV
    P_new = build_scaled_table(
        E_grid, P_ref, E_ref_raw, gap_new, eh_new, args.ref_gap_eV, args.ref_eh_eV
    )
    out = args.out or default_out_path(gap_new, eh_new)
    write_p100k_csv(
        out,
        E_grid,
        P_new,
        gap_new=gap_new,
        eh_new=eh_new,
        gap_ref=args.ref_gap_eV,
        eh_ref=args.ref_eh_eV,
        scenario=args.scenario,
        ref_path=str(ref_path),
    )

    if args.plot:
        label = args.scenario or f"gap={gap_new} eh={eh_new}"
        plot_pdf = Path(
            f"outplots/band_gap_pheno/step2_p100K/Pne_gap{fmt_ev_tag(gap_new)}_eh{fmt_ev_tag(eh_new)}.pdf"
        )
        plot_p100k_compare(
            E_ref_raw,
            P_ref,
            [(label, E_grid, P_new)],
            plot_pdf,
        )


if __name__ == "__main__":
    main()
