#!/usr/bin/env python3
"""
Extend the reference p100K table (Si, E_gap=1.2 eV, eh=3.8 eV) from 50 eV to
a user-specified maximum energy using a Gaussian approximation calibrated on the
existing high-energy rows.

Physics rationale
-----------------
At large n_e (many electron-hole pairs), the central limit theorem guarantees the
ionization distribution converges to a Gaussian:

    P(n_e | E) ~ Normal(mu(E), sigma(E))
    mu(E)    = E / eh_eV - offset   (threshold offset ~0.69 e- stabilized from data)
    sigma(E) = sqrt(Fano * mu(E))   (Fano factor ~0.122 from last unsaturated rows)

Both Fano and offset are fit from the last ~10 rows of the existing table where
the distribution is fully converged (not clipped by the column limit). The
extension is stitched directly onto the existing table with the n_e column range
extended to cover ceil(E_max / eh_eV) + a few sigma.

The extended table can then be fed to the existing scaling script to regenerate
p100K_gap0p34_eh1p6.csv (or any other material) with the full energy range.

Usage
-----
  python3 utils/extend_p100K_reference.py \\
      --input  data/p100K_table.csv \\
      --output data/p100K_table_extended.csv \\
      --E_max  200.0 \\
      --eh_eV  3.8
"""
from __future__ import annotations

import argparse
import math
import numpy as np
from pathlib import Path


def read_table(path: str) -> tuple[np.ndarray, np.ndarray]:
    data = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            data.append([float(x) for x in line.split(",")])
    arr = np.array(data)
    return arr[:, 0], arr[:, 1:]   # (E, P)


def fit_gaussian_params(E: np.ndarray, P: np.ndarray,
                        eh_eV: float, n_fit_rows: int = 15
                        ) -> tuple[float, float]:
    """
    Fit Fano factor and mean offset from the last n_fit_rows of the existing
    table, where the Gaussian approximation is fully valid.
    Returns (fano, offset) where mu(E) = E/eh_eV - offset.
    """
    ne_vals = np.arange(P.shape[1])
    fanos, offsets = [], []
    for i in range(-n_fit_rows, 0):
        row = P[i]
        s   = row.sum()
        if s < 0.99:   # skip rows where distribution is clipped
            continue
        mean = (ne_vals * row).sum() / s
        var  = (ne_vals**2 * row).sum() / s - mean**2
        n_mean = E[i] / eh_eV
        fanos.append(var / mean if mean > 0 else np.nan)
        offsets.append(n_mean - mean)
    fano   = float(np.nanmedian(fanos))
    offset = float(np.nanmedian(offsets))
    print(f"[extend] Fitted Fano  = {fano:.4f}  (from {len(fanos)} rows)")
    print(f"[extend] Fitted offset = {offset:.4f} e-")
    return fano, offset


def gaussian_row(E_val: float, eh_eV: float, fano: float,
                 offset: float, ne_max: int) -> np.ndarray:
    """
    Gaussian P(n_e | E) for n_e in 0..ne_max-1, normalised to sum=1.
    """
    mu  = E_val / eh_eV - offset
    sig = math.sqrt(fano * mu) if mu > 0 else 1e-6
    ne  = np.arange(ne_max, dtype=float)
    row = np.exp(-0.5 * ((ne - mu) / sig)**2)
    row[row < 0] = 0.0
    s = row.sum()
    if s > 0:
        row /= s
    return row


def extend_table(input_path: str, output_path: str,
                 E_max: float, eh_eV: float, dE: float = 0.05,
                 n_sigma_extra: int = 5):
    E, P = read_table(input_path)
    E_last = E[-1]
    ne_orig = P.shape[1]

    if E_max <= E_last:
        print(f"[extend] E_max={E_max} <= table max {E_last}; nothing to do.")
        return

    fano, offset = fit_gaussian_params(E, P, eh_eV)

    # Determine new ne_max: cover E_max/eh_eV + n_sigma_extra * sigma
    mu_max  = E_max / eh_eV - offset
    sig_max = math.sqrt(fano * mu_max)
    ne_max  = int(math.ceil(mu_max + n_sigma_extra * sig_max)) + 1
    ne_max  = max(ne_max, ne_orig)
    print(f"[extend] Extending n_e columns: {ne_orig} -> {ne_max}")
    print(f"[extend] Extending E rows: {E_last:.2f} -> {E_max:.2f} eV")

    # Pad existing rows with zeros for new n_e columns
    P_pad = np.hstack([P, np.zeros((len(E), ne_max - ne_orig))])

    # Build new energy rows
    E_new   = np.arange(E_last + dE, E_max + dE/2, dE)
    P_new   = np.zeros((len(E_new), ne_max))
    for i, e in enumerate(E_new):
        P_new[i] = gaussian_row(e, eh_eV, fano, offset, ne_max)

    # Stitch
    E_out = np.concatenate([E, E_new])
    P_out = np.concatenate([P_pad, P_new], axis=0)

    # Write
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        f.write(f"# Extended p100K reference table\n")
        f.write(f"# Original: {input_path}  (E up to {E_last:.2f} eV, {ne_orig} n_e columns)\n")
        f.write(f"# Extension: Gaussian approx with Fano={fano:.4f}, offset={offset:.4f} e-\n")
        f.write(f"# E_gap_ref=1.2 eV  eh_ref={eh_eV} eV\n")
        f.write(f"# Extended to E_max={E_max:.1f} eV, ne_max={ne_max}\n")
        for e_val, row in zip(E_out, P_out):
            vals = [f"{e_val:.2f}"] + [f"{v:.5f}" for v in row]
            f.write(",".join(vals) + "\n")

    print(f"[extend] Wrote {out}  ({len(E_out)} rows x {ne_max+1} columns)")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input",  default="data/p100K_table.csv")
    ap.add_argument("--output", default="data/p100K_table_extended.csv")
    ap.add_argument("--E_max",  type=float, default=200.0,
                    help="Maximum energy to extend to (eV). Default 200.")
    ap.add_argument("--eh_eV",  type=float, default=3.8,
                    help="eh pair energy of reference table (eV). Default 3.8.")
    ap.add_argument("--dE",     type=float, default=0.05,
                    help="Energy step for new rows (eV). Default 0.05.")
    args = ap.parse_args()
    extend_table(args.input, args.output, args.E_max, args.eh_eV, args.dE)


if __name__ == "__main__":
    main()
