#!/usr/bin/env python3
import sys, glob
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import re

def smart_read_two_col_csv(path):
    """
    Robustly read a CCDarkSens QEDark CSV:
      - Skips all comment lines (# ...)
      - Detects delimiter (comma or whitespace)
      - Assumes the first non-comment line is a header with 2 names
      - Returns (E, R) as float arrays
    """
    with open(path, "r", encoding="utf-8-sig") as f:
        lines = f.readlines()

    # Strip and drop empty lines
    lines = [ln.strip() for ln in lines if ln.strip()]

    # Keep only non-comment lines
    data_lines = [ln for ln in lines if not ln.lstrip().startswith("#")]
    if not data_lines:
        raise ValueError("No data lines found (only comments?)")

    # First non-comment is the header row
    header = data_lines[0]
    rest   = data_lines[1:]

    # Detect delimiter from header
    if "," in header:
        delim = ","
        cols = [c.strip() for c in header.split(",")]
    else:
        # whitespace-split fallback
        delim = None
        cols = header.split()

    if len(cols) < 2:
        raise ValueError(f"Header has <2 columns: {header!r}")

    # Parse numeric rows
    E_vals, R_vals = [], []
    for ln in rest:
        if delim is None:
            parts = ln.split()
        else:
            parts = [p.strip() for p in ln.split(delim)]
        # tolerate extra columns by taking the first two numeric ones
        if len(parts) < 2:
            # skip malformed lines quietly
            continue
        try:
            E = float(parts[0])
            R = float(parts[1])
        except ValueError:
            # skip lines that aren't numeric (defensive)
            continue
        E_vals.append(E)
        R_vals.append(R)

    if not E_vals:
        raise ValueError("No numeric data parsed from file.")

    return np.asarray(E_vals, float), np.asarray(R_vals, float), cols[0], cols[1]

def main(mediator: str):
    base_dir = Path(f"data/qedark_rates/Si/{mediator}")
    if not base_dir.exists():
        print(f"Error: directory {base_dir} not found.")
        sys.exit(1)

    # pick all matching CSVs
    csv_files = sorted(glob.glob(str(base_dir / "dRdE_*.csv")))
    if not csv_files:
        print(f"No CSV files found in {base_dir}")
        sys.exit(1)

    plt.figure(figsize=(7.5, 5.2))
    plotted = 0
    for fpath in csv_files:
        try:
            E, R, Ename, Rname = smart_read_two_col_csv(fpath)
        except Exception as e:
            print(f"[warn] failed to read {fpath}: {e}")
            continue

        # Build a compact label from filename
        # fname = Path(fpath).name
        # parts = fname.replace(".csv", "").split("_")
        # mpart = next((p for p in parts if p.startswith("m")), None)
        # spart = next((p for p in parts if p.startswith("s")), None)
        # label = fname
        # if mpart and spart:
        #     label = f"mχ={mpart[1:]} MeV, σₑ={spart[1:]} cm²"

        fname = Path(fpath).name
        # Match e.g. dRdE_Si_heavy_m30.000000_s1e-37.csv
        m = re.search(
            r"_m(?P<mchi>[0-9.eE+-]+)_s(?P<sigma>[0-9.eE+-]+)\.csv$", fname
        )
        if m:
            mchi = m.group("mchi")
            sigma = m.group("sigma")
            label = f"mχ={mchi} MeV, σₑ={sigma} cm²"
        else:
            label = fname  # fallback

        plt.plot(E, R, lw=1.6, label=label)
        plotted += 1

    if plotted == 0:
        print("No valid curves were plotted.")
        sys.exit(1)

    plt.xlabel("Recoil energy $E_e$ [eV]", fontsize=13)
    plt.ylabel(r"$dR/dE$ [events / kg / year / eV]", fontsize=13)
    plt.title(f"QEDark rates: Si ({mediator} mediator)", fontsize=14)
    plt.xlim(0, 20)
    plt.yscale("log")
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend(fontsize=8, ncol=1)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: python3 utils/plot_qedark_rates.py <mediator>")
        print("example: python3 utils/plot_qedark_rates.py heavy")
        print("example: python3 utils/plot_qedark_rates.py massless")
        sys.exit(1)
    main(sys.argv[1])
