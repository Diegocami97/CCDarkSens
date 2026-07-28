#!/usr/bin/env python3
"""
Build a lower-envelope CSV for the DM-electron direct detection limits.

Reads all relevant experiment curves for heavy and light mediator cases,
computes the minimum sigma at each sampled mass (log-spaced), and writes:
  data/previous_limits/heavy_mediator/direct_detection_envelope.csv
  data/previous_limits/light_mediator/direct_detection_envelope.csv

The DAMIC-M 2025 line is included in the envelope (it defines the solid black
line AND the lower boundary of the shaded region).

Usage:
  python3 utils/build_dme_direct_detection_envelope.py
"""
from __future__ import annotations

import csv
import math
import numpy as np
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "data" / "previous_limits"

# DAMIC-M 2025 paper-export files (masses in eV — converted to MeV on load)
DAMIC_2025_BASE = (REPO / "collab_frameworks/pydme/analysis/DailyModulation"
                   "/LBC-Sep2024/paper_figures/data/LBC_results"
                   "/ScienceRun2024_results-Pattern")

# Mass range matching the DM-e plot window (MeV)
MASS_MIN = 1e-1
MASS_MAX = 1e3
N_SAMPLES = 2000


def read_curve(path: Path, mass_scale: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
    """Read a two-column file, skipping comment/header lines.
    mass_scale converts the stored mass unit to MeV (e.g. 1e-6 if file is in eV)."""
    masses, sigmas = [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(",") if "," in line else line.split()
            if len(parts) < 2:
                continue
            try:
                m = float(parts[0].strip()) * mass_scale
                s = float(parts[1].strip())
            except ValueError:
                continue
            if m > 0 and s > 0:
                masses.append(m)
                sigmas.append(s)
    return np.array(masses), np.array(sigmas)


def lower_envelope(curves: list[tuple[np.ndarray, np.ndarray]],
                   mass_min: float, mass_max: float,
                   n: int = N_SAMPLES) -> tuple[np.ndarray, np.ndarray]:
    """Compute lower envelope using log-log interpolation to avoid linear-space
    artifacts on curves that span many decades."""
    log_masses = np.linspace(math.log10(mass_min), math.log10(mass_max), n)
    masses_out, sigmas_out = [], []
    # Pre-compute log versions of each curve for log-log interpolation
    log_curves = []
    for ms, ss in curves:
        log_curves.append((np.log10(ms), np.log10(ss)))
    for logm in log_masses:
        m = 10.0 ** logm
        y_min = np.inf
        for (ms, ss), (lms, lss) in zip(curves, log_curves):
            if m < ms.min() or m > ms.max():
                continue
            log_y = np.interp(logm, lms, lss)
            y = 10.0 ** log_y
            if y > 0 and y < y_min:
                y_min = y
        if y_min < np.inf:
            masses_out.append(m)
            sigmas_out.append(y_min)
    return np.array(masses_out), np.array(sigmas_out)


def write_csv(path: Path, masses: np.ndarray, sigmas: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write("# Direct detection lower envelope (DM-electron)\n")
        f.write("# mass_MeV, sigma_cm2\n")
        for m, s in zip(masses, sigmas):
            f.write(f"{m:.6e},{s:.6e}\n")
    print(f"Wrote {path}  ({len(masses)} points)")


def build(mediator: str) -> None:
    d = DATA / mediator
    damic_key = ("DAMIC-M_2025_QEDark_DMe_heavymediator.txt"
                 if mediator == "heavy_mediator"
                 else "DAMIC-M_2025_QEDark_DMe_ulightmediator.txt")
    damic_collab = DAMIC_2025_BASE / damic_key  # masses in eV

    curves = []

    # DAMIC-M 2025 from the collab framework file (same source as the solid black line)
    if damic_collab.exists():
        ms, ss = read_curve(damic_collab, mass_scale=1e-6)  # eV → MeV
        print(f"  loaded {damic_collab.name} (collab)  ({len(ms)} pts, mass {ms.min():.3g}–{ms.max():.3g} MeV)")
        curves.append((ms, ss))
    else:
        print(f"  [skip] collab DAMIC-M 2025 not found: {damic_collab}")

    # Heavy mediator only: add Panda4T for high-mass coverage
    if mediator == "heavy_mediator":
        p = d / "Panda4T.csv"
        if p.exists():
            ms, ss = read_curve(p)
            print(f"  loaded {p.name}  ({len(ms)} pts, mass {ms.min():.3g}–{ms.max():.3g} MeV)")
            curves.append((ms, ss))

    masses, sigmas = lower_envelope(curves, MASS_MIN, MASS_MAX)
    out = d / "direct_detection_envelope.csv"
    write_csv(out, masses, sigmas)


if __name__ == "__main__":
    for med in ("heavy_mediator", "light_mediator"):
        print(f"\n=== {med} ===")
        build(med)
