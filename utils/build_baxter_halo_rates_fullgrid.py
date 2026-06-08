#!/usr/bin/env python3
# ============================================================================
#  CCDarkSens — build_baxter_halo_rates_fullgrid
#  Build full-grid Baxter halo dR/dE library by per-mass scaling
#
#  Author: Diego Venegas-Vargas
# ============================================================================
"""Build full-grid Baxter halo rate library by per-mass scaling of long_scan."""
from __future__ import annotations

import argparse
import json
import re
import shutil
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

SRC = ROOT / "data/qedark_rates/Si/heavy/long_scan"
DST = ROOT / "data/qedark_rates/Si/heavy/long_scan_baxter_vE253p7"
RATIO_CACHE = ROOT / "data/qedark_rates/Si/heavy/baxter_vE253p7_mass_ratios.json"

V0 = np.array([0.0, 238.0, 0.0])
V_SUN = np.array([11.1, 12.2, 7.3])
V_EARTH_MAR9 = np.array([29.2, -0.1, 5.9])
BAXTER_VE = float(np.linalg.norm(V0 + V_SUN + V_EARTH_MAR9))

BASE = {"v0_kms": 238.0, "vE_kms": 263.0, "vesc_kms": 544.0}
BAXTER = {"v0_kms": 238.0, "vE_kms": BAXTER_VE, "vesc_kms": 544.0}


def integrated_rate(m_mev: float, halo: dict) -> float:
    from ccdarkphys.qedark.entry import compute_dRdE

    out = compute_dRdE("Si", "heavy", m_mev * 1e6, 1e-36, halo, binsize_eV=0.1)
    E = out["E_eV"]
    mask = E >= 2.0
    return float(np.trapz(out["dRdE_kg_year_eV"][mask], E[mask]) / 365.25 / 1000)


def mass_ratios(masses: list[float]) -> dict[str, float]:
    out = {}
    for m in masses:
        rb = integrated_rate(m, BASE)
        rx = integrated_rate(m, BAXTER)
        out[f"{m:.6f}"] = rx / rb if rb > 0 else 1.0
    return out


def list_masses(src: Path) -> list[float]:
    pat = re.compile(r"_m([0-9.]+)_s")
    return sorted({float(pat.search(p.name).group(1)) for p in src.glob("dRdE_Si_heavy_m*_s*.csv") if pat.search(p.name)})


def scale_one(src: Path, dst: Path, scale: float) -> None:
    lines = []
    for ln in src.read_text().splitlines():
        if ln.startswith("#") or not ln.strip():
            lines.append(ln)
            continue
        parts = ln.split(",")
        if len(parts) < 2:
            lines.append(ln)
            continue
        try:
            lines.append(f"{parts[0]},{float(parts[1]) * scale:.16e}")
        except ValueError:
            lines.append(ln)
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()

    if DST.exists() and not args.force:
        n = len(list(DST.glob("*.csv")))
        if n > 1000:
            print(f"skip build: {DST} already has {n} files (use --force)")
            return

    masses = list_masses(SRC)
    print(f"Baxter v_E = {BAXTER_VE:.2f} km/s; scaling {len(masses)} masses from {SRC}")

    if RATIO_CACHE.exists() and not args.force:
        ratios = json.loads(RATIO_CACHE.read_text())
    else:
        print("computing per-mass rate ratios (this may take a few minutes)...")
        ratios = mass_ratios(masses)
        RATIO_CACHE.write_text(json.dumps(ratios, indent=2) + "\n")
        print(f"wrote {RATIO_CACHE}")

    if DST.exists() and args.force:
        shutil.rmtree(DST)
    DST.mkdir(parents=True, exist_ok=True)

    tasks = []
    pat = re.compile(r"_m([0-9.]+)_s")
    for src in SRC.glob("dRdE_Si_heavy_m*_s*.csv"):
        m = pat.search(src.name)
        if not m:
            continue
        scale = ratios.get(m.group(1)) or ratios.get(f"{float(m.group(1)):.6f}")
        if scale is None:
            scale = ratios[f"{float(m.group(1)):.6f}"]
        tasks.append((src, DST / src.name, float(scale)))

    print(f"scaling {len(tasks)} CSV files with {args.workers} workers...")
    done = 0
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(scale_one, s, d, sc): s for s, d, sc in tasks}
        for fut in as_completed(futs):
            fut.result()
            done += 1
            if done % 5000 == 0:
                print(f"  {done}/{len(tasks)}")
    print(f"done → {DST} ({done} files)")


if __name__ == "__main__":
    main()
