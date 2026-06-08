#!/usr/bin/env python3
# ============================================================================
#  CCDarkSens — baxter_halo_paper_smoke
#  Smoke test comparing Baxter 2021 halo model at v_E=253.7 vs 263 km/s
#
#  Author: Diego Venegas-Vargas
# ============================================================================
"""Baxter 2021 halo smoke: paper-export masses at v_E=253.7 vs 263 km/s."""
from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

_pydme_ref_dir = os.environ.get("PYDME_REF_DIR", "")
PAPER_EXPORT = (
    Path(_pydme_ref_dir)
    / "ScienceRun2024-figures/data/ScienceRun2024_results-1/"
      "DAMIC-M_2025_QEDark_DMe_heavymediator.txt"
) if _pydme_ref_dir else None
SRC_RATES = ROOT / "data/qedark_rates/Si/heavy/long_scan"

# Baxter 2021 Table 1 + March 9 static v_⊕ (EPJC 81:907, Sec. 3.1)
V0 = np.array([0.0, 238.0, 0.0])
V_SUN = np.array([11.1, 12.2, 7.3])
V_EARTH_MAR9 = np.array([29.2, -0.1, 5.9])
BAXTER_VE = float(np.linalg.norm(V0 + V_SUN + V_EARTH_MAR9))

# Exact paper-export mass points (MeV); no 5.0 MeV in export — use 5.5
PAPER_SMOKE_MASSES = [1.0, 2.0, 5.5]

BASE_HALO = {"v0_kms": 238.0, "vE_kms": 263.0, "vesc_kms": 544.0}
BAXTER_HALO = {"v0_kms": 238.0, "vE_kms": BAXTER_VE, "vesc_kms": 544.0}


def load_paper_export(path: Path):
    rows = []
    for ln in path.read_text().splitlines():
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        parts = ln.replace("\t", ",").split(",")
        if parts[0].lower().startswith("mass"):
            continue
        m_raw, s = float(parts[0]), float(parts[1])
        rows.append((m_raw / 1e6 if m_raw > 1e4 else m_raw, s))
    arr = np.asarray(rows)
    return arr[np.argsort(arr[:, 0])]


def fmt_mchi(m: float) -> str:
    return f"{m:.6f}"


def integrated_rate(m_mev: float, halo: dict) -> float:
    from ccdarkphys.qedark.entry import compute_dRdE

    out = compute_dRdE(
        "Si", "heavy", m_mev * 1e6, 1e-36, halo, binsize_eV=0.1
    )
    E = out["E_eV"]
    mask = E >= 2.0
    return float(np.trapz(out["dRdE_kg_year_eV"][mask], E[mask]) / 365.25 / 1000)


def rate_ratio(m_mev: float) -> float:
    return integrated_rate(m_mev, BAXTER_HALO) / integrated_rate(m_mev, BASE_HALO)


def scale_rates(masses: list[float], ratios: dict[float, float], dst_dir: Path):
    pat = re.compile(r"(_m)([0-9.]+)(_s)")
    avail = sorted(
        {
            float(re.search(r"_m([0-9.]+)_s", p.name).group(1))
            for p in SRC_RATES.glob("dRdE_Si_heavy_m*_s*.csv")
            if re.search(r"_m([0-9.]+)_s", p.name)
        }
    )
    if dst_dir.exists():
        shutil.rmtree(dst_dir)
    dst_dir.mkdir(parents=True)
    n = 0
    for m in masses:
        m_src = min(avail, key=lambda x: abs(x - m))
        scale = ratios[m]
        for src in SRC_RATES.glob(f"dRdE_Si_heavy_m{fmt_mchi(m_src)}_s*.csv"):
            dst = dst_dir / pat.sub(rf"\g<1>{fmt_mchi(m)}\g<3>", src.name)
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
            dst.write_text("\n".join(lines) + "\n")
            n += 1
    print(f"scaled {n} files → {dst_dir}")


def write_config(path: Path, label: str, outdir: str, rates_dir: str, masses: list[float]):
    base = json.loads(
        (ROOT / "configs/scan_dmelectron_pattern_data_qedark_fewmass_asis.json").read_text()
    )
    base["_comment"] = f"Baxter halo smoke ({label})"
    base["run"]["label"] = label
    base["run"]["outdir"] = outdir
    base["model"]["rates_dir"] = rates_dir
    base["model"]["grid"]["mchi_MeV"]["values"] = masses
    path.write_text(json.dumps(base, indent=2) + "\n")


def run_scan(config: Path) -> int:
    exe = ROOT / "build/ccdarksens_scan_dmelectron_pattern"
    cmd = [str(exe), str(config)]
    print("run:", " ".join(cmd))
    return subprocess.call(cmd, cwd=ROOT)


def load_ul_by_index(root_path: Path) -> np.ndarray:
    import uproot

    v = uproot.open(root_path)["upper_limit_sigma_e_mchi"].values()
    ok = (v > 0) & (v < 0.9e-26)
    return v[ok]


def main():
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--run-scan", action="store_true")
    ap.add_argument("--skip-rates", action="store_true")
    args = ap.parse_args()

    paper = load_paper_export(PAPER_EXPORT)
    pm, ps = paper[:, 0], paper[:, 1]

    print("=== Baxter 2021 static halo (March 9 v_⊕) ===")
    print(f"v_lab = v0 + v_sun + v_earth = {V0 + V_SUN + V_EARTH_MAR9}")
    print(f"v_E = |v_lab| = {BAXTER_VE:.2f} km/s  (compare baseline vE=263, DIM 253.7)")
    print(f"paper smoke masses [MeV]: {PAPER_SMOKE_MASSES}")
    print()

    ratios = {m: rate_ratio(m) for m in PAPER_SMOKE_MASSES}
    print(f"{'m [MeV]':>8}  R(Baxter)/R(263)  expected UL ratio")
    for m in PAPER_SMOKE_MASSES:
        print(f"{m:8.1f}  {ratios[m]:13.4f}  {1.0/ratios[m]:13.4f}")

    rates263 = ROOT / "data/qedark_rates/Si/heavy/baxter_smoke_vE263"
    rates254 = ROOT / "data/qedark_rates/Si/heavy/baxter_smoke_vE253p7"
    cfg263 = ROOT / "configs/scan_baxter_halo_smoke_vE263.json"
    cfg254 = ROOT / "configs/scan_baxter_halo_smoke_vE253p7.json"

    if not args.skip_rates:
        # vE=263: copy nearest long_scan rates to exact paper mass filenames
        scale_rates(PAPER_SMOKE_MASSES, {m: 1.0 for m in PAPER_SMOKE_MASSES}, rates263)
        scale_rates(PAPER_SMOKE_MASSES, ratios, rates254)

    write_config(cfg263, "baxter_smoke_vE263", "outputs/scan_baxter_halo_smoke_vE263",
                 "data/qedark_rates/Si/heavy/baxter_smoke_vE263", PAPER_SMOKE_MASSES)
    write_config(cfg254, "baxter_smoke_vE253p7", "outputs/scan_baxter_halo_smoke_vE253p7",
                 "data/qedark_rates/Si/heavy/baxter_smoke_vE253p7", PAPER_SMOKE_MASSES)

    roots = {}
    if args.run_scan:
        for tag, cfg in [("vE263", cfg263), ("Baxter253.7", cfg254)]:
            rc = run_scan(cfg)
            if rc != 0:
                sys.exit(rc)
            roots[tag] = ROOT / f"outputs/scan_baxter_halo_smoke_{'vE263' if tag=='vE263' else 'vE253p7'}" / "scan_dmelectron_pattern.root"
    else:
        for tag, sub in [("vE263", "vE263"), ("Baxter253.7", "vE253p7")]:
            p = ROOT / f"outputs/scan_baxter_halo_smoke_{sub}/scan_dmelectron_pattern.root"
            if p.exists():
                roots[tag] = p

    fullgrid = ROOT / "outputs/scan_pattern_data_qedark_fullgrid/scan_dmelectron_pattern.root"
    fg_ul = {}
    if fullgrid.exists():
        import uproot

        fg = uproot.open(fullgrid)["upper_limit_sigma_e_mchi"]
        fxc = 0.5 * (fg.axis().edges()[:-1] + fg.axis().edges()[1:])
        fv = fg.values()
        for m in PAPER_SMOKE_MASSES:
            i = int(np.argmin(np.abs(fxc - m)))
            fg_ul[m] = float(fv[i])

    print("\n=== UL vs paper export ===")
    hdr = f"{'m':>6} {'paper':>12}"
    for tag in roots:
        hdr += f" {tag+'/pap':>14}"
    if fg_ul:
        hdr += f" {'fullgrid/pap':>14}"
    print(hdr)
    print("-" * len(hdr))

    rows = []
    for m in PAPER_SMOKE_MASSES:
        pap = float(np.interp(m, pm, ps))
        line = f"{m:6.1f} {pap:12.3e}"
        row = {"m_mev": m, "paper_ul": pap}
        for tag, rp in roots.items():
            ul = load_ul_by_index(rp)
            # scan order matches PAPER_SMOKE_MASSES
            ix = PAPER_SMOKE_MASSES.index(m)
            val = float(ul[ix]) if ix < len(ul) else float("nan")
            ratio = val / pap
            line += f" {ratio:14.3f}"
            row[f"{tag}_over_paper"] = ratio
            row[f"{tag}_ul"] = val
        if fg_ul and m in fg_ul:
            r = fg_ul[m] / pap
            line += f" {r:14.3f}"
            row["fullgrid_over_paper"] = r
        print(line)
        rows.append(row)

    out = ROOT / "outplots/qedark_repro/baxter_halo_paper_smoke_summary.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    keys = sorted({k for r in rows for k in r})
    with out.open("w") as f:
        f.write(",".join(keys) + "\n")
        for r in rows:
            f.write(",".join(str(r.get(k, "")) for k in keys) + "\n")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
