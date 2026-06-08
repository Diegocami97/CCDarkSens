#!/usr/bin/env python3
# ============================================================================
#  CCDarkSens — halo_smoke_test_qedark
#  Halo sensitivity smoke test comparing v_E=253.7 km/s vs v_E=263 km/s baseline
#
#  Author: Diego Venegas-Vargas
# ============================================================================
"""Halo sensitivity smoke test: few-mass scan at v_E=253.7 km/s vs baseline v_E=263."""
from typing import Optional
import argparse
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

BASE_HALO = {"v0_kms": 238.0, "vE_kms": 263.0, "vesc_kms": 544.0}
ALT_HALO = {"v0_kms": 238.0, "vE_kms": 253.7, "vesc_kms": 544.0}
SMOKE_MASSES = [1.000194, 1.99987, 5.001944]
EMIN_EV = 2.0


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


def integrated_rate_mev(m_mev: float, halo: dict, sigma: float = 1e-36) -> float:
    from ccdarkphys.qedark.entry import compute_dRdE

    out = compute_dRdE(
        material="Si",
        mediator="heavy",
        mchi_eV=m_mev * 1e6,
        sigma_e_cm2=sigma,
        halo=halo,
        binsize_eV=0.1,
    )
    E = out["E_eV"]
    mask = E >= EMIN_EV
    return float(np.trapz(out["dRdE_kg_year_eV"][mask], E[mask]) / 365.25 / 1000)


def rate_ratio(m_mev: float) -> float:
    r_base = integrated_rate_mev(m_mev, BASE_HALO)
    r_alt = integrated_rate_mev(m_mev, ALT_HALO)
    return r_alt / r_base


def scale_rates_dir(src_dir: Path, dst_dir: Path, masses: list[float], ratios: dict[float, float]):
    """Copy sigma grid for each mass, scaling dRdE by halo ratio (uniform in sigma)."""
    dst_dir.mkdir(parents=True, exist_ok=True)
    pat = re.compile(r"_m([0-9.]+)_s")
    n_copy = 0
    for src in src_dir.glob("dRdE_Si_heavy_m*_s*.csv"):
        m = pat.search(src.name)
        if not m:
            continue
        mval = float(m.group(1))
        nearest = min(masses, key=lambda x: abs(x - mval))
        if abs(nearest - mval) > 1e-3:
            continue
        scale = ratios[nearest]
        dst = dst_dir / src.name
        lines_out = []
        for ln in src.read_text().splitlines():
            if ln.startswith("#") or not ln.strip():
                lines_out.append(ln)
                continue
            parts = ln.split(",")
            if len(parts) < 2:
                lines_out.append(ln)
                continue
            try:
                e, rate = parts[0], float(parts[1]) * scale
            except ValueError:
                lines_out.append(ln)
                continue
            lines_out.append(f"{e},{rate:.16e}")
        dst.write_text("\n".join(lines_out) + "\n")
        n_copy += 1
    return n_copy


def ul_at_mass(ulmap: dict, m: float, tol: float = 0.35) -> Optional[float]:
    if not ulmap:
        return None
    nearest = min(ulmap.keys(), key=lambda x: abs(x - m))
    if abs(nearest - m) > tol:
        return None
    return ulmap[nearest]


def load_scan_ul(root_path: Path):
    import uproot

    f = uproot.open(root_path)
    ul = f["upper_limit_sigma_e_mchi"]
    xe = ul.axis().edges()
    xc = 0.5 * (xe[:-1] + xe[1:])
    v = ul.values()
    ok = (v > 0) & (v < 0.9e-26)
    return xc[ok], v[ok]


def make_config(out_path: Path, rates_dir: str, label: str, outdir: str):
    base = json.loads(
        (ROOT / "configs/scan_dmelectron_pattern_data_qedark_fewmass_asis.json").read_text()
    )
    base["run"]["label"] = label
    base["run"]["outdir"] = outdir
    base["model"]["rates_dir"] = rates_dir
    base["model"]["grid"]["mchi_MeV"]["values"] = SMOKE_MASSES
    out_path.write_text(json.dumps(base, indent=2) + "\n")


def run_scan(config_path: Path, build_dir: Path) -> int:
    exe = build_dir / "ccdarksens_scan_dmelectron_pattern"
    if not exe.exists():
        print(f"skip scan: {exe} not found")
        return 1
    cmd = [str(exe), str(config_path)]
    print("run:", " ".join(cmd))
    return subprocess.call(cmd, cwd=ROOT)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build-dir", default=str(ROOT / "build"))
    ap.add_argument("--run-scan", action="store_true", help="Run CCDarkSens few-mass scans")
    ap.add_argument("--skip-rate-copy", action="store_true")
    args = ap.parse_args()

    paper = load_paper_export(PAPER_EXPORT)
    pm, ps = paper[:, 0], paper[:, 1]

    print("=== Halo rate ratios (integrated E >= 2 eV, σ=1e-36) ===")
    print(f"{'m [MeV]':>10} {'R(253.7)/R(263)':>16} {'linear UL shift':>16}")
    ratios = {}
    for m in SMOKE_MASSES:
        rr = rate_ratio(m)
        ratios[m] = rr
        print(f"{m:10.4f} {rr:16.4f} {1.0 / rr:16.4f}")

    src_rates = ROOT / "data/qedark_rates/Si/heavy/long_scan"
    alt_rates = ROOT / "data/qedark_rates/Si/heavy/halo_smoke_vE253p7"
    if not args.skip_rate_copy:
        if alt_rates.exists():
            shutil.rmtree(alt_rates)
        n = scale_rates_dir(src_rates, alt_rates, SMOKE_MASSES, ratios)
        print(f"\nScaled {n} rate files → {alt_rates}")

    cfg_base = ROOT / "configs/scan_qedark_halo_smoke_vE263_baseline.json"
    cfg_alt = ROOT / "configs/scan_qedark_halo_smoke_vE253p7.json"
    make_config(
        cfg_base,
        "data/qedark_rates/Si/heavy/long_scan",
        "halo_smoke_vE263",
        "outputs/scan_qedark_halo_smoke_vE263",
    )
    make_config(
        cfg_alt,
        "data/qedark_rates/Si/heavy/halo_smoke_vE253p7",
        "halo_smoke_vE253p7",
        "outputs/scan_qedark_halo_smoke_vE253p7",
    )
    print(f"wrote {cfg_base}")
    print(f"wrote {cfg_alt}")

    results = {}
    scan_outputs = [
        ("vE263", ROOT / "outputs/scan_qedark_halo_smoke_vE263/scan_dmelectron_pattern.root"),
        ("vE253.7", ROOT / "outputs/scan_qedark_halo_smoke_vE253p7/scan_dmelectron_pattern.root"),
    ]
    if args.run_scan:
        for tag, cfg, _ in [
            ("vE263", cfg_base, scan_outputs[0][1]),
            ("vE253.7", cfg_alt, scan_outputs[1][1]),
        ]:
            rc = run_scan(cfg, Path(args.build_dir))
            if rc != 0:
                print(f"scan failed for {tag} (rc={rc})")

    for tag, out in scan_outputs:
        if out.exists():
            mc, sg = load_scan_ul(out)
            results[tag] = {float(m): float(s) for m, s in zip(mc, sg)}
            print(f"loaded {tag} smoke scan: masses={list(results[tag].keys())}")

    # Compare stored fullgrid at same masses if available
    fullgrid = ROOT / "outputs/scan_pattern_data_qedark_fullgrid/scan_dmelectron_pattern.root"
    if fullgrid.exists():
        mc, sg = load_scan_ul(fullgrid)
        results["fullgrid_vE263"] = {}
        for m in SMOKE_MASSES:
            i = int(np.argmin(np.abs(mc - m)))
            results["fullgrid_vE263"][m] = sg[i]

    # Masses to report: smoke-scan bin centers + nominal targets
    report_masses = list(SMOKE_MASSES)
    for tag in ("vE263", "vE253.7"):
        if tag in results:
            for m in results[tag]:
                if all(abs(m - x) > 0.05 for x in report_masses):
                    report_masses.append(m)
    report_masses = sorted(set(report_masses))

    print("\n=== UL comparison vs paper export ===")
    print(f"{'m [MeV]':>10} {'paper':>12} ", end="")
    for tag in results:
        print(f"{tag + '/pap':>18}", end="")
    print()
    print("-" * (22 + 18 * len(results)))

    rows_csv = []
    for m in report_masses:
        pap = float(np.interp(m, pm, ps))
        line = f"{m:10.4f} {pap:12.3e} "
        row = {"m_mev": m, "paper_ul": pap}
        for tag, ulmap in results.items():
            val = ul_at_mass(ulmap, m) if isinstance(ulmap, dict) else None
            if val is not None:
                ratio = val / pap
                line += f"{ratio:18.3f} "
                row[f"{tag}_ratio"] = ratio
            else:
                line += f"{'n/a':>18} "
        print(line)
        rows_csv.append(row)

    out_csv = ROOT / "outplots/qedark_repro/halo_smoke_test_summary.csv"
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    if rows_csv:
        keys = sorted({k for r in rows_csv for k in r})
        with out_csv.open("w") as f:
            f.write(",".join(keys) + "\n")
            for r in rows_csv:
                f.write(",".join(str(r.get(k, "")) for k in keys) + "\n")
        print(f"\nwrote {out_csv}")

    print("\nNote: v_E=253.7 weakens signal ~11% at 2 MeV → UL ~12% higher (scan/paper closer to 1).")
    print("Remaining low-m gap (stored/paper ~0.83 at 2 MeV) is not explained by halo alone.")


if __name__ == "__main__":
    main()
