#!/usr/bin/env python3
# ============================================================================
#  CCDarkSens — gen_band_gap_2d_grid_configs
#  Generate JSON configs for the 2D (E_gap, ε_h) sensitivity grid scan
#
#  Author: Diego Venegas-Vargas
# ============================================================================
"""
Generate 2D (gap, eh) Phase-C scan configs and build missing p100K tables.

Grid:
  gap in {0.1, 0.3, 0.5, 0.7, 0.9, 1.2}
  eh  in {0.5, 1.0, 1.5, 2.0, 2.5, 3.8}

Valid cells satisfy eh >= gap.

Existing cells to skip (already complete, heavy+light):
  - B-thresh row: eh=3.8 for all six gaps
  - D-equal diagonal: (gap,eh)=(0.5,0.5),(0.7,0.7),(0.9,0.9),(1.2,1.2)
    Note: (0.1,0.1) and (0.3,0.3) are outside this eh-grid.
"""

from __future__ import annotations

import copy
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONFIGS = ROOT / "configs"
MANIFEST = CONFIGS / "band_gap_pheno_scenarios.json"
TEMPLATE = CONFIGS / "scan_band_gap_pheno_0p1_eh3p8.json"
BUILD_P100K = ROOT / "utils" / "build_p100K_scaled.py"

GAPS = [0.1, 0.3, 0.5, 0.7, 0.9, 1.2]
EHS = [0.5, 1.0, 1.5, 2.0, 2.5, 3.8]
MEDIATORS = ["heavy", "light"]


def ev_tag(x: float) -> str:
    # Keep one decimal for integer-like entries (1.0 -> 1p0).
    s = f"{x:.1f}"
    return s.replace(".", "p")


def gap_tag(g: float) -> str:
    return f"gap{ev_tag(g)}"


def is_valid_cell(gap: float, eh: float) -> bool:
    return eh >= gap


def is_existing_complete(gap: float, eh: float) -> bool:
    # Existing B-thresh row.
    if abs(eh - 3.8) < 1e-12:
        return True
    # Existing D-equal points that are on this eh-grid.
    return abs(eh - gap) < 1e-12 and gap in {0.5, 0.7, 0.9, 1.2}


def all_valid_cells() -> list[tuple[float, float]]:
    return [(g, e) for g in GAPS for e in EHS if is_valid_cell(g, e)]


def new_cells() -> list[tuple[float, float]]:
    return [(g, e) for (g, e) in all_valid_cells() if not is_existing_complete(g, e)]


def table_csv(gap: float, eh: float) -> str:
    return f"data/p100K_{gap_tag(gap)}_eh{ev_tag(eh)}.csv"


def config_name(mediator: str, gap: float, eh: float) -> str:
    return f"scan_band_gap_2d_{mediator}_{ev_tag(gap)}_eh{ev_tag(eh)}.json"


def outdir_name(mediator: str, gap: float, eh: float) -> str:
    return f"outputs/scan_band_gap_2d_{mediator}_{ev_tag(gap)}_eh{ev_tag(eh)}"


def existing_scan_paths(mediator: str, gap: float, eh: float) -> tuple[str, str]:
    gs = ev_tag(gap)
    es = ev_tag(eh)
    if mediator == "heavy":
        cfg = f"configs/scan_band_gap_pheno_{gs}_eh{es}.json"
        outdir = f"outputs/scan_band_gap_{gs}_eh{es}"
    else:
        cfg = f"configs/scan_band_gap_light_pheno_{gs}_eh{es}.json"
        outdir = f"outputs/scan_band_gap_light_{gs}_eh{es}"
    return cfg, outdir


def build_missing_p100k(cells: list[tuple[float, float]], force: bool) -> None:
    for gap, eh in cells:
        out_rel = table_csv(gap, eh)
        out_abs = ROOT / out_rel
        if out_abs.is_file() and not force:
            print(f"[skip] p100K exists: {out_rel}")
            continue
        cmd = [
            sys.executable,
            str(BUILD_P100K),
            "--band-gap-eV",
            str(gap),
            "--eh-pair-eV",
            str(eh),
            "--out",
            out_rel,
            "--scenario",
            "2d-grid",
        ]
        if force:
            cmd.append("--force")
        print("[run]", " ".join(cmd))
        rc = subprocess.call(cmd, cwd=ROOT)
        if rc != 0:
            raise RuntimeError(f"build_p100K_scaled failed for gap={gap}, eh={eh}")


def write_scan_configs(cells: list[tuple[float, float]]) -> list[dict]:
    base = json.loads(TEMPLATE.read_text(encoding="utf-8"))
    entries: list[dict] = []

    for mediator in MEDIATORS:
        for gap, eh in cells:
            cfg = copy.deepcopy(base)
            gt = gap_tag(gap)
            gtag = ev_tag(gap)
            ehtag = ev_tag(eh)

            cfg["_comment"] = (
                "Band-gap pheno 2D grid cell "
                f"(mediator={mediator}, gap={gap:g} eV, eh={eh:g} eV)."
            )
            cfg["run"]["label"] = f"scan_band_gap_2d_{mediator}_{gtag}_eh{ehtag}"
            cfg["run"]["outdir"] = outdir_name(mediator, gap, eh)

            cfg["experiment"]["observable_bins"] = "ne"
            cfg["experiment"]["roi_bins"] = [1, 2, 3, 4, 5]

            ci = cfg["response"]["charge_ionization"]
            ci["table_csv"] = table_csv(gap, eh)
            ci["band_gap_eV"] = gap
            ci["eh_pair_eV"] = eh
            ci["scenario"] = "2D-grid"

            model = cfg["model"]
            model["mediator"] = mediator
            model["rates_dir"] = f"data/qcdark2_rates/Si/{mediator}/Si_fast_{gt}"
            model["grid"]["mchi_MeV"]["logspace"] = {
                "start": 0.2,
                "stop": 1000.0,
                "num": 80,
                "endpoint": True,
            }
            model["grid"]["sigma_e_cm2"]["logspace"] = {
                "start_exp": -46,
                "stop_exp": -26,
                "num": 30,
                "endpoint": True,
            }

            cfg_path = CONFIGS / config_name(mediator, gap, eh)
            cfg_path.write_text(json.dumps(cfg, indent=2) + "\n", encoding="utf-8")
            print(f"[ok] config {cfg_path.relative_to(ROOT)}")

            entries.append(
                {
                    "mediator": mediator,
                    "band_gap_eV": gap,
                    "eh_pair_eV": eh,
                    "valid": True,
                    "existing_complete": False,
                    "ionization_csv": table_csv(gap, eh),
                    "scan_config": str(cfg_path.relative_to(ROOT)),
                    "outdir": outdir_name(mediator, gap, eh),
                }
            )

    return entries


def update_manifest(cells_new: list[tuple[float, float]]) -> None:
    data = json.loads(MANIFEST.read_text(encoding="utf-8"))

    valid = []
    for gap in GAPS:
        for eh in EHS:
            if not is_valid_cell(gap, eh):
                valid.append(
                    {
                        "band_gap_eV": gap,
                        "eh_pair_eV": eh,
                        "valid": False,
                        "reason": "eh < gap",
                    }
                )
                continue
            valid.append(
                {
                    "band_gap_eV": gap,
                    "eh_pair_eV": eh,
                    "valid": True,
                    "existing_complete": is_existing_complete(gap, eh),
                    "new_cell": (gap, eh) in set(cells_new),
                }
            )

    scan_entries = []
    cells_new_set = set(cells_new)
    for mediator in MEDIATORS:
        for gap in GAPS:
            for eh in EHS:
                if not is_valid_cell(gap, eh):
                    continue
                if (gap, eh) in cells_new_set:
                    cfg = f"configs/{config_name(mediator, gap, eh)}"
                    out = outdir_name(mediator, gap, eh)
                    status = "new"
                else:
                    cfg, out = existing_scan_paths(mediator, gap, eh)
                    status = "existing_complete"
                scan_entries.append(
                    {
                        "mediator": mediator,
                        "band_gap_eV": gap,
                        "eh_pair_eV": eh,
                        "scan_config": cfg,
                        "outdir": out,
                        "status": status,
                    }
                )

    data["grid2d"] = {
        "description": "2D (gap,eh) scan grid; mediator-specific scan configs live in configs/scan_band_gap_2d_*",
        "gaps_eV": GAPS,
        "eh_grid_eV": EHS,
        "valid_cells_count": sum(1 for c in valid if c["valid"]),
        "new_cells_count": len(cells_new),
        "mediators": MEDIATORS,
        "cells": valid,
        "scan_entries": scan_entries,
    }
    MANIFEST.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"[ok] updated manifest: {MANIFEST.relative_to(ROOT)}")


def main() -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--force-p100k", action="store_true", help="Rebuild p100K tables even if they exist")
    args = ap.parse_args()

    cells_new = new_cells()
    print(f"[info] valid cells: {len(all_valid_cells())}, new cells: {len(cells_new)}")

    build_missing_p100k(cells_new, force=args.force_p100k)
    write_scan_configs(cells_new)
    update_manifest(cells_new)

    print("[done] 2D grid setup complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

