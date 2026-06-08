#!/usr/bin/env python3
# ============================================================================
#  CCDarkSens — gen_band_gap_pheno_scan_configs
#  Generate scan configs for Phase C band-gap pheno scan (6 gaps × B-thresh + D-equal)
#
#  Author: Diego Venegas-Vargas
# ============================================================================
"""Generate scan_band_gap_pheno_*.json for Phase C (6 gaps x D-equal + B-thresh)."""

from __future__ import annotations

import copy
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "configs" / "scan_band_gap_pheno_0p1_eh3p8.json"
GAPS = [0.1, 0.3, 0.5, 0.7, 0.9, 1.2]
EH_B_THRESH = 3.8


def gap_tag(g: float) -> str:
    return "gap1p2" if abs(g - 1.2) < 1e-9 else f"gap{g:.1f}".replace(".", "p")


def eh_tag(eh: float) -> str:
    return f"{eh:g}".replace(".", "p")


def config_stem(gap: float, scenario: str) -> str:
    eh = gap if scenario == "D-equal" else EH_B_THRESH
    gs = gap_tag(gap)[3:]  # gap0p1 -> 0p1
    return f"scan_band_gap_pheno_{gs}_eh{eh_tag(eh)}"


def build_case(gap: float, scenario: str, base: dict) -> dict:
    eh = gap if scenario == "D-equal" else EH_B_THRESH
    gt = gap_tag(gap)
    stem = config_stem(gap, scenario)
    cfg = copy.deepcopy(base)
    cfg["_comment"] = (
        f"Band-gap pheno Phase C: scissor {gap:g} eV, {scenario} "
        f"(E_gap={gap:g} eV, e_h={eh:g} eV). Rates: Si_fast_{gt}."
    )
    cfg["run"]["label"] = stem
    gap_short = gt[3:] if gt.startswith("gap") else gt  # gap0p1 -> 0p1
    cfg["run"]["outdir"] = f"outputs/scan_band_gap_{gap_short}_eh{eh_tag(eh)}"
    ci = cfg["response"]["charge_ionization"]
    ci["table_csv"] = f"data/p100K_{gt}_eh{eh_tag(eh)}.csv"
    ci["band_gap_eV"] = gap
    ci["eh_pair_eV"] = eh
    ci["scenario"] = scenario
    cfg["model"]["rates_dir"] = f"data/qcdark2_rates/Si/heavy/Si_fast_{gt}"
    return cfg


def main() -> int:
    base = json.loads(TEMPLATE.read_text(encoding="utf-8"))
    written = []
    for gap in GAPS:
        for scenario in ("D-equal", "B-thresh"):
            stem = config_stem(gap, scenario)
            path = ROOT / "configs" / f"{stem}.json"
            cfg = build_case(gap, scenario, base)
            path.write_text(json.dumps(cfg, indent=2) + "\n", encoding="utf-8")
            written.append(path.relative_to(ROOT).as_posix())
            print(f"[ok] {path.name}")
    print(f"\nWrote {len(written)} configs.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
