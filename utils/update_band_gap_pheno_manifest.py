#!/usr/bin/env python3
# ============================================================================
#  CCDarkSens — update_band_gap_pheno_manifest
#  Refresh configs/band_gap_pheno_scenarios.json Phase C entries
#
#  Author: Diego Venegas-Vargas
# ============================================================================
"""Refresh configs/band_gap_pheno_scenarios.json Phase C entries."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "configs" / "band_gap_pheno_scenarios.json"
GAPS = [0.1, 0.3, 0.5, 0.7, 0.9, 1.2]
EH_B = 3.8


def gap_tag(g: float) -> str:
    return "gap1p2" if abs(g - 1.2) < 1e-9 else f"gap{g:.1f}".replace(".", "p")


def eh_tag(eh: float) -> str:
    return f"{eh:g}".replace(".", "p")


def scenario_entry(gap: float, scenario: str) -> dict:
    eh = gap if scenario == "D-equal" else EH_B
    gt = gap_tag(gap)
    gs = gt[3:] if gt.startswith("gap") else gt
    return {
        "label": f"{gap:g} eV gap, {eh:g} eV e_h ({scenario})",
        "band_gap_eV": gap,
        "eh_pair_eV": eh,
        "scenario": scenario,
        "epsilon_h5": f"data/qcdark2_epsilon/Si/Si_fast_{gt}.h5",
        "rates_dir": f"data/qcdark2_rates/Si/heavy/Si_fast_{gt}",
        "ionization_csv": f"data/p100K_{gt}_eh{eh_tag(eh)}.csv",
        "scan_config": f"configs/scan_band_gap_pheno_{gs}_eh{eh_tag(eh)}.json",
        "outdir": f"outputs/scan_band_gap_{gs}_eh{eh_tag(eh)}",
        "phase_c": "pending",
    }


def main() -> int:
    man = json.loads(MANIFEST.read_text(encoding="utf-8"))
    for g in GAPS:
        man["gap_sweep_status"][f"{g:g}"] = {
            "epsilon": "done",
            "rates": "done",
            "scissor_scan": "done" if g not in (0.7, 0.9) else "config_ready",
            "p100K": "done",
        }

    scenarios = []
    for g in GAPS:
        for scen in ("D-equal", "B-thresh"):
            scenarios.append(scenario_entry(g, scen))

    # Keep optional A-ratio at 0.3
    scenarios.append(
        {
            "label": "0.3 eV gap, Si eh ratio (A-ratio)",
            "band_gap_eV": 0.3,
            "eh_pair_eV": 0.95,
            "scenario": "A-ratio",
            "epsilon_h5": "data/qcdark2_epsilon/Si/Si_fast_gap0p3.h5",
            "rates_dir": "data/qcdark2_rates/Si/heavy/Si_fast_gap0p3",
            "ionization_csv": "data/p100K_gap0p3_eh0p95.csv",
            "scan_config": "configs/scan_band_gap_pheno_0p3_eh0p95.json",
            "outdir": "outputs/scan_band_gap_0p3_eh0p95",
            "phase_c": "optional",
        }
    )

    man["scenarios"] = scenarios
    man["phase_c"] = {
        "description": "12 coupled scans: 6 gaps x (D-equal, B-thresh); ROI n_e 1-5",
        "runner": "python3 utils/run_band_gap_phase_c.py",
        "status": "in_progress",
    }
    MANIFEST.write_text(json.dumps(man, indent=2) + "\n", encoding="utf-8")
    print(f"[ok] updated {MANIFEST}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
