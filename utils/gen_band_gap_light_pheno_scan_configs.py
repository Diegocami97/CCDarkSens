#!/usr/bin/env python3
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: gen_band_gap_light_pheno_scan_configs.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  gen_band_gap_light_pheno_scan_configs.py -- Generate scan configs for
#  band-gap pheno scans with the light mediator (QCDark2)
# ============================================================================
"""
Generate Phase C scan JSONs for band-gap pheno study with light mediator.

Creates 12 configs:
  configs/scan_band_gap_light_pheno_{gap_short}_eh{eh_tag}.json

where:
  gap_short in {0p1,0p3,0p5,0p7,0p9,1p2}
  eh_tag is either {gap_short} (D-equal) or 3p8 (B-thresh)
"""

from __future__ import annotations

import copy
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "configs"

TEMPLATE = ROOT / "configs" / "scan_band_gap_pheno_0p1_eh3p8.json"

GAPS = [0.1, 0.3, 0.5, 0.7, 0.9, 1.2]
EH_B = 3.8


# ----------------------------------------------------------------------------
# gap_short
#   Band gap without the "gap" prefix, e.g. "0p7" (1.2 eV gives "1p2").
# ----------------------------------------------------------------------------
def gap_short(g: float) -> str:
    return "1p2" if abs(g - 1.2) < 1e-9 else f"{g:.1f}".replace(".", "p")


# ----------------------------------------------------------------------------
# gap_tag
#   File-name tag of a band gap, e.g. "gap0p7" (1.2 eV gives "gap1p2").
# ----------------------------------------------------------------------------
def gap_tag(g: float) -> str:
    return "gap1p2" if abs(g - 1.2) < 1e-9 else f"gap{g:.1f}".replace(".", "p")


# ----------------------------------------------------------------------------
# eh_tag
#   Electron-hole pair energy formatted for file names with '.' replaced by 'p'.
# ----------------------------------------------------------------------------
def eh_tag(eh: float) -> str:
    return f"{eh:g}".replace(".", "p")


# ----------------------------------------------------------------------------
# main
#   Write the light-mediator Phase C pheno scan configs: one per band gap and scenario (D-equal, B-thresh), all derived from the template.
# ----------------------------------------------------------------------------
def main() -> int:
    base = json.loads(TEMPLATE.read_text(encoding="utf-8"))

    written = 0
    for g in GAPS:
        for scenario in ("D-equal", "B-thresh"):
            eh = g if scenario == "D-equal" else EH_B

            gs = gap_short(g)  # 0p1
            gt = gap_tag(g)  # gap0p1
            et = eh_tag(eh)  # 0p1 or 3p8

            cfg = copy.deepcopy(base)

            cfg["_comment"] = (
                f"Band-gap pheno Phase C (light mediator): scissor {g:g} eV, {scenario} "
                f"(E_gap={g:g} eV, e_h={eh:g} eV). Rates: Si_fast_{gt} (light)."
            )
            cfg["run"]["label"] = f"scan_band_gap_light_pheno_{gs}_eh{et}"
            cfg["run"]["outdir"] = f"outputs/scan_band_gap_light_{gs}_eh{et}"

            # Ionization: only table changes
            ci = cfg["response"]["charge_ionization"]
            ci["table_csv"] = f"data/p100K_{gt}_eh{eh_tag(eh)}.csv"
            ci["band_gap_eV"] = float(g)
            ci["eh_pair_eV"] = float(eh)
            ci["scenario"] = scenario

            # Rates: mediator changes the QCDark2 backend output
            model = cfg["model"]
            model["mediator"] = "light"
            model["rates_dir"] = f"data/qcdark2_rates/Si/light/Si_fast_{gt}"
            # filename_template keeps {mediator} token so it matches the generated files.

            out = OUTDIR / f"scan_band_gap_light_pheno_{gs}_eh{et}.json"
            out.write_text(json.dumps(cfg, indent=2) + "\n", encoding="utf-8")
            written += 1
            print(f"[ok] {out.name}")

    print(f"\nWrote {written} light-mediator scan configs.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

