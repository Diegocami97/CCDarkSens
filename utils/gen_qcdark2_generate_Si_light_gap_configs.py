#!/usr/bin/env python3
# ============================================================================
#  CCDarkSens — gen_qcdark2_generate_Si_light_gap_configs
#  Generate QCDark2 rate configs for Si with varying band gaps (light mediator)
#
#  Author: Diego Venegas-Vargas
# ============================================================================
"""
Generate QCDark2 dR/dE grid configs for Si with light mediator (Phase C input).

Writes:
  configs/qcdark2_generate_Si_light_gap0p1.json
  ...

These configs are compatible with:
  python3 utils/qcdark2_generate_grid.py <config.json>
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "configs"

EH_REF = 3.8
PARALLEL = 1

TEMPLATE = ROOT / "configs" / "qcdark2_generate_Si_heavy_gap0p7.json"

GAPS = {
    "0p1": 0.1,
    "0p3": 0.3,
    "0p5": 0.5,
    "0p7": 0.7,
    "0p9": 0.9,
    "1p2": 1.2,
}


def gap_tag(short: str) -> str:
    return f"gap{short}"


def main() -> int:
    base = json.loads(TEMPLATE.read_text(encoding="utf-8"))

    written = 0
    for short, gap_eV in GAPS.items():
        tag = gap_tag(short)  # gap0p7 etc
        cfg = json.loads(json.dumps(base))  # cheap deep copy

        cfg["model_type"] = "dm_electron"
        cfg["material"] = "Si"
        cfg["mediator"] = "light"

        cfg["detector"]["band_gap_eV"] = float(gap_eV)
        cfg["detector"]["eh_pair_eV"] = float(EH_REF)

        cfg["epsilon_h5"] = f"data/qcdark2_epsilon/Si/Si_fast_{tag}.h5"
        cfg["rates_dir"] = f"data/qcdark2_rates/Si/light/Si_fast_{tag}"

        # Keep filename_template with {mediator} token.
        # qcdark2_generate_grid expands it using cfg["mediator"].

        cfg["options"]["parallel"] = PARALLEL

        out = OUTDIR / f"qcdark2_generate_Si_light_{tag}.json"
        out.write_text(json.dumps(cfg, indent=2) + "\n", encoding="utf-8")
        written += 1
        print(f"[ok] {out.name}")

    print(f"\nWrote {written} light rate configs.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

