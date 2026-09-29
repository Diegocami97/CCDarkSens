#!/usr/bin/env python3
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: generate_dc_scaled_limit_configs.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  generate_dc_scaled_limit_configs.py -- Generate limit-scan and ne_imaging
#  configs with scaled dark current (10x, 100x).
# ============================================================================

"""Generate limit-scan and ne_imaging configs with scaled dark current (10x, 100x)."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONFIG_DIR = ROOT / "configs"

BASE_LAMBDA = 0.0365

KLEIN_LIMIT_BASES = sorted(
    p
    for p in CONFIG_DIR.glob("scan_band_gap*klein*.json")
    if "_dc10x" not in p.name and "_dc100x" not in p.name
)
SI_REF_LIMIT_BASES = [
    CONFIG_DIR / "scan_band_gap_pheno_1p2_eh3p8__refix.json",
    CONFIG_DIR / "scan_band_gap_light_pheno_1p2_eh3p8__refix.json",
]

NE_IMAGING_BASES = [
    CONFIG_DIR / "ne_imaging_one_point_si_ref.json",
    CONFIG_DIR / "ne_imaging_one_point_si_ref_light.json",
    CONFIG_DIR / "ne_imaging_one_point_klein_gap0p1.json",
    CONFIG_DIR / "ne_imaging_one_point_klein_gap0p5.json",
    CONFIG_DIR / "ne_imaging_one_point_klein_gap0p1_light.json",
    CONFIG_DIR / "ne_imaging_one_point_klein_gap0p5_light.json",
]

DC_TIERS = {
    "dc10x": 10.0,
    "dc100x": 100.0,
}


# ----------------------------------------------------------------------------
# scale_config
#   Copy of a config with the dark current set to the base value times factor and the label, output directory and comment tagged with the tier.
# ----------------------------------------------------------------------------
def scale_config(cfg: dict, tier: str, factor: float) -> dict:
    out = json.loads(json.dumps(cfg))
    lam = BASE_LAMBDA * factor
    out["backgrounds"]["dark_current"]["lambda_e_per_pix_per_year"] = lam

    run = out["run"]
    label = run["label"]
    outdir = run["outdir"]
    run["label"] = f"{label}_{tier}"
    run["outdir"] = f"{outdir}_{tier}"

    comment = out.get("_comment", "")
    out["_comment"] = (
        f"{comment}  Dark current {tier}: "
        f"lambda_e_per_pix_per_year={lam:g} ({factor:g}x baseline {BASE_LAMBDA})."
    )
    return out


# ----------------------------------------------------------------------------
# write_scaled_family
#   Write the scaled configs (every DC tier) of each existing base config and return the written paths; missing base files are skipped.
# ----------------------------------------------------------------------------
def write_scaled_family(base_paths: list[Path], kind: str) -> list[Path]:
    written: list[Path] = []
    for base_path in base_paths:
        if not base_path.is_file():
            print(f"SKIP missing {base_path}")
            continue
        with base_path.open() as f:
            cfg = json.load(f)
        stem = base_path.stem
        for tier, factor in DC_TIERS.items():
            scaled = scale_config(cfg, tier, factor)
            out_path = CONFIG_DIR / f"{stem}_{tier}.json"
            with out_path.open("w") as f:
                json.dump(scaled, f, indent=2)
                f.write("\n")
            written.append(out_path)
            lam = scaled["backgrounds"]["dark_current"]["lambda_e_per_pix_per_year"]
            print(f"[{kind}] Wrote {out_path.name}  lambda={lam:g}")
    return written


# ----------------------------------------------------------------------------
# main
#   Command line: --kind limit | ne_imaging | all selects which config families get 10x and 100x dark-current versions.
# ----------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--kind",
        choices=["limit", "ne_imaging", "all"],
        default="all",
        help="Which config families to generate (default: all)",
    )
    args = ap.parse_args()

    written: list[Path] = []
    if args.kind in ("limit", "all"):
        written += write_scaled_family(
            KLEIN_LIMIT_BASES + SI_REF_LIMIT_BASES, "limit"
        )
    if args.kind in ("ne_imaging", "all"):
        written += write_scaled_family(NE_IMAGING_BASES, "ne_imaging")

    print(f"\nTotal: {len(written)} configs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
