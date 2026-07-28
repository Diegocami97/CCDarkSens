#!/usr/bin/env python3
"""
Generate scan configs for 1, 10, and 100 g-yr exposures.
Copies existing 1 kg-yr configs and changes mass_kg only.

Produces:
  configs/Sr2Cb2Sd/scan_srcd_gap0p34_{heavy|light}_dc{1e2|1e3|1e5}_{1|10|100}gyr.json
  configs/darkphoton_scan_hypmat_unscreened_ne{|_dc1e3|_dc1e2}_{1|10|100}gyr.json

Usage:
  python3 utils/generate_exposure_configs.py
"""
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

EXPOSURES = {
    "1gyr":   0.001,
    "10gyr":  0.01,
    "100gyr": 0.1,
}

# --- DM-electron configs ---
DME_SOURCES = {
    "heavy": {
        "dc1e5": "configs/Sr2Cb2Sd/scan_srcd_gap0p34_heavy_dc1e5.json",
        "dc1e3": "configs/Sr2Cb2Sd/scan_srcd_gap0p34_heavy_dc1e3.json",
        "dc1e2": "configs/Sr2Cb2Sd/scan_srcd_gap0p34_heavy_dc1e2.json",
    },
    "light": {
        "dc1e5": "configs/Sr2Cb2Sd/scan_srcd_gap0p34_light_dc1e5.json",
        "dc1e3": "configs/Sr2Cb2Sd/scan_srcd_gap0p34_light_dc1e3.json",
        "dc1e2": "configs/Sr2Cb2Sd/scan_srcd_gap0p34_light_dc1e2.json",
    },
}

# --- Dark photon unscreened configs ---
DP_SOURCES = {
    "dc1e5": "configs/darkphoton_scan_hypmat_unscreened_ne.json",
    "dc1e3": "configs/darkphoton_scan_hypmat_unscreened_ne_dc1e3.json",
    "dc1e2": "configs/darkphoton_scan_hypmat_unscreened_ne_dc1e2.json",
}


def write_config(cfg: dict, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(cfg, f, indent=2)
    print(f"  wrote {out_path.relative_to(REPO)}")


def make_dme(mediator: str, dc_tag: str, exp_tag: str, mass_kg: float, src_path: str) -> None:
    with open(REPO / src_path) as f:
        cfg = json.load(f)

    base_label = f"srcd_gap0p34_{mediator}_{dc_tag}_{exp_tag}"
    cfg["run"]["label"]  = base_label
    cfg["run"]["outdir"] = f"outputs/Sr2Cb2Sd/{base_label}"
    cfg["detector"]["mass_kg"] = mass_kg

    # update comment
    cfg["_comment"] = (
        f"SrCd2Sb2 gap=0.34 eV, {mediator} mediator, {exp_tag}, "
        f"Klein eh=1.7778, DC={dc_tag.replace('dc','').replace('e','e-')} e-/pix/day, ROI n_e 1-40."
    )

    out = REPO / f"configs/Sr2Cb2Sd/scan_{base_label}.json"
    write_config(cfg, out)


def make_dp(dc_tag: str, exp_tag: str, mass_kg: float, src_path: str) -> None:
    with open(REPO / src_path) as f:
        cfg = json.load(f)

    dc_suffix = "" if dc_tag == "dc1e5" else f"_{dc_tag}"
    base_label = f"darkphoton_hypmat_unscreened_ne{dc_suffix}_{exp_tag}"
    cfg["run"]["label"]  = base_label
    cfg["run"]["outdir"] = f"outputs/darkphoton/hypmat_unscreened_ne{dc_suffix}_{exp_tag}"
    cfg["detector"]["mass_kg"] = mass_kg

    out_name = f"darkphoton_scan_hypmat_unscreened_ne{dc_suffix}_{exp_tag}.json"
    out = REPO / "configs" / out_name
    write_config(cfg, out)


if __name__ == "__main__":
    print("=== DM-electron configs ===")
    for mediator, dcs in DME_SOURCES.items():
        for dc_tag, src in dcs.items():
            for exp_tag, mass_kg in EXPOSURES.items():
                make_dme(mediator, dc_tag, exp_tag, mass_kg, src)

    print("\n=== Dark photon unscreened configs ===")
    for dc_tag, src in DP_SOURCES.items():
        for exp_tag, mass_kg in EXPOSURES.items():
            make_dp(dc_tag, exp_tag, mass_kg, src)

    print("\nDone. Run scans with ccdarksens_scan_dmelectron_pattern for each new config.")
