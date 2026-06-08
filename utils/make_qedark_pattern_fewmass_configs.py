#!/usr/bin/env python3
# ============================================================================
#  CCDarkSens — make_qedark_pattern_fewmass_configs
#  Generate few-mass pattern-scan configs for QEdark reproduction studies
#
#  Author: Diego Venegas-Vargas
# ============================================================================
"""Generate few-mass pattern-scan configs for qedark reproduction matrix."""
from __future__ import annotations

import copy
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE_HEAVY = ROOT / "configs/scan_dmelectron_pattern_data_qedark_fullgrid.json"
OUT_DIR = ROOT / "configs"  # keep configs at repo configs/ so relative data/ paths resolve

MASSES = [
    1.000194, 1.999870, 5.001944, 10.001296, 19.997409,
    50.016198, 100.006478, 199.961135, 500.129579, 1000.000000,
]

# From build/data_pattern.root (TParameter)
DATA_MASS_KG = 0.00014622166557101568
DATA_EXPOSURE_KG_YEAR = 3.4238720356355986e-05
DATA_LIVETIME_DAYS = DATA_EXPOSURE_KG_YEAR * 365.25 / DATA_MASS_KG  # ~85.5 d

VARIANTS = [
    {
        "name": "fewmass_tau_config_exp",
        "label": "fewmass_tau_config_exp",
        "outdir": "outputs/qedark_repro_matrix/fewmass_tau_config_exp",
        "constrain_use_tau_weighted": True,
        "constrain_n_bins": 1,
        "mass_kg": None,
        "livetime_days": None,
        "note": "Baseline: tau constraint, config exposure (current JSON).",
    },
    {
        "name": "fewmass_tau_data_exp",
        "label": "fewmass_tau_data_exp",
        "outdir": "outputs/qedark_repro_matrix/fewmass_tau_data_exp",
        "constrain_use_tau_weighted": True,
        "constrain_n_bins": 1,
        "mass_kg": DATA_MASS_KG,
        "livetime_days": DATA_LIVETIME_DAYS,
        "note": "Tau constraint, exposure aligned to data_pattern.root.",
    },
    {
        "name": "fewmass_pydme_multibin_config_exp",
        "label": "fewmass_pydme_multibin_config_exp",
        "outdir": "outputs/qedark_repro_matrix/fewmass_pydme_multibin_config_exp",
        "constrain_use_tau_weighted": False,
        "constrain_n_bins": 4450,
        "mass_kg": None,
        "livetime_days": None,
        "note": "pydme multibin constraint, config exposure (Bp/Br scale).",
    },
]


def make_config(base: dict, variant: dict, mediator: str, rates_dir: str, filename_template: str) -> dict:
    j = copy.deepcopy(base)
    j["_comment"] = variant["note"]
    j["run"]["label"] = variant["label"]
    j["run"]["outdir"] = variant["outdir"]
    j["run"]["constrain_use_tau_weighted"] = variant["constrain_use_tau_weighted"]
    j["run"]["constrain_n_bins"] = variant["constrain_n_bins"]
    if variant["mass_kg"] is not None:
        j["detector"]["mass_kg"] = variant["mass_kg"]
    if variant["livetime_days"] is not None:
        j["experiment"]["livetime_days"] = variant["livetime_days"]
    j["model"]["mediator"] = mediator
    j["model"]["rates_dir"] = rates_dir
    j["model"]["filename_template"] = filename_template
    j["model"]["grid"]["mchi_MeV"] = {"values": MASSES}
    return j


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    base = json.loads(BASE_HEAVY.read_text())

    manifest = []
    for v in VARIANTS:
        path = OUT_DIR / f"scan_qedark_repro_{v['name']}.json"
        cfg = make_config(
            base, v,
            mediator="heavy",
            rates_dir="data/qedark_rates/Si/heavy/long_scan",
            filename_template="dRdE_{material}_{mediator}_m{mchi_MeV}_s{sigma_e_cm2}.csv",
        )
        path.write_text(json.dumps(cfg, indent=2) + "\n")
        manifest.append({"variant": v["name"], "config": str(path.relative_to(ROOT)), "note": v["note"]})
        print("wrote", path)

    # Light mediator few-mass (config exposure — same Bp/Br scale as heavy reference workflow)
    light = {
        "name": "fewmass_light_tau_config_exp",
        "label": "fewmass_light_tau_config_exp",
        "outdir": "outputs/qedark_repro_matrix/fewmass_light_tau_config_exp",
        "constrain_use_tau_weighted": True,
        "constrain_n_bins": 1,
        "mass_kg": None,
        "livetime_days": None,
        "note": "Light (massless) mediator QEdark, tau constraint, config exposure.",
    }
    path = OUT_DIR / f"scan_qedark_repro_{light['name']}.json"
    cfg = make_config(
        base, light,
        mediator="massless",
        rates_dir="data/qedark_rates/Si/ultralight/long_scan",
        filename_template="dRdE_{material}_massless_m{mchi_MeV}_s{sigma_e_cm2}.csv",
    )
    path.write_text(json.dumps(cfg, indent=2) + "\n")
    manifest.append({"variant": light["name"], "config": str(path.relative_to(ROOT)), "note": light["note"]})
    print("wrote", path)

    for tag, med, rdir, tmpl, constr in [
        ("fullgrid_heavy_multibin_config_exp", "heavy", "data/qedark_rates/Si/heavy/long_scan",
         "dRdE_{material}_{mediator}_m{mchi_MeV}_s{sigma_e_cm2}.csv",
         {"tau": False, "n_bins": 4450, "mass": None, "lt": None}),
        ("fullgrid_heavy_tau_config_exp", "heavy", "data/qedark_rates/Si/heavy/long_scan",
         "dRdE_{material}_{mediator}_m{mchi_MeV}_s{sigma_e_cm2}.csv",
         {"tau": True, "n_bins": 1, "mass": None, "lt": None}),
        ("fullgrid_light_tau_config_exp", "massless", "data/qedark_rates/Si/ultralight/long_scan",
         "dRdE_{material}_massless_m{mchi_MeV}_s{sigma_e_cm2}.csv",
         {"tau": True, "n_bins": 1, "mass": None, "lt": None}),
    ]:
        v = {
            "name": tag,
            "label": tag,
            "outdir": f"outputs/qedark_repro_matrix/{tag}",
            "constrain_use_tau_weighted": constr["tau"],
            "constrain_n_bins": constr["n_bins"],
            "mass_kg": constr["mass"],
            "livetime_days": constr["lt"],
            "note": f"Full grid {med}, config exposure.",
        }
        path = OUT_DIR / f"scan_qedark_repro_{tag}.json"
        cfg = make_config(base, v, med, rdir, tmpl)
        # restore full grid
        cfg["model"]["grid"]["mchi_MeV"] = base["model"]["grid"]["mchi_MeV"]
        path.write_text(json.dumps(cfg, indent=2) + "\n")
        manifest.append({"variant": tag, "config": str(path.relative_to(ROOT)), "note": v["note"]})
        print("wrote", path)

    (ROOT / "configs/qedark_repro_matrix_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"\nData exposure: {DATA_EXPOSURE_KG_YEAR:.6e} kg·year")
    print(f"Data mass: {DATA_MASS_KG:.6e} kg, livetime: {DATA_LIVETIME_DAYS:.3f} days")


if __name__ == "__main__":
    main()
