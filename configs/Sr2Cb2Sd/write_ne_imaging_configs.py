#!/usr/bin/env python3
"""Write Sr2Cb2Sd phase-2 n_e imaging one-point scan JSON configs."""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CFG = ROOT / "configs" / "Sr2Cb2Sd"
OUT_ROOT = "outputs/Sr2Cb2Sd/ne_imaging"
RATES_ROOT = "data/Sr2Cb2Sd/rates"
EFF_CSV = "../data/Efficiencies_patterns_Nsims1000000_DCTrue_alpha1.csv"

MASS_KG = 1.0
MCHI_MEV = 1.000194  # nearest Sr2Cb2Sd rate-grid mass to 1 MeV
SIGMA_CM2 = 1.1e-35

DC_BASELINE_YR = 0.00365  # 1e-5 e-/pix/day
DC_TIERS = {
    "1x": {"lambda_yr": DC_BASELINE_YR, "suffix": "1x", "dc_tag": "1e-5 e-/pix/day"},
    "dc100x": {"lambda_yr": 0.365, "suffix": "dc100x", "dc_tag": "1e-3 e-/pix/day (100x)"},
    "dc1000x": {"lambda_yr": 3.65, "suffix": "dc1000x", "dc_tag": "1e-2 e-/pix/day (1000x)"},
}

SCENARIOS = {
    "damic": {
        "tag": "si_ref_gap1p2",
        "gap_tag": "1p2",
        "gap_eV": 1.2,
        "eh_eV": 3.8,
        "p100k": "data/p100K_gap1p2_eh3p8.csv",
        "ion_scenario": "B-thresh",
        "label": "DAMIC-M Si ref",
    },
    "srcd_indirect": {
        "tag": "srcd_indirect_gap0p556",
        "gap_tag": "0p556",
        "gap_eV": 0.556,
        "eh_eV": 2.06,
        "p100k": "data/p100K_gap0p556_eh2p06.csv",
        "ion_scenario": "Klein",
        "label": "SrCd2Sb2 indirect",
    },
    "srcd_direct": {
        "tag": "srcd_direct_gap0p603",
        "gap_tag": "0p603",
        "gap_eV": 0.603,
        "eh_eV": 2.19,
        "p100k": "data/p100K_gap0p603_eh2p19.csv",
        "ion_scenario": "Klein",
        "label": "SrCd2Sb2 direct",
    },
}


def ne_imaging_config(
    scenario_key: str,
    mediator: str,
    *,
    dc_tier: str,
) -> dict:
    sc = SCENARIOS[scenario_key]
    dc = DC_TIERS[dc_tier]
    run_tag = f"Sr2Cb2Sd_ne_{scenario_key}_{mediator}_{dc['suffix']}"
    outdir = f"{OUT_ROOT}/{scenario_key}_{mediator}_{dc['suffix']}"
    return {
        "_comment": (
            f"Sr2Cb2Sd phase 2 n_e imaging: {sc['label']}, {mediator} mediator, "
            f"1.0 kg·yr, DAMIC pattern efficiency, DC={dc['dc_tag']}, flat 1 d.r.u., "
            f"observable_bins=ne, dump S_true/S_obs/B_tot at mchi={MCHI_MEV}, sigma={SIGMA_CM2}."
        ),
        "run": {
            "label": run_tag,
            "outdir": outdir,
            "n_toys": 0,
            "rng_seed": 12345,
            "verbosity": 1,
            "cl": 0.9,
            "test_stat": "PLR",
            "use_profile_likelihood": False,
            "data_path": "",
            "background_source": "dc_flat",
            "background_model": "scale",
            "profile_minimizer": "brent",
            "pydme_style_ul": False,
            "dump_point_spectra_root": True,
        },
        "detector": {
            "rows": 1300,
            "cols": 6300,
            "pixel_size_um": 15.0,
            "thickness_mm": 0.67,
            "active_fraction": 0.98,
            "target_element": "Si",
            "density_g_cm3": 2.329,
            "mass_kg": MASS_KG,
        },
        "experiment": {
            "mode": "asimov",
            "livetime_days": 365.25,
            "duty_cycle": 1.0,
            "binning": {"ne_min": 1, "ne_max": 20},
            "roi_bins": [1, 2, 3, 4, 5],
            "pattern_roi": [11, 21, 111, 31, 22, 211],
            "observable_bins": "ne",
        },
        "response": {
            "mode": "pattern",
            "analysis_space": "pattern",
            "charge_ionization": {
                "table_csv": sc["p100k"],
                "band_gap_eV": sc["gap_eV"],
                "eh_pair_eV": sc["eh_eV"],
                "scenario": sc["ion_scenario"],
            },
            "efficiency_mc": {
                "efficiency_csv": EFF_CSV,
                "n_events_per_ne": 1000000,
                "sigma_readout_e": 0.16,
                "Qmin_e": 0.6,
                "Qmax_e": 20.0,
                "A_um2": 803.25,
                "b_umInv": 0.00065,
                "alpha": 1.0,
                "beta_per_keV": 0.0,
                "half_window_pix": 2,
                "row_length": 32,
                "rows_bin": 1,
                "cols_bin": 1,
                "pileup_with_dc": False,
                "enable_MN": True,
                "enable_MNL": True,
                "use_2d_image_efficiency": False,
                "rng_seed": 987654321,
            },
            "pattern_image": {
                "nrows_binned": 3,
                "ncols": 50,
                "row_binning": 100,
                "col_binning": 1,
                "pixel_size_um": 15.0,
                "sigma_readout_e": 0.16,
                "lambda_dc": dc["lambda_yr"],
                "rng_seed": 987654321,
            },
            "pattern_classifier": {
                "max_e_per_pixel": 5,
                "Qmin_e": 0.6,
                "neighbor_Qmax_e": 0.6,
                "thr_M": 3.5,
                "thr_MN": 4.0,
                "thr_MNL": 5.5,
            },
            "pcd": {"q_min": 0.0, "q_max": 20.0, "nbins": 200, "mc_trials": 50000},
        },
        "background_source": "dc_flat",
        "backgrounds": {
            "dark_current": {
                "lambda_e_per_pix_per_year": dc["lambda_yr"],
                "norm_scale": 1.0,
            },
            "pattern_efficiency": {"type": "flat", "epsilon": 0.95},
            "timing": {"exposure_time_s": 1800, "n_exposures": None},
            "flat_background": {
                "norm_per_kg_year_keV": 1.0,
                "Emin_eV": 2,
                "Emax_eV": 20,
                "nbins": 200,
            },
            "background_efficiency_csv": "",
        },
        "model": {
            "type": "dm_electron",
            "material": "Si",
            "mediator": mediator,
            "rates_dir": f"{RATES_ROOT}/Si/{mediator}/Si_fast_gap{sc['gap_tag']}",
            "filename_template": "dRdE_{material}_{mediator}_m{mchi_MeV}_s{sigma_e_cm2}.csv",
            "Emin_eV": 2,
            "Emax_eV": 20,
            "nbins": 200,
            "example_point": {"mchi_MeV": MCHI_MEV, "sigma_e_cm2": SIGMA_CM2},
            "grid": {
                "mchi_MeV": {"values": [MCHI_MEV]},
                "sigma_e_cm2": {"values": [SIGMA_CM2]},
                "format": {"mchi": ".6f", "sigma": ".1e"},
            },
        },
    }


def main() -> None:
    CFG.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for med in ("heavy", "light"):
        for dc_tier in ("1x", "dc100x", "dc1000x"):
            p = CFG / f"ne_imaging_damic_{med}_{dc_tier}.json"
            p.write_text(json.dumps(ne_imaging_config("damic", med, dc_tier=dc_tier), indent=2) + "\n")
            paths.append(p)
        for scen in ("srcd_indirect", "srcd_direct"):
            p = CFG / f"ne_imaging_{scen}_{med}_1x.json"
            p.write_text(json.dumps(ne_imaging_config(scen, med, dc_tier="1x"), indent=2) + "\n")
            paths.append(p)
    print(f"Wrote {len(paths)} ne_imaging configs under {CFG}/")
    for p in paths:
        print(f"  {p.name}")


if __name__ == "__main__":
    main()
