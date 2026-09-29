#!/usr/bin/env python3
"""Write Sr2Cb2Sd phase-1 rate-grid and limit-scan JSON configs."""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CFG = ROOT / "configs" / "Sr2Cb2Sd"
RATES_ROOT = "data/Sr2Cb2Sd/rates"
OUT_ROOT = "outputs/Sr2Cb2Sd"
# Scan JSONs live in configs/Sr2Cb2Sd/; app resolves efficiency_csv as config_dir/../path.
EFF_CSV_DAMIC = "../data/Efficiencies_patterns_Nsims1000000_DCTrue_alpha1.csv"

SCENARIOS = [
    {
        "tag": "si_ref_gap1p2",
        "gap_tag": "1p2",
        "gap_eV": 1.2,
        "p100k": "data/p100K_gap1p2_eh3p8.csv",
        "eh_eV": 3.8,
        "ion_scenario": "B-thresh",
        "label_suffix": "Si ref (1.2 eV)",
    },
    {
        "tag": "srcd_indirect_gap0p556",
        "gap_tag": "0p556",
        "gap_eV": 0.556,
        "p100k": "data/p100K_gap0p556_eh2p06.csv",
        "eh_eV": 2.06,
        "ion_scenario": "Klein",
        "label_suffix": "SrCd2Sb2 indirect (0.556 eV)",
        "dc_sensitivity": True,
    },
    {
        "tag": "srcd_direct_gap0p603",
        "gap_tag": "0p603",
        "gap_eV": 0.603,
        "p100k": "data/p100K_gap0p603_eh2p19.csv",
        "eh_eV": 2.19,
        "ion_scenario": "Klein",
        "label_suffix": "SrCd2Sb2 direct (0.603 eV)",
        "dc_sensitivity": True,
    },
]

MEDIATORS = ("heavy", "light")
DC_BASELINE_YR = 0.00365  # 1e-5 e-/pix/day
# Sr2Cb2Sd DC ladder (same rates; background-only change)
DC_TIERS = {
    "1x": {"lambda_yr": DC_BASELINE_YR, "suffix": "", "label": "1e-5 e-/pix/day"},
    "dc100x": {"lambda_yr": 0.365, "suffix": "_dc100x", "label": "1e-3 e-/pix/day (100x)"},
    "dc1000x": {"lambda_yr": 3.65, "suffix": "_dc1000x", "label": "1e-2 e-/pix/day (1000x)"},
}
MASS_KG = 1.0
ROI = list(range(1, 11))

# OSCURA: same Si_fast gap 1.2 rates as si_ref; different exposure / backgrounds / ROI / efficiency.
OSCURA_MASS_KG = 30.0
OSCURA_DC_YR = 3.6525e-4  # 1e-6 e-/pix/day
OSCURA_FLAT_DRU = 0.01
OSCURA_ROI = list(range(2, 11))  # drop n_e = 1
OSCURA_EFF_CSV = "../data/Sr2Cb2Sd/efficiency_unity_single_pixel_ne2to10.csv"
SI_REF_GAP = SCENARIOS[0]


def rate_grid_config(mediator: str, gap_tag: str, gap_eV: float) -> dict:
    return {
        "_comment": f"Sr2Cb2Sd phase 1: Si_fast epsilon, {mediator} mediator, gap {gap_eV} eV.",
        "backend": "qcdark2",
        "material": "Si",
        "mediator": mediator,
        "halo": {"v0_kms": 238.0, "vE_kms": 263.0, "vesc_kms": 544.0},
        "detector": {
            "band_gap_eV": gap_eV,
            "eh_pair_eV": 3.8,
            "binsize_eV": 0.1,
        },
        "epsilon_h5": f"data/qcdark2_epsilon/Si/Si_fast_gap{gap_tag}.h5",
        "rates_dir": f"{RATES_ROOT}/Si/{mediator}/Si_fast_gap{gap_tag}",
        "filename_template": "dRdE_{material}_{mediator}_m{mchi_MeV}_s{sigma_e_cm2}.csv",
        "grid": {
            "mchi_MeV": {
                "logspace": {"start": 0.2, "stop": 1000.0, "num": 800, "endpoint": True}
            },
            "sigma_e_cm2": {
                "logspace": {
                    "start_exp": -46,
                    "stop_exp": -26,
                    "num": 300,
                    "endpoint": True,
                }
            },
            "format": {"mchi": ".6f", "sigma": ".1e"},
        },
        "options": {
            "skip_existing": True,
            "overwrite": False,
            "parallel": 8,
            "progress": True,
            "use_sigma_scaling": True,
            "sigma_ref": 1e-37,
            "format": {"mchi": ".6f", "sigma": ".1e"},
        },
    }


def oscura_scan_config(mediator: str) -> dict:
    """OSCURA projection: Si_fast 1.2 eV rates (reuse si_ref grid), 30 kg·yr, low DC / flat bkg."""
    sc = SI_REF_GAP
    gap_tag = sc["gap_tag"]
    gap_eV = sc["gap_eV"]
    run_tag = f"Sr2Cb2Sd_oscura_gap1p2_{mediator}_30kgy"
    return {
        "_comment": (
            "Sr2Cb2Sd OSCURA: Si_fast 1.2 eV (same rates as si_ref), 30 kg·yr, "
            "DC=1e-6 e-/pix/day, flat 0.01 d.r.u., ROI n_e 2-10 (ignore 1e), "
            "100% single-pixel efficiency (unity CSV + row_length=1)."
        ),
        "run": {
            "label": run_tag,
            "outdir": f"{OUT_ROOT}/{run_tag}",
            "cl": 0.9,
            "test_stat": "PLR",
            "n_toys": 0,
            "rng_seed": 12345,
            "verbosity": 1,
            "use_profile_likelihood": True,
            "constrain_prior_strength": 98,
            "constrain_use_gamma_sign": False,
            "constrain_use_tau_weighted": True,
            "constrain_n_bins": 1,
            "theta_lo": 0.5,
            "theta_hi": 10,
            "profile_minimizer": "pydme",
            "pydme_style_ul": True,
            "background_source": "dc_flat_migration",
            "background_model": "scale",
        },
        "detector": {
            "rows": 1300,
            "cols": 6300,
            "pixel_size_um": 15.0,
            "thickness_mm": 0.67,
            "active_fraction": 0.98,
            "target_element": "Si",
            "density_g_cm3": 2.329,
            "mass_kg": OSCURA_MASS_KG,
        },
        "experiment": {
            "mode": "asimov",
            "livetime_days": 365.25,
            "duty_cycle": 1.0,
            "binning": {"ne_min": 1, "ne_max": 20},
            "roi_bins": OSCURA_ROI,
            "pattern_roi": [11, 21, 111, 31, 22, 211],
            "observable_bins": "n_e",
        },
        "response": {
            "mode": "pattern",
            "analysis_space": "pattern",
            "charge_ionization": {
                "table_csv": sc["p100k"],
                "band_gap_eV": gap_eV,
                "eh_pair_eV": sc["eh_eV"],
                "scenario": sc["ion_scenario"],
            },
            "efficiency_mc": {
                "efficiency_csv": OSCURA_EFF_CSV,
                "n_events_per_ne": 1000,
                "sigma_readout_e": 0.16,
                "Qmin_e": 0.6,
                "Qmax_e": 20.0,
                "A_um2": 803.25,
                "b_umInv": 0.00065,
                "alpha": 1.0,
                "beta_per_keV": 0.0,
                "half_window_pix": 2,
                "row_length": 1,
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
                "lambda_dc": OSCURA_DC_YR,
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
        "backgrounds": {
            "dark_current": {
                "lambda_e_per_pix_per_year": OSCURA_DC_YR,
                "norm_scale": 1.0,
            },
            "pattern_efficiency": {"type": "flat", "epsilon": 1.0},
            "timing": {"exposure_time_s": 1800, "n_exposures": None},
            "flat_background": {
                "norm_per_kg_year_keV": OSCURA_FLAT_DRU,
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
            "rates_dir": f"{RATES_ROOT}/Si/{mediator}/Si_fast_gap{gap_tag}",
            "filename_template": "dRdE_{material}_{mediator}_m{mchi_MeV}_s{sigma_e_cm2}.csv",
            "Emin_eV": 2,
            "Emax_eV": 20,
            "nbins": 200,
            "example_point": {"mchi_MeV": 0.5, "sigma_e_cm2": 1e-27},
            "grid": {
                "mchi_MeV": {
                    "logspace": {
                        "start": 0.2,
                        "stop": 1000.0,
                        "num": 800,
                        "endpoint": True,
                    }
                },
                "sigma_e_cm2": {
                    "logspace": {
                        "start_exp": -46,
                        "stop_exp": -26,
                        "num": 300,
                        "endpoint": True,
                    }
                },
                "format": {"mchi": ".6f", "sigma": ".1e"},
            },
        },
    }


def scan_config(scenario: dict, mediator: str, *, dc_tier: str = "1x") -> dict:
    tag = scenario["tag"]
    gap_tag = scenario["gap_tag"]
    gap_eV = scenario["gap_eV"]
    dc = DC_TIERS[dc_tier]
    lambda_yr = dc["lambda_yr"]
    run_tag = f"Sr2Cb2Sd_{tag}_{mediator}_1kgy{dc['suffix']}"
    return {
        "_comment": (
            f"Sr2Cb2Sd phase 1: {scenario['label_suffix']}, {mediator} mediator, "
            f"1.0 kg·yr, Si_fast dielectric, DC={dc['label']}, flat 1 d.r.u., ROI n_e 1-10."
        ),
        "run": {
            "label": run_tag,
            "outdir": f"{OUT_ROOT}/{run_tag}",
            "cl": 0.9,
            "test_stat": "PLR",
            "n_toys": 0,
            "rng_seed": 12345,
            "verbosity": 1,
            "use_profile_likelihood": True,
            "constrain_prior_strength": 98,
            "constrain_use_gamma_sign": False,
            "constrain_use_tau_weighted": True,
            "constrain_n_bins": 1,
            "theta_lo": 0.5,
            "theta_hi": 10,
            "profile_minimizer": "pydme",
            "pydme_style_ul": True,
            "background_source": "dc_flat_migration",
            "background_model": "scale",
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
            "roi_bins": ROI,
            "pattern_roi": [11, 21, 111, 31, 22, 211],
            "observable_bins": "n_e",
        },
        "response": {
            "mode": "pattern",
            "analysis_space": "pattern",
            "charge_ionization": {
                "table_csv": scenario["p100k"],
                "band_gap_eV": gap_eV,
                "eh_pair_eV": scenario["eh_eV"],
                "scenario": scenario["ion_scenario"],
            },
            "efficiency_mc": {
                "efficiency_csv": EFF_CSV_DAMIC,
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
                "lambda_dc": lambda_yr,
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
        "backgrounds": {
            "dark_current": {
                "lambda_e_per_pix_per_year": lambda_yr,
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
            "rates_dir": f"{RATES_ROOT}/Si/{mediator}/Si_fast_gap{gap_tag}",
            "filename_template": "dRdE_{material}_{mediator}_m{mchi_MeV}_s{sigma_e_cm2}.csv",
            "Emin_eV": 2,
            "Emax_eV": 20,
            "nbins": 200,
            "example_point": {"mchi_MeV": 0.5, "sigma_e_cm2": 1e-27},
            "grid": {
                "mchi_MeV": {
                    "logspace": {
                        "start": 0.2,
                        "stop": 1000.0,
                        "num": 800,
                        "endpoint": True,
                    }
                },
                "sigma_e_cm2": {
                    "logspace": {
                        "start_exp": -46,
                        "stop_exp": -26,
                        "num": 300,
                        "endpoint": True,
                    }
                },
                "format": {"mchi": ".6f", "sigma": ".1e"},
            },
        },
    }


def main() -> None:
    CFG.mkdir(parents=True, exist_ok=True)
    for sc in SCENARIOS:
        for med in MEDIATORS:
            p = CFG / f"qcdark2_rates_{sc['tag']}_{med}.json"
            p.write_text(json.dumps(rate_grid_config(med, sc["gap_tag"], sc["gap_eV"]), indent=2) + "\n")
            print(f"wrote {p.relative_to(ROOT)}")
            tier_keys = list(DC_TIERS) if sc.get("dc_sensitivity") else ["1x"]
            for tier_key in tier_keys:
                tier = DC_TIERS[tier_key]
                scan_name = f"scan_{sc['tag']}_{med}_1kgy{tier['suffix']}.json"
                p = CFG / scan_name
                p.write_text(json.dumps(scan_config(sc, med, dc_tier=tier_key), indent=2) + "\n")
                print(f"wrote {p.relative_to(ROOT)}")
    for med in MEDIATORS:
        p = CFG / f"scan_oscura_gap1p2_{med}_30kgy.json"
        p.write_text(json.dumps(oscura_scan_config(med), indent=2) + "\n")
        print(f"wrote {p.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
