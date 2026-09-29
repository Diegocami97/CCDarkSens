#!/usr/bin/env python3
"""Generate band_gap_one_point_gap*.json configs.

Two scenarios per gap value:
  B-thresh  — eh_pair_eV = 3.8 (silicon, standard)
  D-equal   — eh_pair_eV = band_gap_eV (eh pair energy equals band gap)

Usage:
    python3 configs/generate_band_gap_configs.py
"""

import json, pathlib

GAPS = [0.1, 0.3, 0.5, 0.7, 0.9, 1.2]

SCENARIOS = {
    "B-thresh": {"eh_pair_eV": 3.8},
    "D-equal":  {"eh_pair_eV": None},  # filled per gap: eh_pair_eV = band_gap_eV
}

HERE = pathlib.Path(__file__).parent


def gap_str(g: float) -> str:
    return f"{g:.1f}".replace(".", "p")


def make_config(gap: float, scenario: str, eh: float) -> dict:
    gs = gap_str(gap)
    label = f"gap{gs}_{scenario}"
    return {
        "_comment": f"One-point spectra dump (n_e space): gap {gap} eV, eh={eh} eV ({scenario}). "
                    "dRdE + S_true(n_e); no pattern-space observable.",
        "run": {
            "label": label,
            "outdir": f"outputs/band_gap_one_point_spectra/{label}",
            "n_toys": 0,
            "rng_seed": 12345,
            "verbosity": 1,
            "cl": 0.9,
            "test_stat": "PLR",
            "dump_point_spectra_root": True,
            "use_profile_likelihood": False,
            "data_path": "",
            "background_source": "dc_flat",
            "background_model": "scale",
            "profile_minimizer": "brent",
            "pydme_style_ul": False,
        },
        "detector": {
            "rows": 1300,
            "cols": 6300,
            "pixel_size_um": 15.0,
            "thickness_mm": 0.67,
            "active_fraction": 0.98,
            "target_element": "Si",
            "density_g_cm3": 2.329,
            "mass_kg": 0.5,
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
            "efficiency_mc": {
                "efficiency_csv": "data/efficiencies_paolo.csv",
                "efficiency_csv_reference": "data/efficiencies_paolo.csv",
                "n_events_per_ne": 1000000,
                "sigma_readout_e": 0.16,
                "Qmin_e": 0.6,
                "Qmax_e": 20.0,
                "A_um2": 803.25,
                "b_umInv": 0.00065,
                "alpha": 1.0,
                "beta_per_keV": 0.0,
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
                "lambda_dc": 0.00041,
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
            "charge_ionization": {
                "table_csv": f"data/p100K_gap{gs}_eh{str(eh).replace('.', 'p')}.csv",
                "band_gap_eV": gap,
                "eh_pair_eV": eh,
                "scenario": scenario,
            },
        },
        "background_source": "dc_flat",
        "backgrounds": {
            "dark_current": {
                "lambda_e_per_pix_per_year": 0.0365,
                "norm_scale": 1.0,
            },
            "pattern_efficiency": {"type": "flat", "epsilon": 0.95},
            "timing": {"exposure_time_s": 1800, "n_exposures": None},
            "flat_background": {
                "norm_per_kg_year_keV": 1.0,
                "Emin_eV": 0,
                "Emax_eV": 20,
                "nbins": 200,
            },
            "background_efficiency_csv": "",
        },
        "model": {
            "type": "dm_electron",
            "material": "Si",
            "mediator": "heavy",
            "rates_dir": f"data/qcdark2_rates/Si/heavy/Si_fast_gap{gs}",
            "filename_template": "dRdE_{material}_{mediator}_m{mchi_MeV}_s{sigma_e_cm2}.csv",
            "Emin_eV": 0,
            "Emax_eV": 20,
            "nbins": 200,
            "example_point": {"mchi_MeV": 0.5, "sigma_e_cm2": 1e-27},
            "grid": {
                "mchi_MeV": {"values": [10.800876]},
                "sigma_e_cm2": {"values": [1.1e-35]},
                "format": {"mchi": ".6f", "sigma": ".1e"},
            },
        },
    }


if __name__ == "__main__":
    for gap in GAPS:
        gs = gap_str(gap)
        for scenario, opts in SCENARIOS.items():
            eh = opts["eh_pair_eV"] if opts["eh_pair_eV"] is not None else gap
            cfg = make_config(gap, scenario, eh)
            out = HERE / f"band_gap_one_point_gap{gs}_{scenario}.json"
            out.write_text(json.dumps(cfg, indent=2) + "\n")
            print(f"wrote {out.name}")
