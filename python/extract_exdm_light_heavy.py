# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: extract_exdm_light_heavy.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  extract_exdm_light_heavy.py -- Extract EXCEED-DM light and heavy (mA'=100
#  keV) rate CSVs from the Stratman & Trickle 2026 Zenodo HDF5, for use as
#  reference light/heavy curves from the same DFT dataset as the intermediate
#  mediator scans.
# ============================================================================

"""
Extract EXCEED-DM light and heavy (mA'=100 keV) rate CSVs from the
Stratman & Trickle 2026 Zenodo HDF5, for use as reference light/heavy curves
from the same DFT dataset as the intermediate mediator scans.

Usage:
    python python/extract_exdm_light_heavy.py
"""
import sys
import os
import numpy as np

sys.path.insert(0, os.path.expanduser(
    "~/Documents/Software/intermediate_mass_mediators/output_parser"))
sys.path.insert(0, os.path.dirname(__file__))

from EXDMDataHandler import EXDMData
from ccdarkphys.common import io as CIO

HDF5 = os.path.expanduser(
    "~/Documents/Software/intermediate_mass_mediators/data/EXDM_out_mX_range_Si.hdf5")

SIGMA_REF  = 1e-36   # cm^2 — reference for scaling
SIGMA_LIST = [1e-35, 1e-36, 1e-37, 1e-38, 1e-39, 1e-40, 1e-41, 1e-42]

OUT_BASE = "data/exdm_rates/Si"
LIGHT_DIR = os.path.join(OUT_BASE, "light_exdm")
HEAVY_DIR = os.path.join(OUT_BASE, "heavy_exdm")   # mA'=100 keV

MA_HEAVY_EV = 100e3  # 100 keV in eV — effectively heavy limit for Si (>>  alpha*me=3.73 keV)


# ----------------------------------------------------------------------------
# write_csvs
#   For one mass: pad the spectrum with zeros below the band gap, rescale the reference rate to every cross section
#   in sigma_list (the rate is proportional to sigma) and write one CSV per cross section.
# ----------------------------------------------------------------------------
def write_csvs(E_bins, dRdE_ref, band_gap, mchi_MeV, sigma_list, sigma_ref, out_dir, label):
    pad_E = np.array([0.0, band_gap - 0.01])
    pad_R = np.zeros(2)
    E_arr = np.concatenate([pad_E, np.array(E_bins)])
    R_ref = np.concatenate([pad_R, dRdE_ref])

    mchi_str = f"{mchi_MeV:.6f}"
    for sigma in sigma_list:
        scale = sigma / sigma_ref
        R_arr = R_ref * scale
        fname = f"dRdE_mchi_{mchi_str}MeV_sigma_e_{sigma:.1e}.csv"
        header = [
            f"# EXCEED-DM {label} mediator, Si",
            f"# hdf5 = {os.path.abspath(HDF5)}",
            f"# mX (eV) = {mchi_MeV*1e6:.6g}",
            f"# sigma_e (cm^2) = {sigma:.6e}",
            f"# Output units: dR/dE in events / kg / year / eV",
            f"# Columns: E (eV), dRdE (events/kg/year/eV)",
        ]
        CIO.write_csv_generic(os.path.join(out_dir, fname), E_arr, R_arr, header)


# ----------------------------------------------------------------------------
# main
#   Read every mass of the EXCEED-DM HDF5 and write the light (mA' = 0) and heavy (mA' = 100 keV) rate CSVs for all cross sections.
# ----------------------------------------------------------------------------
def main():
    d = EXDMData(HDF5)
    masses_MeV = d.get_masses_MeV()
    band_gap   = d.get_material_band_gap()
    dE         = d.get_numerics_scatter_binned_rate_E_bin_width()

    os.makedirs(LIGHT_DIR, exist_ok=True)
    os.makedirs(HEAVY_DIR, exist_ok=True)

    print(f"Masses (MeV): {len(masses_MeV)} values  {masses_MeV[0]:.3f} – {masses_MeV[-1]:.3f}")
    print(f"Band gap: {band_gap} eV,  E bin width: {dE} eV")
    print(f"Extracting {len(masses_MeV)} masses × {len(SIGMA_LIST)} sigmas each for light and heavy")
    print()

    for mchi in masses_MeV:
        # --- light (med_FF=2, mA=0) ---
        E_bins, rate_ref = d.get_binned_scatter_rate_E(
            mass_MeV=mchi, med_FF=2, mass_A=0,
            sigma_cm2=SIGMA_REF, E_bin_width=dE)
        dRdE_ref = np.array(rate_ref) / dE
        write_csvs(E_bins, dRdE_ref, band_gap, mchi, SIGMA_LIST, SIGMA_REF, LIGHT_DIR, "light")

        # --- heavy (mA'=100 keV >> alpha*me=3.73 keV) ---
        E_bins, rate_ref = d.get_binned_scatter_rate_E(
            mass_MeV=mchi, med_FF=2, mass_A=MA_HEAVY_EV,
            sigma_cm2=SIGMA_REF, E_bin_width=dE)
        dRdE_ref = np.array(rate_ref) / dE
        write_csvs(E_bins, dRdE_ref, band_gap, mchi, SIGMA_LIST, SIGMA_REF, HEAVY_DIR, "heavy_100keV")

    n = len(masses_MeV) * len(SIGMA_LIST)
    print(f"  Light -> {LIGHT_DIR}  ({n} files)")
    print(f"  Heavy -> {HEAVY_DIR}  ({n} files)")
    print("Done.")


if __name__ == "__main__":
    main()
