# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: entry.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  entry.py -- EXCEED-DM HDF5 -> CCDarkSens CSV converter for intermediate-
#  mass mediators.
# ============================================================================

"""
EXCEED-DM HDF5 -> CCDarkSens CSV converter for intermediate-mass mediators.

Reads an EXCEED-DM binned_scatter_rate HDF5 output (produced with explicit mA
values) and writes one CSV per (mA', mchi, sigma) combination in the standard
CCDarkSens rate-table format:
    data/exdm_rates/Si/intermediate/mA_{X}keV/dRdE_mchi_{M}MeV_sigma_e_{S}.csv

The rate at each sigma is obtained by linearly scaling the reference rate
(computed at sigma_ref by EXCEED-DM) by sigma/sigma_ref.  This avoids
re-running EXCEED-DM for each cross section value.

Usage:
    python -m ccdarkphys.exdm.entry \\
        --hdf5   /path/to/EXDM_out.hdf5 \\
        --outdir data/exdm_rates/Si/intermediate \\
        [--sigma-ref  1e-36] \\
        [--sigma-list 1e-35 1e-36 1e-37 1e-38 1e-39 1e-40 1e-41 1e-42] \\
        [--exdm-utils /path/to/EXCEED-DM/utilities/output_parser]
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

from ccdarkphys.common import io as CIO

DEFAULT_SIGMA_LIST = [1e-35, 1e-36, 1e-37, 1e-38, 1e-39, 1e-40, 1e-41, 1e-42]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load_exdm_data(hdf5_path: str, exdm_utils_dir: str | None):
    """Import EXDMData from EXCEED-DM utilities and open the HDF5 file."""
    if exdm_utils_dir:
        if exdm_utils_dir not in sys.path:
            sys.path.insert(0, exdm_utils_dir)
    try:
        from EXDMDataHandler import EXDMData
    except ImportError:
        raise ImportError(
            "EXDMDataHandler not found. Pass --exdm-utils pointing to "
            "EXCEED-DM/utilities/output_parser/, or add it to PYTHONPATH."
        )
    return EXDMData(hdf5_path)


def _mA_label(mA_eV: float) -> str:
    """Format mA' in eV as a compact keV label for directory names."""
    keV = mA_eV / 1e3
    if keV == int(keV):
        return f"{int(keV)}keV"
    return f"{keV:g}keV"


def _mchi_str(mchi_MeV: float) -> str:
    """Format mchi in MeV with 6 decimal places — matches C++ format_mchi_6f."""
    return f"{mchi_MeV:.6f}"


def _sigma_str(sigma_cm2: float) -> str:
    """Format sigma with 1 decimal place scientific — matches scan app .1e default."""
    return f"{sigma_cm2:.1e}"


def _prepend_subgap_zeros(E_bins: list, rate: np.ndarray, band_gap_eV: float
                          ) -> tuple[np.ndarray, np.ndarray]:
    """
    Prepend two zero-rate points below the band gap so that RateTable's
    front-clamp (Ec <= E_eV_.front() -> first CSV value) returns 0 rather
    than leaking sub-gap rate into the signal integral.
    """
    pad_E = np.array([0.0, band_gap_eV - 0.01])
    pad_R = np.zeros(2)
    E_arr = np.concatenate([pad_E, np.array(E_bins)])
    R_arr = np.concatenate([pad_R, rate])
    return E_arr, R_arr


# ---------------------------------------------------------------------------
# Core conversion
# ---------------------------------------------------------------------------

def convert(hdf5_path: str,
            out_dir: str,
            sigma_ref_cm2: float = 1e-36,
            sigma_list: list[float] | None = None,
            exdm_utils_dir: str | None = None,
            verbose: bool = True) -> None:
    """
    Read EXCEED-DM HDF5 and write one CSV per (mA', mchi, sigma).

    Output tree:
        out_dir/mA_{label}/dRdE_mchi_{M}MeV_sigma_e_{S}.csv

    The rate at each sigma is sigma/sigma_ref * rate_ref, where rate_ref is
    extracted from the HDF5 at sigma_ref.
    """
    if sigma_list is None:
        sigma_list = DEFAULT_SIGMA_LIST

    d = _load_exdm_data(hdf5_path, exdm_utils_dir)

    mA_list_eV  = d.get_mediator_masses_eV()
    mX_list_MeV = d.get_masses_MeV()
    E_bin_width = d.get_numerics_scatter_binned_rate_E_bin_width()
    band_gap    = d.get_material_band_gap()

    if len(mA_list_eV) == 0:
        raise ValueError("HDF5 contains no mA values — was it run with mA specified?")

    n_files = len(mA_list_eV) * len(mX_list_MeV) * len(sigma_list)

    if verbose:
        print(f"EXCEED-DM HDF5  : {hdf5_path}")
        print(f"mA' values (keV): {mA_list_eV / 1e3}")
        print(f"mchi range (MeV): [{mX_list_MeV.min():.3g}, {mX_list_MeV.max():.3g}]  ({len(mX_list_MeV)} masses)")
        print(f"sigma values    : {[f'{s:.0e}' for s in sigma_list]}")
        print(f"sigma_ref       : {sigma_ref_cm2:.2e} cm^2")
        print(f"Total CSVs      : {n_files}")
        print()

    for mA_eV in mA_list_eV:
        mA_label = _mA_label(mA_eV)
        mA_out_dir = os.path.join(out_dir, f"mA_{mA_label}")
        os.makedirs(mA_out_dir, exist_ok=True)

        for mchi_MeV in mX_list_MeV:
            # Fetch reference rate once per (mA', mchi)
            E_bins, rate_ref = d.get_binned_scatter_rate_E(
                mass_MeV    = mchi_MeV,
                med_FF      = 2,        # ignored when mass_A != 0
                mass_A      = mA_eV,
                sigma_cm2   = sigma_ref_cm2,
                E_bin_width = E_bin_width,
            )

            # dR/dE at reference sigma [events/kg/yr/eV]
            dRdE_ref = np.array(rate_ref) / E_bin_width

            # Prepend sub-gap zeros (done once; shape is shared across sigmas)
            E_arr, R_ref_arr = _prepend_subgap_zeros(E_bins, dRdE_ref, band_gap)

            mchi_label = _mchi_str(mchi_MeV)

            for sigma in sigma_list:
                scale = sigma / sigma_ref_cm2
                R_arr = R_ref_arr * scale

                filename = f"dRdE_mchi_{mchi_label}MeV_sigma_e_{_sigma_str(sigma)}.csv"
                out_path = os.path.join(mA_out_dir, filename)

                header = [
                    f"# Differential Rates computed with CCDarkSens (EXCEED-DM entry)",
                    f"# material = Si, mediator = intermediate, mA_eV = {mA_eV:.6g}",
                    f"# hdf5 = {os.path.abspath(hdf5_path)}",
                    f"# mX (eV) = {mchi_MeV * 1e6:.6g}",
                    f"# sigma_e (cm^2) = {sigma:.6e}",
                    f"# Output units: dR/dE in events / kg / year / eV",
                    f"# Columns: E (eV), dRdE (events/kg/year/eV)",
                ]
                CIO.write_csv_generic(out_path, E_arr, R_arr, header)

        if verbose:
            print(f"  mA' = {mA_label:8s}  -> {mA_out_dir}  ({len(mX_list_MeV) * len(sigma_list)} files)")

    if verbose:
        print("\nDone.")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# main
#   Command-line front end of the EXCEED-DM converter: reads the HDF5, the output directory, the reference and target cross sections, and calls convert().
# ----------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--hdf5",       required=True, help="Path to EXCEED-DM output HDF5")
    p.add_argument("--outdir",     required=True, help="Root output directory for CSV files")
    p.add_argument("--sigma-ref",  type=float, default=1e-36,
                   help="Reference sigma used in EXCEED-DM run (default: 1e-36)")
    p.add_argument("--sigma-list", type=float, nargs="+",
                   default=DEFAULT_SIGMA_LIST,
                   help="Sigma values to generate CSVs for (default: 1e-35 to 1e-42)")
    p.add_argument("--exdm-utils", default=None,
                   help="Path to EXCEED-DM/utilities/output_parser/ if not on PYTHONPATH")
    p.add_argument("--quiet",      action="store_true")
    args = p.parse_args()

    convert(
        hdf5_path      = args.hdf5,
        out_dir        = args.outdir,
        sigma_ref_cm2  = args.sigma_ref,
        sigma_list     = args.sigma_list,
        exdm_utils_dir = args.exdm_utils,
        verbose        = not args.quiet,
    )


if __name__ == "__main__":
    main()
