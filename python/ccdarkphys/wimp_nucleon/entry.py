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
#  entry.py -- CSV entry point for spin-independent WIMP-nucleus elastic
#  scattering rates.
# ============================================================================

"""
CSV entry point for spin-independent WIMP-nucleus elastic scattering rates.

Two outputs:
  - compute_dRdE: raw rate on the NUCLEAR RECOIL energy axis (E_R, not
    electron-equivalent). Unchanged since Phase 1.
  - compute_dRdE_ee: the same rate, quenched to electron-equivalent energy
    (E_ee) via a NuclearQuenching model (see ccdarkphys.wimp_nucleon.quenching)
    and re-histogrammed. This is what downstream detector-response work
    needs; compute_dRdE stays available on its own for debugging/audits.

Both use the column convention (E in eV, dRdE in events/kg/year/eV) that
matches every other CCDarkSens rate CSV, but the physical meaning of the
energy axis differs between the two -- documented explicitly in each CSV's
header (E_R vs E_ee).

Public:
    compute_dRdE(A, mchi_MeV, sigma_n_cm2, target_nucleus="Si28", mediator="heavy",
                 nr_Emin_keV=0.001, nr_Emax_keV=30.0, nr_nbins=3000,
                 rho_chi_gev_cm3=None, v0_kms=None, vE_kms=None, vesc_kms=None)
      -> dict(E_eV, dRdE_kg_year_eV, meta)      # E_R axis

    compute_dRdE_ee(..., ee_Emin_eV=0.0, ee_Emax_eV=8000.0, ee_nbins=800,
                     quenching_model="lindhard" | "chavarria_table" | "julian_table")
      -> dict(E_eV, dRdE_kg_year_eV, meta)      # E_ee axis

CLI examples:
    # raw E_R
    PYTHONPATH=python python3 -m ccdarkphys.wimp_nucleon.entry \\
      --target_nucleus Si28 --A 28 --mchi_MeV 5000 --sigma_n_cm2 1e-40 \\
      --out_csv data/wimp_nucleon_rates/Si/heavy/dRdE_Si28_heavy_m5000.000000_s1.0e-40.csv

    # quenched E_ee
    PYTHONPATH=python python3 -m ccdarkphys.wimp_nucleon.entry \\
      --target_nucleus Si28 --A 28 --mchi_MeV 5000 --sigma_n_cm2 1e-40 \\
      --quenching_model chavarria_table \\
      --out_csv data/wimp_nucleon_rates_ee/Si/heavy/dRdE_ee_Si28_heavy_m5000.000000_s1.0e-40.csv
"""
from __future__ import annotations

import argparse

import numpy as np

from ccdarkphys.common import io as CIO
from ccdarkphys.wimp_nucleon.rate import dRdE_nr_kg_day_keV
from ccdarkphys.wimp_nucleon.halo import (
    RHO_CHI_GEV_CM3_DEFAULT,
    V0_KMS_DEFAULT,
    VE_KMS_DEFAULT,
    VESC_KMS_DEFAULT,
)
from ccdarkphys.wimp_nucleon.quenching import QUENCHING_MODELS

_DAY_TO_YEAR = 365.25
_KEV_TO_EV = 1000.0


# ----------------------------------------------------------------------------
# _make_recoil_grid_keV
#   Bin centres of a uniform nuclear-recoil grid in keV.
# ----------------------------------------------------------------------------
def _make_recoil_grid_keV(Emin_keV: float, Emax_keV: float, nbins: int):
    edges = np.linspace(Emin_keV, Emax_keV, nbins + 1)
    return 0.5 * (edges[:-1] + edges[1:])


def compute_dRdE(
    A: float,
    mchi_MeV: float,
    sigma_n_cm2: float,
    target_nucleus: str = "Si28",
    mediator: str = "heavy",
    nr_Emin_keV: float = 0.001,
    nr_Emax_keV: float = 30.0,
    nr_nbins: int = 3000,
    *,
    rho_chi_gev_cm3: float | None = None,
    v0_kms: float | None = None,
    vE_kms: float | None = None,
    vesc_kms: float | None = None,
) -> dict:
    """
    Same outward contract as ccdarkphys.migdal.entry.compute_dRdE, but the
    energy axis is nuclear recoil energy E_R, not electronic omega. Halo
    parameters default to Baxter et al. 2021 (matching ccdarkphys.migdal's
    darkelf calls) when not explicitly overridden.
    """
    rho_chi_gev_cm3 = RHO_CHI_GEV_CM3_DEFAULT if rho_chi_gev_cm3 is None else rho_chi_gev_cm3
    v0_kms = V0_KMS_DEFAULT if v0_kms is None else v0_kms
    vE_kms = VE_KMS_DEFAULT if vE_kms is None else vE_kms
    vesc_kms = VESC_KMS_DEFAULT if vesc_kms is None else vesc_kms

    mchi_GeV = mchi_MeV * 1.0e-3
    E_R_keV = _make_recoil_grid_keV(nr_Emin_keV, nr_Emax_keV, nr_nbins)

    dRdE_kg_day_keV = dRdE_nr_kg_day_keV(
        E_R_keV, mchi_GeV, sigma_n_cm2, A,
        rho_chi_gev_cm3=rho_chi_gev_cm3, v0_kms=v0_kms, vE_kms=vE_kms, vesc_kms=vesc_kms,
    )
    dRdE_kg_year_eV = dRdE_kg_day_keV * _DAY_TO_YEAR / _KEV_TO_EV
    E_eV = E_R_keV * _KEV_TO_EV

    meta = {
        "target_nucleus": target_nucleus,
        "A": float(A),
        "mediator": mediator,
        "mchi_MeV": float(mchi_MeV),
        "sigma_n_cm2": float(sigma_n_cm2),
        "nr_Emin_keV": float(nr_Emin_keV),
        "nr_Emax_keV": float(nr_Emax_keV),
        "nr_nbins": int(nr_nbins),
        "rho_chi_gev_cm3": float(rho_chi_gev_cm3),
        "v0_kms": float(v0_kms),
        "vE_kms": float(vE_kms),
        "vesc_kms": float(vesc_kms),
    }
    return {
        "E_eV": np.asarray(E_eV, dtype=float),
        "dRdE_kg_year_eV": np.asarray(dRdE_kg_year_eV, dtype=float),
        "meta": meta,
    }


def compute_dRdE_ee(
    A: float,
    mchi_MeV: float,
    sigma_n_cm2: float,
    target_nucleus: str = "Si28",
    mediator: str = "heavy",
    nr_Emin_keV: float = 0.001,
    nr_Emax_keV: float = 30.0,
    nr_nbins: int = 3000,
    ee_Emin_eV: float = 0.0,
    ee_Emax_eV: float = 8000.0,
    ee_nbins: int = 800,
    quenching_model: str = "lindhard",
    *,
    rho_chi_gev_cm3: float | None = None,
    v0_kms: float | None = None,
    vE_kms: float | None = None,
    vesc_kms: float | None = None,
) -> dict:
    """
    dR/dE_ee -- the nuclear recoil rate (compute_dRdE, unchanged) remapped
    through a quenching model into electron-equivalent energy.

    Uses the same weighted-histogram technique WIMPyCCD's own notebook uses
    (rate_with_NR.ipynb): bin the E_R rate into counts, map each bin's
    energy through the yield function, re-histogram in E_ee. No Jacobian
    needed -- this is a change of variables on binned counts, not a density
    remap.

    quenching_model: "lindhard" | "chavarria_table" | "julian_table" -- see
    ccdarkphys.wimp_nucleon.quenching for what each one is.
    """
    if quenching_model not in QUENCHING_MODELS:
        raise ValueError(f"unknown quenching_model {quenching_model!r}; choose from {list(QUENCHING_MODELS)}")
    yield_fn = QUENCHING_MODELS[quenching_model]

    res_nr = compute_dRdE(
        A=A, mchi_MeV=mchi_MeV, sigma_n_cm2=sigma_n_cm2,
        target_nucleus=target_nucleus, mediator=mediator,
        nr_Emin_keV=nr_Emin_keV, nr_Emax_keV=nr_Emax_keV, nr_nbins=nr_nbins,
        rho_chi_gev_cm3=rho_chi_gev_cm3, v0_kms=v0_kms, vE_kms=vE_kms, vesc_kms=vesc_kms,
    )
    E_r_keV = res_nr["E_eV"] / _KEV_TO_EV
    bin_width_eV = (nr_Emax_keV - nr_Emin_keV) * _KEV_TO_EV / nr_nbins
    counts_kg_year = res_nr["dRdE_kg_year_eV"] * bin_width_eV  # density -> per-bin counts

    E_ee_keV = yield_fn(E_r_keV)
    E_ee_eV = E_ee_keV * _KEV_TO_EV

    edges = np.linspace(ee_Emin_eV, ee_Emax_eV, ee_nbins + 1)
    hist_counts, _ = np.histogram(E_ee_eV, bins=edges, weights=counts_kg_year)
    ee_bin_width_eV = (ee_Emax_eV - ee_Emin_eV) / ee_nbins
    dRdE_ee_kg_year_eV = hist_counts / ee_bin_width_eV
    E_ee_centers = 0.5 * (edges[:-1] + edges[1:])

    meta = dict(res_nr["meta"])
    meta.update({
        "quenching_model": quenching_model,
        "ee_Emin_eV": float(ee_Emin_eV),
        "ee_Emax_eV": float(ee_Emax_eV),
        "ee_nbins": int(ee_nbins),
    })
    return {
        "E_eV": np.asarray(E_ee_centers, dtype=float),
        "dRdE_kg_year_eV": np.asarray(dRdE_ee_kg_year_eV, dtype=float),
        "meta": meta,
    }


# ----------------------------------------------------------------------------
# _header_lines_ee
#   Comment header of a WIMP-nucleon CSV quenched to electron-equivalent energy: quenching model, nucleus, mass, cross section, both grids, halo and units.
# ----------------------------------------------------------------------------
def _header_lines_ee(meta: dict) -> list:
    return [
        "# Differential Rates computed with CCDarkSens (wimp_nucleon entry) -- "
        "spin-independent WIMP-nucleus elastic scattering, quenched to electron-equivalent energy",
        f"# quenching_model = {meta['quenching_model']} -- see "
        "python/ccdarkphys/wimp_nucleon/quenching.py for the model and its extrapolation conventions.",
        f"# target_nucleus = {meta['target_nucleus']}, A = {meta['A']}, mediator = {meta['mediator']} "
        "(informational only -- SI is a contact interaction, no light/heavy distinction here)",
        f"# mchi (MeV) = {meta['mchi_MeV']}",
        f"# sigma_n (cm^2, DM-nucleon-labeled -- see rate.py docstring re: reduced-mass convention) "
        f"= {meta['sigma_n_cm2']}",
        f"# nr_Emin_keV = {meta['nr_Emin_keV']}, nr_Emax_keV = {meta['nr_Emax_keV']}, "
        f"nr_nbins = {meta['nr_nbins']}  (raw recoil-energy binning, pre-quenching)",
        f"# ee_Emin_eV = {meta['ee_Emin_eV']}, ee_Emax_eV = {meta['ee_Emax_eV']}, "
        f"ee_nbins = {meta['ee_nbins']}  (output electron-equivalent binning)",
        f"# halo: rho_chi = {meta['rho_chi_gev_cm3']} GeV/cm^3, v0 = {meta['v0_kms']} km/s, "
        f"vE = {meta['vE_kms']} km/s, vesc = {meta['vesc_kms']} km/s",
        "# Rate formula ported from WIMPyCCD (analysis/rate_models.py); see "
        "python/ccdarkphys/wimp_nucleon/rate.py docstring for open parity/convention caveats.",
        "# Output units: dR/dE_ee in events / kg / year / eV (electron-equivalent axis)",
        "# Columns: E_ee (eV), dRdE (events/kg/year/eV)",
    ]


# ----------------------------------------------------------------------------
# _header_lines
#   Comment header of a WIMP-nucleon CSV on the raw nuclear-recoil axis E_R (no quenching applied), with a warning about that.
# ----------------------------------------------------------------------------
def _header_lines(meta: dict) -> list:
    return [
        "# Differential Rates computed with CCDarkSens (wimp_nucleon entry) -- "
        "spin-independent WIMP-nucleus elastic scattering",
        "# WARNING: energy axis is E_R (nuclear recoil energy), NOT electron-equivalent -- "
        "quenching to E_ee is a separate downstream step (NuclearQuenching), not applied here.",
        f"# target_nucleus = {meta['target_nucleus']}, A = {meta['A']}, mediator = {meta['mediator']} "
        "(informational only -- SI is a contact interaction, no light/heavy distinction here)",
        f"# mchi (MeV) = {meta['mchi_MeV']}",
        f"# sigma_n (cm^2, DM-nucleon-labeled -- see rate.py docstring re: reduced-mass convention) "
        f"= {meta['sigma_n_cm2']}",
        f"# nr_Emin_keV = {meta['nr_Emin_keV']}, nr_Emax_keV = {meta['nr_Emax_keV']}, "
        f"nr_nbins = {meta['nr_nbins']}",
        f"# halo: rho_chi = {meta['rho_chi_gev_cm3']} GeV/cm^3, v0 = {meta['v0_kms']} km/s, "
        f"vE = {meta['vE_kms']} km/s, vesc = {meta['vesc_kms']} km/s",
        "# Rate formula ported from WIMPyCCD (analysis/rate_models.py); see "
        "python/ccdarkphys/wimp_nucleon/rate.py docstring for open parity/convention caveats.",
        "# Output units: dR/dE_R in events / kg / year / eV (E_R axis)",
        "# Columns: E_R (eV), dRdE (events/kg/year/eV)",
    ]


# ----------------------------------------------------------------------------
# _cli
#   Command-line front end: compute one spin-independent WIMP-nucleus rate table, quenched to E_ee when --quenching_model is given and on the raw E_R axis otherwise, and write it to --out_csv.
# ----------------------------------------------------------------------------
def _cli():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--target_nucleus", default="Si28")
    ap.add_argument("--A", type=float, required=True)
    ap.add_argument("--mchi_MeV", type=float, required=True)
    ap.add_argument("--sigma_n_cm2", type=float, required=True)
    ap.add_argument("--mediator", default="heavy")
    ap.add_argument("--nr_Emin_keV", type=float, default=0.001)
    ap.add_argument("--nr_Emax_keV", type=float, default=30.0)
    ap.add_argument("--nr_nbins", type=int, default=3000)
    ap.add_argument("--rho_chi_gev_cm3", type=float, default=None)
    ap.add_argument("--v0_kms", type=float, default=None)
    ap.add_argument("--vE_kms", type=float, default=None)
    ap.add_argument("--vesc_kms", type=float, default=None)
    ap.add_argument("--quenching_model", default=None, choices=list(QUENCHING_MODELS),
                     help="If set, output is quenched to E_ee instead of raw E_R.")
    ap.add_argument("--ee_Emin_eV", type=float, default=0.0)
    ap.add_argument("--ee_Emax_eV", type=float, default=8000.0)
    ap.add_argument("--ee_nbins", type=int, default=800)
    ap.add_argument("--out_csv", required=True)
    args = ap.parse_args()

    common_kwargs = dict(
        A=args.A,
        mchi_MeV=args.mchi_MeV,
        sigma_n_cm2=args.sigma_n_cm2,
        target_nucleus=args.target_nucleus,
        mediator=args.mediator,
        nr_Emin_keV=args.nr_Emin_keV,
        nr_Emax_keV=args.nr_Emax_keV,
        nr_nbins=args.nr_nbins,
        rho_chi_gev_cm3=args.rho_chi_gev_cm3,
        v0_kms=args.v0_kms,
        vE_kms=args.vE_kms,
        vesc_kms=args.vesc_kms,
    )

    if args.quenching_model:
        res = compute_dRdE_ee(
            **common_kwargs,
            ee_Emin_eV=args.ee_Emin_eV,
            ee_Emax_eV=args.ee_Emax_eV,
            ee_nbins=args.ee_nbins,
            quenching_model=args.quenching_model,
        )
        CIO.write_csv_generic(args.out_csv, res["E_eV"], res["dRdE_kg_year_eV"], _header_lines_ee(res["meta"]))
    else:
        res = compute_dRdE(**common_kwargs)
        CIO.write_csv_generic(args.out_csv, res["E_eV"], res["dRdE_kg_year_eV"], _header_lines(res["meta"]))
    print(f"[wimp_nucleon] wrote {args.out_csv}  (N={len(res['E_eV'])})")


if __name__ == "__main__":
    _cli()
