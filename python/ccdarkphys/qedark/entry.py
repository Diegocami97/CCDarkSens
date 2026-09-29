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
#  entry.py -- QEDark entry point for DM-electron scattering (silicon). It is
#  built only from the QEDark source files: the constants, the SHM halo and
#  Si_f2.txt.
# ============================================================================

"""
QEDark entry point for DM–electron scattering (Silicon).
Built ONLY from your provided source: constants, SHM halo, Si_f2.txt.

Public:
    compute_dRdE(material, mediator, mchi_eV, sigma_e_cm2, halo,
                 band_gap_eV=1.2, eh_pair_eV=3.8, binsize_eV=0.1)
      -> dict(E_eV, dRdE_kg_year_eV)

CLI example:
    PYTHONPATH=python python3 -m ccdarkphys.qedark.entry \
      --material Si --mediator heavy \
      --mchi_MeV 10 --sigma_e_cm2 1e-37 \
      --v0_kms 220 --vE_kms 232 --vesc_kms 544 \
      --out_csv data/qedark_rates/Si/heavy/dRdE_Si_heavy_m10.000000_s1e-37.csv
"""

from __future__ import annotations
import argparse
import numpy as np

from ccdarkphys.common import constants as QEC
from ccdarkphys.common import halo as HALO
from ccdarkphys.common import io as CIO
from ccdarkphys.common.mediator_map import MEDIATOR_TO_FDM_INDEX

_MEDIATOR_TO_INDEX = MEDIATOR_TO_FDM_INDEX

def _load_si_table_as_notebook(nE: int, nq: int) -> np.ndarray:
    """
    Load Si_f2.txt exactly like the notebook:
      fcrys_Si = transpose( resize(loadtxt(..., skiprows=1), (nE, nq)) )
    Returns array of shape (nq, nE).
    """
    path = CIO.data_path(__file__, "data", "Si_f2.txt")
    raw = np.loadtxt(path, skiprows=1)
    fcrys = np.transpose(np.resize(raw, (nE, nq)))
    return fcrys

def compute_dRdE(material: str,
                 mediator: str,
                 mchi_eV: float,
                 sigma_e_cm2: float,
                 halo: dict,
                 band_gap_eV: float = 1.2,
                 eh_pair_eV: float = 3.8,
                 binsize_eV: float = 0.1) -> dict:
    """
    Returns:
      dict with 'E_eV' and 'dRdE_kg_year_eV' (events / kg / year / eV)
    """
    if material.lower() not in ("si", "silicon"):
        raise NotImplementedError("Only Silicon is supported (material='Si').")
    nFDM = _MEDIATOR_TO_INDEX.get(mediator, None)
    if nFDM is None:
        raise ValueError("mediator: heavy/massive/0 or light/massless/2")

    # Halo velocities → cm/s
    if "v0_cm_s" in halo:
        v0_cm_s = float(halo["v0_cm_s"]); vE_cm_s = float(halo["vE_cm_s"]); vesc_cm_s = float(halo["vesc_cm_s"])
    else:
        v0_cm_s = HALO.kms_to_cms(halo["v0_kms"])
        vE_cm_s = HALO.kms_to_cms(halo["vE_kms"])
        vesc_cm_s = HALO.kms_to_cms(halo["vesc_kms"])

    # ---------------- NOTEBOOK-CORRECT RATE KERNEL ----------------
    # Notebook global parameters
    nE = 500
    nq = 900

    # Load/reshape fcrys exactly like the notebook → shape (nq, nE)
    fcrys_Si = _load_si_table_as_notebook(nE=nE, nq=nq)

    # Notebook constants
    dQ = 0.02 * QEC.alpha * QEC.me_eV   # eV
    dE = binsize_eV                             # eV (native energy step)
    wk = 2.0 / 137.0                     # ~ 2*alpha

    # materials[mat] = [Mcell(kg), Eprefactor, Egap(eV), epsilon(eV), fcrys]
    Mcell_Si    = 2.0 * 28.0855 * QEC.amu_kg
    Eprefactor  = 2.0
    Egap_Si     = band_gap_eV
    epsilon_Si  = eh_pair_eV
    fcrys_Si    = (wk / 4.0) * fcrys_Si  # NOTEBOOK: remove wk/4 only if fcrys was regenerated locally

    materials = {
        "Si": [Mcell_Si, Eprefactor, Egap_Si, epsilon_Si, fcrys_Si],
    }

    def FDM(q_eV: float, n: int) -> float:
        """
        DM form factor:
          n = 0: 1
          n = 1: ~(alpha*me/q)^1 (not used here but supported)
          n = 2: ~(alpha*me/q)^2
        """
        if n == 0:
            return 1.0
        qsafe = max(q_eV, 1e-12)
        return (QEC.alpha * QEC.me_eV / qsafe) ** n

    def mu_Xe(mX_eV: float) -> float:
        """DM–electron reduced mass in eV."""
        return (mX_eV * QEC.me_eV) / (mX_eV + QEC.me_eV)

    # Halo params vector as in the notebook helpers
    vparams = [v0_cm_s, vE_cm_s, vesc_cm_s]  # not passed explicitly; kept for reference

    # Precompute η(vmin) on a 1D grid once per spectrum (avoids ~nE*nq nquad calls)
    vmax_cut = (vesc_cm_s + vE_cm_s) * 1.1
    vmin_min_cm_s = 1.0e5
    vmin_max_cm_s = vmax_cut * 1.01
    n_eta_grid = 500  # larger → less interpolation error in η(vmin)
    vmin_grid = np.logspace(
        np.log10(max(vmin_min_cm_s, 1.0)),
        np.log10(max(vmin_max_cm_s, vmin_min_cm_s * 1.1)),
        num=n_eta_grid,
        dtype=float,
    )
    eta_grid = HALO.eta_shm_numeric(vmin_grid, v0_cm_s, vE_cm_s, vesc_cm_s)

    # q grid (1..nq) in eV
    q_arr = np.arange(1, nq + 1, dtype=float) * dQ
    q_safe = np.maximum(q_arr, 1e-12)
    if nFDM == 0:
        FDM_sq_arr = np.ones_like(q_arr)
    else:
        FDM_sq_arr = (QEC.alpha * QEC.me_eV / q_safe) ** (2 * nFDM)

    Mcell, Epref, Egap, epsilon, f_arr = materials["Si"]
    prefactor = (QEC.ccms**2) * QEC.sec_per_year * (QEC.rho_X_eVcm3 / mchi_eV) * (1.0 / Mcell) \
                * QEC.alpha * (QEC.me_eV**2) / (mu_Xe(mchi_eV)**2)

    def rate_at_E(Ee: float) -> float:
        """Vectorized over q: one η interpolation and one sum per E."""
        if Ee < Egap:
            return 0.0
        Ei = int(np.floor(Ee * 10.0))
        if Ei < 1 or Ei > nE:
            return 0.0
        vmin_arr = (q_arr / (2.0 * mchi_eV) + Ee / q_safe) * QEC.ccms
        mask_skip = vmin_arr > vmax_cut
        eta_arr = np.interp(vmin_arr, vmin_grid, eta_grid)
        eta_arr[mask_skip] = 0.0
        fcrys_col = f_arr[:, Ei - 1]
        integrand = Epref * (1.0 / q_safe) * eta_arr * FDM_sq_arr * fcrys_col
        integrand[mask_skip] = 0.0
        return float(prefactor * np.sum(integrand))

    # Native energy grid (0.1 eV) starting at Egap
    E_min = materials["Si"][2]
    E_eV = E_min + np.arange(nE, dtype=float) * dE

    # Notebook dRdE returns integrated rate per dE bin [events/(kg·year)]. True dR/dE = value/dE → events/(kg·year·eV).
    dRdE_kg_year = np.array([rate_at_E(E) for E in E_eV], dtype=float)
    dRdE_kg_year *= float(sigma_e_cm2)
    dRdE_kg_year_eV = dRdE_kg_year / dE

    # Derived units for g/day if needed elsewhere
    dRdE_g_day_eV = dRdE_kg_year_eV / 1000.0 / 365.25
    dRdE_g_day_eV[~np.isfinite(dRdE_g_day_eV)] = 0.0

    # ---------------------------------------------------------------

    # Minimal metadata for CSV header
    si_path = CIO.data_path(__file__, "data", "Si_f2.txt")
    return {"E_eV": E_eV, "dRdE_kg_year_eV": dRdE_kg_year_eV, "meta": {
        "material": material, "mediator": mediator,
        "table_path": si_path, "table_sha1": CIO.sha1sum(si_path),
        "v0_cm_s": v0_cm_s, "vE_cm_s": vE_cm_s, "vesc_cm_s": vesc_cm_s,
        "mchi_eV": float(mchi_eV), "sigma_e_cm2": float(sigma_e_cm2)
    }}

# ----------------------------------------------------------------------------
# _cli
#   Command-line front end: compute one QEDark rate table for the given mediator, mass and cross section (default halo v0 = 220, vE = 232, vesc = 544 km/s) and write it to --out_csv.
# ----------------------------------------------------------------------------
def _cli():
    ap = argparse.ArgumentParser()
    ap.add_argument("--material", default="Si")
    ap.add_argument("--mediator", required=True, choices=["heavy","massive","light","massless"])
    ap.add_argument("--mchi_MeV", type=float, required=True)
    ap.add_argument("--sigma_e_cm2", type=float, required=True)
    ap.add_argument("--v0_kms", type=float, default=220.0)
    ap.add_argument("--vE_kms", type=float, default=232.0)
    ap.add_argument("--vesc_kms", type=float, default=544.0)
    ap.add_argument("--band_gap_eV", type=float, default=1.2)
    ap.add_argument("--eh_pair_eV", type=float, default=3.8)
    ap.add_argument("--binsize_eV", type=float, default=0.1)
    ap.add_argument("--out_csv", required=True)
    args = ap.parse_args()

    res = compute_dRdE(
        material=args.material,
        mediator=args.mediator,
        mchi_eV=args.mchi_MeV*1.0e6,
        sigma_e_cm2=args.sigma_e_cm2,
        halo={"v0_kms": args.v0_kms, "vE_kms": args.vE_kms, "vesc_kms": args.vesc_kms},
        band_gap_eV=args.band_gap_eV,
        eh_pair_eV=args.eh_pair_eV,
        binsize_eV=args.binsize_eV,
    )
    E, R, meta = res["E_eV"], res["dRdE_kg_year_eV"], res["meta"]
    CIO.write_csv(args.out_csv, E, R, meta)
    print(f"[qedark] wrote {args.out_csv}  (N={len(E)})")

if __name__ == "__main__":
    _cli()
