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

_MEDIATOR_TO_INDEX = {
    "heavy": 0, "massive": 0, "0": 0, 0: 0,
    "light": 2, "massless": 2, "2": 2, 2: 2,
}

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
    fcrys_Si    = (wk / 4.0) * fcrys_Si  # NOTEBOOK: remove wk/4 only if you regenerated fcrys yourself

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

    # dRdE(material, mX, Ee, FDMn, 'shm', params) → events / kg / year at that Ee
    def dRdE_notebook(material: str, mX: float, Ee: float, FDMn: int) -> float:
        if Ee < materials[material][2]:   # Egap
            return 0.0

        qunit = dQ
        Mcell, Epref, Egap, epsilon, f_arr = materials[material]

        # Energy bin index in 0.1 eV steps (Ei-1 used for array indexing)
        Ei = int(np.floor(Ee * 10.0))
        if Ei < 1 or Ei > nE:
            return 0.0

        # Prefactor [(kg-year)^-1] when sigma_e = 1 cm^2; we multiply sigma_e outside
        prefactor = (QEC.ccms**2) * QEC.sec_per_year * (QEC.rho_X_eVcm3 / mX) * (1.0 / Mcell) \
                    * QEC.alpha * (QEC.me_eV**2) / (mu_Xe(mX)**2)

        # Sum over q grid
        acc = 0.0
        for qi in range(1, nq + 1):
            q = qi * qunit
            # vmin = (q/(2 mX) + Ee/q) * c   (in cm/s)
            qsafe = max(q, 1e-12)
            vmin = (q / (2.0 * mX) + Ee / qsafe) * QEC.ccms

            # rough kinematic cutoff as notebook
            if vmin > (vesc_cm_s + vE_cm_s) * 1.1:
                continue

            # SHM halo η (cm/s)^-1
            eta = HALO.eta_shm_analytic(np.array([vmin]), v0_cm_s, vE_cm_s, vesc_cm_s)[0]
            # print(f"[debug] vmin={vmin:.2e} cm/s → eta={eta:.2e} (q={q:.2e} eV, Ei={Ei})")
            # eta = HALO.eta_shm_numeric(np.array([vmin]), v0_cm_s, vE_cm_s, vesc_cm_s)[0]

            # array_[qi-1] = Eprefactor * (1/q) * eta * FDM(q,n)^2 * fcrys[qi-1, Ei-1]
            acc += Epref * (1.0 / qsafe) * eta * (FDM(q, FDMn) ** 2) * f_arr[qi - 1, Ei - 1]

        # This is for sigma_e = 1 cm^2; scale by sigma_e outside
        return prefactor * acc  # [(kg-year)^-1 at Ee]

    # Native energy grid (0.1 eV) starting at Egap
    E_min = materials["Si"][2]
    E_eV = E_min + np.arange(nE, dtype=float) * dE

    # Evaluate notebook rate (sigma_e = 1), then scale by sigma_e
    dRdE_kg_year = np.array([dRdE_notebook("Si", mchi_eV, E, nFDM) for E in E_eV], dtype=float)
    # print(f"[qedark] computed dRdE for mchi={mchi_eV/1.0e6:.6f} MeV, sigma_e=1 cm²")
    # print(f"dRdE_kg_year: {dRdE_kg_year}")
    dRdE_kg_year *= float(sigma_e_cm2)
    # print(f"dRdE_kg_year: {dRdE_kg_year}")

    dRdE_kg_year_eV = dRdE_kg_year / dE  # Convert to events / kg / year / eV


    # Convert to events / g / day / eV
    dRdE_g_day_eV = dRdE_kg_year / 1000.0 / 365.25 / dE
    dRdE_g_day_eV[~np.isfinite(dRdE_g_day_eV)] = 0.0
    # print(f"dRdE_g_day_eV: {dRdE_kg_year_eV}")

    # ---------------------------------------------------------------

    # Minimal metadata for CSV header
    si_path = CIO.data_path(__file__, "data", "Si_f2.txt")
    return {"E_eV": E_eV, "dRdE_kg_year_eV": dRdE_kg_year_eV, "meta": {
        "material": material, "mediator": mediator,
        "table_path": si_path, "table_sha1": CIO.sha1sum(si_path),
        "v0_cm_s": v0_cm_s, "vE_cm_s": vE_cm_s, "vesc_cm_s": vesc_cm_s,
        "mchi_eV": float(mchi_eV), "sigma_e_cm2": float(sigma_e_cm2)
    }}

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
