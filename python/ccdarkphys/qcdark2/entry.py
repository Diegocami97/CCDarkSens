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
#  entry.py -- QCDark2 entry point for DM-electron scattering using
#  dielectric-function inputs.
# ============================================================================

"""
QCDark2 entry point for DM-electron scattering using dielectric-function inputs.

Expected input file is a QCDark2 dielectric HDF5 (e.g. ``Si_comp.h5``) with
datasets ``epsilon``, ``q``, ``E`` and attrs including ``M_cell``, ``V_cell``,
``dE``. Rates are computed through ``qcdark2.dark_matter_rates.get_dR_dE``.

Public:
    compute_dRdE(..., epsilon_h5=None, screening="Lindhard")
"""

from __future__ import annotations

import argparse
import os

import numpy as np

from ccdarkphys.common import constants as QEC
from ccdarkphys.common import halo as HALO
from ccdarkphys.common import io as CIO
from ccdarkphys.common.mediator_map import MEDIATOR_TO_FDM_INDEX

_ENV_EPSILON = "CCDARK_SENS_QCDARK2_EPSILON"


# ----------------------------------------------------------------------------
# _resolve_epsilon_path
#   Path of the QCDark2 dielectric-function HDF5: the explicit argument, else the environment variable; raises FileNotFoundError if neither is given.
# ----------------------------------------------------------------------------
def _resolve_epsilon_path(explicit: str | None) -> str:
    if explicit:
        return os.path.abspath(os.path.expanduser(explicit))
    env = os.environ.get(_ENV_EPSILON)
    if env:
        return os.path.abspath(os.path.expanduser(env.strip()))
    raise FileNotFoundError(
        f"QCDark2 epsilon file is required. Pass epsilon_h5=... or set {_ENV_EPSILON}."
    )


def compute_dRdE(
    material: str,
    mediator: str,
    mchi_eV: float,
    sigma_e_cm2: float,
    halo: dict,
    band_gap_eV: float = 1.2,
    eh_pair_eV: float = 3.8,
    binsize_eV: float = 0.1,
    *,
    epsilon_h5: str | None = None,
    screening: str = "Lindhard",
) -> dict:
    """
    Same outward contract as other CCDarkSens backends.

    Returns:
      dict with ``E_eV``, ``dRdE_kg_year_eV`` (events / kg / year / eV), ``meta``.
    """
    if material.lower() not in ("si", "silicon"):
        raise NotImplementedError("Only Silicon is currently supported for QCDark2 backend.")
    # mediator can be a named string ("heavy"/"light") or a float/numeric string mA' in eV
    try:
        mA_eV = float(mediator)
        q2_mediator: str | float = mA_eV   # pass numeric mA' directly to QCDark2
    except (TypeError, ValueError):
        if mediator not in MEDIATOR_TO_FDM_INDEX:
            raise ValueError("mediator: 'heavy'/'light' or a float mA' in eV (e.g. 5000.0)")
        q2_mediator = "heavy" if MEDIATOR_TO_FDM_INDEX[mediator] == 0 else "light"

    path = _resolve_epsilon_path(epsilon_h5)
    if not os.path.isfile(path):
        raise FileNotFoundError(path)

    try:
        import qcdark2.dark_matter_rates as dm
    except ImportError as exc:
        raise ImportError(
            "qcdark2 package is required for backend='qcdark2' "
            "(install in this Python env with `pip install qcdark2`)."
        ) from exc

    if "v0_cm_s" in halo:
        v0_cm_s = float(halo["v0_cm_s"])
        v_e_cm_s = float(halo["vE_cm_s"])
        v_esc_cm_s = float(halo["vesc_cm_s"])
    else:
        v0_cm_s = HALO.kms_to_cms(halo["v0_kms"])
        v_e_cm_s = HALO.kms_to_cms(halo["vE_kms"])
        v_esc_cm_s = HALO.kms_to_cms(halo["vesc_kms"])

    astro_model = {
        "v0": v0_cm_s / 1.0e5,
        "vEarth": v_e_cm_s / 1.0e5,
        "vEscape": v_esc_cm_s / 1.0e5,
        "rhoX": float(QEC.rho_X_eVcm3),
        "sigma_e": float(sigma_e_cm2),
    }

    eps = dm.load_epsilon(path)
    d_rd_e, e_ev = dm.get_dR_dE(
        eps,
        m_X=float(mchi_eV),
        mediator=q2_mediator,
        astro_model=astro_model,
        screening=str(screening),
        velocity_dist="MB",
    )

    # Keep binsize consistency check for downstream assumptions.
    tol = max(1e-12, abs(float(eps.dE)) * 1e-6)
    if abs(float(binsize_eV) - float(eps.dE)) > tol:
        raise ValueError(
            f"binsize_eV ({binsize_eV}) must match epsilon file dE={float(eps.dE)} eV."
        )

    scissor_attr = None
    try:
        import h5py

        with h5py.File(path, "r") as h5:
            if "scissor_bandgap_eV" in h5.attrs:
                scissor_attr = float(h5.attrs["scissor_bandgap_eV"])
    except Exception:
        pass

    meta = {
        "material": material,
        "mediator": mediator,
        **({"mA_eV": float(mediator)} if isinstance(q2_mediator, float) else {}),
        "table_path": path,
        "table_sha1": CIO.sha1sum(path),
        "v0_cm_s": v0_cm_s,
        "vE_cm_s": v_e_cm_s,
        "vesc_cm_s": v_esc_cm_s,
        "mchi_eV": float(mchi_eV),
        "sigma_e_cm2": float(sigma_e_cm2),
        "eh_pair_eV": float(eh_pair_eV),
        "band_gap_eV": float(band_gap_eV),
        "scissor_bandgap_eV": scissor_attr,
        "screening": str(screening),
        "epsilon_dE_eV": float(eps.dE),
    }
    return {
        "E_eV": np.asarray(e_ev, dtype=float),
        "dRdE_kg_year_eV": np.asarray(d_rd_e, dtype=float),
        "meta": meta,
    }


# ----------------------------------------------------------------------------
# _cli
#   Command-line front end: compute one QCDark2 rate table (mediator 'heavy', 'light' or a mediator mass in eV; default halo v0 = 238, vE = 263, vesc = 544 km/s) and write it to --out_csv.
# ----------------------------------------------------------------------------
def _cli():
    ap = argparse.ArgumentParser()
    ap.add_argument("--material", default="Si")
    ap.add_argument(
        "--mediator", required=True,
        help="'heavy', 'light', or a float mA' in eV (e.g. 5000.0 for 5 keV)"
    )
    ap.add_argument("--mchi_MeV", type=float, required=True)
    ap.add_argument("--sigma_e_cm2", type=float, required=True)
    ap.add_argument("--v0_kms", type=float, default=238.0)
    ap.add_argument("--vE_kms", type=float, default=263.0)
    ap.add_argument("--vesc_kms", type=float, default=544.0)
    ap.add_argument("--band_gap_eV", type=float, default=1.2)
    ap.add_argument("--eh_pair_eV", type=float, default=3.8)
    ap.add_argument("--binsize_eV", type=float, default=0.1)
    ap.add_argument("--epsilon_h5", required=True)
    ap.add_argument("--screening", default="Lindhard")
    ap.add_argument("--out_csv", required=True)
    args = ap.parse_args()

    res = compute_dRdE(
        material=args.material,
        mediator=args.mediator,
        mchi_eV=args.mchi_MeV * 1.0e6,
        sigma_e_cm2=args.sigma_e_cm2,
        halo={"v0_kms": args.v0_kms, "vE_kms": args.vE_kms, "vesc_kms": args.vesc_kms},
        band_gap_eV=args.band_gap_eV,
        eh_pair_eV=args.eh_pair_eV,
        binsize_eV=args.binsize_eV,
        epsilon_h5=args.epsilon_h5,
        screening=args.screening,
    )
    CIO.write_csv(args.out_csv, res["E_eV"], res["dRdE_kg_year_eV"], res["meta"], entry="QCDark2")
    print(f"[qcdark2] wrote {args.out_csv}  (N={len(res['E_eV'])})")


if __name__ == "__main__":
    _cli()
